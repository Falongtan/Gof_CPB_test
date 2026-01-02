

###generation nonlinear
xpand_features <- function(X) {
  # make sure the input is a data.frame
  if (!is.data.frame(X)) {
    X <- as.data.frame(X)
  }
  
  # extract nocl and colnames 
  col_names <- names(X)
  n_cols <- ncol(X)
  
  # create square terms
  squares <- lapply(col_names, function(col) {
    X[[col]]^2
  })
  names(squares) <- paste0(col_names, "_sq")
  squares_df <- as.data.frame(squares)
  
  # create interaction terms
  cross_products <- list()
  cp_names <- character()
  
  for (i in 1:(n_cols-1)) {
    for (j in (i+1):n_cols) {
      col_i <- col_names[i]
      col_j <- col_names[j]
      cp_name <- paste0(col_i, "_x_", col_j)
      
      cross_products[[cp_name]] <- X[[col_i]] * X[[col_j]]
      cp_names <- c(cp_names, cp_name)
    }
  }
  
  cross_products_df <- as.data.frame(cross_products)
  names(cross_products_df) <- cp_names
  
  # combine X with square terms and interaction terms
  result <- cbind(X, squares_df, cross_products_df)
  return(result)
}


library(xlsx)
library('Hmisc')
library('lmtest')
library('glmnet')
library(psych)
library(MASS)
library(matrixcalc)
library(LassoSIR)
library(VariableScreening)
library(msda)
library(energy)
library(expm)
library(gee)
library(gsl)
library(car)
library(pracma)
library(caret)


linear_data <- read.xlsx("E:/桌面/CommViolPredUnnormalizedData.xlsx",sheetIndex=1,header=T)
linear_data_nona <- linear_data[which(rowSums(is.na(linear_data))==0),]

linear_data_final <- matrix(nrow=0,ncol=120)
for (i in 1:2197){
  if(any(linear_data_nona[i,]=="?")){
    linear_data_final <- linear_data_final
  }else{
    linear_data_final <- rbind(linear_data_final,linear_data_nona[i,])
  }
}
Y_matrix <- linear_data_final[,103:120]
colnames(Y_matrix) <- c("Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y19","Y10","Y11","Y12","Y13","Y14","Y15","Y16","Y17","Y18")
for(i in 1:8){
  o <- c(3,4,9,10,15,16,17,18)[i]
  Y_matrix[,o] <- as.numeric(Y_matrix[,o])
}
###X,Y
Y <- Y_matrix$Y17
Y <- as.numeric(Y)
Y <- scale(Y)
X_0 <- linear_data_final[,-c(103:120)]
X <- scale(X_0)


n <- nrow(X)
p <- ncol(X)
set.seed(222)
index_x<-sample(1:n, floor(n/2),replace = FALSE)  
index_x<-sort(index_x)
p <- ncol(X)
x1 <- X[-index_x,] 
y1 <- as.matrix(Y[-index_x]) 

x2 <- X[index_x,]
y2 <- as.matrix(Y[index_x])

n1 <- nrow(x1)
n2 <- nrow(x2)

# residuals and projections based on x1 y1 and dimension reduction
# beta_hat projection and residuals based on lasso and data x1 y1
lasso_model1 <- cv.glmnet(as.matrix(x1), y1, family="gaussian", intercept = F)
lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
index_beta1_non0 <- seq(1:p)[as.numeric(lasso_beta1) != 0]
len_beta1_non0 <- length(index_beta1_non0)
if(len_beta1_non0 == 0){
  U1 <- y1 - mean(y1)                                                        # residual based on x1 y1
}else{
  x1_sec <- x1[,index_beta1_non0]                                            # second estimation
  # cor_matrix1 <- cor(x1_sec, use = "pairwise.complete.obs")
  # high_cor_vars1 <- findCorrelation(cor_matrix1, cutoff = 0.8, verbose = TRUE)
  # vars_to_remove1 <- colnames(cor_matrix1)[high_cor_vars1]
  # x1_sec <- x1_sec[,-which(colnames(x1_sec)%in%vars_to_remove1)]
  sec_model1 <- glm(y1~x1_sec-1, family = gaussian)
  pred1 = predict(sec_model1, newx = x1_sec, type="response")
  U1 <- y1-pred1   # residual based on x1 y1
}

# projections based on screening and lassosir or seas for x1 U1
screen_num1 <- floor(n1/log(n1))
if(p <= screen_num1){
  # projection based on lassosir without screening
  sir_U1 <- LassoSIR(x1, U1, H=10, choosing.d="automatic", solution.path=FALSE, 
                     categorical=FALSE, nfolds=5, screening=FALSE)
  sir_Ubeta1<- sir_U1$beta
  sdr_Upro1 <- sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
}else{
  #screening
  rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                     # screening method based on distance correlation
  index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                       # chosen index of X
  x1_Uscreen <- x1[,index_U1]
  
  # projection based on lassosir after screening
  sir_U1 <- LassoSIR(x1_Uscreen, U1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_Ubeta1<- sir_U1$beta
  sdr_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
  sdr_Upro1[index_U1,] <- sir_Ubeta1
  sdr_Upro1 <- sdr_Upro1/sqrt(colSums(sdr_Upro1^2))
}

# projections based on screening and lassosir or seas for x1 y1
if(p <= screen_num1){
  # laosssir projection without screening
  sir_y1 <- LassoSIR(x1, y1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_ybeta1<- sir_y1$beta
  sdr_ypro1 <- sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
}else{
  #screening
  rank_y1 <- screenIID(X=x1, Y=y1, method = "DC-SIS")                     # screening method 
  index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                       # chosen index of X
  x1_yscreen <- x1[,index_y1]
  
  # projection based on lassosir after screening
  sir_y1 <- LassoSIR(x1_yscreen, y1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_ybeta1<- sir_y1$beta
  sdr_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
  sdr_ypro1[index_y1,] <- sir_ybeta1
  sdr_ypro1 <- sdr_ypro1/sqrt(colSums(sdr_ypro1^2))
}

# residuals and projections based on x2 y2 and dimension reduction
# residuals based on lasso and data x2 y2
lasso_model2 <- cv.glmnet(as.matrix(x2), y2, family="gaussian", intercept = F)
lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]                   # estimated beta0 without intercept 
index_beta2_non0 <- seq(1:p)[as.numeric(lasso_beta2)!=0]                  # index of beta_n with nonzero components
len_beta2_non0 <- length(index_beta2_non0)
if(len_beta2_non0 == 0){
  U2 <- y2-mean(y2)                                                       # residual based on x2 y2
}else{
  x2_sec <- x2[,index_beta2_non0]
  # cor_matrix2 <- cor(x2_sec, use = "pairwise.complete.obs")
  # high_cor_vars2 <- findCorrelation(cor_matrix2, cutoff = 0.8, verbose = TRUE)
  # vars_to_remove2 <- colnames(cor_matrix2)[high_cor_vars2]
  # x2_sec <- x2_sec[,-which(colnames(x2_sec)%in%vars_to_remove2)]
  sec_model2 = glm(y2~x2_sec-1,family = gaussian)                           # second estimation
  pred2 = predict(sec_model2, newx = x2_sec, type="response")
  U2 <- y2-pred2   # residual based on x1 y1                                                  # residual based on x2 y2
}

#sir projections based on screening and x2 U2
screen_num2 <- floor(n2/log(n2))
if(p <= screen_num2){
  # lassosir projection without screen
  sir_U2 <- LassoSIR(x2, U2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_Ubeta2<- sir_U2$beta
  sdr_Upro2 <- sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))                     # sdr projection based on lassosir and x2 U2
}else{
  #screening
  rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                     # screening method based on distance correlation
  index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                       # chosen index of X
  x2_Uscreen <- x2[,index_U2]
  
  # projection based on lassosir after screening
  sir_U2 <- LassoSIR(x2_Uscreen, U2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_Ubeta2 <- sir_U2$beta
  sdr_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
  sdr_Upro2[index_U2,] <- sir_Ubeta2
  sdr_Upro2 <- sdr_Upro2/sqrt(colSums(sdr_Upro2^2))                       # sdr projection based on lassosir and x2 U2
}

#sir projections based on screening and x2 y2
if(p <= screen_num2){
  # Lassosir projection without screen
  sir_y2 <- LassoSIR(x2, y2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_ybeta2<- sir_y2$beta
  sdr_ypro2 <- sir_ybeta2/sqrt(colSums(sir_ybeta2^2))                     # sdr projection based on lassosir and x2 y2
}else{ 
  #screening
  rank_y2 <- screenIID(X=x2, Y=y2, method = "DC-SIS")                     # screening method based on distance correlation
  index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                       # chosen index of X
  x2_yscreen <- x2[,index_y2]
  
  # projection based on lassosir after screening
  sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                          nfolds=5, screening=FALSE)
  sir_ybeta2 <- sir_ylasso2$beta
  sdr_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
  sdr_ypro2[index_y2,] <- sir_ybeta2
  sdr_ypro2 <- sdr_ypro2/sqrt(colSums(sdr_ypro2^2))                       # sdr projection based on lassosir and x2 y2
}

#construct test statistic based on x1 y1 and projections based on x2 y2
pro2 <- cbind(sdr_Upro2, sdr_ypro2)              # projections based on x1 y1
pro2_num <- ncol(pro2)                                                    # the number of projections

h1 <- (n1^(-2/9))  #bandwidth 

pval_matrix1 <- matrix(nrow=1, ncol=pro2_num) # p-value matrix
errormat1 <- U1%*%t(U1)  #residual matrix based on x1 y1
for(q in 1:pro2_num){
  x_pro1 <- x1%*%pro2[,q]
  x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1  # kernel function matrix
  # kermat1 <-(1/sqrt(2*pi))*exp(-(x_pro1_mat^2)/2)                           # Gaussian kernel
  indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
  kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                # Epanechnikov kernel
  #test statistics
  Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
  pval1 <- 1-pnorm(Tn1)
  pval_matrix1[,q] <- pval1
}

#construct test statistics based on x2 y2 and projections based on x1 y1
pro1 <- cbind(sdr_Upro1, sdr_ypro1)                           # projections based on x1 y1
pro1_num <- ncol(pro1)                                                    # the number of projections

h2 <- (n2^(-2/9))                                                      # bandwidth 

pval_matrix2 <- matrix(nrow=1, ncol=pro1_num) # p-value matrix
errormat2 <- U2%*%t(U2)  #residual matrix based on x2 y2
for(l in 1:pro1_num){
  x_pro2 <- x2%*%pro1[,l]
  x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2  #kernel function matrix
  #kermat2 <-(1/sqrt(2*pi))*exp(-(x_pro2_mat^2)/2)                            # Gaussian kernel
  indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0) 
  kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                # Epanechnikov kernel 
  #test statistics
  Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
  pval2 <- 1-pnorm(Tn2)
  pval_matrix2[,l] <- pval2
}

pval_matrix <- cbind(pval_matrix1, pval_matrix2)
pval_num <- pro1_num + pro2_num

pval_fisher1 <- 1- pchisq(-2*sum(log(pval_matrix1[1:pro2_num])),df=2*pro2_num)
pval_fisher2 <- 1- pchisq(-2*sum(log(pval_matrix2[1:pro1_num])),df=2*pro1_num)
pval_fisher_mat <- cbind(pval_fisher1,pval_fisher2)
pval_cauchy_fisher <- 1- pcauchy(mean(tan((0.5-pval_fisher_mat)*pi))) 


pval_min1 <- 1 - (1 - min(pval_matrix1))^pro2_num
pval_min2 <- 1 - (1 - min(pval_matrix2))^pro1_num
pval_min_mat <- cbind(1 - (1 - min(pval_matrix1))^pro2_num, 1 - (1 - min(pval_matrix2))^pro1_num) 
pval_cauchy_min <- 1- pcauchy(mean(tan((0.5-pval_min_mat)*pi))) 
pval_matrix;pval_fisher1;pval_fisher2;pval_cauchy_fisher;pval_min1;pval_min2;pval_cauchy_min

# pval_cauchy1 <- 1- pcauchy(mean(tan((0.5-pval_matrix1[1:pro2_num])*pi)))  # cauchy combination based on x1 y1
# pval_cauchy2 <- 1- pcauchy(mean(tan((0.5-pval_matrix2[1:pro1_num])*pi)))  # cauchy combination based on x1 y1
# pval_cauchy  <- 1- pcauchy(mean(tan((0.5-pval_matrix[1:pval_num])*pi)))   # cauchy combination based on x y
# pval_cauchy1;pval_cauchy2;pval_cauchy


######scatterplot_1
# Y <- scale(Y)
p <- ncol(X)
lasso_model <- cv.glmnet(as.matrix(X), Y, family="gaussian", intercept = F)
beta_lasso <- coef(lasso_model, s = "lambda.min")[-1]
index_beta_non0 <- seq(1:p)[as.numeric(beta_lasso) != 0]
x_sec <- X[,index_beta_non0]                                          
sec_model = glm(Y~x_sec-1,family = gaussian)                           

pred = predict(sec_model, newx = x_sec, type="response")
U <- Y-pred  
xbeta <- x_sec%*%sec_model$coefficients   #[-1]+sec_model$coefficients[1]
#####plot_1

plot_range <- range(c(Y, xbeta))
plot(xbeta, Y, 
     xlab = "",  
     ylab = "Y",
     pch = 16,
     ylim = plot_range,
     cex.lab = 2.5,
     cex.axis = 2,
     mgp = c(2, 0.8, 0)  
)
mtext(expression(hat(beta)^T*X), side = 1,  
      line = 3,  
      cex = 2.5)   




###############Polynomial Regression


set.seed(222)
p <- ncol(X)
lasso_model <- cv.glmnet(as.matrix(X), Y, family="gaussian", intercept = F)
beta_lasso <- coef(lasso_model, s = "lambda.min")[-1]
index_beta_non0 <- seq(1:p)[as.numeric(beta_lasso) != 0]
x_sec <- X[,index_beta_non0]   

X_quad_inter <- xpand_features(x_sec)
p <- ncol(X_quad_inter)
new_names <- paste("col", 1:p, sep = "_")
colnames(X_quad_inter) <- new_names
X_quad_inter[,(ncol(x_sec)+1):ncol(X_quad_inter)] <- scale(X_quad_inter[,(ncol(x_sec)+1):ncol(X_quad_inter)])

x1 <- X_quad_inter[-index_x,] 
x2 <- X_quad_inter[index_x,]

n1 <- nrow(x1)
n2 <- nrow(x2)
# residuals and projections based on x1 y1 and dimension reduction
# beta_hat projection and residuals based on lasso and data x1 y1
lasso_model1 <- cv.glmnet(as.matrix(x1), y1, family="gaussian", intercept = F)
lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
index_beta1_non0 <- seq(1:p)[as.numeric(lasso_beta1) != 0]
len_beta1_non0 <- length(index_beta1_non0)
if(len_beta1_non0 == 0){
  U1 <- y1 - mean(y1)                                                        # residual based on x1 y1
}else{
  x1_sec <- as.matrix(x1[,index_beta1_non0])                                            # second estimation
  # cor_matrix1 <- cor(x1_sec, use = "pairwise.complete.obs")
  # high_cor_vars1 <- findCorrelation(cor_matrix1, cutoff = 0.9, verbose = TRUE)
  # vars_to_remove1 <- colnames(cor_matrix1)[high_cor_vars1]
  # x1_sec <- x1_sec[,-which(colnames(x1_sec)%in%vars_to_remove1)]
  sec_model1 <- glm(y1~x1_sec-1, family = gaussian)
  pred1 = predict(sec_model1, newx = x1_sec, type="response")
  U1 <- y1-pred1
}

# projections based on screening and lassosir or seas for x1 U1
screen_num1 <- floor(n1/log(n1))
if(p <= screen_num1){
  # projection based on lassosir without screening
  sir_U1 <- LassoSIR(x1, U1, H=10, choosing.d="automatic", solution.path=FALSE, 
                     categorical=FALSE, nfolds=5, screening=FALSE)
  sir_Ubeta1<- sir_U1$beta
  sdr_Upro1 <- sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
}else{
  #screening
  rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                     # screening method based on distance correlation
  index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                       # chosen index of X
  x1_Uscreen <- x1[,index_U1]
  
  # projection based on lassosir after screening
  sir_U1 <- LassoSIR(x1_Uscreen, U1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_Ubeta1<- sir_U1$beta
  sdr_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
  sdr_Upro1[index_U1,] <- sir_Ubeta1
  sdr_Upro1 <- sdr_Upro1/sqrt(colSums(sdr_Upro1^2))
}

# projections based on screening and lassosir or seas for x1 y1
if(p <= screen_num1){
  # laosssir projection without screening
  sir_y1 <- LassoSIR(x1, y1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_ybeta1<- sir_y1$beta
  sdr_ypro1 <- sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
}else{
  #screening
  rank_y1 <- screenIID(X=x1, Y=y1, method = "DC-SIS")                     # screening method 
  index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                       # chosen index of X
  x1_yscreen <- x1[,index_y1]
  
  # projection based on lassosir after screening
  sir_y1 <- LassoSIR(x1_yscreen, y1, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_ybeta1<- sir_y1$beta
  sdr_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
  sdr_ypro1[index_y1,] <- sir_ybeta1
  sdr_ypro1 <- sdr_ypro1/sqrt(colSums(sdr_ypro1^2))
}

# residuals and projections based on x2 y2 and dimension reduction
# residuals based on lasso and data x2 y2
lasso_model2 <- cv.glmnet(as.matrix(x2), y2, family="gaussian", intercept = F)
lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]                   # estimated beta0 without intercept 
index_beta2_non0 <- seq(1:p)[as.numeric(lasso_beta2)!=0]                  # index of beta_n with nonzero components
len_beta2_non0 <- length(index_beta2_non0)
if(len_beta2_non0 == 0){
  U2 <- y2-mean(y2)                                                       # residual based on x2 y2
}else{
  x2_sec <- as.matrix(x2[,index_beta2_non0])
  # cor_matrix2 <- cor(x2_sec, use = "pairwise.complete.obs")
  # high_cor_vars2 <- findCorrelation(cor_matrix2, cutoff = 0.9, verbose = TRUE)
  # vars_to_remove2 <- colnames(cor_matrix2)[high_cor_vars2]
  # x2_sec <- x2_sec[,-which(colnames(x2_sec)%in%vars_to_remove2)]
  sec_model2 = glm(y2~x2_sec-1,family = gaussian)                           # second estimation
  # sec_beta2 = unname(sec_model2$coefficients)[-1]
  # U2 <- y2-x2_sec%*%sec_beta2-sec_model2$coefficients[1]
  pred2 = predict(sec_model2, newx = x2_sec, type="response")
  U2 <- y2-pred2
}

#sir projections based on screening and x2 U2
screen_num2 <- floor(n2/log(n2))
if(p <= screen_num2){
  # lassosir projection without screen
  sir_U2 <- LassoSIR(x2, U2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_Ubeta2<- sir_U2$beta
  sdr_Upro2 <- sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))                     # sdr projection based on lassosir and x2 U2
}else{
  #screening
  rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                     # screening method based on distance correlation
  index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                       # chosen index of X
  x2_Uscreen <- x2[,index_U2]
  
  # projection based on lassosir after screening
  sir_U2 <- LassoSIR(x2_Uscreen, U2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_Ubeta2 <- sir_U2$beta
  sdr_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
  sdr_Upro2[index_U2,] <- sir_Ubeta2
  sdr_Upro2 <- sdr_Upro2/sqrt(colSums(sdr_Upro2^2))                       # sdr projection based on lassosir and x2 U2
}

#sir projections based on screening and x2 y2
if(p <= screen_num2){
  # Lassosir projection without screen
  sir_y2 <- LassoSIR(x2, y2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                     nfolds=5, screening=FALSE)
  sir_ybeta2<- sir_y2$beta
  sdr_ypro2 <- sir_ybeta2/sqrt(colSums(sir_ybeta2^2))                     # sdr projection based on lassosir and x2 y2
}else{ 
  #screening
  rank_y2 <- screenIID(X=x2, Y=y2, method = "DC-SIS")                     # screening method based on distance correlation
  index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                       # chosen index of X
  x2_yscreen <- x2[,index_y2]
  
  # projection based on lassosir after screening
  sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                          nfolds=5, screening=FALSE)
  sir_ybeta2 <- sir_ylasso2$beta
  sdr_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
  sdr_ypro2[index_y2,] <- sir_ybeta2
  sdr_ypro2 <- sdr_ypro2/sqrt(colSums(sdr_ypro2^2))                       # sdr projection based on lassosir and x2 y2
}

#construct test statistic based on x1 y1 and projections based on x2 y2
pro2 <- cbind(sdr_Upro2, sdr_ypro2)              # projections based on x1 y1
pro2_num <- ncol(pro2)                                                    # the number of projections

h1 <- (n1^(-2/9))  #bandwidth 

pval_matrix1 <- matrix(nrow=1, ncol=pro2_num) # p-value matrix
errormat1 <- U1%*%t(U1)  #residual matrix based on x1 y1
for(q in 1:pro2_num){
  x_pro1 <- as.matrix(x1)%*%pro2[,q]
  x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1  # kernel function matrix
  #kermat1 <-(1/sqrt(2*pi))*exp(-(x_pro1_mat^2)/2)                           # Gaussian kernel
  indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0) 
  kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                # Epanechnikov kernel 
  #test statistics
  Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
  pval1 <- 1-pnorm(Tn1)
  pval_matrix1[,q] <- pval1
}

#construct test statistics based on x2 y2 and projections based on x1 y1
pro1 <- cbind(sdr_Upro1, sdr_ypro1)                           # projections based on x1 y1
pro1_num <- ncol(pro1)                                                    # the number of projections

h2 <- (n2^(-2/9))                                                      # bandwidth 

pval_matrix2 <- matrix(nrow=1, ncol=pro1_num) # p-value matrix
errormat2 <- U2%*%t(U2)  #residual matrix based on x2 y2
for(l in 1:pro1_num){
  x_pro2 <- as.matrix(x2)%*%pro1[,l]
  x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2  #kernel function matrix
  #kermat2 <-(1/sqrt(2*pi))*exp(-(x_pro2_mat^2)/2)                            # Gaussian kernel
  indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0) 
  kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                # Epanechnikov kernel 
  #test statistics
  Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
  pval2 <- 1-pnorm(Tn2)
  pval_matrix2[,l] <- pval2
}

pval_matrix <- cbind(pval_matrix1, pval_matrix2)
pval_num <- pro1_num + pro2_num

pval_fisher1 <- 1- pchisq(-2*sum(log(pval_matrix1[1:pro2_num])),df=2*pro2_num)
pval_fisher2 <- 1- pchisq(-2*sum(log(pval_matrix2[1:pro1_num])),df=2*pro1_num)
pval_fisher_mat <- cbind(pval_fisher1,pval_fisher2)
pval_cauchy_fisher <- 1- pcauchy(mean(tan((0.5-pval_fisher_mat)*pi))) 


pval_min1 <- 1 - (1 - min(pval_matrix1))^pro2_num
pval_min2 <- 1 - (1 - min(pval_matrix2))^pro1_num
pval_min_mat <- cbind(1 - (1 - min(pval_matrix1))^pro2_num, 1 - (1 - min(pval_matrix2))^pro1_num) 
pval_cauchy_min <- 1- pcauchy(mean(tan((0.5-pval_min_mat)*pi))) 
pval_matrix;pval_fisher1;pval_fisher2;pval_cauchy_fisher;pval_min1;pval_min2;pval_cauchy_min
# pval_cauchy1 <- 1- pcauchy(mean(tan((0.5-pval_matrix1[1:pro2_num])*pi)))  # cauchy combination based on x1 y1
# pval_cauchy2 <- 1- pcauchy(mean(tan((0.5-pval_matrix2[1:pro1_num])*pi)))  # cauchy combination based on x1 y1
# pval_cauchy  <- 1- pcauchy(mean(tan((0.5-pval_matrix[1:pval_num])*pi)))   # cauchy combination based on x y
# pval_cauchy1;pval_cauchy2;pval_cauchy



######scatterplot_2
p <- ncol(X_quad_inter)
lasso_model <- cv.glmnet(as.matrix(X_quad_inter), Y, family="gaussian", intercept = F)
beta_lasso <- coef(lasso_model, s = "lambda.min")[-1]
index_beta_non0 <- seq(1:p)[as.numeric(beta_lasso) != 0]
x_sec <- as.matrix(X_quad_inter[,index_beta_non0] )                                         
sec_model = glm(Y~x_sec-1,family = gaussian)                           # second estimation
pred = predict(sec_model, newx = x_sec, type="response")
U <- Y-pred
xbeta <- x_sec%*%sec_model$coefficients   
###plot
plot_range <- range(c(U, xbeta))
plot(xbeta, U, 
     xlab = "",  
     ylab = "Residuals",
     pch = 16,
     ylim = plot_range,
     cex.lab = 2.5,
     cex.axis = 2,
     mgp = c(2, 0.8, 0)  
)
mtext(expression(hat(Y)), side = 1, 
      line = 3,  
      cex = 2.5)   




#### prediction：MSE MAE
# mtkl_pre <- matrix(nrow=50,ncol=4)
# for(i in 1:50){
#   n <- nrow(X)
#   p <- ncol(X)
#   new_names <- paste("col", 1:p, sep = "_")
#   colnames(X) <- new_names
#   index_x<-sample(1:n, floor(0.75*n),replace = FALSE)  
#   index_x<-sort(index_x)
#   x_train <- X[index_x,] 
#   y_train <- as.matrix(Y[index_x]) 
#   
#   x_test <- X[-index_x,]
#   y_test <- as.matrix(Y[-index_x])
#   
#   
#   x_train_quad_inter <- X_quad_inter[index_x,]
#   x_test_quad_inter <- X_quad_inter[-index_x,]
#   
#   ####linear
#   p1 <- ncol(x_train)
#   lasso_model_linear <- cv.glmnet(as.matrix(x_train), y_train, family="gaussian", intercept = F)
#   beta_lasso_linear <- coef(lasso_model_linear, s = "lambda.min")[-1]
#   index_beta_non0_linear <- seq(1:p1)[as.numeric(beta_lasso_linear) != 0]
#   x_sec_train <- x_train[,index_beta_non0_linear]                                          
#   sec_model_linear = glm(y_train ~ . - 1, data = as.data.frame(x_sec_train), family = gaussian)                           
#   
#   pred_linear = predict(sec_model_linear, newdata = as.data.frame(x_test[, index_beta_non0_linear]), type="response")
#   U_linear <- y_test-pred_linear
#   MSE_linear <- mean(U_linear^2)
#   MAE_linear <- mean(abs(U_linear))
#   
#   
#   ###nonlinear
#   p2 <- ncol(x_train_quad_inter)
#   lasso_model_nonlinear <- cv.glmnet(as.matrix(x_train_quad_inter), y_train, family="gaussian", intercept = F)
#   beta_lasso_nonlinear <- coef(lasso_model_nonlinear, s = "lambda.min")[-1]
#   index_beta_non0_nonlinear <- seq(1:p2)[as.numeric(beta_lasso_nonlinear) != 0]
#   x_sec_train_nonlinear <- x_train_quad_inter[,index_beta_non0_nonlinear]                                          
#   sec_model_nonlinear = glm(y_train ~ . - 1, data = as.data.frame(x_sec_train_nonlinear), family = gaussian)                          
#   
#   pred_nonlinear = predict(sec_model_nonlinear, newdata = as.data.frame(x_test_quad_inter[, index_beta_non0_nonlinear]), type="response")
#   U_nonlinear <- y_test-pred_nonlinear
#   MSE_nonlinear <- mean(U_nonlinear^2)
#   MAE_nonlinear <- mean(abs(U_nonlinear))
#   
#   
#   #result
#   mtkl_pre[i,1] <- MSE_linear
#   mtkl_pre[i,2] <- MAE_linear
#   mtkl_pre[i,3] <- MSE_nonlinear
#   mtkl_pre[i,4] <- MAE_nonlinear
#   
# }
# apply(mtkl_pre,2,mean)





