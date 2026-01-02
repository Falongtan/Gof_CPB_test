

Gof_CPB_test <- function(X,Y,fam = c("gaussian", "binomial", "poisson"),penalize){
  #library packages
  library('glmnet')
  library(LassoSIR)
  library(VariableScreening)
  library(psych)
  #check if X , Y and the parameter are set correctly
  if (!is.matrix(X)){
    stop("X should be a matrix with at least one column.")
  }
  np <- dim(X)
  if (is.null(np) | (np[2] < 1L)){
    stop("X should be a matrix with at least one column.")
  }
  n <- as.integer(np[1])
  p <- as.integer(np[2])
  
  Y <- as.numeric(Y)
  if(fam == "binomial"){
    categorical_Y <- TRUE
  }else{
    categorical_Y <- FALSE
  }
  if (length(Y) != n){
    stop("Y must have nrow(X) components.")
  }
  
  if ((p >= n-1) && (penalize == FALSE)){
    stop("When penalize=FALSE we must have ncol(X) < nrow(X)-1; try setting penalize=TRUE.")
  }
  
  if(fam != "gaussian" && fam != "binomial" && fam != "poisson"){
    stop("family must be gaussian, binomial or poisson; try changing parameter fam.")
  }
  #data splitting
  index_x<-sample(1:n, floor(n/2),replace = FALSE)  
  index_x<-sort(index_x)
  #x1 y1 n1
  x1 <- X[-index_x,] 
  y1 <- as.matrix(Y[-index_x]) 
  n1 <- nrow(x1)
  #x2 y2 n2
  x2 <- X[index_x,]
  y2 <- as.matrix(Y[index_x])
  n2 <- nrow(x2)
  if(penalize == TRUE){
    #penalize
    # residuals  based on x1 y1 
    lasso_model1 <- cv.glmnet(as.matrix(x1), y1, family=fam)
    lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
    index_beta1_non0 <- seq(1:p)[as.numeric(lasso_beta1) != 0]
    len_beta1_non0 <- length(index_beta1_non0)
    if(len_beta1_non0 == 0){
      U1 <- as.numeric(y1) - mean(as.numeric(y1) )                                                        
    }else{
      x1_sec <- x1[,index_beta1_non0]                                            
      sec_model1 <- glmnet(as.matrix(x1_sec), y1, family = fam, lambda = 0)     # second estimation
      pred1 = predict(sec_model1, newx = x1_sec, type="response")
      U1 <- as.numeric(y1)-pred1   
    }
    # residuals  based on x2 y2 
    lasso_model2 <- cv.glmnet(as.matrix(x2), y2, family=fam)
    lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]                   
    index_beta2_non0 <- seq(1:p)[as.numeric(lasso_beta2)!=0]                  
    len_beta2_non0 <- length(index_beta2_non0)
    if(len_beta2_non0 == 0){
      U2 <- as.numeric(y2) - mean(as.numeric(y2))                                                       
    }else{
      x2_sec <- x2[,index_beta2_non0]
      sec_model2 = glmnet(as.matrix(x2_sec), y2, family = fam, lambda = 0)      # second estimation
      pred2 = predict(sec_model2, newx = x2_sec, type="response")
      U2 <- as.numeric(y2)-pred2                                                
    }
  }else{
    #non-penalize
    # residuals  based on x1 y1 
    fit1 <- glmnet(as.matrix(x1), y1, family = fam, lambda = 0)
    pred1 <- predict(fit1,newx=x1,type="response")
    U1 <- as.numeric(y1)-pred1
    # residuals  based on x2 y2
    fit2 <- glmnet(as.matrix(x2), y2, family = fam, lambda = 0)
    pred2 <- predict(fit2,newx=x2,type="response")
    U2 <- as.numeric(y2)-pred2
  }
  
  ####projection based x1 y1 U1
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
    sir_y1 <- LassoSIR(x1, as.numeric(y1), H=10, choosing.d="automatic", solution.path=FALSE, categorical=categorical_Y, 
                       nfolds=5, screening=FALSE)
    sir_ybeta1<- sir_y1$beta
    sdr_ypro1 <- sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
  }else{
    #screening
    rank_y1 <- screenIID(X=x1, Y=as.numeric(y1), method = "DC-SIS")                     # screening method 
    index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                       # chosen index of X
    x1_yscreen <- x1[,index_y1]
    
    # projection based on lassosir after screening
    sir_y1 <- LassoSIR(x1_yscreen, as.numeric(y1), H=10, choosing.d="automatic", solution.path=FALSE, categorical=categorical_Y, 
                       nfolds=5, screening=FALSE)
    sir_ybeta1<- sir_y1$beta
    sdr_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
    sdr_ypro1[index_y1,] <- sir_ybeta1
    sdr_ypro1 <- sdr_ypro1/sqrt(colSums(sdr_ypro1^2))
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
    sir_y2 <- LassoSIR(x2, as.numeric(y2), H=10, choosing.d="automatic", solution.path=FALSE, categorical=categorical_Y, 
                       nfolds=5, screening=FALSE)
    sir_ybeta2<- sir_y2$beta
    sdr_ypro2 <- sir_ybeta2/sqrt(colSums(sir_ybeta2^2))                     # sdr projection based on lassosir and x2 y2
  }else{ 
    #screening
    rank_y2 <- screenIID(X=x2, Y=as.numeric(y2), method = "DC-SIS")                     # screening method based on distance correlation
    index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                       # chosen index of X
    x2_yscreen <- x2[,index_y2]
    
    # projection based on lassosir after screening
    sir_ylasso2 <- LassoSIR(x2_yscreen, as.numeric(y2), H=10, choosing.d="automatic", solution.path=FALSE, categorical=categorical_Y, 
                            nfolds=5, screening=FALSE)
    sir_ybeta2 <- sir_ylasso2$beta
    sdr_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
    sdr_ypro2[index_y2,] <- sir_ybeta2
    sdr_ypro2 <- sdr_ypro2/sqrt(colSums(sdr_ypro2^2))                       # sdr projection based on lassosir and x2 y2
  }
  #construct test statistic based on x1 y1 and projections based on x2 y2
  pro2 <- cbind(sdr_Upro2, sdr_ypro2)              # projections based on x1 y1
  pro2_num <- ncol(pro2)                                                    # the number of projections
  
  h1 <- ch*(n1^(-2/9))  #bandwidth 
  
  pval_matrix1 <- matrix(nrow=1, ncol=pro2_num) # p-value matrix
  errormat1 <- U1%*%t(U1)  #residual matrix based on x1 y1
  for(q in 1:pro2_num){
    x_pro1 <- x1%*%pro2[,q]
    x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1  # kernel function matrix
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
  
  h2 <- ch*(n2^(-2/9))                                                      # bandwidth 
  
  pval_matrix2 <- matrix(nrow=1, ncol=pro1_num) # p-value matrix
  errormat2 <- U2%*%t(U2)  #residual matrix based on x2 y2
  for(l in 1:pro1_num){
    x_pro2 <- x2%*%pro1[,l]
    x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2  #kernel function matrix
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
  return(list(pval_cauchy_fisher=pval_cauchy_fisher,
              pval_cauchy_min=pval_cauchy_min))
  
}