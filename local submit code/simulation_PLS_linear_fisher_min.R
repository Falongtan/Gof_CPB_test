####PLS_linear
###test whether the linear model is properly formulated
###use min and fisher combination method to combine the p-values on half of the data, respectively
###use the Cauchy combination method to combine the p-values from two data 
###compare with rp test and grp test
PLS_lm_parallel <- function(p,n,a,example,c_h,pho,rpgrp_switch){
  #library packages
  library(MASS)
  library(glmnet)
  library(LassoSIR) 
  library(VariableScreening)
  library(psych)
  library(foreach)
  library(parallel)
  library(iterators)
  library(doParallel)
  library(RPtests)
  library(GRPtests)
  #the number of Monte Carlo
  s <- 1000
  mu <- rep(0, p)
  if(pho == 0){
    sigma <- diag(rep(1,p))
  }else{
    v <- pho^(0:(p-1))
    sigma <- toeplitz(v)
  }   
  beta0 <- c(rep(1,5),rep(0,p-5))
  beta1 <- c(rep(1,10),rep(0,p-10))
  #parallel
  # tic()
  if (detectCores(logical=F) >= 16){
    cores <- 14
  }else if(detectCores(logical=F)<=6){
    cores <- detectCores(logical=F)
  }else{
    cores <- detectCores(logical=F)-4
  }
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores=cores)
  PLS_power <- foreach(k = 1:s, .combine='rbind', 
                       .packages = c('MASS','glmnet','LassoSIR','VariableScreening','psych',
                                     'foreach','parallel','iterators','doParallel','RPtests','GRPtests')) %dopar%  
    {
      #generate X and Y
      x <- mvrnorm(n, mu, sigma)
      if(example==1){
        y <-  x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)  # H11
      }else if(example==2){
        y <- x %*% beta0 + a * cos(0.6 * pi * x %*% beta0) + rnorm(n)            # H12
      }else{
        y <- x %*% beta0 + a *exp(0.5*x%*%beta1) + rnorm(n)                      # H13
      }
      
      ##rp and grp: not related to c_h,  run only once and run rp and grp if needed
      if(rpgrp_switch && c_h==1){
        pval_rp <- RPtest(x, y, test="nonlin", B=49L, nperms=2, resid_type = "Lasso")
        pval_grp <- GRPtest(x, y, fam = "gaussian", nsplits = 1,RP_function = NULL,output_all=FALSE,penalize=TRUE)
      }else{
        pval_rp <- NA
        pval_grp <- NA
      }
      
                        
      # data splitting  
      index_x<-sample(1:n, floor(n/2),replace = FALSE)  
      index_x<-sort(index_x)
      
      x1 <- x[-index_x,] 
      y1 <- as.matrix(y[-index_x,]) 
      
      x2 <- x[index_x,]
      y2 <- as.matrix(y[index_x,])
      
      n1 <- nrow(x1)
      n2 <- nrow(x2)
      
      # residuals based on lasso and data x1 y1
      lasso_model1 <- cv.glmnet(x1, y1, family="gaussian", intercept = T)
      lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
      index_beta1_non0 <- seq(1:p)[as.numeric(lasso_beta1) != 0]
      len_beta1_non0 <- length(index_beta1_non0)
      if(len_beta1_non0 == 0){
        U1 <- y1 - mean(y1)                                                        # residual based on x1 y1
      }else{
        x1_sec <- x1[,index_beta1_non0]                                            # second estimation
        sec_model1 <- glm(y1~x1_sec, family = gaussian)
        sec_beta1 = unname(sec_model1$coefficients)[-1]
        lasso_beta1[index_beta1_non0] <- sec_beta1
        beta1_hat <- lasso_beta1
        pred1 = predict(sec_model1, newx = x1_sec, type="response")
        pred1 = matrix(unname(pred1), ncol = 1)
        U1 <- y1-pred1                                                             # residual based on x1 y1
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
      
      # residuals based on lasso and data x2 y2
      lasso_model2 <- cv.glmnet(x2, y2, family="gaussian", intercept = T)
      lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]                    # estimated beta0 without intercept 
      index_beta2_non0 <- seq(1:p)[as.numeric(lasso_beta2)!=0]                  # index of beta_n with nonzero components
      len_beta2_non0 <- length(index_beta2_non0)
      if(len_beta2_non0 == 0){
        U2 <- y2-mean(y2)                                                       # residual based on x2 y2
      }else{
        x2_sec <- x2[,index_beta2_non0]
        sec_model2 = glm(y2~x2_sec,family = gaussian)                           # second estimation
        sec_beta2 <- unname(sec_model2$coefficients)[-1]
        lasso_beta2[index_beta2_non0] <- sec_beta2
        beta2_hat <- lasso_beta2
        pred2 = predict(sec_model2, newx = x2_sec, type="response")
        pred2 = matrix(unname(pred2), ncol = 1)
        U2 <- y2-pred2                                                          # residual based on x2 y2
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
      pro2 <- cbind(sdr_Upro2, sdr_ypro2)                                       # projections based on x1 y1
      pro2_num <- ncol(pro2)                                                    # the number of projections
      h1 <- c_h*(n1^(-2/9))                                                     # bandwidth 
      
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
      pro1 <- cbind(sdr_Upro1, sdr_ypro2)                                          # projections based on x1 y1
      pro1_num <- ncol(pro1)                                                       # the number of projections
      h2 <- c_h*(n2^(-2/9))                                                        # bandwidth 
      
      pval_matrix2 <- matrix(nrow=1, ncol=pro1_num) # p-value matrix
      errormat2 <- U2%*%t(U2)  #residual matrix based on x2 y2
      for(l in 1:pro1_num){
        x_pro2 <- x2%*%pro1[,l]
        x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2  # kernel function matrix
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
      pval_fisher_m <- cbind(pval_fisher1, pval_fisher2)
      pval_fisher <- 1- pcauchy(mean(tan((0.5-pval_fisher_m)*pi)))  
      
      pval_min1 <-   min(pval_matrix1)<=1-(0.95)^(1/pro2_num) 
      pval_min2 <-  min(pval_matrix2)<=1-(0.95)^(1/pro1_num) 
      pval_min_m <- cbind(1 - (1 - min(pval_matrix1))^pro2_num, 1 - (1 - min(pval_matrix2))^pro1_num) 
      pval_min <- 1- pcauchy(mean(tan((0.5-pval_min_m)*pi))) 
      
      result <- c(pval_fisher1,pval_fisher2,pval_fisher,
                  pval_min1,pval_min2,pval_min, pval_rp, pval_grp)  
    }
  #end parallel
  stopImplicitCluster()
  stopCluster(cl)
  
  # toc()
  #power
  PLS_pval_fisher1 <- PLS_power[,1]
  PLS_pval_fisher2 <- PLS_power[,2]
  PLS_pval_fisher  <- PLS_power[,3]
  PLS_pval_min1  <- PLS_power[,4]
  PLS_pval_min2  <- PLS_power[,5]
  PLS_pval_min  <- PLS_power[,6]
  ##rp and grp
  rp_pval <- PLS_power[,7]
  grp_pval <- PLS_power[,8]
  
  
  PLS_fisher1_power <- mean(PLS_pval_fisher1<=0.05,na.rm=TRUE)
  PLS_fisher2_power <- mean(PLS_pval_fisher2<=0.05,na.rm=TRUE)
  PLS_fisher_power <- mean(PLS_pval_fisher<=0.05,na.rm=TRUE)
  PLS_min1_power <- mean(PLS_pval_min1,na.rm=TRUE)
  PLS_min2_power <- mean(PLS_pval_min2,na.rm=TRUE)
  PLS_min_power <- mean(PLS_pval_min <=0.05,na.rm=TRUE)
  rp_power <- mean( rp_pval<=0.05,na.rm=TRUE)
  grp_power <- mean(grp_pval<=0.05,na.rm=TRUE)
  
  
  return(list(PLS_fisher1_power=PLS_fisher1_power,
              PLS_fisher2_power=PLS_fisher2_power,
              PLS_fisher_power=PLS_fisher_power,
              PLS_min1_power=PLS_min1_power,
              PLS_min2_power=PLS_min2_power,
              PLS_min_power=PLS_min_power,
              rp_power=rp_power,grp_power=grp_power))

}



#####analysis of Empirical sizes and powers of our tests with different bandwidths in linear models
example <- c(1,2,3) 
a <- c(0,1) #null and alternative
n <- 300                      #sample size
p<-  c(50,100,300,600, 900,1200) #c(1200,1200)  #dimension
c_h <-  c(0.75, 1, 1.25)     # 1                      #bandwidths
pho <- c(0,.4,.8)           #correlation
###don't run rp and grp
rpgrp_run <- FALSE
####if needed to compare with rp and grp :rpgrp_run <- TRUE


r <- 8
compar_lm_result <- matrix(0,nrow = r*length(example)*length(pho)*length(a)*length(c_h),ncol = length(p)+5)
for(ll in 1:length(example)){
  for (kk in 1:length(pho)){
    for (ii in 1:length(a)){
      for (jj in 1:length(p)){
        for (cc in 1:length(c_h)){ 
          t3 <- Sys.time()
          result_pls <- PLS_lm_parallel(p[jj], n, a[ii],example[ll],c_h[cc],pho[kk],rpgrp_switch = rpgrp_run) 
          
          ##index
          
          base_index <- r * ((ll-1)*length(pho)*length(a)*length(c_h) +  (kk-1)*length(a)*length(c_h) + (ii-1)*length(c_h) + (cc-1) )
          
          # ---- result ----
          compar_lm_result[base_index + 1, jj]  <- result_pls$PLS_fisher1_power
          compar_lm_result[base_index + 2, jj]  <- result_pls$PLS_fisher2_power
          compar_lm_result[base_index + 3, jj]  <- result_pls$PLS_fisher_power
          compar_lm_result[base_index + 4, jj]  <- result_pls$PLS_min1_power
          compar_lm_result[base_index + 5, jj]  <- result_pls$PLS_min2_power
          compar_lm_result[base_index + 6, jj]  <- result_pls$PLS_min_power
          compar_lm_result[base_index + 7, jj]  <- result_pls$rp_power
          compar_lm_result[base_index + 8, jj]  <- result_pls$grp_power
          
          t4 <- Sys.time()
          cat("example = ",example[ll],"p = ", p[jj],", n = ", n, ", a = ", a[ii] ,", c_h = ", c_h[cc], ", pho = ", pho[kk],
              ", PLS_fisher1_power = ", result_pls$PLS_fisher1_power,
              ", PLS_fisher2_power= ", result_pls$PLS_fisher2_power,
              ", PLS_fisher_power = ", result_pls$PLS_fisher_power,
              ", PLS_min1_power = ", result_pls$PLS_min1_power, 
              ", PLS_min2_power = ", result_pls$PLS_min2_power,
              ", PLS_min_power= ", result_pls$PLS_min_power, 
              ", rp_power = ", result_pls$rp_power,
              ", grp_power = ",  result_pls$grp_power,  
              ", time:",t4-t3, "\n")
          
          
          compar_lm_result[base_index + 1:r, (length(p)+1):(length(p)+5)] <- 
            matrix(c(example[ll],n, a[ii], pho[kk], c_h[cc]), nrow = r, ncol = 5, byrow = TRUE)
          
        }
      }
    }
  }
}
outcome.name <- paste("PLS_linear_fisher_min_H1", example, ".txt", sep="")

write.table(compar_lm_result,file = outcome.name ,sep = "  ", row.names = FALSE, col.names = FALSE, eol = "\r\n")
