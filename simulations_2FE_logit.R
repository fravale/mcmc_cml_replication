rm(list=ls())


library("alpaca")
source("PCML.R")
source("MCMC_CML.R")

set.seed(157956)
## Simulate data from a 2-FE Logit model
iters = 1000

BC = SEBC = CC = CCP = SE = SEP = CMC = SEMC  = CMC2 = SEMC2 = matrix(0, iters,2)

TIME = matrix(NA,iters,5)

RES = array(0,c(8,5,2))

b1 = 1
be2 = 2.5

NOBS =  c(25,50)


######## Generate DATA ################


for(Nobs in NOBS){

    nu = T = Nobs


    filename = sprintf("sim_2FE_charb3_Nobs%g",nu)

    ## ID
    idi = seq(1:nu)%x%rep(1,T) ## INDEX BASED ON INDIVIDUALS

    idt = rep(1,nu)%x%seq(1:T) ## INDEX BASED ON TIME OCCASIONS

    
    for(nit in 1:iters){

        
        ## INDIVIDUAL EFFECTS
        alpha = rnorm(nu)
        
        Alpha = alpha%x%rep(1,T)
        
        gamma = rnorm(nu)
        
        GGamma = rep(1,T)%x%gamma


        ## COVARIATES
        
        X = Alpha + GGamma + rnorm(nu*T)
        Z = (runif(nu*T)<=0.5)

        ## DEP. VAR.
        
        pr = exp(Alpha + GGamma + b1*X + be2*Z)/(1+exp(Alpha + GGamma + b1*X+ be2*Z))
        
        y = as.numeric(pr > runif(nu*T))
        
        y[idi==idt] = NA


        
################################################################
        TRUEVAL = c(b1,be2)
        covariates = cbind(X,Z)
        
        
        
        
        ## ML
        DAT = as.data.frame(cbind(idi,idt,y,X,Z))
        t0 = proc.time()
        out_ml <- try(feglm(y ~ X + Z  | idi + idt, DAT, binomial("logit")))
        if(inherits(out_ml,"try-error")){
            CCP[nit,] = NA
            SEP[nit,] = NA
            BC[nit,] = NA
            SEBC[nit,] = NA
            CMC[nit,] = NA
            SEMC[nit,] = NA
        }else{
            

            t1 = proc.time() - t0
            TIME[nit,1] = t1[3]

            
            ## BC
            t0 = proc.time()
            out_bc <- biasCorr(out_ml)
            t2 = proc.time() - t0
            TIME[nit,2] = t2[3]

            ## Store ML
            CCP[nit,] = out_ml$coefficients
            SEP[nit,] = sqrt(diag(solve(out_ml$Hessian)))


            ## Store BC

            BC[nit,] = out_bc$coefficients
            SEBC[nit,] = sqrt(diag(solve(out_bc$Hessian)))
            
            
            
            ## MCMC - 100k replications
            cat("MCMC1 started\n")

            
            initb = out_ml$coefficients[1:2]
            t0 = proc.time()
            outmc = try(MCMC_MLE(idi,y,covariates,initb,100000,subele=5,bin_par=0.2,thin=100))
            if(inherits(outmc,"try-error")){
                CMC[nit,]  = NA
                SEMC[nit,] = NA
            }else{
                t1 = proc.time() - t0
                TIME[nit,4] = t1[3]
                CMC[nit,] = outmc$theta
                SEMC[nit,] = outmc$se_theta
            }

            ## MCMC - 200k replications
            cat("MCMC2 started\n")

            
            initb = out_ml$coefficients[1:2]
            t0 = proc.time()
            outmc2 = try(MCMC_MLE(idi,y,covariates,initb,200000,subele=5,bin_par=0.2,thin=100))
            if(inherits(outmc2,"try-error")){
                CMC2[nit,]  = NA
                SEMC2[nit,] = NA
            }else{
                t1 = proc.time() - t0
                TIME[nit,5] = t1[3]
                CMC2[nit,] = outmc2$theta
                SEMC2[nit,] = outmc2$se_theta
            }


        }
        
        ## PCML
        t0 = proc.time()
        out = try(jochmans_charb(idi,y,covariates,verbose=FALSE))
        if(inherits(out,"try-error")){
            CC[nit,] = NA
            SE[nit,] = NA
        }else{
            t1 = proc.time() - t0
            TIME[nit,3] = t1[3]
            
            CC[nit,] = out$beta
            SE[nit,] = out$ser
        }
        
######## REPORT RESULTS #############
        cat("Iteration:",nit,"\n")    
        
        statnames = c("mean bias","median bias ","sd","iqr","rmse","mae","se/std","size")
        estinames = c("MLE","BC","PCML","MCMC-100k","MCMC-200k")
        
        if(nit>1){
            
            for(cof in 1:2){
                
                ## ML
                RES[1,1,cof] = mean(CCP[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[2,1,cof] = median(CCP[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[3,1,cof] = sd(CCP[1:nit,cof],na.rm=TRUE)
                RES[4,1,cof] = IQR(CCP[1:nit,cof],na.rm=TRUE)
                RES[5,1,cof] = sqrt(mean((CCP[1:nit,cof] - TRUEVAL[cof])^2,na.rm=TRUE))
                RES[6,1,cof] = median(abs(CCP[1:nit,cof] - TRUEVAL[cof]),na.rm=TRUE)
                RES[7,1,cof] = mean(SEP[1:nit,cof],na.rm=TRUE)/sd(CCP[1:nit,cof],na.rm=TRUE)
                RES[8,1,cof] = mean(as.numeric(abs((CCP[1:nit,cof]-TRUEVAL[cof])/SEP[1:nit,cof])>1.96),na.rm=TRUE)
                
                
                ## BC
                RES[1,2,cof] = mean(BC[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[2,2,cof] = median(BC[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[3,2,cof] = sd(BC[1:nit,cof],na.rm=TRUE)
                RES[4,2,cof] = IQR(BC[1:nit,cof],na.rm=TRUE)
                RES[5,2,cof] = sqrt(mean((BC[1:nit,cof] - TRUEVAL[cof])^2,na.rm=TRUE))
                RES[6,2,cof] = median(abs(BC[1:nit,cof] - TRUEVAL[cof]),na.rm=TRUE)
                RES[7,2,cof] = mean(SEBC[1:nit,cof],na.rm=TRUE)/sd(BC[1:nit,cof],na.rm=TRUE)
                RES[8,2,cof] = mean(as.numeric(abs((BC[1:nit,cof]-TRUEVAL[cof])/SEBC[1:nit,cof])>1.96),na.rm=TRUE)
                
                
                ## P-CML
                RES[1,3,cof] = mean(CC[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[2,3,cof] = median(CC[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[3,3,cof] = sd(CC[1:nit,cof],na.rm=TRUE)
                RES[4,3,cof] = IQR(CC[1:nit,cof],na.rm=TRUE)
                RES[5,3,cof] = sqrt(mean((CC[1:nit,cof] - TRUEVAL[cof])^2,na.rm=TRUE))
                RES[6,3,cof] = median(abs(CC[1:nit,cof] - TRUEVAL[cof]),na.rm=TRUE)
                RES[7,3,cof] = mean(SE[1:nit,cof],na.rm=TRUE)/sd(CC[1:nit,cof],na.rm=TRUE)
                RES[8,3,cof] = mean(as.numeric(abs((CC[1:nit,cof]-TRUEVAL[cof])/SE[1:nit,cof])>1.96),na.rm=TRUE)
                
                ## MCMC - 100k
                
                RES[1,4,cof] = mean(CMC[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[2,4,cof] = median(CMC[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[3,4,cof] = sd(CMC[1:nit,cof],na.rm=TRUE)
                RES[4,4,cof] = IQR(CMC[1:nit,cof],na.rm=TRUE)
                RES[5,4,cof] = sqrt(mean((CMC[1:nit,cof] - TRUEVAL[cof])^2,na.rm=TRUE))
                RES[6,4,cof] = median(abs(CMC[1:nit,cof] - TRUEVAL[cof]),na.rm=TRUE)
                RES[7,4,cof] = mean(SEMC[1:nit,cof],na.rm=TRUE)/sd(CMC[1:nit,cof],na.rm=TRUE)
                RES[8,4,cof] = mean(as.numeric(abs((CMC[1:nit,cof]-TRUEVAL[cof])/SEMC[1:nit,cof])>1.96),na.rm=TRUE)

                ## MCMC - 200k
                
                RES[1,5,cof] = mean(CMC2[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[2,5,cof] = median(CMC2[1:nit,cof] - TRUEVAL[cof],na.rm=TRUE)
                RES[3,5,cof] = sd(CMC2[1:nit,cof],na.rm=TRUE)
                RES[4,5,cof] = IQR(CMC2[1:nit,cof],na.rm=TRUE)
                RES[5,5,cof] = sqrt(mean((CMC2[1:nit,cof] - TRUEVAL[cof])^2,na.rm=TRUE))
                RES[6,5,cof] = median(abs(CMC2[1:nit,cof] - TRUEVAL[cof]),na.rm=TRUE)
                RES[7,5,cof] = mean(SEMC2[1:nit,cof],na.rm=TRUE)/sd(CMC2[1:nit,cof],na.rm=TRUE)
                RES[8,5,cof] = mean(as.numeric(abs((CMC2[1:nit,cof]-TRUEVAL[cof])/SEMC2[1:nit,cof])>1.96),na.rm=TRUE)                     
                
                
                resplot = RES[,,cof]
                
                if(cof==1) cat("b1 = 1\n\n") else cat("b2 = 2.5\n\n")
                colnames(resplot) = estinames
                rownames(resplot) = statnames
                print(resplot)
                
            }
            cat("\n")
            
            ## Covergence Failures
            cat("Convergence Failed\n")
            FAIL = cbind(sum(is.na(CCP[1:nit,1])),sum(is.na(BC[1:nit,1])),sum(is.na(CC[1:nit,1])),sum(is.na(CMC[1:nit,1])),sum(is.na(CMC2[1:nit,1])))

            rownames(FAIL) = c("Failures")
            colnames(FAIL) = estinames
            print(FAIL)
            cat("\n")
            
            ## Times
            cat("Computatinal Time (in secs)\n")
            CTIME = t(as.matrix(colMeans(TIME[1:nit,],na.rm=TRUE)))
            rownames(CTIME) = c("Time (secs)")
            colnames(CTIME) = estinames
            print(CTIME)
            cat("\n\n\n")
        }
    }
    
    save.image(filename)

}
