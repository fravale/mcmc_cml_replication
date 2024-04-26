library(cquad)
library(purrr)

prop_rec_and_hex <- function(oss,S,A,pow,ind){

    oss1 = oss
    strct = which(is.na(oss1))
    oss[strct]=0
    
    ind1 = c(t(A%*%oss)%*%pow+1)

    strct = which(is.na(oss1))
    sav = rowSums(as.matrix(S[,strct]==0))==length(strct)
    S = S[sav,]
    ind = ind[sav]
    sele = matrix(S[ind==ind1,],,9)
    
    
    if(length(sele)>0){
        exl = rep(0,nrow(sele))
        for(j in 1:nrow(sele)){
            exl[j] = all(sele[j,]==oss)
        }

        rownew = nrow(sele)-1
        sele = matrix(sele[!exl,],rownew)
           
        if(length(sele)>0){
            newmat = sele[rdunif(1,as.numeric(nrow(sele)),1),]
            newmat[strct] = NA
            error = FALSE
        }else{
            newmat = matrix(0,0,0)
            error = TRUE
        }
    }
    out = list(newmat = newmat, error=error)
    return(out)
}

reduce_table2 <-  function(ycur,par,random = TRUE,warn=FALSE){
    ## par = minimum dimension of the reduced matrix

    x = ycur
    
    nc = ncol(x)
    nr = nrow(x)
    
    rl = seq(1,nr)
    cl = seq(1,nc)
   
    if(random){
        dimr = min(par,nr)
        dimc = min(par,nc)
        rl = sort(sample(rl,dimr))       
        cl = sort(sample(cl,dimc))        
    }
    
    xt = x[rl,cl]
    
    
    rs = rowSums(xt,na.rm=TRUE)
    cs = colSums(xt,na.rm=TRUE)
    
    maxc = colSums(!is.na(xt))
    maxr = rowSums(!is.na(xt))


    ## Drop uniformative columns and rows
    keepc = (cs>0 & cs<maxc)
    keepr = (rs>0 & rs<maxr)
    xt = xt[keepr,keepc]
    cl = cl[keepc]
    rl = rl[keepr]

    if(length(cl)<2 | length(rl)<2){
        err = TRUE
        out = list(error=err)
        return(out)
    }

    ## Scan xt
    n1 = length(rl)
    c1 = length(cl)

    M = (matrix(1:(n1*c1),n1))


    ntab = n1*(n1-1)*c1*(c1-1)/4
    ind1 = ind2 = ind3 = ind4 = rep(0,ntab)
    j = 0
    LABEL = matrix(0,ntab,4)
    for(i1 in 1:(n1-1)){
        for(i2 in (i1+1):n1){                
            for(t1 in 1:(c1-1)){
                if(t1!=i1){
                    for(t2 in (t1+1):c1){
                        if(t2!=i2){
                            
                            j = j+1
                            
                            LABEL[j,] = cbind(rl[i1],rl[i2],cl[t1],cl[t2])
                            ind1[j] = M[i1,t1]; ind2[j] = M[i1,t2]
                            ind3[j] = M[i2,t1]; ind4[j] = M[i2,t2]
                        }
                        
                    }
                }
            }
        }
    }

    IND = cbind(ind1,ind2,ind3,ind4)
    

    sw0 = (xt[ind1]==0) & (xt[ind2]==1) & (xt[ind3]==1) & (xt[ind4]==0)
    sw1 = (xt[ind1]==1) & (xt[ind2]==0) & (xt[ind3]==0) & (xt[ind4]==1)
    
    tab0 = c(1,0,0,1)                                    
    tab1 = c(0,1,1,0)

    if(length(c(which(sw0),which(sw1)))==0){
        err = TRUE
        out = list(error=err)
        return(out)
    }else if(length(c(which(sw0),which(sw1)))==1){
        j = c(which(sw0),which(sw1))
    }else{
        j = sample(c(which(sw0),which(sw1)),1)
    }

    rl = LABEL[j,1:2]
    cl = LABEL[j,3:4]


    err = FALSE
    out = list(rl=rl,cl=cl, error=err)
    

return(out)


    
}





samplernew  <- function(id,y,X,param,reps,subele=subele,bin_par=0.2,thin=100){
    pid = id
    param=as.matrix(param)
    label = unique(pid)
    nn = length(label)

    TT = length(pid)/nn
    
    X = as.matrix(X)
    ycurrent = y
    bin = bin_par*reps
    stored = matrix(NA,(((1-bin_par)*reps)/thin),length(id))


    S = sq(9)
    A = rbind(diag(3)%x%matrix(1,1,3),
              matrix(1,1,3)%x%diag(3))
    Tot = t(A%*%t(S))
    pow = 4^(5:0)
    ind = c(Tot%*%pow+1)
    
    

    i = 0; j = 0
    change = 0
    while(i<=(reps-1)){
        i = i + 1


        Y = matrix(ycurrent,TT)
        

            
        err=TRUE
        while(err){
            rl = sample(1:nrow(Y),3) #tmp$rl
            cl = sample(1:ncol(Y),3) #tmp$cl
            
            oss0 = Y[rl,cl]
            oss = c(oss0)
            
            oss1 = prop_rec_and_hex(oss,S,A,pow,ind)
            err=oss1$error
        }            

     
        YC = Y
        YC[rl,cl] = matrix(c(oss1$newmat),3)
        
        ycand = c(YC)

        
        alpha = min(c(1,(exp((ycand[!is.na(ycand)]-ycurrent[!is.na(ycand)])%*%(X[!is.na(ycand),]%*%param)))))
        
        
        if(alpha>runif(1)){
            change = change +1
            ycurrent = ycand
        }

        if(i>bin & i%%thin==0){
            j = j + 1
            stored[j,] = ycurrent
        }
        
        if(i%%10000==0){print(i)}
    }


    return(stored)

}


    
## Maximize the approximant Conditioanal Likelihood
MC_NR2 <- function(yy,XX,psi,Y_MC){
    error=FALSE
    XX = as.matrix(XX)
    y = yy[!is.na(yy)]
    X = XX[!is.na(yy),]
    
    psi = as.matrix(psi)
    sc = 0
    it = 0; lk = -Inf; lk0 = -Inf
    theta = psi
    
    cat(" |--------------|--------------|--------------|\n")
    cat(" |   iteration  |      lk      |    lk-lko    |\n")
    cat(" |--------------|--------------|--------------|\n")
    
    while(abs(lk-lk0)>10^-6 | it==0){
        
        it = it+1        
        if(it>1){theta= theta - iJ%*%sc} #update beta
        
        db = theta - psi
        xb = X%*%db

       
        cat("", sprintf("%12g", c(it, lk, lk - lk0)), "\n", sep = " | ")

        if(it>1 & (lk-lk0)<0){error=TRUE;cat("Non-concave function");return(list(err=error,theta=theta))}
        
        lk0 = lk

        repl = nrow(Y_MC)
       
        ratios = NULL
        scorenum = h1 = 0
        h2 = c()
        for(j in 1:repl){
            yyi = Y_MC[j,]
            yi = yyi[!is.na(yyi)]
            ratio_i = as.vector(exp(yi%*%xb))

            T_i = yi%*%X
            ratios = c(ratios,ratio_i)
            scorenum = scorenum + ((ratio_i)%*%T_i)/repl
            h1 = h1 + ratio_i*(t(T_i)%*%T_i)
        }


       
        meanMC = mean(ratios)
        lk = y%*%xb - log(meanMC)
        
        sc = t(as.matrix(y%*%X - (scorenum/mean(ratios))))

        J =  - ((h1/repl)/(mean(ratios))) + (t(scorenum/mean(ratios))%*%scorenum/(mean(ratios)))

        iJ = try(solve(J), silent = TRUE)
        if (inherits(iJ, "try-error")) {
            iJ = ginv(J)
            print("Inversion problems of the information matrix")
        }

       


    }
    cat(" |--------------|--------------|--------------|\n")
    
    return (list(theta=theta,J=J,lk = lk,err=error))
}

## The main function for MCMC-MLE of a double-CML estimator
MCMC_MLE <- function(id,y,X,psi,rep_mc,subele=15,bin_par=0.2,thin=100){

    ## Generate samples
    YY = samplernew(id,y,X,psi,rep_mc,subele,bin_par=bin_par,thin=thin)
    
    ## Maximize likelihood
    maxim =  try(MC_NR2(y,X,psi=psi,Y_MC = YY))
    if(inherits(maxim,"try-error")){
        err_mess= "MCMC-ML failed"; out = list(theta=NA,se_theta = NA,err_mess=err_mess);return(out)
    }else{
        
    
            
        theta = maxim$theta
        se = sqrt(diag(solve(-maxim$J)))
        
        out = list(theta=theta,se_theta = se)
        return(out)
    }
}

