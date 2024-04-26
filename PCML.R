jochmans_charb  <-  function(id,y,X,verbose=TRUE){

    rr = length(id)
    label = unique(id)
    nn = length(label)
    TT = rr/nn
    
    X = as.matrix(X)
    ncov = ncol(X)

    
    M = t(matrix(1:(nn*TT),TT))
    ntab = nn*(nn-1)*TT*(TT-1)/4
    ind1 = ind2 = ind3 = ind4 = rr = rep(0,ntab)
    j = 0

    for(i1 in 1:(nn-1))
        for(i2 in (i1+1):nn){
            for(t1 in 1:(TT-1))
                for(t2 in (t1+1):TT){
                    j = j+1
                    ind1[j] = M[i1,t1]; ind2[j] = M[i1,t2]
                    ind3[j] = M[i2,t1]; ind4[j] = M[i2,t2]
                }
        }
    IND = cbind(ind1,ind2,ind3,ind4)
    

    Z = (y[IND[,1]] - y[IND[,2]] - y[IND[,3]] + y[IND[,4]])/2
    
    #INDG = IND[abs(z)==1,]
    
    R = as.matrix(X[IND[,1],] - X[IND[,2],] - X[IND[,3],] + X[IND[,4],])

    r = as.matrix(R[abs(Z)==1 & !is.na(Z),])
    z = as.matrix(Z[abs(Z)==1 & !is.na(Z)])
    IN = as.matrix(IND[abs(Z)==1 & !is.na(Z),])
    #browser()
    if(verbose){
        cat("NR started\n")
        cat("\n")
    }
    
    beta = rep(0,ncov)
    sc = 0
    it = 0; lk = 0; lk0 = -Inf
    
    if(verbose){
    cat(" |--------------|--------------|--------------|\n")
    cat(" |   iteration  |      lk      |    lk-lko    |\n")
    cat(" |--------------|--------------|--------------|\n")
    }
    
    while(abs(lk-lk0)>10^-6 | it==0){

        it = it+1
        lk0 = lk
        rb = r%*%beta
        ex = exp(rb)
        F = ex/(1+ex)
        f = ex/((1+ex)^2)
        #browser()
        lk = sum(((z==1)*log(F) + (z==-1)*log(1-F))) #,na.rm=TRUE)

        fac = ((z==1)*(1-F) - (z==-1)*(F))
        scv = matrix(NA,length(z),ncov)
        for(cc in 1:ncov){
        scv[,cc] = fac * r[,cc]
        }
        
        sc = colSums(scv) #, na.rm=TRUE)
        
        J = 0
        #browser()
        
        for(j in 1:length(z)){
            hj = f[j]*(r[j,]%o%r[j,])
            J = J - hj
        }
                                        #scv = matrix(0,nrow(ISC),ncov)
        #browser()
        
        iJ = solve(J)

        beta = beta - iJ%*%sc
                                        # for idx
    

    #sc = colSums(scv)
    
    
    if(verbose){cat("", sprintf("%12g", c(it, lk, lk - lk0)), "\n", sep = " | ")}
    
    }
    
if(verbose){
    cat(" |--------------|--------------|--------------|\n")
    
    cat("NR is over\n")
}
    
###### Compute standard errors clustered at ij level

    ## Score :
                                        #browser()
    ## ii = 0
    OPG = matrix(0,ncov,ncov)
    for(i in 1:length(y)){

            clus = IN==i #& IN[,3]==j ## CONTROLLARE BENE
       #if(any(clus)){browser()}
        scvi = matrix(scv[rowSums(clus)>0,],,ncov)
        #browser()
        sscvi = colSums(scvi)
        OPG = OPG + sscvi%o%sscvi
       
    } # for i
    
    ## ii = 0
    ## OPG = matrix(0,ncov,ncov)
    ## for(i in unique(IN[,1])){
    ##     #browser()
    ##     clus = IN[,1]==i 
        
    ##     sscvi = colSums(matrix(scv[clus,],,ncov))
        
    ##     OPG = OPG + sscvi%o%sscvi
       
    ## } # for i
    
    iH = solve(J)
    se = sqrt(diag(-iH))
    ser = sqrt(diag(iH%*%OPG%*%iH))

    if(verbose){
        cat("Estimation Process is over\n")
        cat("\n")
    }
        #browser()
    out = list(beta=beta,scv=scv,J=J,se=se,ser=ser)

}


  ## sw0 = (ycurrent[ind1]==0) & (ycurrent[ind2]==1) & (ycurrent[ind3]==1) & (ycurrent[ind4]==0)
  ##   sw1 = (ycurrent[ind1]==1) & (ycurrent[ind2]==0) & (ycurrent[ind3]==0) & (ycurrent[ind4]==1)

  ##   len = length(c(which(sw0),which(sw1)))
    
  ##   
  ##   rr = NULL
  ##   for(j in 1:nrow(INDG)){
  ##       rj = X[INDG[j,1],] - X[INDG[j,2],] - X[INDG[j,3],] + X[INDG[j,4],]
  ##   }
