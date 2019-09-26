
if("glmnet" %in% rownames(installed.packages()) == FALSE) {install.packages("glmnet")}
if("pls" %in% rownames(installed.packages()) == FALSE) {install.packages("pls")}
if("randomForest" %in% rownames(installed.packages()) == FALSE) {install.packages("randomForest")}
if("ncdf4" %in% rownames(installed.packages()) == FALSE) {install.packages("ncdf4")}
if("stringr" %in% rownames(installed.packages()) == FALSE) {install.packages("stringr")}

library(glmnet)
library(pls)
library(randomForest)
library(ncdf4)
library(stringr)

### Function which loads the database ###
load_proxies=function(path){
    
    spl=stringr::str_split(path,'\\.')
    if (spl[[1]][length(spl[[1]])]=='csv'){
        a=read.csv(path,header=T)
        if (ncol(a)==1){
            cnt=1
            seps=c(",",";","\t")
            while(ncol(a)==1){
                a=read.csv(path,sep=seps[cnt],header=T)
                cnt=cnt+1
            }
        }
    }

    if (spl[[1]][length(spl[[1]])]=='txt'){
        a=read.table(path,header=T)
        if (ncol(a)==1){
            cnt=1
            seps=c(',',';','\t')
            while(ncol(a)==1){
                a=read.table(path,sep=seps[cnt],header=T)
                cnt=cnt+1
            }
        }
    }
    datas=a
    return(datas)
}

### Substract the period y1-y2 from the database with all complete records available on this period ### 
subdf=function(df,y1,y2){
    years=df[,1]
    sub=df[,intersect(which(!is.nan(as.numeric(df[years==y1,]))),which(!is.nan(as.numeric(df[years==y2,]))))]
    sub=sub[sub[,1] %in% seq(y1,y2),]
    return(sub)
}

### Correlation test from McCarthy et al. 2015 (Default in this code) ###
test_cor=function(x,y,conf){
    ax=as.numeric(acf(x,plot=F)[[1]][2])
    ay=as.numeric(acf(y,plot=F)[[1]][2])
    r=cor(x,y)
    N=length(x)
    Neff=N*((1-ax*ay)/(1+ax*ay))
    TT=sqrt(Neff)*(r/sqrt(1-r^2))
    if (abs(TT)>qt(conf+(1-conf)/2,Neff)) return(T)
    else return(F)
}

### Correlation test from Bretherton et al. 1999 ###
test_cor2=function(x,y,conf){
    ax=as.numeric(acf(x,plot=F,length(x))[[1]])
    ay=as.numeric(acf(y,plot=F,length(y))[[1]])
    r=cor(x,y)
    N=length(x)
    tt=(-(N-1)):(N-1)
    
    ss=c()
    
    for (tau in tt){
        ss=c(ss,as.numeric((1-abs(tau)/N)*ax[abs(tau)]*ay[abs(tau)]))
    }
    
    ss=sum(ss)
    Neff1=N*((1-ax[2]*ay[2])/(1+ax[2]*ay[2]))
    Neff2=N/ss

    TT1=sqrt(Neff1)*(r/sqrt(1-r^2))
    TT2=sqrt(Neff2)*(r/sqrt(1-r^2))

    if (abs(TT2)>qt(1-((1-conf)/2),Neff2)) return(T)
    else return(F)
}


### Function which reads the file that contains the climate index ###
get_mdv=function(path){
    spl=stringr::str_split(path,'\\.')
    if (spl[[1]][length(spl[[1]])]=='csv'){
        fmov=read.csv(path,header=T)
        if (ncol(fmov)==1){
            cnt=1
            seps=c(',',';','\t')
            while(ncol(fmov)==1){
                fmov=read.csv(path,sep=seps[cnt],header=T)
                cnt=cnt+1
            }
        }
    }
    if (spl[[1]][length(spl[[1]])]=='txt'){
        fmov=read.table(path,header=T)  
        if (ncol(fmov)==1){
            cnt=1
            seps=c(',',';','\t')
            while(ncol(fmov)==1){
                fmov=read.table(path,sep=seps[cnt],header=T)
                cnt=cnt+1
            }
        }
    }
    time=fmov[,1]
    mov=fmov[,2]
    return(list(mov,time))
}

### Function that normalize a given time serie ###
normalize_ser=function(x) return((x-mean(x))/sqrt(var(x)))

### Substract the portion of the mdv that overlaps with the reconstruction period ###
trunc_mdv=function(list_mdv,y2){
    mdv=list_mdv[[1]]
    time=list_mdv[[2]]
    t2=seq(min(time),y2)
    mdv=mdv[time %in% t2]
    return(list(mdv,t2))
}

### Truncate the proxy database to a given period ###
trunc_df=function(df,time_sub){
    sub=df[df[,1] %in% time_sub,]
    return(sub)
}

### K-folds cross validation for the PLS method ###
kf_pls=function(method,x,y,blockstyle_cv=T,kcv=10){
    if (method=='pls'){
        rmses=rep(0,ncol(x))

        if (!blockstyle_cv) {
            ids=split(sample(1:nrow(x),nrow(x)), rep_len(1:kcv, nrow(x)))
        }
        else{
            ids=split(1:nrow(x),cut(1:nrow(x), kcv,labels=F))
        }

        for (idx in 1:length(ids)){
            id=as.numeric(ids[[idx]])

            xtrain=x[-id,]
            xtest=x[id,]
            ytrain=y[-id]
            ytest=y[id]
            

            for (k in 1:ncol(xtrain)){
                xtrain[,k]=(xtrain[,k]-mean(xtrain[,k]))/sqrt(var(xtrain[,k]))
                xtest=(xtest-mean(xtrain[,k]))/sqrt(var(xtrain[,k]))
            }
	    
            dfall=data.frame(ytrain,xtrain)

            names(xtest)=names(dfall)[2:ncol(dfall)]
            plsall=plsr(ytrain~.,data=dfall)
            preds=predict(plsall,xtest)

            for (j in 1:dim(preds)[3]){
                pred=preds[,,j]
                rmses[j]=rmses[j]+mean((pred-ytest)^2)
            }
        }
        
        rmses[rmses==0]=NA 
        rmses[max(which(!is.na(rmses)))]=NA 
        return(which.min(sqrt(rmses/nrow(x)))) 
    }
}

### Calculate the Nash-Sutcliffe Coefficient of Efficiency ###
nse=function(obs,pred){
    s1=sum((obs-pred)^2)
    s2=sum((obs-mean(obs))^2)
    nse=1-s1/s2
    return(nse)
}

### 10-folds cross validation for the PCR method ###
kf_pcr=function(method,x,y,blockstyle_cv=T,kcv=10){
    if (method=='pcr'){
        rmses=rep(0,ncol(x))
        print(dim(x))
        
        if (blockstyle_cv) {
            ids=split(sample(1:nrow(x),nrow(x)), rep_len(1:kcv, nrow(x)))
        }
        else{
            ids=split(1:nrow(x),cut(1:nrow(x), kcv,labels=F))
        }

        for (idx in 1:length(ids)){
            id=as.numeric(ids[[idx]])   
            xtrain=x[-id,]
            xtest=x[id,]
            ytrain=y[-id]
            ytest=y[id]
            for (k in 1:ncol(xtrain)){
                xtrain[,k]=(xtrain[,k]-mean(xtrain[,k]))/sqrt(var(xtrain[,k]))
                xtest[,k]=(xtest[,k]-mean(xtrain[,k]))/sqrt(var(xtrain[,k])) 
            }
            pca=prcomp(xtrain)
            pcs=pca$x

            for (j in 1:ncol(pcs)){
                pcss=pcs[,1:j]
                lmfit=lm(ytrain~.,data.frame(ytrain,pcss))
                coefs=lmfit$coefficients
                new_pcs=predict(pca,xtest)
                pred=coefs[1]
                for (j in 1:(length(coefs)-1)) pred=pred+coefs[j+1]*new_pcs[,j]
                rmses[j]=rmses[j]+mean((pred-ytest)^2)
            }
        }
        rmses[rmses==0]=NA
        print(rmses)
        print(which.min(sqrt(rmses/nrow(x))))
        rmses[max(which(!is.na(rmses)))]=NA 
        return(which.min(sqrt(rmses/nrow(x))))
    }
}

### Calculate the Root Mean Squared Error between x and y ###
RMSE=function(x,y){
    return(sqrt(mean((x-y)^2)))
}

### Main function ###
apply_rec=function(workdir='.',path_db,path_mode,y1,y2,method,R,freq_calib,tests=T,conf=0.9,seed=3,trace=T,blockstyle_cv=T,blockstyle_holdout=T,K_cv=10){
    trace=T
    y1=as.numeric(y1)
    y2=as.numeric(y2)

    freq_calib=as.numeric(freq_calib)
    R=as.numeric(R)

    conf=as.numeric(conf)
    seed=as.numeric(seed)

    setwd(workdir)

    datas=load_proxies(path_db)
    set.seed(seed)
    datas[is.na(datas)]=NaN
    t=seq(y1,y2)

    sub_datas=subdf(datas,y1,y2)

    datas=sub_datas
    
    tdatas=datas[,1]
    datas=datas[,-1]

    datasb=datas

    mdv_list=trunc_mdv(get_mdv(path_mode),y2)
    ty=mdv_list[[2]]
    Y=mdv_list[[1]]

    K_cv=ifelse(K_cv=='n',length(ty),as.numeric(K_cv))
    
    datasb=data.frame(datas)  

    Y=Y[ty %in% t]
    ty=ty[ty %in% t]
    
    ind_s=c()
    
    X=datas[t %in% ty,]

    val_cors=rep(NA,ifelse(blockstyle_holdout,round(nrow(X)*freq_calib)+1,R))
    incerts=rep(NA,ifelse(blockstyle_holdout,round(nrow(X)*freq_calib)+1,R))
    namdata=c()
    namdatas=c()
    nses=rep(NA,ifelse(blockstyle_holdout,round(nrow(X)*freq_calib)+1,R))
    rmses=rep(NA,ifelse(blockstyle_holdout,round(nrow(X)*freq_calib)+1,R))
    npx=rep(NA,ifelse(blockstyle_holdout,round(nrow(X)*freq_calib)+1,R))
    shapiros=rep(NA,ifelse(blockstyle_holdout,round(nrow(X)*freq_calib)+1,R))
    
    preds=array(NA,dim=c(ifelse(blockstyle_holdout,round(nrow(X)*freq_calib)+1,R),length(t)))

    my=mean(Y)
    sdy=sqrt(var(Y))
    ttrains=array(NA,c(ifelse(blockstyle_holdout,round(nrow(X)*freq_calib)+1,R),round(nrow(X)*freq_calib)))
    ttests=array(NA,c(ifelse(blockstyle_holdout,round(nrow(X)*freq_calib)+1,R),nrow(X)-round(nrow(X)*freq_calib)))

    if (method=='pcr'){
### Loop to calculate scores ###
        if (blockstyle_holdout){R=round(nrow(X)*freq_calib)+1}
        for (r in 1:R){

            if (trace) print(paste(round(r/R,3)*100,'% completed',sep=""))
            datas=datasb
            inds=c(F)
            
            
            datas=datasb
            
            if (!blockstyle_holdout){
                samp=sample(1:nrow(X),round((1-freq_calib)*nrow(X)))
            }
            else{
                if (r==1){samp=1:(nrow(X)-round(nrow(X)*freq_calib))}
                else{samp=samp+1}
            }

            
            Xtrain=X[-samp,]
            Xtest=X[samp,]
            Ytrain=Y[-samp]
            Ytest=Y[samp]
            tsamp=ty[-samp]
            tnsamp=ty[samp]
            if(sum(is.na(tsamp))==0 && sum(is.na(tnsamp))==0){
                

		ttrains[r,]=tsamp
		ttests[r,]=tnsamp

### selection of proxy records significantly correlated at the "conf" level with the climate index over the training period ###
		if (tests){
		    inds=c()
  	 	    for (j in 1:ncol(Xtrain)){
    	     	    	ser=Xtrain[,j]
	     		if (sum(!is.na(ser))>0){
                            if (test_cor(ser,Ytrain,conf)){
                                inds=c(inds,T)
                            }
                            else{
                                inds=c(inds,F)
                            }
    	     	    	}
	     		else{ inds=c(inds,F) }
      	 	    }
         	    Xtrain=Xtrain[,inds]
        	    datas=datas[,inds]
        	    Xtest=Xtest[,inds]
	 	}
	 	else {inds=c(T,T)}
                
                if (sum(inds)>0){
                    namdatas=c(namdatas,names(Xtrain))

                    npx[r]=ncol(Xtrain)
                    
                    dfall=data.frame(Ytrain,Xtrain)

### Each series is normalized to the mean and the standard deviation of the training sample ###
                    for (j in 1:ncol(Xtrain)){
                        Xtest[,j]=(Xtest[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                        datas[,j]=(datas[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                        Xtrain[,j]=(Xtrain[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                    }
                    q_opt=kf_pcr(method,Xtrain,Ytrain,kcv=K_cv)
                    
                    pca=prcomp(Xtrain)

                    prcs=pca$x[,1:q_opt]
                    
                    lmfit=lm(Ytrain~.,data=data.frame(prcs,Ytrain))
                    
                    coefs=lmfit$coefficients
                    
                    new_pcs=predict(pca,Xtest)
                    
                    pred=coefs[1]
                    
                    for (i in 2:length(coefs)) pred=pred+coefs[i]*new_pcs[,i-1]
                    
                    val_cors[r]=cor(pred,Ytest)
                    nses[r]=nse(Ytest,pred)
                    
                    pcs_all=predict(pca,datas)
                    rmses[r]=RMSE(pred,Ytest)
                    
                    pred_all=coefs[1]
                    
                    for (i in 2:length(coefs)) pred_all=pred_all+coefs[i]*pcs_all[,i-1]

                    pall=rep(NA,length(ty[-samp]))
                    for (ttt in 1:length(ty[-samp])){
                        pall[ttt]=pred_all[t==ty[-samp][ttt]]
                    }

                    incerts[r]=sqrt(sum((pall-Ytrain)^2)/(length(Ytrain)-2))

                    preds[r,]=pred_all

                    residuals=pall-Ytrain                                   
                    
                    shapiros[r]=shapiro.test(residuals)$p.value 

                }
            }
        } 
### Final reconstruction ###
        datas=datasb

        
        if (tests){
            inds=c()
            for (j in 1:ncol(X)){
                ser=X[,j]
                if (sum(!is.na(ser))>0){
                    if (test_cor(ser,Y,conf)){
                        inds=c(inds,T)
                    }
                    else{
                        inds=c(inds,F)
                    }
                }
                else{ inds=c(inds,F) }
            }
            datas=datas[,inds]
            X=X[,inds]
        }
        else {
            inds=c(T,T)
        }
        

        for (j in 1:ncol(datas)){
            datas[,j]=(datas[,j]-mean(X[,j],na.rm=T))/sqrt(var(X[,j],na.rm=T))
            X[,j]=(X[,j]-mean(X[,j],na.rm=T))/sqrt(var(X[,j],na.rm=T))
        }

        q_opt=kf_pcr(method,X,Y,kcv=K_cv)       

        pca=prcomp(X)

        prcs=pca$x[,1:q_opt]

        lmfit=lm(Y~.,data=data.frame(prcs,Y))

        coefs=lmfit$coefficients

        new_pcs=predict(pca,datas)

        pred_all=coefs[1]

        for (i in 2:length(coefs)) pred_all=pred_all+coefs[i]*new_pcs[,i-1]

        pred_fin=pred_all
        
        namdatas_fin=names(datas)

        r2_fin=nse(Y,pred_fin[t %in% ty])
        
        cor_fin=cor(Y,pred_fin[t %in% ty])
        
        rmse_fin=RMSE(Y,pred_fin[t %in% ty])

        residuals_fin=Y-pred_fin[t %in% ty]

        shap_fin=shapiro.test(residuals_fin)$p.value

        incerts_fin=sqrt(sum((pred_fin[t %in% ty]-Y)^2)/(length(Y)-2))
        
        nms_fin=names(datas)

        npx_fin=length(nms_fin)
    }

    if (method=='pls'){
### Loop to calculate scores ###
        if (blockstyle_holdout){R=round(nrow(X)*freq_calib)+1}
        for (r in 1:R){
            if (trace) print(paste(round(r/R,3)*100,'% completed',sep=""))
            datas=datasb
            inds=c(F)
            
            

            if (!blockstyle_holdout){
                samp=sample(1:nrow(X),round((1-freq_calib)*nrow(X)))
            }
            else{
                if (r==1){samp=1:(nrow(X)-round(nrow(X)*freq_calib))}
                else{samp=samp+1}
            }


            Xtrain=X[-samp,]
            Xtest=X[samp,]
            Ytrain=Y[-samp]
            Ytest=Y[samp]
            tsamp=ty[-samp]
            tnsamp=ty[samp]


            if(sum(is.na(tsamp))==0 && sum(is.na(tnsamp))==0){


                ttrains[r,]=tsamp
                ttests[r,]=tnsamp

                if (tests){
                    inds=c()
                    for (j in 1:ncol(Xtrain)){
                        
                        ser=Xtrain[,j]
                        
                        if (sum(!is.na(ser))>0){
                            
                            if (test_cor(ser,Ytrain,conf)){
                                inds=c(inds,T)
                            }
                            else{
                                inds=c(inds,F)
                            }
                        }
                        else{ inds=c(inds,F) }
                    }
                    Xtrain=Xtrain[,inds]
                    datas=datas[,inds]
                    Xtest=Xtest[,inds]
                }
                else {inds=c(T,T)}
                if (sum(inds)>0){
                    namdatas=c(namdatas,names(Xtrain))
                    
                    for (j in 1:ncol(Xtrain)){
                        Xtest[,j]=(Xtest[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                        datas[,j]=(datas[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                        Xtrain[,j]=(Xtrain[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                    }
                    dfall=data.frame(Ytrain,Xtrain)
                    
                    plsall=plsr(Ytrain~.,data=dfall)
                    
                    q_opt=kf_pls(method,Xtrain,Ytrain,kcv=K_cv)
                    names(Xtest)=names(dfall)[2:ncol(dfall)]
                    
                    pred=predict(plsall,Xtest)[,,q_opt]
                    rmses[r]=RMSE(pred,Ytest)      
                    val_cors[r]=cor(pred,Ytest)
                    incerts[r]=sqrt(sum((pred-Ytest)^2)/(length(Ytest)-2)) 
                    nses[r]=nse(Ytest,pred)
                    
                    pred_all=predict(plsall,datas)[,,q_opt]
                    
                    pall=rep(NA,length(ty[-samp]))
                    for (ttt in 1:length(ty[-samp])){
                        pall[ttt]=pred_all[t==ty[-samp][ttt]]
                    }

                    incerts[r]=sqrt(sum((pall-Ytrain)^2)/(length(Ytrain)-2))

                    preds[r,]=pred_all

                    residuals=pall-Ytrain                                      
                    
                    shapiros[r]=shapiro.test(residuals)$p.value 
                }
            }
        }
### Final reconstruction ###
        datas=datasb
	
        
        if (tests){
            inds=c()
            
            for (j in 1:ncol(X)){
                ser=X[,j]             
                
                if (sum(!is.na(ser))>0){
                    if (test_cor(ser,Y,conf)){
                        inds=c(inds,T)
                    }
                    else{
                        inds=c(inds,F)
                    }
                }
                else{ inds=c(inds,F) }
            }
            datas=datas[,inds]
            X=X[,inds]
        }
        else {
            inds=c(T,T)
        }
        print(dim(datas))
        print(dim(X))
        
        for (j in 1:ncol(datas)){
            datas[,j]=(datas[,j]-mean(X[,j],na.rm=T))/sqrt(var(X[,j],na.rm=T))
            X[,j]=(X[,j]-mean(X[,j],na.rm=T))/sqrt(var(X[,j],na.rm=T))
        }

        print(length(Y))
        print(dim(X))
        dfall=data.frame(Y,X)

        plsall=plsr(Y~.,data=dfall)
        print(dim(X))
        q_opt=kf_pls(method,X,Y,kcv=K_cv)
        names(datas)=names(dfall)[2:ncol(dfall)]

        pred_all=predict(plsall,datas)[,,q_opt]

        pred_fin=pred_all

        namdatas_fin=names(datas)

        r2_fin=nse(Y,pred_fin[t %in% ty])
        
        cor_fin=cor(Y,pred_fin[t %in% ty])
        
        rmse_fin=RMSE(Y,pred_fin[t %in% ty])

        nms_fin=names(datas)

        npx_fin=length(nms_fin)

        residuals_fin=Y-pred_fin[t %in% ty]

        shap_fin=shapiro.test(residuals_fin)$p.value

        incerts_fin=sqrt(sum((pred_fin[t %in% ty]-Y)^2)/(length(Y)-2))
    }


    if (method=='enet'){
### Loop to calculate scores ###
        if (blockstyle_holdout){R=round(nrow(X)*freq_calib)+1}
        for (r in 1:R){


            if (trace) print(paste(round(r/R,3)*100,'% completed',sep=""))
            inds=c(F)


            datas=datasb

            if (!blockstyle_holdout){
                samp=sample(1:nrow(X),round((1-freq_calib)*nrow(X)))
            }
            else{
                if (r==1){samp=1:(nrow(X)-round(nrow(X)*freq_calib))}
                else{samp=samp+1}
            }


            Xtrain=X[-samp,]
            Xtest=X[samp,]
            Ytrain=Y[-samp]
            Ytest=Y[samp]
            tsamp=ty[-samp]
            tnsamp=ty[samp]
            print(tsamp)
            print('?')
            if(sum(is.na(tsamp))==0 && sum(is.na(tnsamp))==0){
                
                ttrains[r,]=tsamp[!is.na(tsamp)] 
                ttests[r,]=tnsamp[!is.na(tnsamp)] 

                
                if (tests){
                    inds=c()
                    for (j in 1:ncol(Xtrain)){
                        ser=Xtrain[,j]
                        if (sum(!is.na(ser))>0){
                            if (test_cor(ser,Ytrain,conf)){
                                inds=c(inds,T)
                            }
                            else{
                                inds=c(inds,F)
                            }
                        }
                        else{ inds=c(inds,F) }
                    }
                    datas=datas[,inds]
                    Xtest=Xtest[,inds]
                    Xtrain=Xtrain[,inds]
                    
                }
                else {
                    inds=c(T,T)
                }
                if (sum(inds)>1){

                    npx[r]=ncol(Xtrain)      
                    namdatas=c(namdatas,names(Xtrain))

                    for (j in 1:ncol(Xtrain)){
                        Xtest[,j]=(Xtest[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                        datas[,j]=(datas[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                        Xtrain[,j]=(Xtrain[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                        
                    }

                    dfall=data.frame(Ytrain,Xtrain)
                    
### K-folds cross validation for the Elastic net method ###
                    a=seq(0,1,0.1)
                    lbds=seq(0,5,0.2)

                    mse=array(0,c(length(a),length(lbds)))
                    
                    if (!blockstyle_cv) {
                        ids=split(sample(1:nrow(Xtrain),nrow(Xtrain)), rep_len(1:K_cv,nrow(X)))
                    }
                    else{
                        ids=split(1:nrow(Xtrain),cut(1:nrow(Xtrain),K_cv,labels=F))
                    }

                    for (idx in 1:length(ids)){
                        id=as.numeric(ids[[idx]])
                        for (ii in 1:length(a)){
                            for (jj in 1:length(lbds)){
                                
                                xloo=Xtrain[id,]
                                yloo=Ytrain[id]
                                xmod=Xtrain[-id,]
                                ymod=Ytrain[-id]

                                for (col in 1:ncol(xmod)){
                                    xloo[,col]=(xloo[,col]-mean(xmod[,col]))/sqrt(var(xmod[,col]))
                                    xmod[,col]=(xmod[,col]-mean(xmod[,col]))/sqrt(var(xmod[,col]))
                                }		  

                                md=glmnet(data.matrix(xmod),data.matrix(ymod),lambda=lbds[jj],alpha=a[ii])
                                pd=predict(md,data.matrix(xloo))
                                mse[ii,jj]=mse[ii,jj]+mean((as.numeric(pd)-yloo)^2)
                            }
                        }
                    }

                    inda=ifelse((which.min(mse) %% length(a))==0,length(a),which.min(mse) %% length(a))
                    aopt=a[inda]
                    lbdopt=lbds[floor((which.min(mse)-1)/length(a))+1]

                    mdf=glmnet(data.matrix(Xtrain),data.matrix(Ytrain),lambda=lbdopt,alpha=aopt)
                    
                    pred=predict(mdf,data.matrix(Xtest))
                    pred=as.numeric(pred)

                    rmses[r]=RMSE(as.numeric(pred),Ytest)
                    val_cors[r]=cor(as.numeric(pred),Ytest)
                    
                    nses[r]=nse(Ytest,pred)   

                    pred_all=predict(mdf,data.matrix(datas))
                    pred_all=as.numeric(pred_all)
                    pall=rep(NA,length(ty[-samp]))
                    for (ttt in 1:length(ty[-samp])){
                        pall[ttt]=pred_all[t==ty[-samp][ttt]]
                    }

                    incerts[r]=sqrt(sum((pall-Ytrain)^2)/(length(Ytrain)-2))

                    preds[r,]=pred_all

                    residuals=pall-Ytrain                                                                      
                    shapiros[r]=shapiro.test(residuals)$p.value 

                }
            }
        }
### Final reconstruction ###
        datas=datasb

        if (tests){
            inds=c()
            for (j in 1:ncol(X)){
                ser=X[,j]
                if (sum(!is.na(ser))>0){
                    if (test_cor(ser,Y,conf)){
                        inds=c(inds,T)
                    }
                    else{
                        inds=c(inds,F)
                    }
                }
                else{ inds=c(inds,F) }
            }
            datas=datas[,inds]
            X=X[,inds]
        }
        else {
            inds=c(T,T)
        }
        
        
        for (j in 1:ncol(datas)){
            datas[,j]=(datas[,j]-mean(X[,j],na.rm=T))/sqrt(var(X[,j],na.rm=T))
            X[,j]=(X[,j]-mean(X[,j],na.rm=T))/sqrt(var(X[,j],na.rm=T))
        }

        dfall=data.frame(Y,X)

### K-folds cross validation for the Elastic net method ###
        a=seq(0,1,0.1)
        lbds=seq(0,5,0.2)

        mse=array(0,c(length(a),length(lbds)))

        if (!blockstyle_cv) {
            ids=split(sample(1:nrow(X),nrow(X)), rep_len(1:K_cv, nrow(X)))
        }
        else{
            ids=split(1:nrow(X),cut(1:nrow(X), K_cv,labels=F))
        }

        for (idx in 1:length(ids)){
            id=as.numeric(ids[[idx]])
            for (ii in 1:length(a)){
                for (jj in 1:length(lbds)){


                    xloo=X[id,]
                    yloo=Y[id]
                    xmod=X[-id,]
                    ymod=Y[-id]

                    for (col in 1:ncol(xmod)){
                        xmod[,col]=(xmod[,col]-mean(xmod[,col]))/sqrt(var(xmod[,col]))
                        xloo[,col]=(xloo[,col]-mean(xmod[,col]))/sqrt(var(xmod[,col]))
                    }

                    md=glmnet(data.matrix(xmod),data.matrix(ymod),lambda=lbds[jj],alpha=a[ii])
                    pd=predict(md,data.matrix(xloo))
                    mse[ii,jj]=mse[ii,jj]+mean((as.numeric(pd)-yloo)^2)
                }
            }
        }

        inda=ifelse((which.min(mse) %% length(a))==0,length(a),which.min(mse) %% length(a))
        aopt=a[inda]
        lbdopt=lbds[floor((which.min(mse)-1)/length(a))+1]

        mdf=glmnet(data.matrix(X),data.matrix(Y),lambda=lbdopt,alpha=aopt)

        pred_all=predict(mdf,data.matrix(datas))
        pred_all=as.numeric(pred_all)
        pred_fin=pred_all      
        namdatas_fin=names(datas)

        r2_fin=nse(Y,pred_fin[t %in% ty])
        
        cor_fin=cor(Y,pred_fin[t %in% ty])
        
        rmse_fin=RMSE(Y,pred_fin[t %in% ty])
        nms_fin=names(datas)

        npx_fin=length(nms_fin)
        residuals_fin=Y-pred_fin[t %in% ty]

        shap_fin=shapiro.test(residuals_fin)$p.value

        incerts_fin=sqrt(sum((pred_fin[t %in% ty]-Y)^2)/(length(Y)-2))}
    
    if (method=='rf'){
### Loop to calculate scores ###
        if (blockstyle_holdout){R=round(nrow(X)*freq_calib)+1}
        for (r in 1:R){
            
            if (trace) print(paste(round(r/R,3)*100,'% completed',sep=""))
            datas=datasb

            inds=c(F)

            datas=datasb


            if (!blockstyle_holdout){
                samp=sample(1:nrow(X),round((1-freq_calib)*nrow(X)))
            }
            else{
                if (r==1){samp=1:(nrow(X)-round(nrow(X)*freq_calib))}
                else{samp=samp+1}
            }


            Xtrain=X[-samp,]
            Xtest=X[samp,]
            Ytrain=Y[-samp]
            Ytest=Y[samp]
            tsamp=ty[-samp]
            tnsamp=ty[samp]

            if(sum(is.na(tsamp))==0 && sum(is.na(tnsamp))==0){      	      
                
                ttrains[r,]=tsamp
                ttests[r,]=tnsamp       

                if (tests){
                    inds=c()
                    for (j in 1:ncol(Xtrain)){
                        ser=Xtrain[,j]
                        if (sum(!is.na(ser))>0){
                            if (test_cor(ser,Ytrain,conf)){
                                inds=c(inds,T)
                            }
                            else{
                                inds=c(inds,F)
                            }
                        }
                        else{ inds=c(inds,F) }
                    }
                    Xtrain=Xtrain[,inds]
                    datas=datas[,inds]
                    Xtest=Xtest[,inds]
                }
                else {inds=c(T,T)}
                
                if (sum(inds)>0){


                    print(ncol(datas))
                    for (j in 1:ncol(Xtrain)){
                        Xtest[,j]=(Xtest[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                        datas[,j]=(datas[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                        Xtrain[,j]=(Xtrain[,j]-mean(Xtrain[,j]))/sqrt(var(Xtrain[,j]))
                    }
                    
                    dfall=data.frame(Ytrain,Xtrain)
                    npx[r]=ncol(Xtrain)
                    
### 10-folds cross validation for Random Forest ###
                    print(ncol(Xtrain))
                    mtrys=seq(round(ncol(Xtrain)*0.1),round(ncol(Xtrain)*0.75))
                    print(mtrys)
                    if (min(mtrys)<2){
                        mtrys=2:mtrys[length(mtrys)]	
                    }
                    scs=rep(0,length(mtrys))
                    
                    if (!blockstyle_cv) {
                        ids=split(sample(1:nrow(Xtrain),nrow(Xtrain)), rep_len(1:K_cv, nrow(Xtrain)))
                    }
                    else{
                        ids=split(1:nrow(Xtrain),cut(1:nrow(Xtrain), K_cv,labels=F))
                    }

                    for (idx in 1:length(ids)){
                        id=as.numeric(ids[[idx]])
                        for (ms in 1:length(mtrys)){
                            xloo=Xtrain[id,]
                            yloo=Ytrain[id]
                            ymod=Ytrain[-id]
                            xmod=Xtrain[-id,]

                            for (col in 1:ncol(xmod)){
                                xmod[,col]=(xmod[,col]-mean(xmod[,col]))/sqrt(var(xmod[,col]))
                                xloo[,col]=(xloo[,col]-mean(xmod[,col]))/sqrt(var(xmod[,col]))
                            }

                            rfst=randomForest(xmod,ymod,mtry=mtrys[ms],ntree=200)
                            pd=predict(rfst,xloo)
                            scs[ms]=scs[ms]+mean((pd-yloo)^2)
                        }
                    }
                    mtryopt=mtrys[which.min(scs)]

                    rfopt=randomForest(Xtrain,Ytrain,mtry=mtryopt,ntree=200)
                    
                    pred=as.numeric(predict(rfopt,Xtest))

                    rmses[r]=RMSE(pred,Ytest)      
                    val_cors[r]=cor(pred,Ytest)

                    nses[r]=nse(Ytest,pred)
                    
                    pred_all=predict(rfopt,datas)
                    
                    pall=rep(NA,length(ty[-samp]))
                    for (ttt in 1:length(ty[-samp])){
                        pall[ttt]=pred_all[t==ty[-samp][ttt]]
                    }

                    incerts[r]=sqrt(sum((pall-Ytrain)^2)/(length(Ytrain)-2))

                    preds[r,]=pred_all

                    residuals=pall-Ytrain
                    shapiros[r]=shapiro.test(residuals)$p.value
                    
                    namdatas=c(namdatas,names(Xtrain))

                }       
            }
        }
        
### Final reconstruction ###
        print('Performing final reconstruction...')
        datas=datasb

        if (tests){
            inds=c()
            for (j in 1:ncol(X)){
                ser=X[,j]
                if (sum(!is.na(ser))>0){
                    if (test_cor(ser,Y,conf)){
                        inds=c(inds,T)
                    }
                    else{
                        inds=c(inds,F)
                    }
                }
                else{ inds=c(inds,F) }
            }
            datas=datas[,inds]
            X=X[,inds]
        }
        else {
            inds=c(T,T)
        }
        print(ncol(datas))    

        for (j in 1:ncol(datas)){
            datas[,j]=(datas[,j]-mean(X[,j],na.rm=T))/sqrt(var(X[,j],na.rm=T))
            X[,j]=(X[,j]-mean(X[,j],na.rm=T))/sqrt(var(X[,j],na.rm=T))
        }

### 10-folds cross validation for Random Forest ###
        mtrys=seq(as.integer(ncol(Xtrain))*0.1,as.integer(ncol(Xtrain)*0.75))
        if (min(mtrys)<2){
            mtrys=2:mtrys[length(mtrys)]
        }
        scs=rep(0,length(mtrys))

        

        if (!blockstyle_cv) {
            ids=split(sample(1:nrow(X),nrow(X)), rep_len(1:K_cv, nrow(X)))
        }
        else{
            ids=split(1:nrow(X),cut(1:nrow(X), K_cv,labels=F))
        }

        for (idx in 1:length(ids)){
            id=as.numeric(ids[[idx]])
            for (ms in 1:length(mtrys)){
                xloo=X[id,]
                yloo=Y[id]
                ymod=Y[-id]
                xmod=X[-id,]

                for (col in 1:ncol(xmod)){
                    xmod[,col]=(xmod[,col]-mean(xmod[,col]))/sqrt(var(xmod[,col]))
                    xloo[,col]=(xloo[,col]-mean(xmod[,col]))/sqrt(var(xmod[,col]))
                }

                rfst=randomForest(xmod,ymod,mtry=mtrys[ms],ntree=200)
                pd=predict(rfst,xloo)
                scs[ms]=scs[ms]+mean((pd-yloo)^2)
            }
        }
        mtryopt=mtrys[which.min(scs)]

        rffin=randomForest(X,Y,mtry=mtryopt,ntree=200)

        pred_fin=predict(rffin,datas)
        
        namdatas_fin=names(datas)

        r2_fin=nse(Y,pred_fin[t %in% ty])
        
        cor_fin=cor(Y,pred_fin[t %in% ty])
        
        rmse_fin=RMSE(Y,pred_fin[t %in% ty])
        nms_fin=names(datas)

        npx_fin=length(nms_fin)
        residuals_fin=Y-pred_fin[t %in% ty]

        shap_fin=shapiro.test(residuals_fin)$p.value

        incerts_fin=sqrt(sum((pred_fin[t %in% ty]-Y)^2)/(length(Y)-2))
    }
    pred_base=pred_fin

### Renormalizing the reconstruction to the mean and the standard deviation of the climate index such as in Ortega et al. (2015, Nature Geosci.)###
### The code provide as output both the renormalized reconstruction and the original reconstruction ###  
    pred_fin[!t %in% ty]=(pred_fin[!t %in% ty]-mean(pred_fin[!t %in% ty]))/sqrt(var(pred_fin[!t %in% ty]))
    pred_fin[!t %in% ty]=pred_fin[!t %in% ty]*sqrt(var(pred_fin[t %in% ty]))+mean(pred_fin[t %in% ty])
    pred_fin=(pred_fin-mean(pred_fin))/sqrt(var(pred_fin))
    pred_fin=pred_fin*sdy+my
    
    namdata=unique(namdatas)
    
    return(list(pred_fin,pred_base,val_cors,rmses,npx,namdatas,preds,nses,incerts,ttrains,ttests,shapiros,namdatas_fin,r2_fin,cor_fin,rmse_fin,incerts_fin,shap_fin,npx_fin,nms_fin))
}
