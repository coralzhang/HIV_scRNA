library(ggplot2)


### This function Threshold.fun2 calculates the threshold for an input vector xx. 
### below are parameters
### @ xx: it is the vector that you want to identify thresholds for
### @ kkt: it is the number of subpopulations/clusters 
### @ muv.init: this is a vector of length kkt (remember this), which is an initial value for the algorithm. You can set it where you think the peaks for each subpopulations are. 
### @ seed: this is a random seed. 
### @ xaxis.text: the title of x axis
### @ saving: an internal argument for better saving the plot
### To use this function, just copy paste this function in the consol for one time. 
### As mentioned, the algorithm is a random algorithm, it **could** be influenced by two things: change the ""seed"" or ""muv.init"". These are included so that the output is reproducible if different people run the function on different platforms. 

### Output
## @ threshold: the chosen thresholds; 
## @ means: the means of the idenfitied subpopulations/clusters 
## @ sds: the standard deviations of the identified subpopulations/clusters 
## @ likelihood: model selection metric; See details below
## @ BIC=BIC: model selection metric; See details below
## @ xxemp.v and lmdv1: only for debug purpose. 

Threshold.fun2=function(xx,kkt,seed=NULL,muv.init=NULL,xaxis.text="xx",saving=FALSE){
    if (is.null(seed)){seed=1234}
    set.seed(seed)
    
    if (is.null(muv.init)){
        if (kkt==2) {muv.init=as.vector(quantile(xx,prob=c(0.3,0.7)))}
        if (kkt==3) {muv.init=as.vector(quantile(xx,prob=c(0.05,0.5,0.9)))}
        if (kkt==4) {muv.init=as.vector(quantile(xx,prob=c(0.05,0.5,0.7,0.95)))}
        if (kkt>4) {muv.init=as.vector(quantile(xx,prob=c(0.03,seq(0.05,0.9,length.out=kkt-2),0.95)))}
    }
    tmp1=NULL
    out1=out2=1
    min1=max1=10
    ttmean1v=0.001 ## the predictive probability of classes
    ttmean33v=0.2 ## difference_between_means/max_sd
    tmp1v=c(NA) ## This is the cutpoint
    
    #xx0=xx; xx=scale(xx)

    xx=sort(xx)
    #while ( ( any(is.na(tmp1v))) | ( any(abs(ttmean1v)>0.99) | any(ttmean33v<0.3) ) ){ 
    while ( ( any(is.na(tmp1v))) | ( any(max(ttmean1v)>0.99) | any(ttmean33v<0.3) ) ){  ## this is added on Dec 28, 2022
        ## the following function is fit the Gaussian mixture model using the EM algorithm. 
        mixmdl1 <- mixtools::normalmixEM(xx,k=kkt,maxit=8000,mu=muv.init,maxrestarts=1000);
        #summary(mixmdl1)
        muv=mixmdl1$mu
        sdv=mixmdl1$sigma
        lmdv=mixmdl1$lambda
        postmat=mixmdl1$posterior
        post.as=apply(postmat,1,which.max)

        ord1=order(muv)
        muv1=muv[ord1]
        sdv1=sdv[ord1]
        lmdv1=lmdv[ord1]
        postmat1=postmat[,ord1]
        
        out1v=out2v=min1v=max1v=ttmean1v=ttmean33v=tmp1v=c()
        
        for (ii in 1:(kkt-1)){
            mean1=muv1[ii];mean2=muv1[ii+1]
            sd1=sdv1[ii];sd2=sdv1[ii+1]
            lambda1=lmdv1[ii];lambda2=lmdv1[ii+1];
    
            ttmean11=mean(postmat1[,ii]==apply(postmat1,1,max))
            ttmean33=abs(mean1-mean2)/max(sd1,sd2)
            mean3=mean(xx)
            sd3=sd(xx)
            
            A = -1/(sd1^2)+1/(sd2^2)
            B = 2*(-1*mean2/(sd2^2)+mean1/(sd1^2))
            C = mean2^2/(sd2^2)- mean1^2/(sd1^2)+log(sd2^2/(sd1^2)) + 2*log(lambda1)-2*log(lambda2)
            
            if (abs((sd1-sd2)/min(c(mean1,mean2)))<0.01){
                out1=out2=-C/B
            } else{ 
                Delta=B^2-4*A*C
                if (Delta>=0){
                    out1=(-1*B+sqrt(Delta))/(2*A);out2=(-1*B-sqrt(Delta))/(2*A)
                    a1=A*out1^2+B*out1+C
                    a2=A*out2^2+B*out2+C
                    if ((abs(a1)>0.001)|(abs(a2)>0.001)) stop("formula wrong")
                }
                else {
                    out1=out2=NA
                }
            }
            
            min1=min(mean1,mean2)
            max1=max(mean1,mean2)
            
            if (!is.na(out1) & !is.na(out2)){
                if (((out1>min1)&(out1<max1))){
                    tmp1=out1
                } else if (((out2>min1)&(out2<max1))) {
                    tmp1=out2
                } else {tmp1=NA}
            } else (tmp1=NA)
            
            tmp1v=c(tmp1v,tmp1)
            
            out1v=c(out1v,out1)
            out2v=c(out2v,out2)
            min1v=c(min1v,min1)
            max1v=c(max1v,max1)
            
            ttmean1v=c(ttmean1v,ttmean11)
            ttmean33v=c(ttmean33v,ttmean33)

        }
    }
    
    #plot(xx,post.as,type="l")
    
    xxemp.v=rep(0,length=kkt)
    for (jj in 1:kkt){
        xxemp.v[jj]=xx[max(which(post.as==jj))]
    }
        
    xxnew=seq(min(xx),max(xx),length.out=1000)
    ff=rep(0,length(xxnew))
    for (jj in 1:kkt){
        ffjj=dnorm(xxnew,mean=muv1[jj],sd=sdv1[jj])*lmdv1[jj]
        ff=ff+ffjj
    } 
    ylims=c(0,max(ff)*1.2)

#    if (save_figure==TRUE){
#        pdf(paste0(dir_to_figure,".pdf"))
#    }
    
    if (saving==FALSE){
        par(mfrow=c(1,2))
        hist(xx,nclass=100,freq=FALSE,ylim=ylims,main="", xlab=xaxis.text)
        for (jj in 1:kkt){
            par(new=TRUE)
            lines(x=xxnew,y=dnorm(xxnew,mean=muv1[jj],sd=sdv1[jj])*lmdv1[jj],col=jj+1,lwd=2)
        }
    
        main.title=paste0("Loglike=",round(mixmdl1$loglik,2))
    } else {
        main.title=""
        par(mfrow=c(1,1))
    }
    hist(xx,nclass=100,freq=FALSE,ylim=ylims,main=main.title, xlab=xaxis.text)
    par(new=TRUE)
    lines(x=xxnew,y=ff,col="grey",lwd=2)

    abline(v=muv1,lty=2,lwd=2)
    abline(v=tmp1v,lty=2,col=2,lwd=2)
    
#    if (save_figure==TRUE){
#        dev.off()
#    }
    
    #BIC=kkt*3*log(length(xx))-2*(mixmdl1$loglik)
    BIC=(kkt*3-1)*log(length(xx))-2*(mixmdl1$loglik)

    return(list(threshold=tmp1v,means=muv1,sds=sdv1,likelihood=mixmdl1$loglik,BIC=BIC, xxemp.v=xxemp.v,lmdv1=lmdv1)) ## the lmdv1 is added on Dec 21
}


remove.iqr=function(xx,nn=1.5,xaxis.text=""){
    Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)
    iqr=IQR(xx)
    up <-  Q[2]+nn*iqr # Upper Range  
    low<- Q[1]-nn*iqr # Lower Range
    out=c(low,up)
    par(mfrow=c(1,1))
    hist(xx,nclass=100,freq=FALSE,xlab=xaxis.text,main="")
    abline(v=up,lty=2,col=2,lwd=2)
    return(as.vector(out))
}

