
Plot_PCA<-function(data,samples,names_col,group_size,no_of_batches,size_of_batches,col_vector,title){
  
  symbols<-c()
  for (i in 1:no_of_batches){
    symbols<-c(symbols,c(rep(i,each=size_of_batches[i])))
  }
  
  pc<-prcomp(t(data))
  pc1<-round((pc$sdev[1])^2 / sum(pc$sdev^2),digits=2)
  pc2<-round((pc$sdev[2])^2 / sum(pc$sdev^2),digits=2) 
  pc3<-round((pc$sdev[3])^2 / sum(pc$sdev^2),digits=2) 
  cols<-rep(sample(col_vector,length(names_col)),each=group_size)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot( pc$x[,1:2],pch=symbols, cex=3,col=cols, main=title,xlab=paste('PC1=',pc1),ylab=paste('PC2=',pc2))
  #plot( pc$x[,2:3],pch=symbols, cex=3,col=cols, main=title,xlab=paste('PC2=',pc2),ylab=paste('PC3=',pc3))
  legend("topright", inset=c(-0.3,0),legend=names_col,pch=rep(unique(symbols),size_of_batches/group_size),col=unique(cols),cex=0.7)
}


Plot_violin<-function(data,samples, group_size){
  
  variable<-group_size-1
  for (i in seq(1,length(samples)-variable, by=group_size))
  {
    disp<-paste(samples[i],'_disp_raw',sep="")
    disp_v<-apply(data[,i:(i+variable)],1,function(x) (
      if (sum(x)>0){
        (sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))*100}
      else{200}
    ))
    assign(paste(disp),disp_v)
  }
  
  df_dispersion_raw<-NULL
  for (i in (1:length(names_col)))
  {
    df_dispersion_raw<-cbind(df_dispersion_raw,get(paste(names_col[i],'_disp_raw',sep='')))
  }
  colnames(df_dispersion_raw)<-paste(names_col,'_disp_raw',sep='')
  df_dispersion_raw<-as.data.frame(df_dispersion_raw)
  
  pr_below_10_raw<-c()
  pr_below_20_raw<-c()
  for (i in seq(1,length(df_dispersion_raw)))
  {
    pr_below_10_raw[i]<-(length(which(df_dispersion_raw[,i]<10))/dim(df_dispersion_raw)[1])*100
    pr_below_20_raw[i]<-(length(which(df_dispersion_raw[,i]<20))/dim(df_dispersion_raw)[1])*100
  }
  
  df.m <- reshape2::melt(df_dispersion_raw, id.vars = NULL)
  pl<-ggplot(df.m, aes(x = variable, y = value)) + geom_violin()
  pl + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}