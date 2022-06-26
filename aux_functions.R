#
yearly_mean<-function(data,ext=FALSE){
ny<-length(data)/12
if(ext) ny<-round(length(data)/12)
ret<-array(0,dim=c(ny))
for (i in 1:ny){
 i1<-(i-1)*12+1
 i2<-(i)*12
 ret[i]<-mean(data[i1:i2],na.rm=TRUE)
}
return(ret)
}
#
get_fld<-function(f,v,start=NULL,count=NULL){
    library(ncdf4)
    nc<-nc_open(f)
    if(is.null(start)||is.null(count)){
    var<-ncvar_get(nc,v)
    } else {
    var<-ncvar_get(nc,v,start=start,count=count)
    }
    nc_close(nc)
    return(var)
}
#
norm<-function(dt){ {dt-mean(dt)}/sd(dt) }
#
anom<-function(dt){ {dt-mean(dt)}}
#
ohct<-function(ohc){
  if( length(dim(ohc)) == 0 ) {
	  nl<-1
	  d1<-TRUE
  } else {
   	nl<-dim(ohc)[2]
	  d1<-FALSE
  }
  ar<-4.42e+14
  ar<- 5.102225e+14
  res<-array(0,dim=nl)
  for (i in 1:nl){
    if(d1) {
      tmp<-ohc
    } else {
      tmp<-ohc[,i]
    }
    res[i]<-lm(y~x,data=data.frame(y=tmp,x=1:length(tmp)))$coefficients[2]
  }
  res<-1.E+21*res/(365.2425*86400*ar)
  return(res)
}
#
mcat<-function(str) cat(paste(str,"\n",sep=""))

do_polisd<-function(x,m,s,col="grey",lty=0,den=NULL){
nt<-length(x)

polygon(c(x,x[seq(nt,1,by=-1)]),c(m-s,m[seq(nt,1,by=-1)]),col=col,lty=lty)
polygon(c(x,x[seq(nt,1,by=-1)]),c(m+s,m[seq(nt,1,by=-1)]),col=col,lty=lty)

}


add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
                     function(x)
                       rgb(x[1], x[2], x[3], alpha=alpha))
}

mylist<-function(chr){
  tmp<-system(paste("ls",chr),intern=TRUE)
  return(as.character(tmp))
}

