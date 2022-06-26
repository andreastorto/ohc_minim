#
source("aux_functions.R")
#
ll<-c(300,700,2000,6000)

levl<-c("0-300","300-700","700-2000","2000-6000")

# List of reanalyses taken
nen<-c("cglo","oras","ecco")
lat_bnd<-"90S-90N"
ne<-length(nen)
nl<-length(ll)
# Reanalysis years
yy<-(1993:2018)+0.5
nt<-length(yy)
n_col<-2
test<-0

data<-array(0,dim=c(nl,ne,nt))
datam<-array(0,dim=c(nl,nt))

#  Read ensemble products
for (e in 1:ne){
 for (l in 1:nl){
   ex<-nen[e]
   lv<-ll[l]
   ff<-sprintf("data_export/ohc_%s_0-%d_%s.dat",ex,lv,lat_bnd)
   print(ff)
   hc<-as.matrix(read.table(ff))[1:nt,n_col]
   data[l,e,]<-hc
 }
}

# Calculate ensemble mean anomalies
for (l in 1:nl){
   for (t in 1:nt){
	   datam[l,t]<-mean(data[l,,t])
   }
   datam[l,]<-datam[l,]-mean(datam[l,])
}

# Decumulate vertically (reanalyses provide OHC from surface to level z)
datan<-data
datamf<-datam
for (l in 2:nl){
	datan[l,,]<-datan[l,,]-data[l-1,,]
	datamf[l,]<-datamf[l,]-datam[l-1,]
}

for (e in 1:ne){
	for (l in c(nl-1,nl)){
	tmp<-datan[l,e,]
	print(paste(l,nen[e],ohct(tmp)))
	}
}

# Calculate and plot Anomalies
dataa<-datan
png("bg_ohc.png",width=800,height=800)
lc<-c("red","orange","blue","grey","cyan","green3","green","purple")
par(mfrow=c(2,2),cex=1.2)
for (l in 1:nl){
	for (e in 1:ne){
		dataa[l,e,]<-dataa[l,e,]-mean( dataa[l,e,] )
	}
	plot(0,type="n",xlim=range(yy),ylim=range(dataa[l,,]),main=paste("Background OHC",levl[l]),ylab="ZJ",xlab="")
	for (e in 1:ne) lines(yy,dataa[l,e,],col=lc[e],lwd=2)
	emean_000<-array(0,dim=dim(data)[3])
	for (t in 1:dim(data)[3]) emean_000[t]<-mean( dataa[l,,t] )
	lines(yy,emean_000,lwd=3)
	if(l==1) legend("bottomright",c(nen,"EnsMean"),fill=c(lc,"black"),cex=1.0,bty="n")
	box(lwd=2)
}
dev.off()

# Ensemble anomalies
for (l in 1:nl){
	for (t in 1:nt){
		datan[l,,t]<-datan[l,,t]-mean( datan[l,,t] )
	}
}

# B-matrix: vertical covariances
B<-array(0,dim=c(nl,nl))
Bsub<-array(0,dim=c(nl,nl))
Bsub2<-array(0,dim=c(nl,nl))

for (l in 1:nl){
	 print(ll[l])
	 print(sd( datan[l,,] ))
	 print(100*sd( datan[l,,]/sd(data[l,,]) ))
	 print(cor( c(datan[l,,]),  c(datan[1,,]) ))
	 for (k in 1:nl){
		 B[k,l]<-cov( c(datan[k,,]),  c(datan[l,,]) )*25
		 Bsub[k,l]<-cov( c(datan[k,,seq(1,nt,by=2)]),  c(datan[l,,seq(1,nt,by=2)]) )*25
		 Bsub2[k,l]<-cov( c(datan[k,,seq(1,nt,by=3)]),  c(datan[l,,seq(1,nt,by=3)]) )*25
	 }
}

dim(datan)
datan2e<-array(0,c(4,26*3))
for (l in 1:nl){
	kp<-0
	for (t in 1:26){
	for (e in 1:3){
		kp<-kp+1
		datan2e[l,kp]<-datan[l,e,t]
	}
	}
}

# Calculate EOFs
pca<-La.svd( t(datan2e) )
neofs<-4
nc<-26*3
d<-pca$d
eigenv<-d[1:neofs]/sqrt(nc-1)
eigenv.perc<-100*d[1:neofs]/sum(d[d>0])
eigenv.perc
cumsum(eigenv.perc)
eigenv.perc<-100*d[1:neofs]^2/sum(d[d>0]^2)
eigenv.perc
cumsum(eigenv.perc)

# Write some outputs
write.table(B,"Bmatrix.dat", row.names = FALSE, col.names = FALSE)
write.table(Bsub,"Bsubmatrix.dat", row.names = FALSE, col.names = FALSE)
write.table(Bsub2,"Bsub2matrix.dat", row.names = FALSE, col.names = FALSE)
B
Bsub
Bsub2

range(B)

m<-200000
sc<-seq(-m,m,length.out=41)
colmm<-colorRampPalette(c("blue","green","white","yellow","red"))
lc<-colmm(length(sc)-1)
lc[length(lc)/2]<-"white"
lc[length(lc)/2+1]<-"white"
png("B.png")
par(cex=1.3)
image(1:nl,1:nl,B,col=lc,breaks=sc)
box(lwd=2)
dev.off()

# Temporal autocorrelations
for (l in 1:nl){
	 print(ll[l])
	 for (e in 1:ne){
		 acf.res<-acf( datan[l,e,], plot=FALSE )
		 if(e==1&&l==1) {
			 res<-acf.res$acf/(ne*nl)
			 lag<-acf.res$lag
		 } else {
			 res<-res+acf.res$acf/(ne*nl)
		 }
                 if(e==1&&l==1) rest<-array(0,dim=c(length(res),nl))
                 rest[,l]<-rest[,l]+acf.res$acf/ne
		}
}
print(res)
print(rest)
resp<-rest
resp[resp<0]<-0
lag
dt<-seq(-max(lag),max(lag))
np<-dim(resp)[1]
np
lag
dim(resp)
length( c(resp[seq(np,2,by=-1),l],resp[,l]) )
resp2<-array(0,dim=c(np*2-1,nl))

for (l in 1:nl) resp2[,l]<-c(resp[seq(np,2,by=-1),l],resp[,l])

lc<-c("orange","red","blue","green3")
png("AC.png")
par(cex=1.3)
plot(dt,resp2[,1],lwd=2,bty="n",type="n")
for (l in 1:nl) lines(dt,resp2[,l],lwd=3,col=lc[l])
legend("topright",levl,fill=lc,bty="n")
axis(1)
axis(2)
dev.off()

# Find OHC vs steric relationship (using reanalyses data)

ny<-2018-1993+1

steric.cm<-array(dim=c(ny,ne))
for (i in 1:ne){
   fscm<-mylist(paste("data_export/steric_",nen[i],"_0-6000_",lat_bnd,".dat",sep=""))
   tmp<-read.table(fscm)
   steric.cm[,i]<-tmp[1:ny,n_col]*1000
}

x1<-NULL 
x2<-NULL 
x3<-NULL 
x4<-NULL 
yy<-NULL 

stregl<-c(1,3)

for (i in stregl){

   x1<-c(x1,c(dataa[1,i,1:ny]))
   x2<-c(x2,c(dataa[2,i,1:ny]))
   x3<-c(x3,c(dataa[3,i,1:ny]))
   x4<-c(x4,c(dataa[4,i,1:ny]))
   yy<-c(yy,c(steric.cm[,i]))

}

# Linear model fitting and calculation of residual errors

summary ( lm ( y ~ x1+x2+x3+x4  , data=data.frame( y=yy, x1=x1, x2=x2, x3=x3, x4=x4) ) )
summary ( lm ( y ~ x1+x2+x3+x4-1, data=data.frame( y=yy, x1=x1, x2=x2, x3=x3, x4=x4) ) )
alm<-lm ( y ~ x1+x2+x3+x4, data=data.frame( y=yy, x1=x1, x2=x2, x3=x3, x4=x4) )$coefficients
ypred<-alm[1]+alm[2]*x1+alm[3]*x2+alm[4]*x3+alm[5]*x4
sqrt(mean( (ypred-yy)^2))
rerr<-sd(ypred-yy)
rerr

a.lm<-lm ( y ~ x1+x2+x3+x4-1, data=data.frame( y=yy, x1=x1, x2=x2, x3=x3, x4=x4) )$coefficients
a.lm
c.lm<-a.lm
write(c.lm,file="steric_regression.dat")
write(rerr,file="steric_regression.dat",append=TRUE)
dim(a.lm)<-c(1,nl)

# Read observed steric (altimetry minus gravimetry)

datin<-as.matrix(read.table("data_export/steric_tot.dat"))
obs<-datin[,2]
obs<-obs-mean(obs)
lm ( y ~ x, data=data.frame( y=obs, x=1:length(obs)) )$coefficients

#  sigma(altim. 1y)  + bias(altim.) + GRACE error + represent.
steric_err<-sqrt( 1.3^2 +0.5^2 + 0.15^2 + rerr^2 + 1.3^2)
sqrt(1.3^2+0.5^2)
steric_err
rerr

# EnsMea Anomalies on the common period
dtb<-datamf[,11:26]
write.table(datamf,file="ooc_bg.dat")
for (l in 1:nl) dtb[l,]<-dtb[l,]-mean(dtb[l,])
nt<-dim(dtb)[2]
dim(obs)<-c(1,nt)
nv<-nt*nl

dx<-array(0,dim=c(nv,1))
xb<-array(0,dim=c(nv,1))

# Background: ensemble mean anomalies
k1<-0
for (t in 1:nt){
        for (l in 1:nl) {
                k1<-k1+1
                xb[k1,1]<-dtb[l,t]
        }
}

# Full background-error covariance matrix
Bt<-array(0,dim=c(nl*nt,nl*nt))
k1<-0
for (t in 1:nt){
         for (l in 1:nl) {
                k1<-k1+1
                k2<-0
                for (u in 1:nt){
                      for (k in 1:nl) {
                        k2<-k2+1
                        dt<-abs(u-t)+1

                        if(dt<10){
                                crt<-max(0,res[dt])
				if(crt<0.1) crt <-0
                                Bt[k1,k2]<-B[k,l]*crt
				print(paste(t,l,u,k,crt,B[k,l],Bt[k1,k2]))
                        }
                      }
                }
        }
}

max(Bt)
isSymmetric(Bt)

write.table(Bt,"Btmatrix.dat", row.names = FALSE, col.names = FALSE)
write.table(res,"autocor.dat", row.names = FALSE, col.names = FALSE)

sc<-seq(-m,m,length.out=41)
colmm<-colorRampPalette(c("blue","green","white","yellow","red"))
lc<-colmm(length(sc)-1)
lc[length(lc)/2]<-"white"
lc[length(lc)/2+1]<-"white"
png("Bt.png")
par(cex=1.3)
image(1:nv,1:nv,Bt,col=lc,breaks=sc)
box(lwd=2)
dev.off()

# Inverse
Btm1<-solve(Bt)

# Steric observation operator
myH<-array(0,dim=c(nt,nv))
for (i in 1:nt){ 
	k1<-(i-1)*nl+1
	k2<-(i  )*nl
        myH[i,k1:k2]<-c.lm[1:nl]
}

# CERES observation operator
myC<-array(0,dim=c(nt,nv))
for (i in 1:nt){ 
	k1<-(i-1)*nl+1
	k2<-(i  )*nl
        myC[i,k1:k2]<-1
}

# OHC observation operator (0-700)
myO7<-array(0,dim=c(nt,nv))
for (i in 1:nt){
        k1<-(i-1)*nl+1
        k2<-(i-1)*nl+2
        myO7[i,k1:k2]<-1
}

# OHC observation operator (700-2000)
myO2<-array(0,dim=c(nt,nv))
for (i in 1:nt){
        k1<-(i-1)*nl+3
        k2<-(i-1)*nl+3
        myO2[i,k1:k2]<-1
}
myC
myO7
myO2

tr12<-function(xx){
  dd<-array(0,dim=c(nt,nl))
  k1<-0
  for (t in 1:nt){
       for (l in 1:nl) {
                k1<-k1+1
                dd[t,l]<-xx[k1]
        }
  }
  return(dd)
}

# OHC data
o_y<-2003:2018
tmp<-get_fld("data_export/GCOS_all_heat_content_2003-2018_ZJ_v22062020.nc","ohc_0-700m")
o_7<-anom(tmp)
e_7<-get_fld("data_export/GCOS_all_heat_content_2003-2018_ZJ_v22062020.nc","ohc_0-700m_uncertainty")
range(e_7)
e_7<-sqrt(e_7^2+2^2)
range(e_7)
tmp<-get_fld("data_export/GCOS_all_heat_content_2003-2018_ZJ_v22062020.nc","ohc_700-2000m")
o_2<-anom(tmp)
e_2<-get_fld("data_export/GCOS_all_heat_content_2003-2018_ZJ_v22062020.nc","ohc_700-2000m_uncertainty")
range(e_2)
e_2<-sqrt(e_2^2+0.3^2)
range(e_2)

tmp<-get_fld("data_export/GCOS_all_heat_content_2003-2018_ZJ_v22062020.nc","ohc_below_2000m")
ohcd_kvs<-anom(tmp)
tmp<-get_fld("data_export/GCOS_all_heat_content_2003-2018_ZJ_v22062020.nc","ohc_below_2000m_uncertainty")
ohcd_kvse<-tmp

dt_7<-anom(dtb[1,]+dtb[2,])
dt_2<-anom(dtb[3,])

sd(dt_7)
sd(o_7)
sd(dt_2)
sd(o_2)

d_7<-o_7-dt_7
d_2<-o_2-dt_2

dim(d_7)<-c(nt,1)
dim(d_2)<-c(nt,1)

e_7

R_7<-diag(e_7^2)
R_7
Rm_7<-solve(R_7)
R_2<-diag(e_2^2)
Rm_2<-solve(R_2)

# INTEGRATED OHC
dts<-array(0,dim=nt)
for (t in 1:nt){
	dts[t]<-sum( dtb[,t] )
}

# CERES Observations AND MISFITS
norm<-function(dt){ {dt-mean(dt)}/sd(dt) }
anom<-function(dt){ {dt-mean(dt)}}
cer<-cumsum(get_fld("data_export/CERES_TOA_2003-2018_ym.nc","gtoa_net_all_mon")*1.1)
mean(get_fld("data_export/CERES_TOA_2003-2018_ym.nc","gtoa_net_all_mon")*1.0)
cconv<-0.07169398
beta<-0.90
obs2<-beta*cer[1:nt]/cconv

for (t in 1:nt){
write(paste(t+2002,dtb[t],obs2[t]),file="oocc.dat",append=TRUE)
}

d2<-anom(obs2)-anom(dts)

source("~/R_library/v1/runmea.R")
rep<-sd(runmea( anom(dts),1,nolim=TRUE)-anom(dts),na.rm=TRUE)*cconv
rep

rep
cer_err<-sqrt( 0.43^2 +0.3^2 +rep^2)/cconv
#             instrum + reprB + repOHC
cer_err
cer_err*cconv

d2
cer0<-norm(cer)
dts0<-norm(dts)
pdf("ceres.pdf")
par(cex=1.2)
plot(cer0,lwd=4,type="l",ylim=c(range(dts,cer)))
lines(dts0,lwd=2,col="red")
box(lwd=2)
dev.off()

obs2
dts
d2
sqrt(sum(d2^2))

dim(dtb)<-c(nv,1)
dim(obs)<-c(nt,1)

hxb<- myH %*% dtb
hxb
obs
dim(hxb)
length(obs)
d<-obs-hxb

obs
hxb

dtb

dim(d)

dtb
myH

d
# Error for steric
obs_err<-steric_err
obs_err2<-cer_err
R<-diag(obs_err^2, nrow=nt,ncol=nt)
Rm1<-solve(R)
Rm1

R2<-diag(obs_err2^2, nrow=nt,ncol=nt)
Rm2<-solve(R2)

Jb<- t(dx) %*% ( Btm1 %*% dx )
Jb

d2b<-d2
dim(d2b)<-c(nt,1)

# Here below the functions called by the minimizer
# to calculate the cot function and its gradient

fcost4<-function(x){
  dx0<-x
  dim(dx0)<-c(nv,1)
  Jb<- t(dx0) %*% ( Btm1 %*% dx0 )

  hdx<- myH %*% dx0
  hdxd<-hdx - d

  hdx2<- myC %*% dx0
  hdxd2<-hdx2 - d2b

  hdx3<- myO7 %*% dx0
  hdxd3<-hdx3 - d_7

  hdx4<- myO2 %*% dx0
  hdxd4<-hdx4 - d_2

  Jo<-  t(hdxd) %*% ( Rm1 %*% hdxd )
  Jo2<-  t(hdxd2) %*% ( Rm2 %*% hdxd2 )
  Jo3<-  t(hdxd3) %*% ( Rm_7 %*% hdxd3 )
  Jo4<-  t(hdxd4) %*% ( Rm_2 %*% hdxd4 )
  mcat(paste("Jb =",0.5*c(Jb),"   Jo1 =",0.5*c(Jo),"   Jo2 =",0.5*c(Jo2),"   Jtot =",0.5*c(Jb)+0.5*c(Jo+Jo2+Jo3+Jo4)))
  0.5*c(Jb)+0.5*c(Jo)+0.5*c(Jo2)+0.5*c(Jo3)+0.5*c(Jo4)
}

gcost4<-function(x){
  dx0<-x
  dim(dx0)<-c(nv,1)
  gb<-( Btm1 %*% dx0 )

  hdx<- myH %*% dx0
  hdxd<-hdx - d
  go<-t(myH) %*% ( Rm1 %*% hdxd )

  hdx2<- myC %*% dx0
  hdxd2<-hdx2 - d2b
  go2<-t(myC) %*% ( Rm2 %*% hdxd2 )

  hdx3<- myO7 %*% dx0
  hdxd3<-hdx3 - d_7
  go3<-t(myO7) %*% ( Rm_7 %*% hdxd3 )

  hdx4<- myO2 %*% dx0
  hdxd4<-hdx4 - d_2
  go4<-t(myO2) %*% ( Rm_2 %*% hdxd4 )

  gt<-gb+go+go2+go3+go4
  c(gt)
}

### 

fcost4(rep(0,nv))
fcost4(rep(1,nv))
gcost4(rep(1,nv))

# Call the minimizer
res<-optim( rep(0,nv), fn=fcost4, gr=gcost4,method="L-BFGS",control=list(factr=1e1,maxit=1000 ))
res$counts

warnings()
d
sqrt(mean(d^2))
nv
dxa<-res$par
dxa
dim(dxa)<-c(nv,1)
hdx<- myH %*% dxa
hxb<- myH %*% dtb
hdxd<-hdx - d
hdxa<- myH %*% ( dxa + dtb )
sqrt(mean(c(hdxd)^2))

png("steric.png")
par(cex=1.2)
plot(2003:2018,anom(obs),lwd=3,type="l",ylab="Steric Sea Level (mm)",xlab="",
     ylim=range(anom(obs),anom(hxb),anom(hdxa)))
lines(2003:2018,anom(hxb),lwd=2,col="red")
lines(2003:2018,anom(hdxa),lwd=2,col="blue")
legend("topleft",c("BG","AN","OBS"),fill=c("red","blue","black"),bty="n",cex=1.5)
box(lwd=2)
dev.off()

xtmp0s<-anom(hdxa)

sqrt(mean(d2^2))
dxa
hdxs<- c(myC %*% dxa)
sqrt(mean(c(hdxs-d2)^2))

sqrt(mean(d_7^2))
hdtmp<-c(myO7 %*% dxa)
sqrt(mean(c(hdtmp-d_7)^2))

sqrt(mean(d_2^2))
hdtmp<-c(myO2 %*% dxa)
sqrt(mean(c(hdtmp-d_2)^2))

png("ceres.png")
par(cex=1.2)
plot(2003:2018,anom(obs2),lwd=3,type="l",ylab="CERES cumulated EEI (ZJ)",xlab="",
     ylim=range(anom(obs2),anom(dts),anom(dts+hdxs)))
lines(2003:2018,anom(dts),lwd=2,col="red")
lines(2003:2018,anom(dts+hdxs),lwd=2,col="blue")
legend("topleft",c("BG","AN","OBS"),fill=c("red","blue","black"),bty="n",cex=1.5)
box(lwd=2)
dev.off()
xtmp0c<-anom(dts+hdxs)

xa<-xb+dxa
xb
xa
xa2<-tr12(xa)
xb2<-tr12(xb)

png("ohc.png")
par(mfrow=c(2,1))
par(cex=1.2)
par(mar=c(5,5,2,1)+0.1)

m<-40
sc<-seq(-m,m,length.out=41)
colmm<-colorRampPalette(c("blue","green","white","yellow","red"))
lc<-colmm(length(sc)-1)
lc[length(lc)/2]<-"white"
lc[length(lc)/2+1]<-"white"

range(xa2,xb2)
dim(xb2)

image(1:nt,1:4,xb2,col=lc,breaks=sc,ylim=c(4.5,0.5),main="Background",xaxt="n",yaxt="n",ylab="",xlab="")
axis(1,at=1:nt,labels=seq(2003,2018))
axis(2,at=1:4,labels=levl,las=1)
box(lwd=2)

image(1:nt,1:4,xa2,col=lc,breaks=sc,ylim=c(4.5,0.5),main="Analysis",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,at=1:nt,labels=seq(2003,2018))
axis(2,at=1:4,labels=levl,las=1)
box(lwd=2)

xb2_s<-xb2
xa2_s<-xa2

dev.off()

png("ohc_aux.png")
par(cex=1.2)

m<-40
sc<-seq(-m,m,length.out=41)
colmm<-colorRampPalette(c("blue","green","white","yellow","red"))
lc<-colmm(length(sc)-1)
lc[length(lc)/2]<-"white"
lc[length(lc)/2+1]<-"white"

filled.contour(1:nt,1:4,xb2,col=lc,levels=sc,ylim=c(4.5,0.5),main="Background",xaxt="n",yaxt="n",ylab="",xlab="")

dev.off()


png("ohct.png")
par(mfrow=c(1,2))
par(cex=1.2,las=1)
par(mar=c(5,5,2,1)+0.1)
plot(0,type="n",ylim=c(4,1),xlim=range( 0, ohct(xa2),ohct(xb2) ),yaxt="n",ylab="",xlab="W m-2",main="Ocean Warming")
axis(2,at=1:4,labels=levl)
lines(c(0,0),c(0,5),lty=2)
lines(ohct(xb2),1:4,lwd=3,col="red")
lines(ohct(xa2),1:4,lwd=3,col="blue")
box(lwd=2)
plot(0,type="n",ylim=c(4,1),xlim=range( 0, cumsum(ohct(xa2)),cumsum(ohct(xb2)) ),yaxt="n",ylab="",xlab="W m-2",main="Cumulative\nOcean Warming")
axis(2,at=1:4,labels=levl)
lines(c(0,0),c(0,5),lty=2)
lines(cumsum(ohct(xb2)),1:4,lwd=3,col="red",lty=1)
lines(cumsum(ohct(xa2)),1:4,lwd=3,col="blue",lty=1)
box(lwd=2)
dev.off()

sum( ohct(xb2) )
ohct(xb2)
sum( ohct(xa2) )
ohct(xa2)


# analysis error


HBHR<- R + myH %*% ( Bt %*% t(myH) )
dim(HBHR)
HBHR2<- R2 + myC %*% ( Bt %*% t(myC) )
dim(HBHR2)
HBHR3<- R_7 + myO7 %*% ( Bt %*% t(myO7) )
dim(HBHR3)
HBHR4<- R_2 + myO2 %*% ( Bt %*% t(myO2) )
dim(HBHR4)

nv
nl
nt

dim(Bt)
dim( ( Bt %*% t(myH) ) )
dim( t(myH) ) 
dim( (myH) ) 


HH<-array(0,dim=c(nt*4,nv))
k<-1
HH[((k-1)*nt+1):(nt*k),]<-myH
k<-2
HH[((k-1)*nt+1):(nt*k),]<-myC
k<-3
HH[((k-1)*nt+1):(nt*k),]<-myO7
k<-4
HH[((k-1)*nt+1):(nt*k),]<-myO2

obst<-c("STERIC","CERES","OHC0-700","OHC700-2000")

RR<-diag( c(diag(R),diag(R2),diag(R_7),diag(R_2)) )

HBHR <- RR + HH %*% ( Bt %*% t(HH) )
HBHRm1<-solve ( HBHR )
K<- Bt %*% ( t(HH) %*% HBHRm1)
KH<- K %*% HH
n<-dim(KH)[1]
In<-diag(n)
A<-(In-KH) %*% Bt
sqrt(diag(Bt))
sqrt(diag(A))

KTHT<-t(K) %*% t(HH)
dim(KTHT)
dim(HH)
dim(K)
d<-diag(KTHT)
a<-sum(d[1:nt])
b<-sum(d[(nt+1):(2*nt)])
c<-sum(d[(2*nt+1):(3*nt)])
d<-sum(d[(3*nt+1):(4*nt)])

e<-(nt*4)-(a+b+c+d)
png("dfs.png",width=640,height=640)
par(cex=1.3,las=1,mar=c(5,7,4,2)+0.1)
bpy<-100*c(a,b,c,d,e)/(4*nt)
bpx<-barplot(100*c(a,b,c,d,e)/(4*nt),names.arg=c(obst,"Background"),horiz=TRUE,xlim=c(0,max(bpy)+10),col=c("lightblue", "mistyrose",
                     "lightcyan", "lavender",add.alpha("salmon",0.6)),xlab="DFS (Degrees of Freedom for Signal)")
bpx
for (i in 1:length(bpx)) text(bpy[i],bpx[i],sprintf("%.1f%%",bpy[i]),pos=4)
dev.off()

100*c(a,b,c,d,e)/(4*nt)

KTHT<-t(K[seq(4,nv,by=4),]) %*% t(HH[,seq(4,nv,by=4)])
dim(KTHT)
dim(HH)
dim(K)
d<-diag(KTHT)


a<-sum(d[1:nt])
b<-sum(d[(nt+1):(2*nt)])
c<-sum(d[(2*nt+1):(3*nt)])
d<-sum(d[(3*nt+1):(4*nt)])
e<-(nt*4)-(a+b+c+d)

png("dfsd.png",width=400,height=400)
par(cex=1.2,las=1,mar=c(5,7,4,2)+0.1)
barplot(100*c(a,b,c,d,e)/(4*nt),names.arg=c(obst,"Background"),horiz=TRUE,col=c("lightblue", "mistyrose",
                     "lightcyan", "lavender",add.alpha("salmon",0.6)),xlab="DFS (Degrees of Freedom for Signal)")
dev.off()

100*c(a,b,c,d,e)/(4*nt)

KTHT<-t(K[seq(1,nv,by=4),]) %*% t(HH[,seq(1,nv,by=4)])
dim(KTHT)
dim(HH)
dim(K)
d<-diag(KTHT)

a<-sum(d[1:nt])
b<-sum(d[(nt+1):(2*nt)])
c<-sum(d[(2*nt+1):(3*nt)])
d<-sum(d[(3*nt+1):(4*nt)])
e<-(nt*4)-(a+b+c+d)


png("dfsu.png",width=400,height=400)
par(cex=1.2,las=1,mar=c(5,7,4,2)+0.1)
barplot(100*c(a,b,c,d,e)/(4*nt),names.arg=c(obst,"Background"),horiz=TRUE,col=c("lightblue", "mistyrose",
                     "lightcyan", "lavender",add.alpha("salmon",0.6)),xlab="DFS (Degrees of Freedom for Signal)")
dev.off()

100*c(a,b,c,d,e)/nt


for (k in 1:4){
H1<-HH[((k-1)*nt+1):(nt*k),]
dim(t(K))
dim(t(H1))
S<-t(K) %*% t(H1)
print(dim(S))
print(paste(k,obst[k],sum(diag(S))))
}

HAHT <- HH %*% ( A %*% t(HH) )
dim(HAHT)
sqrt(diag(HAHT))
sqrt(diag(RR))

png("obs_err.png",width = 640, height = 640)
par(mfrow=c(2,2))
par(cex=1.2)
par(las=1)
par(mar=c(3,4,2,1)+0.1)
k<-1 ; ooe<-sqrt(diag(RR))[((k-1)*nt+1):(nt*k)] ; aae<-sqrt(diag(HAHT))[((k-1)*nt+1):(nt*k)]
plot(2003:2018,ooe,type="l",ylim=range(ooe,aae),lwd=2,main="Steric",ylab="Error (mm)",xlab="")
lines(2003:2018,aae,col="blue",lwd=2)
legend("right",c("OBS","AN"),fill=c("black","blue"),bty="n")
box(lwd=2)

xtmpes<-aae
system("rm -f steric_an.out")
for (y in 2003:2018){ write(sprintf("%d %f %f",y,xtmp0s[y-2002],xtmpes[y-2002]),file="steric_an.out",append=TRUE) }

a<-100*(mean(ooe)-mean(aae))/mean(ooe)
text(2005,2.5,sprintf("%.1f%%",a),pos=1,cex=1.2)

k<-2 ; ooe<-sqrt(diag(RR))[((k-1)*nt+1):(nt*k)] ; aae<-sqrt(diag(HAHT))[((k-1)*nt+1):(nt*k)]
plot(2003:2018,ooe,type="l",ylim=range(ooe,aae),lwd=2,main="CERES",ylab="Error (ZJ)",xlab="")
lines(2003:2018,aae,col="blue",lwd=2)
legend("right",c("OBS","AN"),fill=c("black","blue"),bty="n")
box(lwd=2)

xtmpec<-aae
system("rm -f ceres_an.out")
for (y in 2003:2018){ write(sprintf("%d %f %f",y,xtmp0c[y-2002],xtmpec[y-2002]),file="ceres_an.out",append=TRUE) }


a<-100*(mean(ooe)-mean(aae))/mean(ooe)
text(2005,7.65,sprintf("%.1f%%",a),pos=1,cex=1.2)

k<-3 ; ooe<-sqrt(diag(RR))[((k-1)*nt+1):(nt*k)] ; aae<-sqrt(diag(HAHT))[((k-1)*nt+1):(nt*k)]
plot(2003:2018,ooe,type="l",ylim=range(ooe,aae),lwd=2,main="OHC 0-700m",ylab="Error (ZJ)",xlab="")
lines(2003:2018,aae,col="blue",lwd=2)
legend("topright",c("OBS","AN"),fill=c("black","blue"),bty="n")
box(lwd=2)

xtmpeh1<-aae

a<-100*(mean(ooe)-mean(aae))/mean(ooe)
text(2005,min(ooe,aae),sprintf("%.1f%%",a),pos=3,cex=1.2)

k<-4 ; ooe<-sqrt(diag(RR))[((k-1)*nt+1):(nt*k)] ; aae<-sqrt(diag(HAHT))[((k-1)*nt+1):(nt*k)]
plot(2003:2018,ooe,type="l",ylim=range(ooe,aae),lwd=2,main="OHC 700-2000m",ylab="Error (ZJ)",xlab="")
lines(2003:2018,aae,col="blue",lwd=2)
legend("topright",c("OBS","AN"),fill=c("black","blue"),bty="n")
box(lwd=2)

xtmpeh2<-aae

a<-100*(mean(ooe)-mean(aae))/mean(ooe)
text(2005,min(ooe,aae),sprintf("%.1f%%",a),pos=3,cex=1.2)


png("obs_ohc0-700.png")
par(cex=1.2)
bb<-dt_7 ; aa<-dt_7+c(myO7 %*% dxa) ; oo<-o_7
plot(0,type="n",xlim=c(2003,2018),ylim=range( 0, bb,aa,oo ),xlab="",ylab="Ocean heat content(ZJ)")
#lines(c(2003,2016),c(0,0),lty=2)
aa
xtmpeh1
do_polisd(2003:2018,aa,xtmpeh1,col=add.alpha("cyan",0.5))
write.table(cbind(aa,xtmpeh1),file="ohc07.dat",row.names =FALSE, col.names =FALSE)
lines(2003:2018,oo,lwd=3,col="black")
lines(2003:2018,bb,lwd=3,col="red")
lines(2003:2018,aa,lwd=3,col="blue")
legend("topleft",c("BG","AN","OBS"),fill=c("red","blue","black"),bty="n",cex=1.5)
box(lwd=2)
bb<-dt_2 ; aa<-dt_2+c(myO2 %*% dxa) ; oo<-o_2
png("obs_ohc700-2000.png")
par(cex=1.2)
plot(0,type="n",xlim=c(2003,2018),ylim=range( 0, bb,aa,oo ),xlab="",ylab="Ocean heat content(ZJ)")
#lines(c(2003,2016),c(0,0),lty=2)
do_polisd(2003:2018,aa,xtmpeh2,col=add.alpha("cyan",0.5))
write.table(cbind(aa,xtmpeh2),file="ohc7-2.dat",row.names =FALSE, col.names =FALSE)
lines(2003:2018,oo,lwd=3,col="black")
lines(2003:2018,bb,lwd=3,col="red")
lines(2003:2018,aa,lwd=3,col="blue")
legend("topleft",c("BG","AN","OBS"),fill=c("red","blue","black"),bty="n",cex=1.5)
box(lwd=2)
dev.off()

png("ceres.png")
par(cex=1.2)
plot(2003:2018,anom(obs2),lwd=3,type="n",ylab="CERES cumulated EEI (ZJ)",xlab="",
     ylim=range(anom(obs2),anom(dts),anom(dts+hdxs)))
do_polisd(2003:2018,anom(dts+hdxs),xtmpec,col=add.alpha("cyan",0.5))
write.table(cbind(anom(dts+hdxs),xtmpec),file="aceres.dat",row.names =FALSE, col.names =FALSE)
lines(2003:2018,anom(obs2),lwd=3)
lines(2003:2018,anom(dts),lwd=2,col="red")
lines(2003:2018,anom(dts+hdxs),lwd=2,col="blue")
legend("topleft",c("BG","AN","OBS"),fill=c("red","blue","black"),bty="n",cex=1.5)
box(lwd=2)
dev.off()

png("steric.png")
par(cex=1.2)
plot(2003:2018,anom(obs),lwd=3,type="n",ylab="Steric Sea Level (mm)",xlab="",
     ylim=range(anom(obs),anom(hxb),anom(hdxa)))
do_polisd(2003:2018,anom(hdxa),xtmpes,col=add.alpha("cyan",0.5))
write.table(cbind(anom(hdxa),xtmpes),file="asteric.dat",row.names =FALSE, col.names =FALSE)
lines(2003:2018,anom(obs),lwd=3)
lines(2003:2018,anom(hxb),lwd=2,col="red")
lines(2003:2018,anom(hdxa),lwd=2,col="blue")
legend("topleft",c("BG","AN","OBS"),fill=c("red","blue","black"),bty="n",cex=1.5)
box(lwd=2)
dev.off()

png("ohc_err.png")
par(mfrow=c(2,1))
par(cex=1.2)
par(mar=c(5,5,2,1)+0.1)

m<-150
sc<-seq(0,m,length.out=150)
colmm<-colorRampPalette(c("white","blue","green","yellow","red"))
lc<-colmm(length(sc)-1)
#lc[length(lc)/2]<-"white"
#lc[length(lc)/2+1]<-"white"
xb2<-array(0,dim=c(nt,nl))
xa2<-array(0,dim=c(nt,nl))
k1<-0
for (t in 1:nt){
         for (l in 1:nl) {
                k1<-k1+1
                xb2[t,l]<-sqrt(diag(Bt))[k1]/3
                xa2[t,l]<-sqrt(diag(A))[k1]
    }
}
range(xa2)
range(xb2)
dim(xb2)
image(1:nt,1:4,xb2,col=lc,breaks=sc,ylim=c(4.5,0.5),main="Background error",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,at=1:nt,labels=seq(2003,2018))
axis(2,at=1:4,labels=levl,las=1)
box(lwd=2)
image(1:nt,1:4,xa2,col=lc,breaks=sc,ylim=c(4.5,0.5),main="Analysis error",xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,at=1:nt,labels=seq(2003,2018))
axis(2,at=1:4,labels=levl,las=1)
box(lwd=2)
dev.off()

png("dohc.png")
par(cex=1.2)
db<-xb2_s[,4]
da<-xa2_s[,4]
eb<-xb2[,4]
ea<-xa2[,4]
plot(2003:2018,da,type="n",ylim=range(db-eb/4,db+eb/4,da-ea,da+ea),lwd=2,main="Deep ocean (2000m-bottom) OHC",xlab="",ylab="Ocean heat content (ZJ)",xaxs="i")
do_polisd(2003:2018,db,eb/5,col=add.alpha("orange",0.5))
do_polisd(2003:2018,da,ea,col=add.alpha("blue",0.2))
mean(ea)
detrend<-function(dt){
        a<-lm(y~x,data=data.frame(y=dt,x=1:length(dt)))$coefficients
        yf<-a[1]+a[2]*(1:length(dt))
        return(dt-yf)
}
sd(detrend(da))
write.table(cbind(anom(da),ea),file="depoh.dat",row.names =FALSE, col.names =FALSE)
lines(2003:2018,db,col="red",lwd=3)
lines(2003:2018,da,col="blue",lwd=3)
system("rm -f dohc.dat")
for (t in 1:length(da)){
  write(paste(t+2002.5,db[t],eb[t],da[t],ea[t]),file="dohc.dat",append=TRUE)
}
box(lwd=2)
lines(o_y,anom(ohcd_kvs),lwd=3,lty=1)
lines(o_y,anom(ohcd_kvs)-ohcd_kvse*1.1,lwd=1,lty=2,col="black")
lines(o_y,anom(ohcd_kvs)+ohcd_kvse*1.1,lwd=1,lty=2,col="black")
ohcd_kvs
ohcd_kvse
dev.off()

system("rm -f ohc_all.dat")
for (t in 1:length(da)){
  write(paste(t+2002.5,xb2_s[t,1],xb2[t,1],xa2_s[t,1],xa2[t,1],
	      xb2_s[t,2],xb2[t,2],xa2_s[t,2],xa2[t,2],
	      xb2_s[t,3],xb2[t,3],xa2_s[t,3],xa2[t,3],
	      xb2_s[t,4],xb2[t,4],xa2_s[t,4],xa2[t,4]),file="ohc_all.dat",append=TRUE)
}



q()
