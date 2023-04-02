library (plyr)
library(mice)
library(seewave)
library(corrplot)
library(gstat)
library(raster)
library(rgdal)
library(lattice)
library(RColorBrewer)


ExtFech<-function(nmrs,Fecha){
  #Extrae datos de df's para alguna fecha especifica
  
  df<-data.frame()
  for(i in nmrs){
    df<-rbind(df,cbind(i,get(i)[get(i)$Date==as.Date(Fecha),][,c(2,3)]))
  }  
  colnames(df)<-c("Estacion","Lon","Lat")
  rownames(df)<-NULL
  return(df)
}

CreaDf<-function(lista,Fech1,Fech2){
  # Crea el df de las variables regionalizadas consideradas para este estudio
  
  df1<-ExtFech(lista,Fech1)
  df2<-ExtFech(lista,Fech2)
  Coord<-(df1[,c(2,3)]+df2[,c(2,3)])/2
  DF<-data.frame(lista,sqrt((df1$Lon-df2$Lon)^2+(df1$Lat-df2$Lat)^2)*((2*pi*6371000)/360),abs(df1$Lon-df2$Lon)*((2*pi*6371000)/360),abs(df1$Lat-df2$Lat)*((2*pi*6371000)/360))
  DF<-cbind(DF,Coord)
  colnames(DF)<-c("Estacion","Desplazamiento","Deslon","Deslat","Lon","Lat")
  DF<-DF[,c(1,5,6,2,3,4)]
  return(DF)
}

PMiss<-function(lista,op){
  # Estimación de valores faltantes para ciertas series de tiempo
  # op=1 => Lon y op=2 => Lat
  
  if(op==1){
    aux<-data.frame(nrow=1095)
    for(i in lista){
      aux<-cbind(aux,get(i)$Lon)
    }
    aux<-aux[,-c(1)]
    aux<-scale(aux, center = TRUE, scale = TRUE)
    colnames(aux)<-lista
    tempData <- mice(aux,m=1,maxit=100,meth='lasso.norm',seed=12,remove.collinear=F,printFlag = F)
    completedData <- complete(tempData,1)
    a<-completedData[,c(1)]
    #a <- rmnoise(a,f=10E100)
    a<-(a*sd(get(lista[1])[,c(2)],na.rm = T))+mean(get(lista[1])[,c(2)],na.rm = T)
    return(a)
  }else{
    aux<-data.frame(nrow=1095)
    for(i in lista){
      aux<-cbind(aux,get(i)$Lat)
    }
    aux<-aux[,-c(1)]
    aux<-scale(aux, center = TRUE, scale = TRUE)
    colnames(aux)<-lista
    tempData <- mice(aux,m=1,maxit=100,meth='lasso.norm',seed=12,remove.collinear=F,printFlag = F)
    completedData <- complete(tempData,1)
    a<-completedData[,c(1)]
    #a <- rmnoise(a,f=10E100)
    a<-(a*sd(get(lista[1])[,c(3)],na.rm = T))+mean(get(lista[1])[,c(3)],na.rm = T)
    return(a)
  }
}

chgmonths<-function(df){ 
  # Estandariza formato de fechas}
  
  vec<-df[,2]
  vec<-gsub("JAN" ,"-01-",vec)
  vec<-gsub("FEB" ,"-02-",vec)
  vec<-gsub("MAR" ,"-03-",vec)
  vec<-gsub("APR" ,"-04-",vec)
  vec<-gsub("MAY" ,"-05-",vec)
  vec<-gsub("JUN" ,"-06-",vec)
  vec<-gsub("JUL" ,"-07-",vec)
  vec<-gsub("AUG" ,"-08-",vec)
  vec<-gsub("SEP" ,"-09-",vec)
  vec<-gsub("OCT" ,"-10-",vec)
  vec<-gsub("NOV" ,"-11-",vec)
  vec<-gsub("DEC" ,"-12-",vec)
  vec<-as.Date(vec,"%y-%m-%d")
  df[,2]<-vec
  return(df)
}

ecef2lla<-function(x,y,z){ 
  # Transformacion de coordenadas dado que estan en UTM
  
  a = 6378137
  e = 8.1819190842622e-2
  
  b   = sqrt(a^2*(1-e^2))
  ep  = sqrt((a^2-b^2)/b^2)
  p   = sqrt(x^2+y^2)
  th  = atan2(a*z,b*p)
  lon = atan2(y,x)
  
  lat = atan2((z+ep^2*b*sin(th)^3),(p-e^2*a*cos(th)^3))
  N   = a/sqrt(1-e^2*sin(lat)^2)
  alt = p/cos(lat)-N
  
  
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  Lon<-rad2deg(lon)
  Lat<-rad2deg(lat)
  return(data.frame(Lon,Lat))
}

CandCorr<-function(M,nm){ 
  # Identifica las mejores correlaciones de una matriz
  NM<-names(M[1,])
  n<-dim(M)[1]
  for(i in 1:n){
    cat("\n",NM[i],"\n",names(M[i,])[M[i,]>=nm][names(M[i,])[M[i,]>=nm]!=NM[i]],"\n")
  }
}


CrsVal<-function(df,ct=F,CL="O"){
  # Cross -Validación
  # ct ==> cutoff y CL ==> "O" = Ordinary "U" = Universal
  
  nms<-names(df)[2:4]
  if(CL=="O"){
    cat("\n Ordinary Kriging \n")
    if(ct==FALSE){
      for(i in 1:3){
        vrgm <- variogram(get(nms[i])~1, df) # calculates sample variogram values 
        desm<-fit.variogram(vrgm , model=vgm("Mat")) # fit model
        desn<-fit.variogram(vrgm , model=vgm("Nug")) 
        dese<-fit.variogram(vrgm , model=vgm("Exp")) 
        dess<-fit.variogram(vrgm , model=vgm("Sph")) 
        desg<-fit.variogram(vrgm , model=vgm("Gau")) 
        
        x1 <- krige.cv(get(nms[i])~1,df,desm, nfold=nrow(df),verbose=F)#Validacion cruzada
        x2 <- krige.cv(get(nms[i])~1,df,desn, nfold=nrow(df),verbose=F)#
        x3 <- krige.cv(get(nms[i])~1,df,dese, nfold=nrow(df),verbose=F)#
        x4 <- krige.cv(get(nms[i])~1,df,dess, nfold=nrow(df),verbose=F)#
        x5 <- krige.cv(get(nms[i])~1,df,desg, nfold=nrow(df),verbose=F)#
        
        rmsep1di1<-sqrt(sum(x1$residual^2))/length(x1$residual)
        rmsep1di2<-sqrt(sum(x2$residual^2))/length(x2$residual)
        rmsep1di3<-sqrt(sum(x3$residual^2))/length(x3$residual)
        rmsep1di4<-sqrt(sum(x4$residual^2))/length(x4$residual)
        rmsep1di5<-sqrt(sum(x5$residual^2))/length(x5$residual)
        
        assign(paste0("v",i),c(rmsep1di1,rmsep1di2,rmsep1di3,rmsep1di4,rmsep1di5))
      }
    }else{
      for(i in 1:3){
        vrgm <- variogram(get(nms[i])~1, df,cutoff=ct) # calculates sample variogram values 
        desm<-fit.variogram(vrgm , model=vgm("Mat")) # fit model
        desn<-fit.variogram(vrgm , model=vgm("Nug")) 
        dese<-fit.variogram(vrgm , model=vgm("Exp")) 
        dess<-fit.variogram(vrgm , model=vgm("Sph")) 
        desg<-fit.variogram(vrgm , model=vgm("Gau")) 
        
        x1 <- krige.cv(get(nms[i])~1,df,desm, nfold=nrow(df),verbose=F)#Validacion cruzada
        x2 <- krige.cv(get(nms[i])~1,df,desn, nfold=nrow(df),verbose=F)#
        x3 <- krige.cv(get(nms[i])~1,df,dese, nfold=nrow(df),verbose=F)#
        x4 <- krige.cv(get(nms[i])~1,df,dess, nfold=nrow(df),verbose=F)#
        x5 <- krige.cv(get(nms[i])~1,df,desg, nfold=nrow(df),verbose=F)#
        
        rmsep1di1<-sqrt(sum(x1$residual^2))/length(x1$residual)
        rmsep1di2<-sqrt(sum(x2$residual^2))/length(x2$residual)
        rmsep1di3<-sqrt(sum(x3$residual^2))/length(x3$residual)
        rmsep1di4<-sqrt(sum(x4$residual^2))/length(x4$residual)
        rmsep1di5<-sqrt(sum(x5$residual^2))/length(x5$residual)
        
        assign(paste0("v",i),c(rmsep1di1,rmsep1di2,rmsep1di3,rmsep1di4,rmsep1di5))
      }
    }
    a<-data.frame(rbind(v1,v2,v3))
    colnames(a)<-c("Mat","Nug","Exp","Sph","Gau")
    rownames(a)<-c("Des","Deslon","Deslat")
    print(a)
  }else if(CL=="U"){
    cat("\n Universal Kriging \n")
    if(ct==FALSE){
      for(i in 1:3){
        vrgm <- variogram(get(nms[i])~Lon+Lat, df) # calculates sample variogram values 
        desm<-fit.variogram(vrgm , model=vgm("Mat")) # fit model
        desn<-fit.variogram(vrgm , model=vgm("Nug")) 
        dese<-fit.variogram(vrgm , model=vgm("Exp")) 
        dess<-fit.variogram(vrgm , model=vgm("Sph")) 
        desg<-fit.variogram(vrgm , model=vgm("Gau")) 
        
        x1 <- krige.cv(get(nms[i])~Lon+Lat,df,desm, nfold=nrow(df),verbose=F)#Validacion cruzada
        x2 <- krige.cv(get(nms[i])~Lon+Lat,df,desn, nfold=nrow(df),verbose=F)#
        x3 <- krige.cv(get(nms[i])~Lon+Lat,df,dese, nfold=nrow(df),verbose=F)#
        x4 <- krige.cv(get(nms[i])~Lon+Lat,df,dess, nfold=nrow(df),verbose=F)#
        x5 <- krige.cv(get(nms[i])~Lon+Lat,df,desg, nfold=nrow(df),verbose=F)#
        
        rmsep1di1<-sqrt(sum(x1$residual^2))/length(x1$residual)
        rmsep1di2<-sqrt(sum(x2$residual^2))/length(x2$residual)
        rmsep1di3<-sqrt(sum(x3$residual^2))/length(x3$residual)
        rmsep1di4<-sqrt(sum(x4$residual^2))/length(x4$residual)
        rmsep1di5<-sqrt(sum(x5$residual^2))/length(x5$residual)
        
        assign(paste0("v",i),c(rmsep1di1,rmsep1di2,rmsep1di3,rmsep1di4,rmsep1di5))
      }
    }else{
      for(i in 1:3){
        vrgm <- variogram(get(nms[i])~Lon+Lat, df,cutoff=ct) # calculates sample variogram values 
        desm<-fit.variogram(vrgm , model=vgm("Mat")) # fit model
        desn<-fit.variogram(vrgm , model=vgm("Nug")) 
        dese<-fit.variogram(vrgm , model=vgm("Exp")) 
        dess<-fit.variogram(vrgm , model=vgm("Sph")) 
        desg<-fit.variogram(vrgm , model=vgm("Gau")) 
        
        x1 <- krige.cv(get(nms[i])~Lon+Lat,df,desm, nfold=nrow(df),verbose=F)#Validacion cruzada
        x2 <- krige.cv(get(nms[i])~Lon+Lat,df,desn, nfold=nrow(df),verbose=F)#
        x3 <- krige.cv(get(nms[i])~Lon+Lat,df,dese, nfold=nrow(df),verbose=F)#
        x4 <- krige.cv(get(nms[i])~Lon+Lat,df,dess, nfold=nrow(df),verbose=F)#
        x5 <- krige.cv(get(nms[i])~Lon+Lat,df,desg, nfold=nrow(df),verbose=F)#
        
        rmsep1di1<-sqrt(sum(x1$residual^2))/length(x1$residual)
        rmsep1di2<-sqrt(sum(x2$residual^2))/length(x2$residual)
        rmsep1di3<-sqrt(sum(x3$residual^2))/length(x3$residual)
        rmsep1di4<-sqrt(sum(x4$residual^2))/length(x4$residual)
        rmsep1di5<-sqrt(sum(x5$residual^2))/length(x5$residual)
        
        assign(paste0("v",i),c(rmsep1di1,rmsep1di2,rmsep1di3,rmsep1di4,rmsep1di5))
      }
    }
    a<-data.frame(rbind(v1,v2,v3))
    colnames(a)<-c("Mat","Nug","Exp","Sph","Gau")
    rownames(a)<-c("Des","Deslon","Deslat")
    print(a)
  }
}

PlotKrg<-function(df,shp,vari,CL="O",ct=F){
  # Traza Krigin en funcion de la variable regionalizada
  # ct ==> cutoff y CL ==> "O" = Ordinary "U" = Universal
  # vari ==> variable a trabajar
  
  if(CL=="O"){
    if(ct==FALSE){
      vrgm <- variogram(get(vari)~1, df)  
      des<-fit.variogram(vrgm , model=vgm("Mat"))
      lzn.kriged <- krige(get(vari)~ 1, df, shp, model=des)#Kriging
      return(lzn.kriged)  
    }else{
      vrgm <- variogram(get(vari)~1, df,cutoff=ct)  
      des<-fit.variogram(vrgm , model=vgm("Mat"))
      lzn.kriged<- krige(get(vari)~ 1, df, shp, model=des)#Kriging
      return(lzn.kriged)
    }
  }else if(CL=="U"){
    if(ct==FALSE){
      vrgm <- variogram(get(vari)~Lon+Lat, df)  
      des<-fit.variogram(vrgm , model=vgm("Mat"))
      lzn.kriged <- krige(get(vari)~ Lon+Lat, df, shp, model=des)#Kriging
      return(lzn.kriged)  
    }else{
      vrgm <- variogram(get(vari)~Lon+Lat, df,cutoff=ct)  
      des<-fit.variogram(vrgm , model=vgm("Mat"))
      lzn.kriged<- krige(get(vari)~ Lon+Lat, df, shp, model=des)#Kriging
      return(lzn.kriged)
    }
  }
}

Pltvrgm<-function(df,mod,vari,ct=F,CL="O"){
  # Traza variograma 
  # ct ==> cutoff y CL ==> "O" = Ordinary "U" = Universal
  # vari ==> variable a trabajar
  # mod ==> modelo de variograma ajustar
  
  if(CL=="O"){
    cat("\n Ordinary \n")
    if(ct==FALSE){
      vrgm <- variogram(get(vari)~1, df) # Variograma empirico 
      des<-fit.variogram(vrgm , model=vgm(mod))# ajuste del modelo
      plot(vrgm, des,ylab="Semivarianza",xlab="Distancia")
    }else{
      vrgm <- variogram(get(vari)~1, df,cutoff=ct) # Variograma empirico 
      des<-fit.variogram(vrgm , model=vgm(mod))# ajuste del modelo
      plot(vrgm, des,ylab="Semivarianza",xlab="Distancia")
    }
  }else if(CL=="U"){
    cat("\n Universal \n")
    if(ct==FALSE){
      vrgm <- variogram(get(vari)~Lon+Lat, df) # Variograma empirico 
      des<-fit.variogram(vrgm , model=vgm(mod))# ajuste del modelo
      plot(vrgm, des,ylab="Semivarianza",xlab="Distancia")
    }else{
      vrgm <- variogram(get(vari)~Lon+Lat, df,cutoff=ct) # Variograma empirico 
      des<-fit.variogram(vrgm , model=vgm(mod))# ajuste del modelo
      plot(vrgm, des,ylab="Semivarianza",xlab="Distancia")
    }
  }
}

PltVarDir<-function(df,vari,mod,ct=F,CL="O"){
  # Traza variograma direccional
  # ct ==> cutoff y CL ==> "O" = Ordinary "U" = Universal
  # vari ==> variable a trabajar
  # mod ==> modelo de variograma ajustar
  
  if(CL=="O"){
    cat("\n Ordinary \n")
    if(ct==FALSE){
      v.dir <- variogram(get(vari) ~ 1,df, alpha = (0:3) *45)
      des<-fit.variogram(v.dir , model=vgm(mod))
      plot(v.dir, des)
    }else{
      v.dir <- variogram(get(vari) ~ 1,df,cutoff=ct ,alpha = (0:3) *45)
      des<-fit.variogram(v.dir , model=vgm(mod))
      plot(v.dir, des)
    }
  }else if(CL=="U"){
    cat("\n Universal \n")
    if(ct==FALSE){
      v.dir <- variogram(get(vari) ~ Lon+Lat,df, alpha = (0:3) *45)
      des<-fit.variogram(v.dir , model=vgm(mod))
      plot(v.dir, des)
    }else{
      v.dir <- variogram(get(vari) ~ Lon+Lat,df,cutoff=ct ,alpha = (0:3) *45)
      des<-fit.variogram(v.dir , model=vgm(mod))
      plot(v.dir, des)
    }
    
  }
}

c.skew<-function(x) {
  m3=mean((x-mean(x))^3)
  skew=m3/(sd(x)^3)
  return(skew)
}

Descrip<-function(df){
  sumdes<-data.frame(t(sapply(df,FUN = function(x) c(summary(x),sd(x),sd(x)/mean(x),c.skew(x)))))
  colnames(sumdes)<-c("Min","Q1","Q2","Media","Q3","Max","sd","CV","ks")
  #sumdes<-round(sumdes,digits = 2)
  print(sumdes)
}