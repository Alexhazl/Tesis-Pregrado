rm(list=ls())
library (plyr)
library(mice)
library(seewave)
library(corrplot)
library(gstat)
library(raster)
library(rgdal)
library(lattice)
library(RColorBrewer)
source("Funciones-Tesis-Pregrado.R")

options(digits=4)
options(warn=-1)

###################################
#### Cargando conjunto de datos ###
###################################

setwd("Regiones")

nms<-list.files()
nms2<-gsub(".txt","",nms) # Nombre de df's

for(i in nms){  # Carga datos y define formato de fechas
  assign(gsub(".txt","",i),read.table(i,header = T))
  assign(gsub(".txt","",i),chgmonths(get(gsub(".txt","",i))))
}


k<-c()
for(i in nms2){
  #seqdt<-seq(min(get(i)$Date),max(get(i)$Date),by=1)
  seqdt<-seq(as.Date("2013-01-01"),as.Date("2015-12-31"),by=1)
  seqdt<-data.frame(seqdt)
  colnames(seqdt)<-"Date"
  k[i]<-list(merge(seqdt,get(i), by="Date",all.x= TRUE))
  rm(seqdt)
}

for(i in 1:length(nms2)){ # Agrega fechas faltantes como NA's
  assign(nms2[i],ldply (k[i], data.frame))
}
rm(k)

###########################################
## Acotando las fechas para 2013 a 2015  ##
###########################################

for(i in nms2){
  assign(i,get(i)[get(i)$Date>=as.Date("2013-01-01") & get(i)$Date<=as.Date("2015-12-31") ,])  
  assign(i,cbind(get(i),ecef2lla(get(i)$Xm,get(i)$Ym,get(i)$Zm)))   # Transforma coordenadas
  assign(i,data.frame(as.Date(get(i)$Date),get(i)$Lon,get(i)$Lat))  # Acota columnas
  assign(i, `colnames<-`(get(i), c("Date","Lon","Lat")))            # Cambia nombres
}
 
for(i in nms2){ # Estado de los conjuntos de datos
  cat(i,"\n")
  print(summary(get(i)))
  cat("\n")
}

# Realiza la estimación de valores faltantes

# EMAT Lon
lst<-c("EMAT","BN17","CMBA","CMPN","CNBA","JUNT","LSCH","LVI1","SLMC")
EMAT$Lon<-PMiss(lst,1)

# EMAT LAT
lst<-c("EMAT","BN17","CMBA","CMPN","CNBA","JUNT","LSCH","LVI1","SLMC")
EMAT$Lat<-PMiss(lst,2)

# SILL Lon
lst<-c("SILL","CMBA","CNBA","CRZL","JUNT","LNDS","PEDR","SLMC","TOLO","ZAPA")
SILL$Lon<-PMiss(lst,1)

# SILL LAT
lst<-c("SILL","CMBA","CNBA","JUNT","LNDS","LSCH","PEDR","SANT","SLMC","TOLO","ZAPA")
SILL$Lat<-PMiss(lst,2)

#excep<-c("COPO","CTPC","SILL","POR2","OVLL","LLCH","CER1")
excep<-c("COPO","CTPC","POR2","LLCH","SANT")
aux<-data.frame(a=0)
for(i in nms2){
  if(!(i%in%excep)){
    aux<-cbind(aux,get(i)$Lon)
  }
}
aux<-aux[,-c(1)]

colnames(aux)<-nms2[!nms2%in%excep]

nms3<-nms2[!nms2%in%excep]

P1<-CreaDf(nms3,"2015-09-01","2015-09-15")
P2<-CreaDf(nms3,"2015-09-15","2015-09-18")
P3<-CreaDf(nms3,"2015-09-18","2015-09-30")

Descrip(P1[,c(4:6)])
Descrip(P2[,c(4:6)])
Descrip(P3[,c(4:6)])

setwd("C:/Users/alex9/Desktop/Carpetas/Carpetas/Codigos/Shapes")
map1<-readOGR("VLP.shp")
map2<-readOGR("CQM.shp")
map3<-readOGR("ATC.shp")

listmp<-list(map1,map2,map3)

mpreg <- bind(listmp[[1]],listmp[[2]],listmp[[3]])

shpmp <-as(mpreg,"SpatialPolygonsDataFrame")

r <- raster(ncol=180, nrow=180)
extent(r) <- extent(mpreg)
gridshp <- rasterize(mpreg,r)
gridshp <- as(gridshp, "SpatialPointsDataFrame")
gridshp<-as.data.frame(gridshp)
colnames(gridshp)<-c("Reg","Lon","Lat")
coordinates(gridshp) <- ~  Lon + Lat
gridshp <- as(gridshp, "SpatialPointsDataFrame")


coordinates(P1) <- ~ Lon + Lat
coordinates(P2) <- ~ Lon + Lat
coordinates(P3) <- ~ Lon + Lat

# Cross-Validacion
CrsVal(P1)  
CrsVal(P2,ct=20,CL="U")  
CrsVal(P3,CL="U")  


# Variogramas-direccionales
# Para cada caso

PltVarDir(P1,"Desplazamiento","Mat")
PltVarDir(P1,"Deslon","Mat")
PltVarDir(P1,"Deslat","Mat")


PltVarDir(P2,"Desplazamiento","Mat",ct=20,CL="U")
PltVarDir(P2,"Deslon","Mat",ct=20,CL="U")
PltVarDir(P2,"Deslat","Mat",ct=20,CL="U")

PltVarDir(P3,"Desplazamiento","Mat",CL="U")
PltVarDir(P3,"Deslon","Mat",CL="U")
PltVarDir(P3,"Deslat","Mat",CL="U")

# Variogramas
# Para cada caso

Pltvrgm(P1,"Mat","Desplazamiento")
Pltvrgm(P1,"Mat","Deslon")
Pltvrgm(P1,"Mat","Deslat")

Pltvrgm(P2,"Mat","Desplazamiento",ct=20,CL="U")
Pltvrgm(P2,"Mat","Deslon",ct=20,CL="U")
Pltvrgm(P2,"Mat","Deslat",ct=20,CL="U")

Pltvrgm(P3,"Mat","Desplazamiento",CL="U")
Pltvrgm(P3,"Mat","Deslon",CL="U")
Pltvrgm(P3,"Mat","Deslat",CL="U")

# Diseños para el grafico

myCols <- adjustcolor(colorRampPalette(brewer.pal(n=9, name = "YlOrRd"))(100), .85)##	Colores
myCols2<-colorRampPalette(c('yellow','orange', 'red'))(100)## Colores
cord<-as.matrix(data.frame(as(P1,"SpatialPoints")))
cord[,2]<-cord[,2]-0.1
sl1 <- list('sp.points', as(P1,"SpatialPoints"), pch=19, cex=.5, col='blue')##Puntos
sl2 <- list('sp.text', cord, nms3,cex=0.3, col='blue')##Nombres
pl1<-list(sl1,sl2)


# Calculo de Kriging para cada caso
P11<-PlotKrg(P1,gridshp,"Desplazamiento")
P12<-PlotKrg(P1,gridshp,"Deslon")
P13<-PlotKrg(P1,gridshp,"Deslat")

P21<-PlotKrg(P2,gridshp,"Desplazamiento",ct=20,CL="U")
P22<-PlotKrg(P2,gridshp,"Deslon",ct=20,CL="U")
P23<-PlotKrg(P2,gridshp,"Deslat",ct=20,CL="U")

P31<-PlotKrg(P3,gridshp,"Desplazamiento",CL="O")
P32<-PlotKrg(P3,gridshp,"Deslon",CL="O")
P33<-PlotKrg(P3,gridshp,"Deslat",CL="O")

# Asignacion de Kriging a cada variable

Krgn11 <- as(P11["var1.pred"], "SpatialPixelsDataFrame")#Kriging Ordinario
VKgn11 <- as(P11["var1.var"], "SpatialPixelsDataFrame")
Krgn12 <- as(P12["var1.pred"], "SpatialPixelsDataFrame")
VKgn12 <- as(P12["var1.var"], "SpatialPixelsDataFrame")
Krgn13 <- as(P13["var1.pred"], "SpatialPixelsDataFrame")
VKgn13 <- as(P13["var1.var"], "SpatialPixelsDataFrame")

Krgn21 <- as(P21["var1.pred"], "SpatialPixelsDataFrame")
VKgn21 <- as(P21["var1.var"], "SpatialPixelsDataFrame")
Krgn22 <- as(P22["var1.pred"], "SpatialPixelsDataFrame")
VKgn22 <- as(P22["var1.var"], "SpatialPixelsDataFrame")
Krgn23 <- as(P23["var1.pred"], "SpatialPixelsDataFrame")
VKgn23 <- as(P23["var1.var"], "SpatialPixelsDataFrame")

Krgn31 <- as(P31["var1.pred"], "SpatialPixelsDataFrame")
VKgn31 <- as(P31["var1.var"], "SpatialPixelsDataFrame")
Krgn32 <- as(P32["var1.pred"], "SpatialPixelsDataFrame")
VKgn32 <- as(P32["var1.var"], "SpatialPixelsDataFrame")
Krgn33 <- as(P33["var1.pred"], "SpatialPixelsDataFrame")
VKgn33 <- as(P33["var1.var"], "SpatialPixelsDataFrame")

# Trazado de Kriging con su respectiva varianza
# Para cada caso

spplot(Krgn11,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(VKgn11,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(Krgn12,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(VKgn12,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(Krgn13,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(VKgn13,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))

spplot(Krgn21,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(VKgn21,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(Krgn22,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(VKgn22,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(Krgn23,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(VKgn23,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))

spplot(Krgn31,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(VKgn31,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(Krgn32,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(VKgn32,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(Krgn33,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))
spplot(VKgn33,par.settings=list(fontsize=list(text=20)),col.regions=myCols ,sp.layout=list(pl1))

#
