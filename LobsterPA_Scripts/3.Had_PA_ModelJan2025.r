
require(sdmTMB)
require(bio.lobster)
require(bio.utilities)
require(lubridate)
require(devtools)
require(dplyr)
require(ggplot2)
require(INLA)
options(stringAsFactors=F)
require(PBSmapping)
require(SpatialHub)
require(sf)
la()
fd=file.path(project.datadirectory('bio.lobster'),'analysis','ClimateModelling')
dir.create(fd,showWarnings=F)
setwd(fd)
survey = readRDS(file='BaseDataForClimateModelAllTemps_Jan2025.rds')
sf_use_s2(FALSE) #needed for cropping

# Project our survey data coordinates:

survey$lZ = log(survey$z)
survey = subset(survey,!is.na(lZ))

#what is the right offset for gear
#km2 for tows is estimated from sensors
#km2 for traps from Watson et al 2009 NJZ MFR 43 1 -- home radius of 17m, bait radius of 11m == 28m 'attraction zone'
# pi*(.014^2) # assuming traps are independent
# 17m is at 14C, so go up or down from there
survey$W = ceiling(yday(survey$DATE)/366*25)
mi= readRDS(file=file.path(project.datadirectory('bio.lobster'),'analysis','ClimateModelling','tempCatchability.rds'))

mi = mi[,c('Temperature','predicted')]
names(mi) = c('HadBTr','OFFSETcorr')
mi$OFFSETcorr = mi$OFFSETcorr * (1/mi$OFFSETcorr[which(mi$HadBT==14)]) # 17m is at 14C, so go up or down from there
# 17m is at 14C, so go up or down from there
survey$HadBTr = round(survey$HadBT*2)/2
survey = dplyr::full_join(survey,mi)
survey$OFFSETcorr[which(survey$HadBTr< -.7)] <- min(mi$OFFSETcorr)
survey$OFFSETcorr[which(survey$HadBTr> 17)] <- max(mi$OFFSETcorr)

i = which(survey$OFFSET_METRIC == 'Number of traps')
survey$OFFSET[i] = survey$OFFSET[i] * pi*(.014^2) * survey$OFFSETcorr[i]

survey$LO = log(survey$OFFSET)
survey = subset(survey,OFFSET>0.00001 & OFFSET< 0.12)
survey$BT = survey$HadBT
survey = subset(survey,!is.na(BT))
survey$pa = ifelse(survey$Berried>0,1,0)
i = which(survey$SOURCE=='ILTS_ITQ' & survey$WEIGHT_KG==0 & is.na(survey$Berried))
survey$pa[i] = 0
survey = subset(survey, !is.na(pa))
ns_coast =readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","CoastSF.rds"))
st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
crs_utm20 <- 32620
ns_coast <- suppressWarnings(suppressMessages(
  st_crop(ns_coast,
          c(xmin = -68, ymin = 41, xmax = -56.5, ymax = 47.5))))

ns_coast <- st_transform(ns_coast, crs_utm20)

st_crs(ns_coast) <- 4326 # 'WGS84'; necessary on some installs
crs_utm20 <- 32620

survey <- survey %>%   
  st_as_sf()
    
surv_utm_coords <- st_coordinates(survey)

survey$X1000 <- surv_utm_coords[,1] 
survey$Y1000 <- surv_utm_coords[,2] 

spde <- make_mesh(as_tibble(survey), xy_cols = c("X1000", "Y1000"),
                   n_knots=600,type = "cutoff_search")
plot(spde)

# Add on the barrier mesh component:
bspde <- add_barrier_mesh(
  spde, ns_coast, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)


#issues with fitting on a biweekly moving to quarters
survey$m = month(survey$DATE) 
survey$Q = ifelse(survey$m %in% c(10,11,12),1,ifelse(survey$m %in% c(1,2,3),2,ifelse(survey$m %in% c(4,5,6),3,4)))
survey$Time = survey$YEAR+survey$Q/4


fitpa = sdmTMB(pa~
               s(lZ,k=3)+s(BT,k=4)+Q,
             data=as_tibble(survey),
            offset = 'LO',
             time='YEAR', 
             mesh=bspde,
             family=binomial(link='logit'),
             spatial='on',
             spatiotemporal='ar1')

saveRDS(list(data=survey,grid=bspde,model=fitpa),file='sdmTMBBerriedpabyQFinal_Had_Jan2025.rds')

x=readRDS(file='sdmTMBBerriedpabyQFinal_Had_Jan2025.rds')

fitpa = x$model
bspde = x$grid
survey= x$data
x = readRDS('HadClimatologiesDiffsbyQ.rds')
x = x[[1]]

x = subset(x,Q==3, select=c(Q,z,CurClim0.5))
x$z = abs(x$z)
x = subset(x,z>0 )

x$lZ = log(x$z)
x$X1000 = st_coordinates(x)[,1]
x$Y1000 = st_coordinates(x)[,2]
x = subset(x,exp(lZ)<400)
x$BT = x$CurClim0.5


rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326
rL = st_transform(rL,32620) 
st_geometry(rL) <- st_geometry(st_as_sf(rL$geometry/1000)) 
st_crs(rL) <- 32620
ff = st_join(x,rL,join=st_within)
x = subset(ff,!is.na(LFA))

x = as_tibble(subset(x,select=c(Q,BT,X1000,Y1000,lZ,LFA)))
x$geometry=NULL
y = unique(survey$YEAR)
o = list()
for(i in 1:length(y)){
    w = x
    w$YEAR=y[i]
    o[[i]] = w
  }
  x = bind_rows(o)

g = predict(fitpa,newdata=x)

  g$pred = fitpa$family$linkinv(g$est)
gC = aggregate(pred~X1000+Y1000+LFA+BT,data=g,FUN=mean)
gsf = st_as_sf(gC,coords = c("X1000","Y1000"),crs=32620,remove=F)
#Maps


mm = c(0.001,max(gsf$pred))
ggplot(subset(gsf))+#,Q==3 &YEAR %in% 2016)) +
  geom_sf(aes(fill=pred,color=pred)) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()
savePlot('BerriedPAMeanSurface2000-2022_Jan2025.png') 

ggplot(survey)+
geom_sf(size=.3,colour='red')+
  geom_sf(data=rL,size=1,colour='black',fill=NA)

with(gsf,plot(BT,pred))
gsfYR = gsf

saveRDS(gsf,file='Had_CurrentClimatologyBerriedBinomialOutput_Jan2025.rds')

##################################################################
###projections

 glo = readRDS('HadClimatologiesDiffsbyQ.rds')
 glo = glo[[1]]
 glo$X1000 = st_coordinates(glo)[,1]
 glo$Y1000 = st_coordinates(glo)[,2]
 glo$lZ = log(glo$z)
glo = subset(glo,z<400 & Q==3)
 glo = subset(glo,!is.na(lZ))

rL = readRDS(file.path( project.datadirectory("bio.lobster"), "data","maps","LFAPolysSF.rds"))
rL = st_as_sf(rL)
st_crs(rL) <- 4326
rL = st_transform(rL,32620) 
st_geometry(rL) <- st_geometry(st_as_sf(rL$geometry/1000)) 
st_crs(rL) <- 32620


glo = st_join(glo,rL,join=st_within)

glo = subset(glo,!is.na(LFA))

 saveRDS(glo,'HadClimatologiesandDiffsbyQ_withinLFA.rds')

################
#Had30

y = unique(survey$YEAR)
cdata = subset(glo,select=c(X1000,Y1000,lZ,ThirClim0.5,Q))
st_geometry(cdata) <- NULL

cdy = as.data.frame(sapply(cdata,rep.int,length(y)))
cdy$YEAR = rep(y,each=dim(cdata)[1])
cdy$BT = cdy$ThirClim0.5
g = predict(fitpa,newdata=cdy)

g$pred = fitpa$family$linkinv(g$est)

#average over spatial domains for years 2000-2022
gHad30 = aggregate(pred~X1000+Y1000+Q+BT,data=g,FUN=mean)

gHad30 = st_as_sf(gHad30,coords = c("X1000","Y1000"),crs=32620,remove=F)


ggplot(subset(gHad30)) +
  geom_sf(aes(fill=pred,color=pred)) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()

 saveRDS(gHad30,file='Had_30sClimatologyBerriedBinomialOutput_Jan2025.rds')


################
#Had50

y = unique(survey$YEAR)
cdata = subset(glo,select=c(X1000,Y1000,lZ,FiftClim0.5,Q))
st_geometry(cdata) <- NULL

cdy = as.data.frame(sapply(cdata,rep.int,length(y)))
cdy$YEAR = rep(y,each=dim(cdata)[1])
cdy$BT = cdy$FiftClim0.5
cdy = subset(cdy,Q==3)
g = predict(fitpa,newdata=cdy)

g$pred = fitpa$family$linkinv(g$est)

gHad50 = aggregate(pred~X1000+Y1000+Q+BT,data=g,FUN=mean)

gHad50 = st_as_sf(gHad50,coords = c("X1000","Y1000"),crs=32620,remove=F)


ggplot(subset(gHad50)) +
  geom_sf(aes(fill=pred,color=pred)) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()

 saveRDS(gHad50,file='Had_50sClimatologyBerriedBinomialOutput_Jan2025.rds')


################
#Had90

y = unique(survey$YEAR)
cdata = subset(glo,select=c(X1000,Y1000,lZ,NineClim0.5,Q))
st_geometry(cdata) <- NULL

cdy = as.data.frame(sapply(cdata,rep.int,length(y)))
cdy$YEAR = rep(y,each=dim(cdata)[1])
cdy$BT = cdy$NineClim0.5
cdy = subset(cdy,Q==3)
g = predict(fitpa,newdata=cdy)

g$pred = fitpa$family$linkinv(g$est)

gHad90 = aggregate(pred~X1000+Y1000+Q+BT,data=g,FUN=mean)

gHad90 = st_as_sf(gHad90,coords = c("X1000","Y1000"),crs=32620,remove=F)


ggplot(subset(gHad90)) +
  geom_sf(aes(fill=pred,color=pred)) + 
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()

 saveRDS(gHad90,file='Had_90sClimatologyBerriedBinomialOutput_Jan2025.rds')

###merge
 gHad90 = readRDS(file='Had_90sClimatologyBerriedBinomialOutput_Jan2025.rds')
 gHad50 = readRDS(file='Had_50sClimatologyBerriedBinomialOutput_Jan2025.rds')
 gHad30 = readRDS(file='Had_30sClimatologyBerriedBinomialOutput_Jan2025.rds')
 gHad   = readRDS(file='Had_CurrentClimatologyBerriedBinomialOutput_Jan2025.rds')
 
 gHad = subset(gHad,select=c(pred,LFA))
 gHad30 = subset(gHad30,select=c(pred))
 gHad50 = subset(gHad50,select=c(pred))
 gHad90 = subset(gHad90,select=c(pred))
 
gHad = bio.utilities::rename.df(gHad,'pred','current')
 gHad30 = bio.utilities::rename.df(gHad30,'pred','thirpre')
 gHad50 = bio.utilities::rename.df(gHad50,'pred','fiftpre')
 gHad90 = bio.utilities::rename.df(gHad90,'pred','ninpre')

g1 = st_join(gHad,gHad30,join=st_intersects)
g1 = st_join(g1,gHad50,join=st_intersects)
g1 = st_join(g1,gHad90,join=st_intersects)

g1$tc = (g1$thirpre - g1$current) /g1$current *100
g1$fc = (g1$fiftpre - g1$current) /g1$current *100
g1$nc = (g1$ninpre - g1$current) /g1$current *100

g1$ft = (g1$fiftpre - g1$thirpre) /g1$thirpre *100
g1$nf = (g1$ninpre - g1$fiftpre) /g1$fiftpre *100


ggplot(subset(g1,current>0.01)) +
  geom_sf(aes(fill=tc,color=tc)) + 
  scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0) +
  scale_color_gradient2(low='blue',mid='white',high='red',midpoint=0) +
#  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf() +
  labs(title='2035 to Current')



ggplot(subset(g1,current>0.01)) +
  geom_sf(aes(fill=fc,color=fc)) + 
  scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0) +
  scale_color_gradient2(low='blue',mid='white',high='red',midpoint=0) +
  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()+
  labs(title='2055 to Current')



ggplot(subset(g1,current>0.01)) +
  geom_sf(aes(fill=nc,color=nc)) + 
  scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0) +
  scale_color_gradient2(low='blue',mid='white',high='red',midpoint=0) +
  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()+
  labs(title='2095 to Current')



ggplot(subset(g1,current>0.01)) +
  geom_sf(aes(fill=ft,color=ft)) + 
  scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0) +
  scale_color_gradient2(low='blue',mid='white',high='red',midpoint=0) +
#  geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()+
  labs(title='2055 to 2035')

ggplot(subset(g1,current>0.01)) +
  geom_sf(aes(fill=nf,color=nf)) + 
  scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0) +
  scale_color_gradient2(low='blue',mid='white',high='red',midpoint=0) +
 # geom_sf(data=rL,size=1,colour='black',fill=NA)+
  theme( axis.ticks.x = element_blank(),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
  ) +
  coord_sf()



#########################################
##keep everything else constant and what does temperature changes result in?

#temp relationship extracted from the model
  ff = plot_smooth(fitpa,select=2,return_data=T)
  ff$Pred = fitpa$family$linkinv(ff$est)
  ff = ff[,c('BT','Pred')]
  ff$BT = as.integer(round(ff$BT,1)*10)/10
  ts = data.frame(BT=seq(min(ff$BT),max(ff$BT),by=.1))
  ts$BT = as.integer(round(ts$BT,1)*10)/10
  ff1 = merge(ff,ts,all.y=T)
  require(zoo)
  ff1$Pred = na.approx(ff1$Pred)
  ff1$PredS = ff1$Pred/max(ff1$Pred)
  #rescale -1 to 1
  ff1$ResP = 2*(ff1$Pred-min(ff1$Pred)) / (max(ff1$Pred)-min(ff1$Pred))-1

ff1$diffs =  c(0,diff(ff1$Pred))

diffMat = matrix(0,nrow=length(ff1$BT),ncol=length(ff1$BT))
for(i in 1:nrow(ff1)){
  for(j in 1:nrow(ff1)){
      b = ff1$Pred[i:j]
      cc = c(0,diff(b))
      diffMat[i,j] = sum(cc)
  }

}
rownames(diffMat) = colnames(diffMat) = ff1$BT

##climate

#base surface
gs =readRDS(file='Had_CurrentClimatologyBerriedBinomialOutput_Jan2025.rds')
glo = readRDS('HadClimatologiesandDiffsbyQ_withinLFA.rds')

gh = subset(glo,Q==3,select=c(CurClim0.5,ThirClim0.5,FiftClim0.5,NineClim0.5))

gh = st_join(gs,gh,join=st_intersects)
gh$t2016 = round(gh$CurClim0.5,1)
gh$t2035 = round(gh$ThirClim0.5,1)
gh$t2055 = round(gh$FiftClim0.5,1)
gh$t2099 = round(gh$NineClim0.5,1)

hp=gh
hp = subset(hp,!is.na(t2016))
hp$diff1635 = NA
for(i in 1:nrow(hp)){
  hpp = hp[i,]
  st_geometry(hpp) <- NULL
  t1=hpp[1,'t2016']
  t2=hpp[1,'t2035']
  if(t2>max(as.numeric(rownames(diffMat)))) t2 = max(as.numeric(rownames(diffMat))) 
  hp$diff1635[i] = diffMat[which(rownames(diffMat)==t1),which(colnames(diffMat)==t2)]
}

hp$diff1655 = NA
for(i in 1:nrow(hp)){
  hpp = hp[i,]
  st_geometry(hpp) <- NULL
  t1=hpp[1,'t2016']
  t2=hpp[1,'t2055']
  if(t2>max(as.numeric(rownames(diffMat)))) t2 = max(as.numeric(rownames(diffMat))) 
  
  hp$diff1655[i] = diffMat[which(rownames(diffMat)==t1),which(colnames(diffMat)==t2)]
}


hp$diff1699 = NA
for(i in 1:nrow(hp)){
  hpp = hp[i,]
  st_geometry(hpp) <- NULL
  t1=hpp[1,'t2016']
  t2=hpp[1,'t2099']
  if(t2>max(as.numeric(rownames(diffMat)))) t2 = max(as.numeric(rownames(diffMat))) 
  
  hp$diff1699[i] = diffMat[which(rownames(diffMat)==t1),which(colnames(diffMat)==t2)]
}


hp$diff3555 = NA
for(i in 1:nrow(hp)){
  hpp = hp[i,]
  st_geometry(hpp) <- NULL
  t1=hpp[1,'t2035']
  t2=hpp[1,'t2055']
  if(t1>max(as.numeric(rownames(diffMat)))) t1 = max(as.numeric(rownames(diffMat))) 
 if(t2>max(as.numeric(rownames(diffMat)))) t2 = max(as.numeric(rownames(diffMat))) 

  hp$diff3555[i] = diffMat[which(rownames(diffMat)==t1),which(colnames(diffMat)==t2)]
}


hp$diff3599 = NA
for(i in 1:nrow(hp)){
  hpp = hp[i,]
  st_geometry(hpp) <- NULL
  t1=hpp[1,'t2035']
  t2=hpp[1,'t2099']
  if(t1>max(as.numeric(rownames(diffMat)))) t1 = max(as.numeric(rownames(diffMat))) 

  if(t2>max(as.numeric(rownames(diffMat)))) t2 = max(as.numeric(rownames(diffMat))) 

  hp$diff3599[i] = diffMat[which(rownames(diffMat)==t1),which(colnames(diffMat)==t2)]
}


hp$diff5599 = NA
for(i in 1:nrow(hp)){
  hpp = hp[i,]
  st_geometry(hpp) <- NULL
  t1=hpp[1,'t2055']
  t2=hpp[1,'t2099']
   if(t1>max(as.numeric(rownames(diffMat)))) t1 = max(as.numeric(rownames(diffMat))) 

    if(t2>max(as.numeric(rownames(diffMat)))) t2 = max(as.numeric(rownames(diffMat))) 

  hp$diff5599[i] = diffMat[which(rownames(diffMat)==t1),which(colnames(diffMat)==t2)]
}

ggplot(subset(hp,pred>.009) )+
   geom_sf(aes(fill=diff5599,color=diff5599), size=2.5) + 
   scale_colour_distiller(palette='RdYlGn') +
   scale_fill_distiller(palette='RdYlGn') + 
   # geom_sf(data=rL,size=1,colour='black',fill=NA)+
   theme( axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()
   ) +
   coord_sf()
 
saveRDS(hp,'HadSpatialPredictionsBerried_Jan2025.rds')