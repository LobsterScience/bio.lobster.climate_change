#7.Can_CV_Model.r



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
#survey = survey[sample(1:nrow(survey),10000),]
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
names(mi) = c('CanBTr','OFFSETcorr')
mi$OFFSETcorr = mi$OFFSETcorr * (1/mi$OFFSETcorr[which(mi$CanBT==14)]) # 17m is at 14C, so go up or down from there
# 17m is at 14C, so go up or down from there
survey$CanBTr = round(survey$CanBT*2)/2
survey = dplyr::full_join(survey,mi)
survey$OFFSETcorr[which(survey$CanBTr< -.7)] <- min(mi$OFFSETcorr)
survey$OFFSETcorr[which(survey$CanBTr> 17)] <- max(mi$OFFSETcorr)

i = which(survey$OFFSET_METRIC == 'Number of traps')
survey$OFFSET[i] = survey$OFFSET[i] * pi*(.014^2) * survey$OFFSETcorr[i]

survey$LO = log(survey$OFFSET)
survey = subset(survey,OFFSET>0.00001 & OFFSET< 0.12)
survey$BT = survey$CanBT
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
#plot(spde)

# Add on the barrier mesh component:
bspde <- sdmTMBextra::add_barrier_mesh(
  spde, ns_coast, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)


#issues with fitting on a biweekly moving to quarters
survey$m = month(survey$DATE) 
survey$Q = ifelse(survey$m %in% c(10,11,12),1,ifelse(survey$m %in% c(1,2,3),2,ifelse(survey$m %in% c(4,5,6),3,4)))
survey$Time = survey$YEAR+survey$Q/4

k_folds <-5

data = as_tibble(survey)
##make the folds for model validation

data$IDS = "I"
data = cv_SpaceTimeFolds(data,idCol = 'IDS',nfolds=k_folds)
source('C:/Users/cooka/Documents/git/bio.lobster.climate_change/LobsterPA_Scripts/setUpCV_Table.r')

models=c('m1','m2','m3')
if ("m1" %in% models) {
  
  mod.label <- "m1" 
  
  
m = sdmTMB(pa~
               s(lZ,k=3)+s(BT,k=4)+Q,
             data=data,
            offset = 'LO',
             time='YEAR', 
             mesh=bspde,
             family=binomial(link='logit'),
             spatial='on',
             spatiotemporal='ar1')

  m_cv <- sdmTMB_cv(
    				pa~
               s(lZ,k=3)+s(BT,k=4)+Q,
             data=data,
            offset = 'LO',
             time='YEAR', 
             mesh=bspde,
             family=binomial(link='logit'),
             spatial='on',
             spatiotemporal='ar1',
             fold_ids = 'fold_id',
    		 k_folds=k_folds
             )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  m1 <-m
  m1_cv <-m_cv
  
}

if ("m2" %in% models) {
  
  mod.label <- "m2" 
  
  
m = sdmTMB(pa~
               s(lZ,k=3)+s(BT,k=4),
             data=data,
            offset = 'LO',
             time='YEAR', 
             mesh=bspde,
             family=binomial(link='logit'),
             spatial='on',
             spatiotemporal='ar1')

  m_cv <- sdmTMB_cv(
    				pa~
               s(lZ,k=3)+s(BT,k=4),
             data=data,
            offset = 'LO',
             time='YEAR', 
             mesh=bspde,
             family=binomial(link='logit'),
             spatial='on',
             spatiotemporal='ar1',
             fold_ids = 'fold_id',
    		 k_folds=k_folds
             )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  m2 <-m
  m2_cv <-m_cv
  
}


if ("m3" %in% models) {
  
  mod.label <- "m3" 
  
  
m = sdmTMB(pa~
               s(BT,k=3)+Q,
             data=data,
            offset = 'LO',
             time='YEAR', 
             mesh=bspde,
             family=binomial(link='logit'),
             spatial='on',
             spatiotemporal='ar1')

  m_cv <- sdmTMB_cv(
    				pa~
               s(BT,k=3)+Q,
             data=data,
            offset = 'LO',
             time='YEAR', 
             mesh=bspde,
             family=binomial(link='logit'),
             spatial='on',
             spatiotemporal='ar1',
             fold_ids = 'fold_id',
    		 k_folds=k_folds
             )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  m3 <-m
  m3_cv <-m_cv
  
}

write.csv(mod.select,'CanModSel.csv')