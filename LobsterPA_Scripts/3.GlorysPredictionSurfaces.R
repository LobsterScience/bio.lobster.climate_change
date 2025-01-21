

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


#prediction grids from glorys by quarter of the year (lobster year line 31)
s = file.path(project.datadirectory('bio.lobster'),'Temperature Data','GLORYS','SummaryFiles')
k = dir(s,full.names=T)
k = k[grep('ShelfBoF',k)]
#k = k[-grep('2020_',k)]

ol = list()
m=0
for(i in 1:length(k)){
			h = readRDS(k[i])
			h$Date = as.Date(h$Date)
			h$m = month(h$Date)
			h$yr = unique(year(h$Date))
			h$Q = ifelse(h$m %in% c(10,11,12),1,ifelse(h$m %in% c(1,2,3),2,ifelse(h$m %in% c(4,5,6),3,4)))
			h = aggregate(bottomT~Q+X+Y+yr,data=h,FUN=mean)
			h = h %>% st_as_sf(coords=c('X','Y'),crs=4326) %>% st_transform(32620)
			st_geometry(h) = st_geometry(h)/1000
			st_crs(h) = 32620
			m=m+1
		ol[[m]] = h       	
			}
		

crs_utm20 <- 32620
da = dplyr::bind_rows(ol)	
da = st_as_sf(da)


ba = readRDS('~/git/bio.lobster.data/mapping_data/bathymetrySF.rds')
ba = ba %>% st_as_sf() 
st_geometry(ba) = st_geometry(ba)/1000
st_crs(ba) = 32620

				 ss = st_nearest_feature(da,ba)
       	 ds = st_distance(da,ba[ss,],by_element=T)
       	 st_geometry(ba) = NULL
       	 da$z = ba$z[ss]
       	 da$z_dist = as.numeric(ds)
#by quarter     	
saveRDS(da,'GlorysPredictSurface.rds')

