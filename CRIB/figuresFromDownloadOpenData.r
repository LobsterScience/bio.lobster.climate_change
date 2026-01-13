##crib data
require(sf)
require(ggplot2)
require(bio.lobster)
require(devtools)
la()


x = read.csv(file.path('C:/Users/cooka/Downloads/EN_ClimateRiskIndex_Spatial-CanEEZ_145Spp_HomAm.csv'))
x = subset(x,longitude< -57 & latitude <49)
xs = st_as_sf(x,coords=c('longitude','latitude'),crs=4326)


ggplot(xs)+geom_sf(size=7.1)

ggplot()+geom_sf(data=xs,size=7.1,aes(fill=Sensitivity,colour=Sensitivity))+scale_fill_viridis_c()+scale_colour_viridis_c()
ggplot()+geom_sf(data=xs,size=7.1,aes(fill=S.Cumulative.impacts,colour=S.Cumulative.impacts))+scale_fill_viridis_c()+scale_colour_viridis_c()

ggplot()+geom_sf(data=xs,size=7.1,aes(fill=AC.Thermal.habitat.availability,colour=AC.Thermal.habitat.availability))+scale_fill_viridis_c()+scale_colour_viridis_c()
ggplot()+geom_sf(data=xs,size=7.1,aes(fill=E.Time.of.climate.emergence,colour=E.Time.of.climate.emergence))+scale_fill_viridis_c()+scale_colour_viridis_c()



b = ggplot()+geom_sf(data=xs,size=2,aes(fill=Overall.Risk,colour=Overall.Risk))+scale_fill_viridis_d()+scale_colour_viridis_d()+theme_test_adam()
b1 = ggplot()+geom_sf(data=xs,size=2,aes(fill=Sensitivity.risk,colour=Sensitivity.risk))+scale_fill_viridis_d()+scale_colour_viridis_d()+theme_test_adam()
b2 = ggplot()+geom_sf(data=xs,size=2,aes(fill=Adaptive.capacity.risk,colour=Adaptive.capacity.risk))+scale_fill_viridis_d()+scale_colour_viridis_d()+theme_test_adam()
b3 = ggplot()+geom_sf(data=xs,size=2,aes(fill=Exposure.risk,colour=Exposure.risk))+scale_fill_viridis_d()+scale_colour_viridis_d()+theme_test_adam()


cowplot::plot_grid(b1,b3,b,nrow = 2)

require(terra)

xss = vect(xs)
r =  rast(ext(xss),resolution=.25)
plot(rasterize(xss,r,field='Sensitivity'))
plot(rasterize(xss,r,field='Adaptive.capacity'))
plot(rasterize(xss,r,field='Exposure'))
plot(rasterize(xss,r,field='Vulnerability'))



