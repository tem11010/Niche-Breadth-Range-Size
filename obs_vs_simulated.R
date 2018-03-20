library(ape)
library(gdata)
library(dplyr)
library(raster)
library(ggplot2)
library(caper)

#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#
# get pelargonium distribution data for analysis
pel<-read.csv('C:\\Users\\Tim\\Desktop\\meeting cindi\\dbase_good_24_Feb_2017.csv', header=T, na.string='.')

pel.areas<- read.csv("C:\\Users\\Tim\\Desktop\\geog_overlap\\species_areas_in_grid.csv")[, -1]


# clean up in same way as for calculating geog area (from: "calculate areas of ahulls.R code")
library(dplyr)
library(gdata)
pel <- pel %>% dplyr::select(TAXNAME, LATITUDE, LONGITUDE)
pel <- unique(pel)
#  subset to secies with < 2 observations
detach("package:plyr")
pel.sub <- drop.levels(
  pel %>%
    group_by(TAXNAME) %>%
    mutate(total=n()) %>%
    filter(total >3 )) %>%
  data.frame() %>%
  drop.levels

mytable <- table(pel.sub$TAXNAME)
pie(mytable)
min(mytable)
levels(pel.sub$TAXNAME)

tree <-read.tree("C:\\Users\\Tim\\Desktop\\geog_overlap\\tree.sub.u.tre")
levels(as.factor(tree$tip.label))

which(tree$tip.label %in% levels(pel.sub$TAXNAME))
#subset to only pels on phylogeny
pel.sub <-  drop.levels(filter(pel.sub, TAXNAME %in% tree$tip.label))
levels(pel.sub$TAXNAME)
tips <- levels(pel.sub$TAXNAME)
tree$tip.label %in% tips
tree.tips <- as.data.frame(tree$tip.label)
colnames(tree.tips) <- "tree.tips"



drop <- filter(tree.tips, !tree.tips %in% tips)
tree.trim <- drop.tip(phy=tree, tip=as.character(drop[, 1]))
tree.tips <- as.data.frame(tree.trim$tip.label)

plot(tree.trim, cex=0.4)


levels(pel.sub$TAXNAME) %in% tree.trim$tip.label

tree$tip.label %in% levels(pel.sub$TAXNAME)

levels(pel.sub$TAXNAME)
levels(as.factor(tree.trim$tip.label))


supp.spp <- read.csv("C:\\Users\\Tim\\Desktop\\geog_overlap\\final_spp_list.csv")

levels(supp.spp$Species) %in% tree$tip.label
levels(pel.sub$TAXNAME)

pel.sub <-filter(pel.sub, TAXNAME !="grandiflorum")


tree <-  drop.tip(tree, "grandiflorum")

#-----------------------------------------------------#
# extract environmental data for observed layers
#use projection for South African Data
crs_use <- "+proj=aea +lat_1=-24 +lat_2=-32 +lat_0=0 +lon_0=24 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
crs_use_obj <- CRS(crs_use, doCheckCRSArgs=TRUE)

bio1  <- raster("D:\\Schulze Data\\grids\\tmean13c\\prj.adf") # used as scaffold for analyses
bio12  <- raster("D:\\Schulze Data\\grids\\gmap\\prj.adf")
bio1 <- raster::resample(bio1, bio1)
bio12 <- raster::resample(bio12, bio1)

res(bio1)
res(bio12)
crs(bio1) <- crs_use_obj
crs(bio12) <- crs_use_obj
bio1.a <- aggregate(bio1, fact = 8, fun = mean)
bio12.a <- aggregate(bio12, fact = 8, fun = mean)

res(bio1.a)
res(bio12.a)
plot(bio1.a)
plot(bio12.a)


#scale data
values(bio1.a) <-as.numeric(scale(values(bio1.a))) 
values(bio12.a) <-as.numeric(scale(values(bio12.a))) 

# check lengths of observed raster layers (will need this for later)
length(values(bio1.a)[!is.na(values(bio1.a))])
length(values(bio12.a)[!is.na(values(bio12.a))])


# extract bio1 and bio12 values for observed env layers
head(pel.sub)
pel.sub.env <- as.data.frame(pel.sub$TAXNAME)
colnames(pel.sub.env) <- "TAXNAME"

pel.sub.env$bio1 <- raster::extract(bio1.a, pel.sub[, c("LONGITUDE", "LATITUDE")])
pel.sub.env$bio12 <- raster::extract(bio12.a, pel.sub[, c("LONGITUDE", "LATITUDE")])


# now make summary values for each species, here I'm only interested in MAD values

spp.env <-  pel.sub.env %>% group_by(TAXNAME) %>%
  
                summarize(bio1.mad = mad(bio1, na.rm = TRUE), 
                          bio12.mad = mad(bio12, na.rm = TRUE), 
                          bio1.min = min(bio1, na.rm = TRUE), 
                          bio1.max = max(bio1, na.rm = TRUE), 
                          bio12.min = min(bio12, na.rm = TRUE), 
                          bio12.max = max(bio12, na.rm = TRUE)) %>% data.frame

spp.env.n <- pel.sub.env %>% group_by(TAXNAME) %>% tally %>% data.frame 


spp.env$n <- spp.env.n$n
spp.env


colnames(spp.env)
spp.env.area <-  left_join(spp.env, pel.areas, by = 'TAXNAME')

spp.env.area
spp.env.area$area_cells <- log(spp.env.area$area_cells)

hist(spp.env.area$area_cells)


plot(area_cells~bio12.mad, data = spp.env.area)
summary(lm(area_cells~bio12.mad, data = spp.env.area))
#-----------------------------------------------------#
# now run caper analysis for observed data

library(caper)
phy <- tree
dat <- spp.env.area


#check
length(phy$tip.label)
length(dat$TAXNAME)
phy$tip.label %in% as.character(dat$TAXNAME)
plot(phy)
hist(spp.env.area$area_cells)

cdat <- comparative.data(data=spp.env.area, phy=phy, names.col="TAXNAME")

coef(summary(pgls(area_cells ~ cdat$data[,"bio1.mad"], cdat, lambda = "ML")))[2,1] 
coef(summary(pgls(area_cells ~ cdat$data[,"bio12.mad"], cdat, lambda = "ML")))[2,1] 

# these are the values to use for the histogram plot (fig 3 in paper)

#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#

###################### null data#######################

# need to:
# read in inla layers
# read in random layers
# extract inla layers
# extract random layer
# calc mad inla layers
# calc mad random layer
# caper inla
# caper random
# plot inla vs random
# produce contours and histograms


#----------------------------------------------------------------------------------------------------------------------------#
# read in inla and random layers
mat.inla <- read.csv("C:\\Users\\Tim\\Desktop\\geog_overlap\\pel_pred_MAT_inla.csv")[-1]
map.inla <- read.csv("C:\\Users\\Tim\\Desktop\\geog_overlap\\pel_pred_MAP_inla.csv")[-1]
mat.rand <- read.csv("C:\\Users\\Tim\\Desktop\\geog_overlap\\pel_pred_MAT_rand.csv")[-1]
map.rand <- read.csv("C:\\Users\\Tim\\Desktop\\geog_overlap\\pel_pred_MAP_rand.csv")[-1]

# these layers have the extracted data for each species for each random layer

# now summarize grouping by species for each layer


# first for inla models
mat.inla.spp <- mat.inla %>% group_by(TAXNAME) %>%
                              summarize_if(is.numeric, mad ,na.rm = TRUE) %>% data.frame
mat.inla.spp <- filter(mat.inla.spp, TAXNAME %in% spp.env.area$TAXNAME)
mat.inla.spp.area <- left_join(mat.inla.spp, spp.env.area, by = "TAXNAME")
library(reshape2)
mat.inla.spp.area.m <- melt(mat.inla.spp.area, id.vars = c(colnames(spp.env.area)))

map.inla.spp <- map.inla %>% group_by(TAXNAME) %>%
  summarize_if(is.numeric, mad ,na.rm = TRUE) %>% data.frame
map.inla.spp <- filter(map.inla.spp, TAXNAME %in% spp.env.area$TAXNAME)
map.inla.spp.area <- left_join(map.inla.spp, spp.env.area, by = "TAXNAME")
map.inla.spp.area.m <- melt(map.inla.spp.area, id.vars = c(colnames(spp.env.area)))

# then for random models
mat.rand.spp <- mat.rand %>% group_by(TAXNAME) %>%
  summarize_if(is.numeric, mad ,na.rm = TRUE) %>% data.frame
mat.rand.spp <- filter(mat.rand.spp, TAXNAME %in% spp.env.area$TAXNAME)
mat.rand.spp.area <- left_join(mat.rand.spp, spp.env.area, by = "TAXNAME")
mat.rand.spp.area.m <- melt(mat.rand.spp.area, id.vars = c(colnames(spp.env.area)))

map.rand.spp <- map.rand %>% group_by(TAXNAME) %>%
  summarize_if(is.numeric, mad ,na.rm = TRUE) %>% data.frame
map.rand.spp <- filter(map.rand.spp, TAXNAME %in% spp.env.area$TAXNAME)
map.rand.spp.area <- left_join(map.rand.spp, spp.env.area, by = "TAXNAME")
map.rand.spp.area.m <- melt(map.rand.spp.area, id.vars = c(colnames(spp.env.area)))


# plot them
ggplot()+
  geom_point(data = mat.inla.spp.area.m, aes(x = value, y = area_cells, group = variable))+
  geom_point(data = mat.rand.spp.area.m, aes(x = value, y = area_cells, group = variable), color = "grey")

ggplot()+
  geom_density_2d(data = mat.inla.spp.area.m, aes(x = value, y = area_cells), color = "darkred")+
  geom_density_2d(data = mat.rand.spp.area.m, aes(x = value, y = area_cells), color = "blue")

ggplot()+
  geom_point(data = map.inla.spp.area.m, aes(x = value, y = area_cells, group = variable))+
  geom_point(data = map.rand.spp.area.m, aes(x = value, y = area_cells, group = variable), color = "grey")

ggplot()+
  geom_density_2d(data = map.inla.spp.area.m, aes(x = value, y = area_cells), color = "darkred")+
  geom_density_2d(data = map.rand.spp.area.m, aes(x = value, y = area_cells), color = "blue")



#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#

# ok. now lets get them ready for the caper analysis. I should just be able to use a for-loop

dat.inla.mat <- mat.inla.spp.area
dat.inla.map <- map.inla.spp.area
dat.rand.mat <- mat.rand.spp.area
dat.rand.map <- map.rand.spp.area


phy <- phy

inla.mat.comp <- comparative.data(phy = phy , data = dat.inla.mat, names.col="TAXNAME")
inla.map.comp <- comparative.data(phy = phy , data = dat.inla.map, names.col="TAXNAME")
rand.mat.comp <- comparative.data(phy = phy , data = dat.inla.mat, names.col="TAXNAME")
rand.map.comp <- comparative.data(phy = phy , data = dat.inla.map, names.col="TAXNAME")

colnames(inla.mat.comp$data)


sims.coefs <- data.frame(mat.inla = 1:1000, 
                         map.inla = 1:1000, 
                         mat.rand = 1:1000, 
                         map.rand = 1:1000)

for (i in 1:1000){

  sims.coefs[i, "mat.inla"]<- coef(summary(pgls(area_cells ~ inla.mat.comp$data[,i], inla.mat.comp, lambda = "ML")))[2,1]
  sims.coefs[i, "map.inla"]<- coef(summary(pgls(area_cells ~ inla.map.comp$data[,i], inla.map.comp, lambda = "ML")))[2,1]
  sims.coefs[i, "mat.rand"]<- coef(summary(pgls(area_cells ~ rand.mat.comp$data[,i], rand.mat.comp, lambda = "ML")))[2,1]
  sims.coefs[i, "map.rand"]<- coef(summary(pgls(area_cells ~ rand.map.comp$data[,i], rand.map.comp, lambda = "ML")))[2,1]
  print(i)
}

library(cowplot)
p.mat.inla <- 
ggplot()+
  geom_histogram(data = sims.coefs, aes(x = mat.inla), fill = "white",color = "black")+
  geom_vline(xintercept = quantile(sims.coefs$mat.inla, c(0.025, 0.975))[1], lty = 2)+
  geom_vline(xintercept = quantile(sims.coefs$mat.inla, c(0.025, 0.975))[2], lty = 2)+
  geom_vline(xintercept = coef(summary(pgls(area_cells ~ cdat$data[,"bio1.mad"], cdat, lambda = "ML")))[2,1], color = "red")

p.map.inla <- 
ggplot()+
  geom_histogram(data = sims.coefs, aes(x = map.inla), fill = "white",color = "black")+
  geom_vline(xintercept = quantile(sims.coefs$map.inla, c(0.025, 0.975))[1], lty = 2)+
  geom_vline(xintercept = quantile(sims.coefs$map.inla, c(0.025, 0.975))[2], lty = 2)+
  geom_vline(xintercept = coef(summary(pgls(area_cells ~ cdat$data[,"bio12.mad"], cdat, lambda = "ML")))[2,1], color = "red")

p.mat.rand <- 
ggplot()+
  geom_histogram(data = sims.coefs, aes(x = mat.rand), fill = "white",color = "black")+
  geom_vline(xintercept = quantile(sims.coefs$mat.rand, c(0.025, 0.975))[1], lty = 2)+
  geom_vline(xintercept = quantile(sims.coefs$mat.rand, c(0.025, 0.975))[2], lty = 2)+
  geom_vline(xintercept = coef(summary(pgls(area_cells ~ cdat$data[,"bio1.mad"], cdat, lambda = "ML")))[2,1], color = "red")

p.map.rand <- 
ggplot()+
  geom_histogram(data = sims.coefs, aes(x = map.rand), fill = "white",color = "black")+
  geom_vline(xintercept = quantile(sims.coefs$map.rand, c(0.025, 0.975))[1], lty = 2)+
  geom_vline(xintercept = quantile(sims.coefs$map.rand, c(0.025, 0.975))[2], lty = 2)+
  geom_vline(xintercept = coef(summary(pgls(area_cells ~ cdat$data[,"bio12.mad"], cdat, lambda = "ML")))[2,1], color = "red")


plot_grid(p.mat.inla, p.map.inla, 
          p.mat.rand, p.map.rand, ncol = 2)

