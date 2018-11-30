library(ape)
library(gdata)
library(dplyr)
library(raster)
library(ggplot2)
library(caper)
library(phytools)
library(ape)

#------------------------------------------------------------------------------#

# Pelargonium distribution data for analysis

#------------------------------------------------------------------------------#

# distribution data
pel<-read.csv("data/dbase_good_24_Feb_2017.csv", 
              header=T, na.string='.')

# range size data (using variable area_cells)
pel.areas <- read.csv("data/species_areas_in_grid.csv")[, -1]


# clean up data as for main code
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

tree <-read.tree("data/tree.sub.u.tre")
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

# check tree
plot(tree.trim, cex=0.4)
levels(pel.sub$TAXNAME) %in% tree.trim$tip.label
tree$tip.label %in% levels(pel.sub$TAXNAME)

# drop one species with range issues
pel.sub <-filter(pel.sub, TAXNAME !="grandiflorum")
tree <-  drop.tip(tree, "grandiflorum")

#-------------------------------------------------------------------------------#

# extract environmental data for observed layers

#-------------------------------------------------------------------------------#

#use projection for South African Data
crs_use <- "+proj=aea +lat_1=-24 +lat_2=-32 +lat_0=0 +lon_0=24 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
crs_use_obj <- CRS(crs_use, doCheckCRSArgs=TRUE)

# bio1 in MAT, bio 12 is MAP

bio1  <- raster("MAT/prj.adf") # used as scaffold for analyses
bio12  <- raster("MAP/prj.adf")
bio1 <- raster::resample(bio1, bio1)
bio12 <- raster::resample(bio12, bio1)

# check resolution and change to same level as simulated data
res(bio1)
res(bio12)
crs(bio1) <- crs_use_obj
crs(bio12) <- crs_use_obj
bio1.a <- aggregate(bio1, fact = 8, fun = mean)
bio12.a <- aggregate(bio12, fact = 8, fun = mean)

res(bio1.a)
res(bio12.a)

#plot them
par(mfrow =c(1, 2))
plot(bio1.a, main = "MAT")
plot(bio12.a, main = "MAP")


#scale data
values(bio1.a) <-as.numeric(scale(values(bio1.a))) 
values(bio12.a) <-as.numeric(scale(values(bio12.a))) 

# check lengths of observed raster layers (will need this for later)
length(values(bio1.a)[!is.na(values(bio1.a))])
length(values(bio12.a)[!is.na(values(bio12.a))])


# extract bio1 (mat) and bio12 (map) values for observed env layers
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

# add count of number of collections
spp.env.n <- pel.sub.env %>% group_by(TAXNAME) %>% tally %>% data.frame 
spp.env$n <- spp.env.n$n
head(spp.env)
colnames(spp.env)

# merge datasets
spp.env.area <-  left_join(spp.env, pel.areas, by = 'TAXNAME')
head(spp.env.area)
spp.env.area$area_cells <- log(spp.env.area$area_cells)

# non-phylogenetic linear model of relationships between range size and niche breadth for
# observed climate data
summary(lm(area_cells~bio1.mad, data = spp.env.area))
summary(lm(area_cells~bio12.mad, data = spp.env.area))

#-------------------------------------------------------------------------------#

# Caper analysis for observed data

#-------------------------------------------------------------------------------#

# set up comparative dataframe
phy <- tree
dat <- spp.env.area

#checks
length(phy$tip.label)
length(dat$TAXNAME)
phy$tip.label %in% as.character(dat$TAXNAME)
plot(phy)
hist(spp.env.area$area_cells)

# make comparative dataframe
cdat <- comparative.data(data=spp.env.area, phy=phy, names.col="TAXNAME")

# run analyses, using ML to estimate values of lambda, and extract coefficients
coef(summary(pgls(area_cells ~ cdat$data[,"bio1.mad"], cdat, lambda = "ML")))[2,1] 
coef(summary(pgls(area_cells ~ cdat$data[,"bio12.mad"], cdat, lambda = "ML")))[2,1] 

# these are the regression coefficients (slopes) for the observed niche breadth-range size
# relationships

#-------------------------------------------------------------------------------#

# Caper analysis for simulated data

#-------------------------------------------------------------------------------#
# need to:
# read in inla layers
# read in random layers
# extract inla layers
# extract random layer
# calc mad for inla layers (i.e., Niche breadths for spatial simulations)
# calc mad for random layer  (i.e., Niche breadths for non-spatial simulations)
# caper inla (PGLS analyses)
# caper random (PGLS analyses)
# plot inla vs random 
# produce contours and histograms


#-------------------------------------------------------------------------------#
# read in inla and random layers generated by main_code
# in this code I am using the rainfall layer
# to compare to bio12 results

mat.inla <- read.csv("results/pel_pred_inla_scaled.csv")[-1]
mat.rand <- read.csv("results/pel_pred_rand_scaled.csv")[-1]
# these layers have the extracted data for each site for each species 
# for each simulated/randomized layer

head(mat.inla[, 990:1001])
head(mat.rand[, 990:1001])

# now summarize grouping by species for each layer
# first for inla models
mat.inla.spp <- mat.inla %>% group_by(TAXNAME) %>%
  summarize_if(is.numeric, mad ,na.rm = TRUE) %>% data.frame

# here the metric of niche breadth I am using is MAD, or median absolute deviation
# any other form could be used here (e.g. SD, range, etc.)
# the code for the observed data should be adjusted accordingly

head(mat.inla.spp[, 1:10])

#make sure I am working only with the species I have in my target dataset
mat.inla.spp <- filter(mat.inla.spp, TAXNAME %in% spp.env.area$TAXNAME)

#add range size data to this dataframe
mat.inla.spp.area <- left_join(mat.inla.spp, spp.env.area, by = "TAXNAME")

# create long format dataframe to facilite subsetting later
library(reshape2)
mat.inla.spp.area.m <- melt(mat.inla.spp.area, id.vars = c(colnames(spp.env.area)))


# then for random models
mat.rand.spp <- mat.rand %>% group_by(TAXNAME) %>%
  summarize_if(is.numeric, mad ,na.rm = TRUE) %>% data.frame
mat.rand.spp <- filter(mat.rand.spp, TAXNAME %in% spp.env.area$TAXNAME)
mat.rand.spp.area <- left_join(mat.rand.spp, spp.env.area, by = "TAXNAME")

# create long format dataframe to facilite subsetting later
mat.rand.spp.area.m <- melt(mat.rand.spp.area, id.vars = c(colnames(spp.env.area)))


# plot them
ggplot()+
  geom_point(data = mat.inla.spp.area.m, aes(x = value, y = area_cells, group = variable))+
  geom_point(data = mat.rand.spp.area.m, aes(x = value, y = area_cells, group = variable), 
             color = "grey")

# each row here represents a species, and each point an estimate of it's niche breadth
# based on simulated climate layers either with (black), or without (grey) spatial structure.
# Note that, for the large ranged species, values of niche breadth are 'pulled'
# towards the MAD value for the whole domain. This an avenue of research that
# I am hoping to persue in the future. 

# we can also produce density plots, like in the paper
ggplot()+
  geom_density_2d(data = mat.inla.spp.area.m, aes(x = value, y = area_cells), 
                  color = "darkred", size = 1.5)+
  geom_density_2d(data = mat.rand.spp.area.m, aes(x = value, y = area_cells), color = "blue", 
                  size = 1.5)


#-------------------------------------------------------------------------------#

# now lets get them ready for the caper (PGLS) analysis
# using a for-loop

# columns should have each of the NB estimates for 1000 layers (rep.1 - rep.1000),  
# species range sizes, and species names
colnames(mat.inla.spp.area)[c(1, 1000:1009)]


dat.inla.mat <- mat.inla.spp.area
dat.rand.mat <- mat.rand.spp.area
phy <- phy

inla.mat.comp <- comparative.data(phy = phy , data = dat.inla.mat, names.col="TAXNAME")
rand.mat.comp <- comparative.data(phy = phy , data = dat.rand.mat, names.col="TAXNAME")

# set up dataframe to store coefficients
sims.coefs <- data.frame(mat.inla = 1:1000, 
                         mat.rand = 1:1000)

for (i in 1:1000){
  
  sims.coefs[i, "mat.inla"]<- coef(summary(pgls(area_cells ~ inla.mat.comp$data[,i], inla.mat.comp, lambda = "ML")))[2,1]
  sims.coefs[i, "mat.rand"]<- coef(summary(pgls(area_cells ~ rand.mat.comp$data[,i], rand.mat.comp, lambda = "ML")))[2,1]
  print(i)
}


# the regression coefficients for the NB-RS relationships for spatial and non-spatial simulations
# are now stored in sims.coefs. 

# We can now produce the key plots: the distributions of regression coefficients, under
# our null models, as histograms, and the observed NB-RS regression coefficients. 

library(cowplot)
p.mat.inla <- 
  ggplot()+
  geom_histogram(data = sims.coefs, aes(x = mat.inla), fill = "white",color = "black")+
  geom_vline(xintercept = quantile(sims.coefs$mat.inla, c(0.025, 0.975))[1], lty = 2)+
  xlim(c(-0.5, 5.5))+
  geom_vline(xintercept = quantile(sims.coefs$mat.inla, c(0.025, 0.975))[2], lty = 2)+
  # add observed coeffiecient
  geom_vline(xintercept = coef(summary(pgls(area_cells ~ cdat$data[,"bio12.mad"], cdat, 
                                            lambda = "ML")))[2,1], color = "red", size = 1.5)

p.mat.rand <- 
  ggplot()+
  geom_histogram(data = sims.coefs, aes(x = mat.rand), fill = "white",color = "black")+
  geom_vline(xintercept = quantile(sims.coefs$mat.rand, c(0.025, 0.975))[1], lty = 2)+
  geom_vline(xintercept = quantile(sims.coefs$mat.rand, c(0.025, 0.975))[2], lty = 2)+
  xlim(c(-0.5, 5.5))+
  # add observed coeffiecient
  geom_vline(xintercept = coef(summary(pgls(area_cells ~ cdat$data[,"bio12.mad"], cdat, 
                                            lambda = "ML")))[2,1], color = "red", size = 1.5)

plot_grid(p.mat.inla,
          p.mat.rand, ncol = 1)
