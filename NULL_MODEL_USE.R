#-------------------------------------------------------------------------------------------------#

# Randomzing environment models 

#-------------------------------------------------------------------------------------------------#


library(geostatsp)
library(maptools)
library(shapefiles)
#library(ecospat)
library(INLA)
library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
rm(list=ls())
set.seed(0409)


# get pelargonium distribution data for analysis
pel<-read.csv('C:\\Users\\Tim\\Desktop\\meeting cindi\\dbase_good_24_Feb_2017.csv', header=T, na.string='.')

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

#--------------------------------------------------------------------------------------------------#

# INLA

#--------------------------------------------------------------------------------------------------#
#use projection for South African Data
crs_use <- "+proj=aea +lat_1=-24 +lat_2=-32 +lat_0=0 +lon_0=24 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
crs_use_obj <- CRS(crs_use, doCheckCRSArgs=TRUE)

mat.temp  <- raster("D:\\Schulze Data\\grids\\tmean13c\\prj.adf") # used as scaffold for analyses
mat  <- raster("D:\\Schulze Data\\grids\\gmap\\prj.adf")
#mat  <- raster("D:\\Schulze Data\\grids\\tmean13c\\prj.adf")
mat <- raster::resample(mat, mat.temp)
res(mat)
crs(mat) <- crs_use_obj
mat.a <- aggregate(mat, fact = 8, fun = mean)

res(mat.a)
plot(mat.a)

#scale data
mat.a.unscaled <- mat.a
values(mat.a) <-as.numeric(scale(values(mat.a))) 
ext <- extent(mat.a)
mat.a


library(usdm)
library(cowplot)
plot(Variogram(mat.a))
#---------------------------------------------------------------------#
# getting maps for presentation

#read in simulated layer
preds <- read.csv("predicted_MAT_inla.csv")
preds <- read.csv("predicted_MAP_inla.csv")


length(preds[, 1])
length(values(mat.a)[!is.na(values(mat.a))])
pred.mat <- mat.a
values(pred.mat)[!is.na(values(pred.mat))] <- preds[, 15]
plot(pred.mat)
plot(Variogram(pred.mat))


var.mat <- (Variogram(mat.a))
var.mat.df <- var.mat@variogram
var.pred  <- (Variogram(pred.mat)) 
var.pred.df <- var.pred@variogram


ggplot()+
  geom_point(data = var.mat.df, aes(x = distance, y = gamma))+
  geom_point(data = var.pred.df, aes(x = distance, y = gamma), color = "blue", alpha = 0.6)

#writeRaster(pred.mat, "MAP_agg_sim10.asc", format = "ascii")


#---------------------------------------------------------------------#



plot(mat.a)
plot(mat.a.unscaled)
hist(mat.a)
nrow(mat.a)
ncol(mat.a)


# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(mat.a)
colnames(dat) <- c("x", "y", "z")

create_grid <- function(x){
  
  if (class(x) == "matrix") {
    
    dat.df <- as.data.frame(dat)[
      with(as.data.frame(dat), order(y, decreasing = TRUE)),
      ]
    # get the rows and columns into the correct order
    xy_grid <-
      dcast( data =  dat.df, 
             formula = y ~ x, value.var = "z") 
    
    xy_grid <- xy_grid[
      with(xy_grid, order(y, decreasing = TRUE)),
      ]
    
    y_vals <- xy_grid[,1]
    x_vals <- as.numeric(colnames(xy_grid)[2:ncol(xy_grid)])
    xy_grid <- as.matrix(xy_grid)
    xy_grid2 <- xy_grid[, -1]
    
    l1 <- list(y_vals = y_vals, x_vals = x_vals, xy_grid2 = xy_grid2)
    return(l1)
  } 
  else { 
    print("x is not a matrix") 
  }
  
}

mat.grid <- create_grid(dat)
mat.grid

nrow(mat.grid$xy_grid2)
ncol(mat.grid$xy_grid2)

#-------------------------------------------------------------------#

# variables needed for the f function and the model
nrow=nrow(mat.grid$xy_grid2) # should be same as alt.xy_grid2
ncol=ncol(mat.grid$xy_grid2) # should be same as alt.xy_grid2
n = nrow*ncol
s.noise = 1
y = inla.matrix2vector(mat.grid$xy_grid2)
node = 1:n

formula= y ~ 1 + f(node, model="matern2d", nu=1, nrow=nrow, ncol=ncol,
                   hyper = list(range = list(param =c(1, 1),
                                             prior = "normal", param = c(0, 3)),
                                prec = list(param=c(1, 1)))) 


# Need to check initial values for optimization
#inla.set.control.inla.default

# get data into format for the model
mat.data=data.frame(y=y, node=node)
class(mat.data)
## fit the model
result=inla(formula, family="gaussian", data=mat.data, verbose=TRUE,
            control.predictor = list(compute = TRUE),
            control.family = list(hyper = list(theta = list(initial = log(1/s.noise^2),
                                                            fixed = FALSE))),
            control.compute=list(cpo=FALSE), 
            control.inla = list(h = 1e-5, tolerance=1e-5),  
            keep=TRUE)

result = inla.rerun(result)
result$mode$mode.status

summ <- summary(result)
summ

## plot the posterior mean for `predictor' and compare with the truth
#dev.new()
INLA:::inla.display.matrix(mat.grid$xy_grid2)
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol))

# get predicted values into matrix
pred_mat <- inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)
colnames(pred_mat) <- mat.grid$x_vals
rownames(pred_mat) <- mat.grid$y_vals

# convert into a format that can be plotted in ggplot
pred_df <- as.data.frame(pred_mat)
colnames(pred_df)<- paste("V", colnames(pred_df), sep="_")

pred_df <- pred_df %>%
  rownames_to_column("y") %>%
  gather(x,z,starts_with("V"))

head(pred_df)
list_x <- unlist(strsplit(pred_df$x, split='_', fixed=TRUE))
pred_df$x <- as.numeric(list_x[c(seq(from=2, to=length(list_x), by = 2))])
pred_df$y <- as.numeric(pred_df$y)

#plot it
ggplot()+
  geom_raster(data = pred_df, aes(x = x, y=y, fill = z))

# get the observed data into the same set up 
obs_mat <- mat.grid$xy_grid2
colnames(obs_mat) <- mat.grid$x_vals
rownames(obs_mat) <- mat.grid$y_vals

obs_df <- as.data.frame(obs_mat)
colnames(obs_df)<- paste("V", colnames(obs_df), sep="_")

obs_df <- obs_df %>%
  rownames_to_column("y") %>%
  gather(x,z,starts_with("V"))
head(obs_df)
list_x_obs <- unlist(strsplit(obs_df$x, split='_', fixed=TRUE))
obs_df$x <- as.numeric(list_x_obs[c(seq(from=2, to=length(list_x_obs), by = 2))])
obs_df$y <- as.numeric(obs_df$y)


# plot it in ggplot
ggplot()+
  geom_raster(data = obs_df, aes(x = x, y=y, fill = z))

pred_df2 <- pred_df[!is.na(obs_df[, 3]),]

ggplot()+
  geom_raster(data = pred_df2, aes(x = x, y=y, fill = z))


# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
#dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)


summary(lm(inla.matrix2vector(pred_mat)~inla.matrix2vector(obs_mat)) )

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
# run with nu = 2 and compare to 'result'

formula2= y ~ 1 + f(node, model="matern2d", nu=2, nrow=nrow, ncol=ncol,
                    hyper = list(range = list(param =c(1, 1),
                                              prior = "normal", param = c(0, 3)),
                                 prec = list(param=c(1, 1)))) 

result2=inla(formula2, family="gaussian", data=mat.data, verbose=TRUE,
             control.predictor = list(compute = TRUE),
             control.family = list(hyper = list(theta = list(initial = log(1/s.noise^2),
                                                             fixed = FALSE))),
             control.compute=list(cpo=FALSE), 
             control.inla = list(h = 1e-5, tolerance=1e-5),  
             keep=TRUE)

result2 = inla.rerun(result2)

summ2 <- summary(result2)
summ2

pred_mat2 <- inla.vector2matrix(result2$summary.linear.predictor$mean,nrow,ncol)
colnames(pred_mat2) <- mat.grid$x_vals
rownames(pred_mat2) <- mat.grid$y_vals



# now plot the relationship between the two models with different nu values
par(mfrow=c(1, 1))
#dev.new()
#dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(pred_mat2)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)


# seems that in both cases nu = 1 has a better fit
# and fit is better when data are scaled

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#

# Simulate data

# make dat raw data values first
dat <- rasterToPoints(mat.a)
colnames(dat) <- c("x", "y", "z")

dmat <- as.matrix(dist(dat))
dim(dmat)
dev.off()
# get parameter estimates from the INLA model
summ
#summ2
beta0 <- 1.0894
sigma2e <- 1 / 3.765e+04# precision for gaussian observations
sigma2x <- 1 / 7.545e-01#precision for node

kappa <- sqrt(8) / 1.900e+01# range for the node
nu <- 1

mcor <- as.matrix(2 ^ (1 - nu) * (kappa * dmat) ^ nu *
                    besselK(dmat * kappa, nu) / gamma(nu))
diag(mcor) <- 1
mcov <- sigma2e * diag(ncol(dmat)) + sigma2x * mcor

#image(mcor)
# convert var covar to sd matrix
gc()
L <- chol(mcor) # convert var covar to sd matrix
gc()

# sort the values of the real gradient in ascending order...
dat <- data.frame(dat)
vals = as.numeric(sort(dat$z))

#set up empty matrix and fill with predicted values
preds <- matrix(nrow = nrow(dat), ncol = 1000)
for (i in 1:1000) {
  preds[, i] <-
    beta0 + drop(rnorm(ncol(dmat)) %*% L) # simulate data
  print(i)
  
  # replace GRF values with real ones in rank order...
  preds[, i] = vals[rank(preds[, i], ties.method = "random")]
}

preds <- as.data.frame(preds)
colnames(preds) <- paste("y", 1:1000, sep = "")

preds$x <- dat$x
preds$y <- dat$y


# check histograms to make sure the data match
hist(preds$y1)
hist(dat$z, add = TRUE, col = "red")


mat.a.df <- rasterToPoints(mat.a) %>% as.data.frame
ggplot(data = mat.a.df, aes(x = x, y = y, fill = prj)) +
  geom_raster() +
  ggtitle("Observed") +
  scale_fill_gradientn(colors = terrain.colors(10))

ggplot(data = preds, aes(x = x, y = y, fill = y36)) +
  geom_raster() +
  ggtitle("Predicted") +
  scale_fill_gradientn(colors = terrain.colors(10))

plot(mat.a)


max(dat[, 2])-min(dat[, 2])

15/0.25

60*30

head(preds[, 999:1002])
length(preds[,1000])

# unscale values of MAP before saving and extracting for species

for (i in 1:1000){
  
  preds[, i] <- preds[, i]*sd(values(mat.a.unscaled), na.rm = TRUE) + mean(values(mat.a.unscaled), na.rm = TRUE)
  
}





write.csv(preds, "predicted_MAT_inla_unscaled.csv")



# check 
#can use this method to make the predicted maps
pred.1 <- mat.a
values(pred.1)[!is.na(values(pred.1))] <- preds[, 20]
plot(pred.1)
# compare
plot(mat.a.unscaled)

#-------------------------------------------------------------------------------------------------#

#create rasters with simulated data

library(rgdal)
sims <- list()
sims_grid <- list()
sims_r <-list()
i=1
for (i in 1:1000){
  
  sims[[i]] <- as.data.frame(preds[, c("x", "y", paste("y", i, sep = ""))])
  sims[[i]] <- sims[[i]][
    with(sims[[i]], order(y, decreasing = TRUE)),
    ]
  sims_grid[[i]] <-
    dcast( data =  sims[[i]], 
           formula = y ~ x, value.var = paste("y", i, sep = "")) 
  
  sims_grid[[i]] <- sims_grid[[i]][
    with(sims_grid[[i]], order(y, decreasing = TRUE)),
    ]
  
  sims_grid[[i]] <- as.matrix(sims_grid[[i]])
  sims_grid[[i]] <- sims_grid[[i]][, -1]
  sims_r[[i]] <- raster::raster(sims_grid[[i]], template = mat.a)
  
}
#sims_r now holds raster layers of each predicted env layer (should be 1000)

plot(sims_r[[1]], col = topo.colors(20))

#check that details are the same (mins and maxes, etc)
mat.a
sims_r[[10]]

range(raster::extract(sims_r[[500]], pel.sub[, c("LONGITUDE", "LATITUDE")]))
plot(sims_r[[500]])
plot(mat.a)

#-------------------------------------------------------------------------------------------------#


# extract data for all pel species
pel.preds.inla <- as.data.frame(matrix(nrow = length(pel.sub$TAXNAME), ncol = 1000))
colnames(pel.preds.inla) <- paste("rep", 1:1000, sep = '.')

for (i in 1:1000) {
  
  pel.preds.inla[, i] <- raster::extract(sims_r[[i]], pel.sub[, c("LONGITUDE", "LATITUDE")])
  
  print(i)
}

pel.preds.inla$TAXNAME <- pel.sub$TAXNAME


# write extracted null data to csv

write.csv(pel.preds.inla, "pel_pred_MAT_inla_unscaled.csv")

#-------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# Completely randomized

#--------------------------------------------------------------------------------------------------#
# create completely randomized layers of climate vars and save them as text files. need to get the
# standard deviations too

# 6858 cells in temp layer
# 6936 cells in precip layer
length(na.omit(values(mat.a)))
length(na.omit(values(bio12)))

mat.a.r <- mat.a.unscaled

random_maps <- list()

for (i in 1:1000){
  
  random_maps[[i]] <- mat.a.r
  values(random_maps[[i]])[!is.na(values(random_maps[[i]]))] <- 
    values(random_maps[[i]])[!is.na(values(random_maps[[i]]))][sample(6858, 6858)] # change numbers!!
}

plot(random_maps[[500]])
hist(random_maps[[1000]])
hist(mat.a)

plot(mat.a)
plot(random_maps[[1]])

# extract data for all pel species
pel.preds.rand <- as.data.frame(matrix(nrow = length(pel.sub$TAXNAME), ncol = 1000))
colnames(pel.preds.rand) <- paste("rep", 1:1000, sep = '.')

for (i in 1:1000) {
  
  pel.preds.rand[, i] <- raster::extract(random_maps[[i]], pel.sub[, c("LONGITUDE", "LATITUDE")])
  
  print(i)
}

pel.preds.rand$TAXNAME <- pel.sub$TAXNAME


write.csv(pel.preds.rand, "pel_pred_MAT_rand_unscaled.csv")



