library(geostatsp)
library(maptools)
library(shapefiles)
library(INLA)
library(reshape2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(usdm)
library(cowplot)
library(rgdal)

rm(list=ls())
set.seed(0409)

#-----------------------------------------------------------------------------#
# create_grid function to make a matrix into a grid
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

#-------------------------------------------------------------------#

# get pelargonium distribution data for analysis
pel<-read.csv("data/dbase_good_24_Feb_2017.csv", header=T, na.string='.')

# clean up species data
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

# CLIMATE DATA

#--------------------------------------------------------------------------------------------------#

#use projection for South African Data
crs_use <- "+proj=aea +lat_1=-24 +lat_2=-32 +lat_0=0 +lon_0=24 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
crs_use_obj <- CRS(crs_use, doCheckCRSArgs=TRUE)

# read in a raster layer 
# NOTE: there were subtle differences between rainfall and temperature layers, and so I use
# the temperaure layer as my 'scaffold' layer. I did this by resampling the rainfall layer
# by the temperature layer. 


# mat.temp  =  scaffold layer
# mat = layer to be simulated

mat.temp  <- raster("MAT/prj.adf") # used as scaffold for analyses
mat  <- raster("MAP/prj.adf") # rainfall
#mat  <- raster("D:\\Schulze Data\\grids\\tmean13c\\prj.adf") # temperature
mat <- raster::resample(mat, mat.temp)

# check resolution of the raster
res(mat)
# set the coordiante reference system for raster layer
crs(mat) <- crs_use_obj

# aggregate cells for lower resolution
mat.a <- aggregate(mat, fact = 8, fun = mean)
res(mat.a)
plot(mat.a)

#scale the climate data 
mat.a.unscaled <- mat.a
values(mat.a) <-as.numeric(scale(values(mat.a))) 
ext <- extent(mat.a)
mat.a

# check the spatial structure in the observed layer
plot(Variogram(mat.a))


# create matrix from the aggregated raster
dat <- rasterToPoints(mat.a)
colnames(dat) <- c("x", "y", "z")
mat.grid <- create_grid(dat)
# mat.grid
# nrow(mat.grid$xy_grid2)
# ncol(mat.grid$xy_grid2)

#--------------------------------------------------------------------------------------------------#

# INLA

#--------------------------------------------------------------------------------------------------#
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


# If you need to check the initial values for optimization
# inla.set.control.inla.default

# get data into format for the model
mat.data=data.frame(y=y, node=node)
class(mat.data)
## fit the model
# NOTE: this may take some time to run if layer is large or at high resolution!

result=inla(formula, family="gaussian", data=mat.data, verbose=TRUE,
            control.predictor = list(compute = TRUE),
            control.family = list(hyper = list(theta = list(initial = log(1/s.noise^2),
                                                            fixed = FALSE))),
            control.compute=list(cpo=FALSE), 
            control.inla = list(h = 1e-5, tolerance=1e-5),  
            keep=TRUE)
# result contains the model output

# Havard Rue suggested checking the mode.status to confirm that the model has been parameterized 
# appropriately. Ideally the mode.status should be 0
result$mode$mode.status 

# if this does not return 0, you can rerun the model using the previously estimated parameters
# as the starting points
result = inla.rerun(result)
# and check the mode.status again
result$mode$mode.status

# store the parameter estimates to be used for simulation later
summ <- summary(result)
summ

## plot the posterior mean for `predictor' and compare with the truth
#dev.new()
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
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
  geom_raster(data = obs_df, aes(x = x, y=y, fill = z)) # observed

# remove the values that are NA in the observed dataset
pred_df2 <- pred_df[!is.na(obs_df[, 3]),]

ggplot()+
  geom_raster(data = pred_df2, aes(x = x, y=y, fill = z)) # modeled


# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
#dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)

#check correlation
summary(lm(inla.matrix2vector(pred_mat)~inla.matrix2vector(obs_mat)) )

#-------------------------------------------------------------------------------------------------#
# this is to rerun the models with a different parameter value, i.e., setting
# nu = 2 and compare to 'result'

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

#--------------------------------------------------------------------------------------------------#

# SIMULATING THE CLIMATE DATA BASED ON INLA MODEL

#--------------------------------------------------------------------------------------------------#

# make dat raw data values first
dat <- rasterToPoints(mat.a)
colnames(dat) <- c("x", "y", "z")

dmat <- as.matrix(dist(dat))
dim(dmat)
#dev.off()
# get parameter estimates from the INLA model
summ

# use these to set values below
beta0 <- -0.2157
sigma2e <- 1 / 35645.437# precision for gaussian observations
sigma2x <- 1 / 1.029#precision for node
kappa <- sqrt(8) / 19.574# range for the node

# use value of nu specified for model
nu <- 1

mcor <- as.matrix(2 ^ (1 - nu) * (kappa * dmat) ^ nu *
                    besselK(dmat * kappa, nu) / gamma(nu))
diag(mcor) <- 1
mcov <- sigma2e * diag(ncol(dmat)) + sigma2x * mcor


# convert varcovar matrix to sd matrix
# NOTE: computationally intensive!!
gc()
L <- chol(mcor) # convert var covar to sd matrix
gc()

# sort the values of the real gradient in ascending order...
dat <- data.frame(dat)
vals = as.numeric(sort(dat$z))

#set up empty matrix and fill with predicted values
preds <- matrix(nrow = nrow(dat), ncol = 1000)
dim(preds)
for (i in 1:1000) {
  preds[, i] <-
    beta0 + drop(rnorm(ncol(dmat)) %*% L) # simulate data
  print(i) # can suppress this
  
  # replace values in simulated layers with real ones in rank order
  # this is so that all layers have the same values in them 
  # (see Chapman 2010, citation in the paper)
  
  preds[, i] = vals[rank(preds[, i], ties.method = "random")]
}

preds <- as.data.frame(preds)
colnames(preds) <- paste("y", 1:1000, sep = "")

# add coordinates
preds$x <- dat$x
preds$y <- dat$y

# check histograms to make sure the data match
hist(preds$y1)
hist(dat$z, add = TRUE, col = "red")

# example of observed vs expected
mat.a.df <- rasterToPoints(mat.a) %>% as.data.frame
ggplot(data = mat.a.df, aes(x = x, y = y, fill = prj)) +
  geom_raster() +
  ggtitle("Observed") +
  scale_fill_gradientn(colors = terrain.colors(10))

# p
ggplot(data = preds, aes(x = x, y = y, fill = y36)) +
  geom_raster() +
  ggtitle("Predicted") +
  scale_fill_gradientn(colors = terrain.colors(10))

max(dat[, 2])-min(dat[, 2])


# you can now write the scaled predicted layers to CSV
write.csv(preds, "results/predicted_inla_scaled.csv")


preds.s <- preds

# unscale values of MAP before saving and extracting for species
for (i in 1:1000){
  
  preds.s[, i] <- preds.s[, i]*sd(values(mat.a.unscaled), na.rm = TRUE) + mean(values(mat.a.unscaled), na.rm = TRUE)
  
}

# pred

# you can now write the unscaled predicted layers to CSV
write.csv(preds.s, "results/predicted_inla_unscaled.csv")

head(preds[, 1:10])

# check 
# you can use this method to make the predicted maps
pred.1 <- mat.a
values(pred.1)[!is.na(values(pred.1))] <- preds[, 1]
plot(pred.1) # first predicted layer
# compare
plot(mat.a) # observed layer
# we'll formalize this below

#--------------------------------------------------------------------------------------------------#

# create raster layers from simulated data

#--------------------------------------------------------------------------------------------------#
sims <- list()
sims_grid <- list()
sims_r <-list()
#i=1

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

# example
plot(sims_r[[1]])

#check that details are the same (mins and maxes, etc)
mat.a
sims_r[[1]]


#-------------------------------------------------------------------------------------------------#


# extract data for all pel species
pel.preds.inla <- as.data.frame(matrix(nrow = length(pel.sub$TAXNAME), ncol = 1000))
colnames(pel.preds.inla) <- paste("rep", 1:1000, sep = '.')

for (i in 1:1000) {
  
  pel.preds.inla[, i] <- raster::extract(sims_r[[i]], pel.sub[, c("LONGITUDE", "LATITUDE")])
  
  print(i)
}

pel.preds.inla$TAXNAME <- pel.sub$TAXNAME
head(pel.preds.inla[, 990:1001])


# write extracted null data to csv
write.csv(pel.preds.inla, "results/pel_pred_inla_scaled.csv")
# this will be used in observed vs simulated NB R code


#--------------------------------------------------------------------------------------------------#

# generate completely randomized climate layers 

#--------------------------------------------------------------------------------------------------#
# create completely randomized layers of climate vars and save them as a csv. 

# 6858 cells in temp layer
# 6936 cells in rainfall layer

n.cell <- length(na.omit(values(mat.a)))


# first randomize data in 1000 maps
mat.a.r <- mat.a
random_maps <- list()

for (i in 1:1000){
  
  random_maps[[i]] <- mat.a.r
  values(random_maps[[i]])[!is.na(values(random_maps[[i]]))] <- 
    values(random_maps[[i]])[!is.na(values(random_maps[[i]]))][sample(n.cell, n.cell)] # change numbers!!
}

# plot example
par(mfrow = c(2, 2))
plot(random_maps[[1]], main = "random")
plot(mat.a, main = "observed")
hist(random_maps[[1]], main = "random") # histograms should match
hist(mat.a, main = "observed") # histograms should match

# extract data for all pel species
pel.preds.rand <- as.data.frame(matrix(nrow = length(pel.sub$TAXNAME), ncol = 1000))
colnames(pel.preds.rand) <- paste("rep", 1:1000, sep = '.')

for (i in 1:1000) {
  
  pel.preds.rand[, i] <- raster::extract(random_maps[[i]], pel.sub[, c("LONGITUDE", "LATITUDE")])
  
  print(i)
}

pel.preds.rand$TAXNAME <- pel.sub$TAXNAME
write.csv(pel.preds.rand, "results/pel_pred_rand_scaled.csv")
# this will be used in observed vs simulated NB R code

