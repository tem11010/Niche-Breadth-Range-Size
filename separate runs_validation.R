#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.1)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.1 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.1]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.2)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.2 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.2]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.3)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.3 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.3]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.4)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.4 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.4]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.5)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.5 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.5]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.6)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.6 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.6]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.7)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.7 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.7]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.8)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.8 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.8]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.9)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.9 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.9]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# convert raster to matrix and then use create_grid function to make it into a grid
dat <- rasterToPoints(test.50.10)
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
INLA:::inla.display.matrix(mat.grid$xy_grid2) # observed
#dev.new()
INLA:::inla.display.matrix(INLA:::inla.vector2matrix(result$summary.linear.predictor$mean,nrow,ncol)) # posterior mean


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

# now plot the relationship between observed and predicted values
par(mfrow=c(1, 1))
#dev.new()
dev.off()
plot(inla.matrix2vector(pred_mat),  inla.matrix2vector(obs_mat)) 
#ylim =c(8, max(dat[, 3])), xlim =c(8, max(dat[, 3])))
abline(b=1, a =0)
#-----------#
#-----------#
# I want to plot predicted values for cells I left out againt values of cells predicted by INLA model


# head(test.50.df) # this has where the holes are
# head(pred_df) # this is all the predicted data
# dim(pred_df)


pred_df.raster <-rasterize(pred_df[, 2:1], field = pred_df[, 3], mat.a)
plot(pred_df.raster)
pred_df.raster <- mask(pred_df.raster, mat.a)
plot(pred_df.raster)

pred.50.10 <- values(pred_df.raster)[!is.na(values(pred_df.raster))][subset.50.10]
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

plot(obs.50.1~pred.50.1)
abline(a = 0, b = 1)
plot(obs.50.2~pred.50.2)
abline(a = 0, b = 1)
plot(obs.50.3~pred.50.3)
abline(a = 0, b = 1)
plot(obs.50.4~pred.50.4)
abline(a = 0, b = 1)
plot(obs.50.5~pred.50.5)
abline(a = 0, b = 1)
plot(obs.50.6~pred.50.6)
abline(a = 0, b = 1)
plot(obs.50.7~pred.50.7)
abline(a = 0, b = 1)
plot(obs.50.8~pred.50.8)
abline(a = 0, b = 1)
plot(obs.50.9~pred.50.9)
abline(a = 0, b = 1)
plot(obs.50.10~pred.50.10)
abline(a = 0, b = 1)


write.csv(obs.50.1, "obs50_1_P.csv", row.names = FALSE)
write.csv(obs.50.2, "obs50_2_P.csv", row.names = FALSE)
write.csv(obs.50.3, "obs50_3_P.csv", row.names = FALSE)
write.csv(obs.50.4, "obs50_4_P.csv", row.names = FALSE)
write.csv(obs.50.5, "obs50_5_P.csv", row.names = FALSE)
write.csv(obs.50.6, "obs50_6_P.csv", row.names = FALSE)
write.csv(obs.50.7, "obs50_7_P.csv", row.names = FALSE)
write.csv(obs.50.8, "obs50_8_P.csv", row.names = FALSE)
write.csv(obs.50.9, "obs50_9_P.csv", row.names = FALSE)
write.csv(obs.50.10, "obs50_10_P.csv", row.names = FALSE)


write.csv(pred.50.1, "pred50_1_P.csv", row.names = FALSE)
write.csv(pred.50.2, "pred50_2_P.csv", row.names = FALSE)
write.csv(pred.50.3, "pred50_3_P.csv", row.names = FALSE)
write.csv(pred.50.4, "pred50_4_P.csv", row.names = FALSE)
write.csv(pred.50.5, "pred50_5_P.csv", row.names = FALSE)
write.csv(pred.50.6, "pred50_6_P.csv", row.names = FALSE)
write.csv(pred.50.7, "pred50_7_P.csv", row.names = FALSE)
write.csv(pred.50.8, "pred50_8_P.csv", row.names = FALSE)
write.csv(pred.50.9, "pred50_9_P.csv", row.names = FALSE)
write.csv(pred.50.10, "pred50_10_P.csv", row.names = FALSE)


writeRaster(test.50.1,"MAPvalidation_50_1.asc",  format = "ascii", overwrite = TRUE)
writeRaster(test.50.2,"MAPvalidation_50_2.asc",  format = "ascii", overwrite = TRUE)
writeRaster(test.50.3,"MAPvalidation_50_3.asc",  format = "ascii", overwrite = TRUE)
writeRaster(test.50.4,"MAPvalidation_50_4.asc",  format = "ascii", overwrite = TRUE)
writeRaster(test.50.5,"MAPvalidation_50_5.asc",  format = "ascii", overwrite = TRUE)
writeRaster(test.50.6,"MAPvalidation_50_6.asc",  format = "ascii", overwrite = TRUE)
writeRaster(test.50.7,"MAPvalidation_50_7.asc",  format = "ascii", overwrite = TRUE)
writeRaster(test.50.8,"MAPvalidation_50_8.asc",  format = "ascii", overwrite = TRUE)
writeRaster(test.50.9,"MAPvalidation_50_9.asc",  format = "ascii", overwrite = TRUE)
writeRaster(test.50.10,"MAPvalidation_50_10.asc",  format = "ascii", overwrite = TRUE)

summary(lm(obs.50.10~pred.50.10))







