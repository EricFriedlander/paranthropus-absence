#####################################
#                                   #
#     Paranthropus EM analyses      #
#                                   #
#####################################

# Author: Andrew Du & Eric Friedlander
# Date: 4-13-20



################### DEFINE FUNCTIONS ######################

# function for the standard reflected beta-binomial pdf (eq. S6)
  # ARGUMENTS:
    # X: number of successes (# Paranthropus specimens)
    # n: number of trials (# total mammalIAN specimens)
    # lambda: shape parameter
stdReflectBBpdf <- function(X, n, lambda){
  
  require(rmutil)
  
  m <- 1 / (2 - lambda)
  s <- 2 - lambda
  
  # where m*s is the first shape parameter of the beta distribution, and s * (1-m) is the second shape parameter
  
  return(dbetabinom(X, n, m, s))
}

# function for posterior probability that Z_i = 1, i.e., tau (eq. S10)
  # ARGUMENTS:
    # X: number of successes (# Paranthropus specimens)
    # n: number of trials (# total mammalian specimens)
    # p: parameter (unconditional probability that Z_i = 1)
    # lambda: standard reflected beta distribution shape parameter
tau <- function(X, n, p, lambda){
  
  dbetabinom <- stdReflectBBpdf(X, n, lambda)
  
  return((dbetabinom * p) / ((X == 0) * (1 - p) + dbetabinom * p))
}

# function for expected log-likelihood (eq. S11)
  # ARGUMENTS:
    # same as in tau() function
Q <- function(X, n, p, lambda){
  
  tau_res <- tau(X, n, p, lambda)
  
  dbetabinom <- stdReflectBBpdf(X, n, lambda)
  
  return(sum((1 - tau_res) * log(1 - p) + tau_res * (log(dbetabinom) + log(p))))
}

# function for estimating lambda in the expected log-likelihood (eq. S11) (to be numerically solved)
  # ARGUMENTS:
    # X: number of successes (# Paranthropus specimens)
    # n: number of trials (# total mammalian specimens)
    # p: parameter (unconditional probability that Z_i = 1)
    # lambda: standard reflected beta distribution shape parameter
    # tau_res: results from using the tau() function
lambda_Mstep <- function(X, n, p, lambda, tau_res){
  
  dbetabinom <- stdReflectBBpdf(X, n, lambda)
  
  return(sum(tau_res * log(dbetabinom)))
}

# primary function for the EM algorithm
  # ARGUMENTS:
    # X: number of successes (# Paranthropus specimens)
    # n: number of trials (# total mammalian specimens)
    # p.init: initial guess for p
    # lambda.init: initial guess for lambda
    # n.step.max: number of maximum steps for the optimazation process. Only used for our simulations.
EM <- function(X, n, p.init, lambda.init, n.step.max = NULL){
  
  # step 1: calculate expected log-likelihood given initial guesses of p and lambda
  Q.res <- Q(X, n, p.init, lambda.init)
  
  p.res <- p.new <- p.init
  lambda.res <- lambda.new <- lambda.init
  
  delta.param <- 1 # place-holder for changes in parameter estimates after each iteration
  
  while(delta.param > 1e-5){
    
    # step 2: calculate tau using new p and lambda (E-step)
    tau.res <- tau(X, n, p.new, lambda.new)
    
    # step 3: estimate new values of p and lambda given new tau (M-step)
    p.new <- mean(tau.res)
    
    lambda.opt <- optim(par = lambda.new, fn = lambda_Mstep, gr = NULL, X = X, n = n, p = p.new, tau_res = tau.res, method = "L-BFGS-B", control = list(fnscale = -1), lower = -Inf, upper = 0)
    
    lambda.new <- lambda.opt$par[1]
    
    # Calculate revised log-likelihood
    Q.res <- c(Q.res, Q(X, n, p.new, lambda.new))
    
    # save new parameter estimates
    p.res <- c(p.res, p.new)
    lambda.res <- c(lambda.res, lambda.new)
    
    # Step 4: calculate how much parameters changed
    delta.p <- p.res[length(p.res)] - p.res[length(p.res) - 1]
    delta.lambda <- lambda.res[length(lambda.res)] - lambda.res[length(lambda.res) - 1]
    
    # Step 5: stop if difference between old and new parameters is <= 1e-5.
    delta.param <- max(abs(delta.p), abs(delta.lambda))
    
    # stop while() loop if n.step.max is specified and exceeded
    if(!is.null(n.step.max)) if(length(Q.res) == n.step.max) break
  }
  
  return(list(p.res = p.res, p_hat = p.res[length(p.res)], lambda.res = lambda.res, lambda_hat = lambda.res[length(lambda.res)], Q.res = Q.res))
}

# function for simulating number of Paranthropus specimens data with known p and lambda
  # ARGUMENTS:
    # n_sites: number of sites
    # n_mammSpec: vector of number of total mammalian specimens from each site
    # p: p parameter (probability that a site truly has Paranthropus)
    # lambda: shape parameter for beta-binomial pdf
simulateData <- function(n_sites, n_mammSpec, p, lambda){
  
  require(rmutil) # for simulating random draws from a beta-binomial distribution
  
  Z <- rbinom(n_sites, 1, p) # Z: whether a site truly has Paranthropus or not
  
  Z_equals_1 <- Z == 1 # index out which sites truly have Paranthropus
  
  X <- Z # save to a new object, X, which is the number of recovered Paranthropus specimens at site i
  
  if(sum(Z_equals_1) > 0) X[Z_equals_1] <- rbetabinom(n = sum(Z_equals_1), size = n_mammSpec[Z_equals_1], m = 1 / (2 - lambda), s = 2 - lambda) # for those sites where Z_i = 1, use a beta-binomial to randomly select a Paranthropus abundance value
  
  res <- data.frame(n_ParanSpec = X, n_mammSpec = n_mammSpec, Z_i = Z) # collate the Paranthropus abundance, total mammal abundance, and Z vectors
  
  return(res)  
}

###########################################################





###################### MAIN ANALYSIS ######################

# Read in data
d <- read.csv("NISP data.csv", header = TRUE)

# Define objects
X <- d$Paran.nisp # Paranthropus abundance
n <- X + d$Non.paran.mamm.nisp # Total large mammalian abundance

# Analyze data using expectation-maximization algorithm
EM.res <- EM(X = X, n = n, p.init = 0.5, lambda.init = -5)

# Get out estimated parameters
p_hat <- EM.res$p_hat
lambda_hat <- EM.res$lambda_hat

# Calculate posterior probabilities (prob. Paranthropus is truly absent given estimated parameters and data)
tau_hat <- 1 - tau(X, n, p_hat, lambda_hat)

tau.df <- data.frame(site = paste(d$Site.Formation, d$Member, sep = "_"), tau_hat)

tau.df <- tau.df[tau.df$tau_hat > 0, ] # remove sites where Paranthropus is present


# Calculated expected frequency of Paranthropus specimens across sites using estimated model parameters
BBpdf_expect <- lapply(n, function(n1) stdReflectBBpdf(seq(0, n1), n1, lambda_hat)) # calculate probabilities for sampling number of Paranthropus specimens from zero all the way to total number of mammalian specimens at site. Do this using the beta-binomial PMF.

BBpdf_expect <- lapply(BBpdf_expect, function(x) c(x, rep(0, (max(n) + 1) - length(x)))) # add zeros to the end of each vector for NISP that exceeds the total found at each site

BBpdf_expect_sum <- Reduce("+", BBpdf_expect) # sum PMF vectors across sites

X_expect <- BBpdf_expect_sum * p_hat # multiply summed probabilities by probability that Paranthropus is at site (p_hat)

X_expect[1] <- X_expect[1] + (1 - p_hat) * length(BBpdf_expect) # add (1 - p_hat) times number of sites to the expected number of sites with zero Paranthropus. (1 - p_hat) is the probability Paranthropus is not at the site.

###########################################################





######################## FIGURES ##########################

### Fig. 2: Histogram of number of Paranthropus and large mammal specimens
par(mfrow = c(3, 1), mar = c(4, 4.5, 3, 2) + 0.1)

hist(X, col = "gray", xlab = expression(italic(Paranthropus) ~ "NISP"), ylab = "Number of sites", main = "", cex.axis = 1.75, cex.lab = 2, breaks = 10)

mtext("A", at = 0, cex = 2.25)

hist(n, col = "gray", xlab = "Total large mammalian NISP", ylab = "Number of sites", main = "", cex.axis = 1.75, cex.lab = 2, breaks = 10)

mtext("B", at = 1, cex = 2.25)

hist(X / n, col = "gray", xlab = expression(italic(Paranthropus) ~ "NISP / mammalian NISP"), ylab = "Number of sites", main = "", cex.axis = 1.75, cex.lab = 2, breaks = 10)

mtext("C", at = 0, cex = 2.25)



### Fig. 3: Probability density function of true Paranthropus relative abundance (lambda parameter)
x <- seq(0, 1, length.out = 100000)


par(mfrow = c(2, 1), mar = c(4, 4.5, 2, 2) + 0.1)

plot(x, dbeta(x, 1, 1 - lambda_hat), type = "l", lwd = 3, xlab = expression("True " * italic(Paranthropus) * " relative abundance (" * italic(pi[i]) * ")"), ylab = "Density", cex.axis = 1.5, cex.lab = 1.5, log = "x", xaxt = "n")

axis(1, at = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0), labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1", "1"), cex.axis = 1.5)

text(0.04, 130, bquote(hat(lambda) == .(round(lambda_hat))), cex = 2, pos = 4)

mtext("A", at = 0.00001, cex = 2.5)

### Expected vs. observed number of Paranthropus specimens across sites
plot(table(X), xlab = expression(italic(Paranthropus) ~ "NISP"), ylab = "Number of sites", cex.lab = 1.5, cex.axis = 1.5, lwd = 3, xaxt = "n", xlim = c(0, 34), type = "n")

points(seq(0, 34), X_expect[seq(1, 35)], pch = 21, bg = "gray85", cex = 1.5)
points(table(X), lwd = 3)

axis(1, at = seq(0, 34, 2), labels = seq(0, 34, 2), cex.axis = 1.5)

legend("topright", legend = c("Observed", "Expected"), lty = c(1, NA), pch = c(NA, 21), pt.bg = c(NA, "gray85"), pt.cex = c(NA, 2), pt.lwd = c(NA, 1), bty = "n", cex = 1.5, lwd = 3)

mtext("B", at = 0, cex = 2.5)


### Fig. 4: Probability Paranthropus absence curve
n1 <- seq(0, 10000)

tau0 <- tau(X = rep(0, length(n1)), n = n1, p = p_hat, lambda = lambda_hat) # get taus for sites where there's no Paranthropus and with increasing number of mammalian specimens

# get out NISP corresponding to probabilities of 0.5, 0.75, 0.9, 0.95, and 0.99
nisp_0.5 <- which(abs((1 - tau0) - 0.5) == min(abs((1 - tau0) - 0.5)))

nisp_0.75 <- which(abs((1 - tau0) - 0.75) == min(abs((1 - tau0) - 0.75)))

nisp_0.9 <- which(abs((1 - tau0) - 0.9) == min(abs((1 - tau0) - 0.9)))

nisp_0.95 <- which(abs((1 - tau0) - 0.95) == min(abs((1 - tau0) - 0.95)))

nisp_0.99 <- which(abs((1 - tau0) - 0.99) == min(abs((1 - tau0) - 0.99)))


par(mar = c(5, 6.5, 4, 2) + 0.1, fig = c(0, 1, 0, 1))

plot(n1, 1 - tau0, xlab = "Number of mammalian specimens (NISP)", ylab = expression(atop("Posterior probability " * italic(Paranthropus), " was truly absent from site")), type = "l", lwd = 3, log = "x", cex.axis = 1.5, cex.lab = 1.5)

abline(v = n1[nisp_0.5], lty = 2, lwd = 1.5)
text(n1[nisp_0.5] - 6, 0.62, paste("0.5 =", n1[nisp_0.5], "NISP"), srt = 90, cex = 1)

abline(v = n1[nisp_0.75], lty = 2, lwd = 1.5)
text(n1[nisp_0.75] - 30, 0.55, paste("0.75 =", n1[nisp_0.75], "NISP"), srt = 90, cex = 1, pos = 3)

abline(v = n1[nisp_0.9], lty = 2, lwd = 1.5)
text(n1[nisp_0.9] - 100, 0.55, paste("0.9 =", n1[nisp_0.9], "NISP"), srt = 90, cex = 1, pos = 3)

abline(v = n1[nisp_0.95], lty = 2, lwd = 1.5)
text(n1[nisp_0.95] - 200, 0.55, paste("0.95 =", n1[nisp_0.95], "NISP"), srt = 90, cex = 1, pos = 3)

abline(v = n1[nisp_0.99], lty = 2, lwd = 1.5)
text(n1[nisp_0.99] - 700, 0.55, paste("0.99 =", n1[nisp_0.99], "NISP"), srt = 90, cex = 1, pos = 3)

# plot the inset (curve on arithmetic axes)
par(fig = c(0.06, 0.5, 0.4, 0.975), new = TRUE) 

plot(n1, 1 - tau0, xlab = "", ylab = "", type = "l", lwd = 2, cex.axis = 0.75)

###########################################################





################ SUPPLEMENTARY MATERIALS ##################

### Simulations ###

# set simulated parameter values
p.sim <- seq(0.1, 0.9, 0.2)
lambda.sim <- c(-200, -100, -50, -25, -10, -5, -3, -1)

# number of iterations for our simulations
n.iter <- 1000 

# need to do parallel computing for the simulations

# load packages
library(doParallel)
library(parallel)

# get number of cores on computer
numCores <- detectCores()

# make the cluster using all cores
cl <- makeCluster(numCores)

# register the parallel backend
registerDoParallel(cl)

# set seed to make simulations replicable
set.seed(100)

sim.res <- foreach(p.iter = p.sim) %:% # iterate through p
  
  foreach(lambda.iter = lambda.sim) %:% # iterate through lambda
  
  foreach(i = seq_len(n.iter), .combine = "rbind") %dopar% { # run through each iteration
    
    X.sim <- simulateData(n_sites = nrow(d), n_mammSpec = n, p = p.iter, lambda = lambda.iter) # simulate Paranthropus abundance data for each of our 51 sites, using the observed number of mammalian specimens at each site (n_i) and the known values of p and lambda
    
    EM.iter <- EM(X.sim$n_ParanSpec, n, p.init = 0.5, lambda.init = -10, n.step.max = 5000) # use EM algorithm to estimate p and lambda, using a maximum of 5000 steps for the optimization algorithm
    
    return(c(EM.iter$p_hat, EM.iter$lambda_hat)) # return vector of estimated parameters
  }

# turn off cluster
stopCluster(cl)


# estimate parameter bias using either the mean or median of estimated parameters from simulations
p.bias.mean <- lambda.bias.mean <- p.bias.median <- lambda.bias.median <- array(data = NA, dim = c(length(p.sim), length(lambda.sim)), dimnames = list(p.sim, lambda.sim))

for(p.index in seq_along(p.sim)){
  
  for(lambda.index in seq_along(lambda.sim)){
    
    sim.p <- sim.res[[p.index]][[lambda.index]][, 1]
    
    p.med <- median(sim.p)
    p.mean <- mean(sim.p)
    
    sim.lambda <- sim.res[[p.index]][[lambda.index]][, 2]
    
    lambda.med <- median(sim.lambda)
    lambda.mean <- mean(sim.lambda)
    
    p.bias.median[p.index, lambda.index] <- (p.med - p.sim[p.index]) / p.sim[p.index]
    p.bias.mean[p.index, lambda.index] <- (p.mean - p.sim[p.index]) / p.sim[p.index]
    
    # the absolute value is needed because lambda is always negative
    lambda.bias.median[p.index, lambda.index] <- (lambda.med - lambda.sim[lambda.index]) / abs(lambda.sim[lambda.index]) 
    lambda.bias.mean[p.index, lambda.index] <- (lambda.mean - lambda.sim[lambda.index]) / abs(lambda.sim[lambda.index])
  }
}



### Assessing data independence assumptions ###

# define objects
rel.abund <- X / n

# calculate distances between sites
RA.dist <- dist(rel.abund) # relative abundance

PA <- as.numeric(X != 0) # presence = 1, absence = 0
PA.dist <- dist(PA) # presence-absence distance matrix
# invert zeros and ones so 1 indicates a match, 0 otherwise
PA.dist1 <- PA.dist 
PA.dist1[PA.dist == 0] <- 1
PA.dist1[PA.dist == 1] <- 0

time.dist <- dist(d$Mean.Age_Ma) # time

spatial.dist <- dist(d[, c("Latitude", "Longitude")]) # space

# fit logistic regression to presence-absence data as a function of temporal & spatial distance
glm.PA <- glm(PA.dist1 ~ scale(c(time.dist)) * scale(c(spatial.dist)), family = "binomial")

# fit OLS to relative abundance data as a function of temporal & spatial distance
lm.rel.abund <- lm(RA.dist ~ scale(c(time.dist)) * scale(c(spatial.dist)))

# put coefficient estimates into a dataframe
lm.df <- data.frame(glm.PA$coefficients, lm.rel.abund$coefficients)

rownames(lm.df) <- c("Intercept", "scale(temporal distance)", "scale(spatial distance)", "Interaction term")
colnames(lm.df) <- c("Presence-absence logistic regression", "Relative abundance OLS")



### Figures ###

# Fig. S1
x <- seq(0, 1, length.out = 100)

lambda <- seq(-10, 0, 1)

plot(1, type = "n", xlim = range(x), ylim = c(0, 11), xlab = expression(italic(pi[i])), ylab = "Density")

library(RColorBrewer)

colors <- brewer.pal(length(lambda), "Spectral")

for(i in seq_along(lambda)) lines(x, dbeta(x, 1, 1 - lambda[i]), col = colors[i])

legend("topright", legend = lambda, title = expression(lambda), lty = rep(1, length(lambda)), col = colors, cex = 0.85)


# Fig. S3
x <- seq(0, 1, length.out = 1000)

plot(1, type = "n", xlim = range(x), ylim = c(0, 11), xlab = expression(italic(pi[i])), ylab = "Density")

colors <- brewer.pal(length(lambda.sim), "Spectral")

for(i in seq_along(lambda.sim)) lines(x, dbeta(x, 1, 1 - lambda.sim[i]), col = colors[i])

legend("topright", legend = lambda.sim, title = expression(lambda), lty = rep(1, length(lambda.sim)), col = colors, cex = 0.85)


# Fig. S4

# load packages
library(ggplot2)
library(reshape2)

# reshape data
melt_p.bias.mean <- melt(p.bias.mean)
melt_lambda.bias.mean <- melt(lambda.bias.mean)

melt_p.bias.median <- melt(p.bias.median)
melt_lambda.bias.median <- melt(lambda.bias.median)

# combine into one dataframe with parameter name as a column for ggplot faceting
melt_res <- data.frame(rbind(melt_p.bias.mean, 
                             melt_lambda.bias.mean, 
                             melt_p.bias.median, 
                             melt_lambda.bias.median), 
                       param = c(
                         rep("p.mean", nrow(melt_p.bias.mean)), 
                         rep("lambda.mean", nrow(melt_lambda.bias.mean)),
                         rep("p.med", nrow(melt_p.bias.median)),
                         rep("lambda.med", nrow(melt_lambda.bias.median))))

# reorder factors so p is plotted on top of lambda
melt_res$param <- ordered(melt_res$param, levels = c("p.mean", "lambda.mean", "p.med", "lambda.med"))

# have to do this, so p is italicized in plot and hats are put on parameters
levels(melt_res$param) <- c("p.mean" = expression(hat(italic(p)) ~ "(mean)"), 
                            "lambda.mean" = expression(hat(lambda) ~ "(mean)"),
                            "p.med" = expression(hat(italic(p)) ~ "(median)"),
                            "lambda.med" = expression(hat(lambda) ~ "(median)"))

# create plot
ggplot(data = melt_res, aes(x = as.factor(Var1), y = as.factor(Var2), fill = value)) + 
  facet_wrap(~param, nrow = 2, labeller = label_parsed) + 
  geom_tile(color = "black") + 
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0, limit = c(-0.2, 0.2), space = "Lab", 
                       name = "Relative bias") +
  theme_minimal() + 
  labs(x = expression("True" ~ italic(p)), y = expression("True" ~ lambda)) + 
  theme(axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 13),
        legend.text.align = 1)


# Fig. S5
unstable.lambda <- sim.res[[1]][[1]][, 2] < -1e6 # index out unstable lambdas

par(mar = c(5, 4.5, 4, 2) + 0.1)

plot(sim.res[[1]][[1]][!unstable.lambda, 1], sim.res[[1]][[1]][!unstable.lambda, 2], xlab = expression(hat(italic(p))), ylab = expression(hat(lambda)), main = expression("True" ~ italic(p) ~ "= 0.1; True" ~ lambda ~ "= -200"))


# Fig. S6
neg.lamb <- sapply(sim.res[[1]], function(x) sum(x[, 2] < -1e6)) # count those iterations where estimated lambda is extremely negative

par(mar = c(5, 4.5, 4, 2) + 0.1)

plot(lambda.sim, neg.lamb, pch = 16, xlab = expression("True" ~ lambda), ylab = expression("Number of" ~ hat(lambda) ~ "that are extremely negative"))


# Fig. S7
hist(scale(spatial.dist), xlab = "scale(spatial distance)", ylab = "Numbers of pairs of sites", col = "gray", main = "")


# Fig. S8
par(mfrow = c(1, 2), mar = c(5, 4, 4, 0) + 0.1)

plot(time.dist, RA.dist, main = "", xlab = "Temporal distance (millions of years)", ylab = expression(italic(Paranthropus) ~ "relative abundance difference"), cex.axis = 0.8, cex.lab = 0.8)

mtext("A", at = 0, cex = 1.5)

cor.time <- round(cor(time.dist, RA.dist, method = "spearman"), 2)

text(1.25, 0.05, bquote(rho == .(cor.time)), pos = 4, cex = 0.8)


plot(spatial.dist, RA.dist, main = "", xlab = "Spatial distance\n(Euclidean distance using lat/long)", ylab = "", cex.axis = 0.8, cex.lab = 0.8)

mtext("B", at = 0, cex = 1.5)

cor.space <- round(cor(spatial.dist, RA.dist, method = "spearman"), 2)

text(15, 0.05, bquote(rho == .(cor.space)), pos = 4, cex = 0.8)