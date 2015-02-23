### R script associated the manuscipt entitled "Habitat features and long-distance dispersal modify the use of
### social information by a long-distance migratory bird"
### Primary author: Clark Rushing
### Contact: crushing05@gmail.com 

rm(list=ls())

devtools::install_github('crushing05/crushingr')
require(crushingr)


###############################################################################################################
#############################    Analysis of settlement data from playback experiment    ######################
###############################################################################################################

### Read playback data ----
# read table
  pb <- read.csv("pb_complete.csv")
  pb$year <- pb$year-2012
  pb$trt <- factor(pb$trt,levels(pb$trt)[c(2,3,1)])
  pb$occ <- ifelse(pb$abun > 0, 1, 0) # Convert abundance to occupancy
  pb$abun.asy <- pb$abun - pb$abun.sy # ASY abundance
  pb$occ.sy <- ifelse(pb$abun > 0 & pb$abun.sy > 0, 1, 0) # SY occupancy
  pb$occ.asy <- ifelse(pb$abun.asy > 0 , 1, 0) # ASY occupancy
  pb$trt1 <- as.numeric(factor(pb$trt, levels=c("control", "ca", "pi")))

### SY settlement ----
## Full model
  full.sy <- glm(occ.sy ~ trt + PC1*trt + PC2*trt + year*trt, data = pb, family = "binomial")
  summary(full.sy)

## First test for interactions; if not signficant, remove to test for main effects
  fit1.sy <- glm(occ.sy ~ trt + PC1*trt + PC2 + year*trt, data = pb, family = "binomial") # Test for PC2*trt
  anova(fit1.sy, full.sy, test="Chisq") # NS

  fit2.sy <- glm(occ.sy ~ trt + PC1 + PC2*trt + year*trt, data = pb, family = "binomial") # Test for PC1*trt 
  anova(fit2.sy, full.sy, test="Chisq") # NS 

  fit3.sy <- glm(occ.sy ~ trt + PC1*trt + PC2*trt + year, data = pb, family = "binomial") # Test for yr*trt
  anova(fit3.sy, full.sy, test="Chisq") # NS 

## No evidence of interactions, so new 'full' model
  full2.sy <- glm(occ.sy ~ trt + PC1 + PC2 + year, data = pb, family = "binomial") # No interactions
  anova(full2.sy, full.sy, test="Chisq") # NS 
  summary(full2.sy)
  
  fit4.sy <- glm(occ.sy ~ trt + PC1 + PC2, data = pb, family = "binomial") # Test for year effect
  anova(fit4.sy, full2.sy, test= "Chisq") # NS
  
  fit5.sy <- glm(occ.sy ~ trt + PC2 + year, data = pb, family = "binomial") # Test for PC1 effect
  anova(fit5.sy, full2.sy, test= "Chisq") # NS
  
  fit6.sy <- glm(occ.sy ~ trt + PC1 + year, data = pb, family = "binomial") # Test for PC2 effect
  anova(fit6.sy, full2.sy, test= "Chisq") # p = 0.015
  
  fit7.sy <- glm(occ.sy ~ PC1 + PC2 + year, data = pb, family = "binomial") # Test for treatment effect
  anova(fit7.sy, full2.sy, test= "Chisq") # p = 0.0033
  
  fit8.sy <- glm(occ.sy ~ trt + PC2, data = pb, family = "binomial") # Top model
  summary(fit8.sy)


### ASY settlement ----
# Full model
  full.asy <- glm(occ.asy ~ trt + PC1*trt+PC2*trt + year*trt, data = pb, family = "binomial")

## Test for interactions
  fit2.asy <- glm(occ.asy ~ trt + PC1*trt + PC2 + year*trt, data = pb, family = "binomial") # Test for PC2*trt
  anova(fit2.asy, full.asy, test="Chisq") # NS                                                          
  
  fit3.asy <- glm(occ.asy ~ trt + PC1 + PC2*trt + year*trt, data = pb, family = "binomial") # Test for PC1*trt
  anova(fit3.asy, full.asy, test="Chisq") # NS                                                          
  
  fit4.asy <- glm(occ.asy ~ trt + PC1*trt + PC2*trt + year, data = pb, family = "binomial") # Test for yr*trt
  anova(fit4.asy, full.asy, test="Chisq") # NS                                                          

## Remove non-significant interactions, new 'full' model
  full2.asy <- glm(occ.asy ~ trt + PC1 + PC2 + year, data = pb, family = "binomial") #
  anova(full2.asy, full.asy, test="Chisq") # NS 
  summary(full2.asy)
  
  fit5.asy <- glm(occ.asy ~ trt + PC1 + PC2, data = pb, family = "binomial") # Test for year effect
  anova(fit5.asy, full2.asy, test= "Chisq") # NS
  
  fit6.asy <- glm(occ.asy ~ trt + PC2 +  year, data = pb, family = "binomial") # Test for PC1 effect
  anova(fit6.asy, full2.asy, test= "Chisq") # 0.0061
  
  fit7.asy <- glm(occ.asy ~ trt + PC1 + year, data = pb, family = "binomial") # Test for PC2 effect
  anova(fit7.asy, full2.asy, test= "Chisq") # p = 0.086
  
  fit8.asy <- glm(occ.asy ~ PC1 + PC2 + year, data = pb, family = "binomial") # Test for trt effect
  anova(fit8.asy, full2.asy, test= "Chisq") # p < 0.001
  
  fit9.asy <- glm(occ.asy ~ trt + PC1*trt, data = pb, family = "binomial") # Top model
  summary(fit9.asy)


### ASY Outlier test occupancy ----
  pb2 <- pb[pb$PC1<4,]
# Full model
  full.asy <- glm(occ.asy ~ trt + PC1*trt+PC2*trt + year*trt, data = pb2, family = "binomial")

## LRT for interactions
  fit2.asy <- glm(occ.asy ~ trt + PC1*trt + PC2 + year*trt, data = pb2, family = "binomial") # Test for PC2*trt
  anova(fit2.asy, full.asy, test="Chisq") # NS                                                          
  
  fit3.asy <- glm(occ.asy ~ trt + PC1 + PC2*trt + year*trt, data = pb2, family = "binomial") # Test for PC1*trt
  anova(fit3.asy, full.asy, test="Chisq") # NS                                                          
  
  fit4.asy <- glm(occ.asy ~ trt + PC1*trt + PC2*trt + year, data = pb2, family = "binomial") # Test for yr*trt
  anova(fit4.asy, full.asy, test="Chisq") # NS                                                          

## Remove non-significant interactions, new 'full' model
  full2.asy <- glm(occ.asy ~ trt + PC1 + PC2 + year, data = pb2, family = "binomial") #
  anova(full2.asy, full.asy, test="Chisq") # NS 
  summary(full2.asy)
  
  fit5.asy <- glm(occ.asy ~ trt + PC1 + PC2, data = pb2, family = "binomial") # Test for year effect
  anova(fit5.asy, full2.asy, test= "Chisq") # NS
  
  fit6.asy <- glm(occ.asy ~ trt + PC2 +  year, data = pb2, family = "binomial") # Test for PC1 effect
  anova(fit6.asy, full2.asy, test= "Chisq") # 0.014
  
  fit7.asy <- glm(occ.asy ~ trt + PC1 + year, data = pb2, family = "binomial") # Test for PC2 effect
  anova(fit7.asy, full2.asy, test= "Chisq") # p = 0.089
  
  fit8.asy <- glm(occ.asy ~ PC1 + PC2 + year, data = pb2, family = "binomial") # Test for trt effect
  anova(fit8.asy, full2.asy, test= "Chisq") # p < 0.001
  
  fit9.asy <- glm(occ.asy ~ trt + PC1, data = pb2, family = "binomial") # Top model
  summary(fit9.asy)

### SY Outlier test occupancy ----
  pb3 <- pb[pb$PC2<4,]
# Full model
  full.sy <- glm(occ.sy ~ trt + PC1*trt+PC2*trt + year*trt, data = pb3, family = "binomial")

## LRT for interactions
  fit2.sy <- glm(occ.sy ~ trt + PC1*trt + PC2 + year*trt, data = pb3, family = "binomial") # Test for PC2*trt
  anova(fit2.sy, full.sy, test="Chisq") # NS                                                          
  
  fit3.sy <- glm(occ.sy ~ trt + PC1 + PC2*trt + year*trt, data = pb3, family = "binomial") # Test for PC1*trt
  anova(fit3.sy, full.sy, test="Chisq") # NS                                                          
  
  fit4.sy <- glm(occ.sy ~ trt + PC1*trt + PC2*trt + year, data = pb3, family = "binomial") # Test for yr*trt
  anova(fit4.sy, full.sy, test="Chisq") # NS                                                          

## Remove non-significant interactions, new 'full' model
  full2.sy <- glm(occ.sy ~ trt + PC1 + PC2 + year, data = pb3, family = "binomial") #
  anova(full2.sy, full.sy, test="Chisq") # NS 
  summary(full2.sy)
  
  fit5.sy <- glm(occ.sy ~ trt + PC1 + PC2, data = pb3, family = "binomial") # Test for year effect
  anova(fit5.sy, full2.sy, test= "Chisq") # NS
  
  fit6.sy <- glm(occ.sy ~ trt + PC2 +  year, data = pb3, family = "binomial") # Test for PC1 effect
  anova(fit6.sy, full2.sy, test= "Chisq") # NS
  
  fit7.sy <- glm(occ.sy ~ trt + PC1 + year, data = pb3, family = "binomial") # Test for PC2 effect
  anova(fit7.sy, full2.sy, test= "Chisq") # p = 0.064
  
  fit8.sy <- glm(occ.sy ~ PC1 + PC2 + year, data = pb3, family = "binomial") # Test for trt effect
  anova(fit8.sy, full2.sy, test= "Chisq") # p = 0.002
  
  fit9.sy <- glm(occ.sy ~ trt + PC2, data = pb3, family = "binomial") # Top model
  summary(fit9.sy)



################################################################################################################
#########################################      Dispersal Analysis      #########################################
################################################################################################################

### Bayesian Binomial model ----
  sink("Disp_binom.txt")
  cat("
      model { # Binomial GLM
      
      # Priors
      for (i in 1:2){                # p.imm[1] = population immigration rate
      p.imm[i] ~ dunif(0,1)        # p.imm[2] = playback immigration rate
      }                            
      # Likelihood
      
      for (i in 1:2) {
      C[i] ~ dbin(p.imm[i], n[i])
      }
      
      # Derived quantities
      p.diff <- p.imm[1] - p.imm[2]   # Difference between population and playback
      }
      ", fill = TRUE)
  sink()

# Read data ----
  disp.df <- read.csv("disp_df.csv")

# JAGS settings
  nc <- 3
  ni <- 25000
  nb <- 10000
  nt <- 2
  
  inits <- function(){list(p.imm=runif(2, .3,.7))}
  parameters<- c("p.imm", "p.diff")

### 4:1 Odds ratio ----
  jags.disp4 <- list(C=disp.df$imm4, n=disp.df$N)
  
  
  jm.disp4 <- jags.model("Disp_binom.txt", data = jags.disp4, inits = inits, 
                         n.chains = nc, n.adapt = nb)
  zm.disp4 <- coda.samples(jm.disp4, variable.names = parameters, n.iter = ni,
                           n.thin=nt)
  summary(zm.disp4)
  #traceplot(zm.disp4)   # Check traceplots for convergence
  #gelman.diag(zm.disp4, transform = TRUE) # Check Gelman-Rubin diagnostics for convergence
  
  df.tot4 <- as.data.frame(zm.disp4[[1]], zm.disp4[[2]], zm.disp4[[3]]) # Save and write MCMC for figure
  write.csv(df.tot4, "disp4_mcmc.csv", row.names = FALSE)

# 4:1 Posterior Predictive Check
  n.sims <- 1500
  y.rep.pop <- array (NA, c(n.sims))
  
  for (sim in 1:n.sims){
    y.rep.pop[sim] <- rbinom (n=1, size = disp.df[1,1], prob = df.tot4[sim,2])}
  par(mfrow=c(1,1))
  hist(y.rep.pop, xlim=c(60,100),xlab="Number of Local Individuals", ylab = "Count", main="")
  abline(v=disp.df[1,2], col="red", lwd=2)
  
  y.rep.pb <- array (NA, c(n.sims))
  for (sim in 1:n.sims){
    y.rep.pb[sim] <- rbinom (n=1, size = disp.df[2,1], prob = (df.tot4[sim,3]))}
  hist(y.rep.pb, xlim=c(0,16),xlab="Number of Local Individuals", ylab = "Count", main="")
  abline(v=(disp.df[2,2]), col="red", lwd=2)

### 9:1 Odds ratio ----
  jags.disp9  <- list(C=disp.df$imm9, n=disp.df$N)
  
  jm.disp9 <- jags.model("Disp_binom.txt", data = jags.disp9, inits = inits, 
                         n.chains = nc, n.adapt = nb)
  zm.disp9 <- coda.samples(jm.disp9, variable.names = parameters, n.iter = ni,
                           n.thin=nt)
  summary(zm.disp9)
  #traceplot(zm.disp9)
  #gelman.diag(zm.disp9, transform = TRUE)
  
  df.tot9 <- as.data.frame(zm.disp9[[1]], zm.disp9[[2]], zm.disp9[[3]]) # Save and write MCMC for figure
  write.csv(df.tot9, "disp9_mcmc.csv", row.names = FALSE)

# 9:1 Posterior Predictive Check
  n.sims <- 1500
  y.rep.pop <- array (NA, c(n.sims))
  
  for (sim in 1:n.sims){
    y.rep.pop[sim] <- rbinom (n=1, size = disp.df[1,1], prob = df.tot9[sim,2])}
  par(mfrow=c(1,1))
  hist(y.rep.pop, xlim=c(60,100),xlab="Number of Local Individuals", ylab = "Count", main="")
  abline(v=disp.df[1,3], col="red", lwd=2)
  
  y.rep.pb <- array (NA, c(n.sims))
  for (sim in 1:n.sims){
    y.rep.pb[sim] <- rbinom (n=1, size = disp.df[2,1], prob = (df.tot9[sim,3]))}
  hist(y.rep.pb, xlim=c(0,16),xlab="Number of Local Individuals", ylab = "Count", main="")
  abline(v=(disp.df[2,3]), col="red", lwd=2)

### 19:1 Odds ratio ----
  jags.disp19 <- list(C=disp.df$imm19, n=disp.df$N)
  
  jm.disp19 <- jags.model("Disp_binom.txt", data = jags.disp19, inits = inits, 
                          n.chains = nc, n.adapt = nb)
  zm.disp19 <- coda.samples(jm.disp19, variable.names = parameters, n.iter = ni,
                            n.thin=nt)
  summary(zm.disp19)
  #traceplot(zm.disp.tot19)
  #gelman.diag(zm.disp.tot19, transform = TRUE)
  
  df.tot19 <- as.data.frame(zm.disp19[[1]], zm.disp19[[2]], zm.disp19[[3]]) # Save and write MCMC for figure
  write.csv(df.tot19, "disp19_mcmc.csv", row.names = FALSE)

# 19:1 Posterior Predictive Check
  n.sims <- 1500
  y.rep.pop <- array (NA, c(n.sims))
  
  for (sim in 1:n.sims){
    y.rep.pop[sim] <- rbinom (n=1, size = disp.df[1,1], prob = df.tot19[sim,2])}
  par(mfrow=c(1,1))
  hist(y.rep.pop, xlim=c(60,100),xlab="Number of Local Individuals", ylab = "Count", main="")
  abline(v=disp.df[1,4], col="red", lwd=2)
  
  y.rep.pb <- array (NA, c(n.sims))
  for (sim in 1:n.sims){
    y.rep.pb[sim] <- rbinom (n=1, size = disp.df[2,1], prob = (df.tot19[sim,3]))}
  hist(y.rep.pb, xlim=c(0,16),xlab="Number of Local Individuals", ylab = "Count", main="")
  abline(v=(disp.df[2,4]), col="red", lwd=2)




##############################################################################################################
##########################################     Plots     #####################################################
##############################################################################################################

#devtools::install_github('crushing05/crushingr)
require(crushingr)


#### Settlement Probability Plot ----
# Fit 'means' parameterization of top models
  fit.sy <- glm(occ.sy ~ trt + PC2 - 1, data = pb, family = "binomial") # Top yearling model
  fit.asy <- glm(occ.asy ~ trt + PC1 - 1, data = pb, family = "binomial") # Top adult model

# Estimate settlement probablities, CI's and store in df
  sy.pred <- invlogit(cbind(coef(fit.sy)[1:3], confint(fit.sy)[1:3,]))
  asy.pred <- invlogit(cbind(coef(fit.asy)[1:3], confint(fit.asy)[1:3,]))
  settlement.pred <- data.frame(rbind(sy.pred, asy.pred))
  names(settlement.pred) <- c("Estimate", "LCI", "UCI")
  settlement.pred$Treatment <- rep(c("Control", "Public information", "Location cues"), 2)
  settlement.pred$Age <- rep(c("Yearling", "Adult"), each = 3)

  rug_breaksx <- data.frame(x= c(0, 1, 2, 1, 3),  # Hack to get tick marks on the inside of figure
                            y=c(0,0.25, 0.5, 0.75, 1.00))

# Figure
  ggplot(data = settlement.pred, aes(x=Treatment, y=Estimate)) + 
    geom_errorbar(aes(ymax=(UCI), ymin=(LCI), group=Age, linetype=Age), 
                  position=position_dodge(width=.5),width=0) +
    geom_point(aes(group=Age, fill=Age),size=6, position=position_dodge(width=.5),shape=21) +
    scale_fill_manual(values = c("black", "white"))+
    scale_linetype_manual(values=c("solid","dashed"))+ 
    scale_y_continuous("Probability of settlement", lim=c(0,1))+
    scale_x_discrete("Treatment", labels=c("Control \n ", "Location \ncues", 
                                           "Public \ninformation"))+
    geom_rug(data=rug_breaksx, aes(x=x, y=y))+
    theme(panel.grid.major = element_line(color="grey95"), 
          panel.grid.major.x = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=28),
          legend.title = element_blank(),
          legend.key.size = unit(1.75, "cm"),
          legend.position=c(.85,.84),
          plot.margin=unit(c(1,1,1,1), "cm"),
          legend.key = element_blank())  


#### Habitat effect plots ----
# Yearling plot
# Estimates of treatment-specific slopes to overlay on figure
  hab.sy <- glm(occ.sy ~ trt*PC2, data = pb, family = "binomial")
  print(coef(hab.sy), digits=1)
  lc.lab.sy <- paste("beta == ", -1.22) # Location cue slope
  pi.lab.sy <- paste("beta == ", -0.73) # Public information slope
  co.lab.sy <- paste("beta == ", -0.11) # Control slope
  rug.sy <- data.frame(x = c(-2, 0, 2, 4, 6), y = c(0, 0.25, 0.5, 0.75, 1))

# Plot
  ggplot(pb, aes(x=PC2, y=occ.sy))+ geom_point(aes(shape=trt))+
    stat_smooth(aes(linetype=trt), method="glm", family="binomial", color="grey75", 
                fullrange=TRUE, se=FALSE, size=0.6) +
    stat_smooth(method="glm",family="binomial", se=FALSE, fullrange=TRUE, color = "black", size=0.8)+
    scale_x_continuous("Habitat PC2")+
    scale_y_continuous("Yearling settlment probability")+
    scale_linetype_manual(values = c("F1", "longdash", "dashed"), 
                          labels = c("Control", "Public information", "Location cues"))+
    scale_shape_manual(values = c(1,19,4),
                       labels = c("Control", "Public information", "Location cues"))+
    guides(shape=guide_legend())+
    annotate("text", x = -0.2, y = 0.75, label = lc.lab.sy, color="grey40", parse=TRUE)+
    annotate("text", x = -1.9, y = 0.25, label = pi.lab.sy, color="grey40", parse=TRUE)+
    annotate("text", x = 4, y = 0.1, label = co.lab.sy, color="grey40", parse=TRUE)+
    theme(panel.grid.major = element_line(color="grey95"), panel.grid.major.x = element_blank(),
          strip.text.x=element_text(size = 18), strip.text.y=element_text(size = 18), legend.position=c(.8,.8),
          legend.text = element_text(size = 18), legend.key = element_blank(), legend.key.size = unit(1.5,"cm"),
          axis.ticks = element_blank())+
    geom_rug(data=rug.sy, aes(x=x, y=y))


# Adult plot
# Estimates of treatment-specific slopes to overlay on figure
  hab.asy <- glm(occ.asy ~ trt*PC1, data = pb, family = "binomial")
  print(coef(hab.asy), digits=2)
  lc.lab.asy <- paste("beta == ", -0.42) # Location cue slope
  co.lab.asy <- paste("beta == ", -0.27) # Control slope (*Public info slope = 0, so not shown*)
  rug.asy <- data.frame(x = c(-2, 0, 2, 4, 6), y = c(1, 0.25, 0.5, 0.75, 1))
  
  ggplot(pb, aes(x=PC1, y=occ.asy))+ geom_point(aes(shape=trt))+
    stat_smooth(aes(linetype=trt), method="glm", family="binomial", color="grey75", 
                fullrange=TRUE, se=FALSE, size=0.6) +
    stat_smooth(method="glm",family="binomial", se=FALSE, fullrange=TRUE, color = "black", size=0.8)+
    scale_x_continuous("Habitat PC1")+
    scale_y_continuous("Adult settlment probability")+
    scale_linetype_manual(values = c("F1", "longdash", "dashed"), 
                          labels = c("Control", "Public information", "Location cues"))+
    scale_shape_manual(values = c(1,19,4),
                       labels = c("Control", "Public information", "Location cues"))+
    guides(shape=guide_legend())+
    annotate("text", x = -0.2, y = 0.75, label = lc.lab.asy, color="grey40", parse=TRUE)+
    annotate("text", x = -2, y = 0.12, label = co.lab.asy, color="grey40", parse=TRUE)+
    theme(panel.grid.major = element_line(color="grey95"), panel.grid.major.x = element_blank(),
          strip.text.x=element_text(size = 18), strip.text.y=element_text(size = 18), legend.position=c(.8,.8),
          legend.text = element_text(size = 18), legend.key = element_blank(), legend.key.size = unit(1.5,"cm"),
          axis.ticks = element_blank())+
    geom_rug(data=rug.asy, aes(x=x, y=y))


#### Plot of posterior distributions from dispersal model ----
  mcmc4 <- read.csv("disp4_mcmc.csv") # Read MCMC output from dispersal model
  mcmc9 <- read.csv("disp9_mcmc.csv")
  mcmc19 <- read.csv("disp19_mcmc.csv")
  
  disp.legend <- data.frame(x=c(0,0,0), y=c(10,10,10), Odds = factor(c("4:1", "9:1", "19:1"),# Hack to add legend
                                                                     levels=c("4:1", "9:1", "19:1")))
  rug.disp <- data.frame(x=c(-.25, 0, .25, .5, .75), y=c(0,2,4,6,8)) # Hack to add tick marks
  
  ggplot(mcmc4, aes(x=p.diff)) + geom_density(adjust=5, size = 1) + 
    geom_segment(aes(x= 0.00694, y= 0, xend = 0.00694, yend= 0.9), # LCI 
                 linetype="dotted", size=0.6) +
    geom_segment(aes(x= 0.462, y= 0, xend = 0.462, yend= 0.7),     # UCI
                 linetype="dotted", size=0.6)+
    geom_density(data= mcmc9, aes(x=p.diff), adjust= 5, size = 1, linetype="longdash") +
    geom_density(data= mcmc19, aes(x=p.diff),  adjust= 5, size = 1, linetype="dashed") + 
    scale_y_continuous("Density", lim=c(0,8)) + scale_x_continuous(lim=c(-.3,.75)) +
    labs(x=bquote(italic(p[diff]))) + 
    geom_errorbar(data=disp.legend, aes(x=-0.3, ymin=y, ymax=y, linetype=Odds), size=1.3) +
    scale_linetype_manual(values=c("solid", "longdash", "dashed")) +
    geom_rug(data=rug.disp, aes(x=x, y=y)) +
    theme(axis.ticks = element_blank(), 
          legend.text = element_text(size=28),
          legend.title = element_blank(),
          legend.key.size = unit(2, "cm"),
          legend.position = c(.65,.75),
          plot.margin = unit(c(1,1,1,1), "cm"),
          legend.key = element_blank())  





