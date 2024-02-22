rm(list=ls())
library(foreign)
library(tidyr)
library(dplyr)
library(ranger)
library(ggplot2)


# Load data ---------------------------------------------------------------

d <- read.spss("STAR_Students.sav", to.data.frame=TRUE)
# These data and their documentation can be found at: 
# https://doi.org/10.7910/DVN/SIWH9F


# Preliminary variable processing -----------------------------------------

d$mob  <- as.numeric(d$birthmonth)
d$race <- as.factor(d$race)
d$yob  <- as.numeric(d$birthyear)
d$flunch <- as.numeric(d$gkfreelunch)
d$specialed <- as.numeric(d$gkspeced)
d$female  <- as.numeric(d$gender=="FEMALE")
d$smallclass <- as.numeric(d$gkclasstype=="SMALL CLASS")
d$aideclass <- as.numeric(d$gkclasstype=="REGULAR + AIDE CLASS")
d$g2ttotalss <- d$g2tmathss + d$g2treadss + d$g2tlistss
d$g1ttotalss <- d$g1tmathss + d$g1treadss + d$g1tlistss
d$gkttotalss <- d$gktmathss + d$gktreadss + d$gktlistss


# Subset to data of interest ----------------------------------------------

vars.keep <- 
  c("stdntid",
    "mob","race","yob","flunch","specialed","female",
    "smallclass","aideclass","g1tmathss","gktmathss",
    "gkttotalss","g1ttotalss","g2ttotalss",
    "gkpresent")
dat <- na.omit(d[,vars.keep])


# Functions for analysis --------------------------------------------------

prep_data <- function(dat, outY2, outY1, prop.perm, this.seed){
  
  #Separate into sets
  set.seed(this.seed)
  dat$set <- sample(c("I","II","III"), nrow(dat), replace = TRUE)
  
  #Designate outcomes
  dat$Y2 <- dat[,outY2]
  dat$Y1 <- dat[,outY1]
  
  #Induce shift
  set.seed(this.seed)
  dat1 <- subset(dat,dat$set == "I")
  dat2 <- subset(dat,dat$set == "II")
  dat3 <- subset(dat,dat$set == "III")
  
  #Permute dataset II outcomes (shift)
  k2 <- sample(1:nrow(dat2), round(prop.perm*nrow(dat2)), replace = F)
  k2perm <- sample(k2)
  dat2[k2,c("Y2","Y1")] <- dat2[k2perm,c("Y2","Y1")]
  rm(k2,k2perm)
  
  #Double permute dataset III outcomes (even more shift)
  k3 <- sample(1:nrow(dat3), round(prop.perm*nrow(dat3)), replace = F)
  k3perm <- sample(k3)
  dat3[k3,c("Y2","Y1")] <- dat3[k3perm,c("Y2","Y1")]
  rm(k3,k3perm)
  k3 <- sample(1:nrow(dat3), round(prop.perm*nrow(dat3)), replace = F)
  k3perm <- sample(k3)
  dat3[k3,c("Y2","Y1")] <- dat3[k3perm,c("Y2","Y1")]
  rm(k3,k3perm)
  
  #Indicate what is observed and what is not observed in practice
  #Set1
  dat1$Y2o <- NA
  dat1$Y1o <- NA
  #Set2
  dat2$Y2o <- NA
  dat2$Y1o <- dat2$Y1
  #Set3
  dat3$Y2o <- dat3$Y2
  dat3$Y1o <- dat3$Y1
  
  return(
    list(dat = dat,
         dat1 = dat1,
         dat2 = dat2,
         dat3 = dat3)
  )
  
}

run_analysis <- function(dat,dat1,dat2,dat3,this.seed){
  
  set.seed(this.seed)
  
  tau <- ranger(Y2 ~ mob + race + yob + flunch + specialed + female + 
                  smallclass + aideclass, dat1)
  taupreds <- tau$predictions
  R2T <- 1 - (
    (sum((dat1$Y2 - taupreds)^2))/(sum((dat1$Y2 - mean(dat1$Y2))^2))
  )
  
  tauinterceptpred <- mean(dat3$Y2o)
  mseintercept <- mean((taupreds - tauinterceptpred)^2)
  
  tauA <- ranger(Y2o ~ mob + race + yob + flunch + specialed + female + 
                   smallclass + aideclass, dat3)
  tauApreds <- predict(tauA, data = dat1)$predictions
  mseA <- mean((taupreds - tauApreds)^2)
  R2A <- 1 - (
    (sum((dat1$Y2 - tauApreds)^2))/(sum((dat1$Y2 - mean(dat1$Y2))^2))
  )
  
  tauB <- ranger(Y1o ~ mob + race + yob + flunch + specialed + female + 
                   smallclass + aideclass, dat2)
  tauBpreds <- predict(tauB, data = dat1)$predictions
  tauBtilde_trans <- lm(Y2o ~ Y1o, dat3)
  tauBtildepreds <- predict(tauBtilde_trans, newdata = data.frame(Y1o = tauBpreds))
  mseB <- mean((taupreds - tauBtildepreds)^2)
  R2B <- 1 - (
    (sum((dat1$Y2 - tauBtildepreds)^2))/(sum((dat1$Y2 - mean(dat1$Y2))^2))
  )
  
  tauC1 <- ranger(Y2o ~ Y1o + mob + race + yob + flunch + specialed + female + 
                    smallclass + aideclass, dat3)
  dat2$Y2o <- predict(tauC1, data = dat2)$predictions
  tauC2 <- ranger(Y2o ~ mob + race + yob + flunch + specialed + female + 
                    smallclass + aideclass, dat2)
  tauCpreds <- predict(tauC2, data = dat1)$predictions
  mseC <- mean((taupreds - tauCpreds)^2)
  R2C <- 1 - (
    (sum((dat1$Y2 - tauCpreds)^2))/(sum((dat1$Y2 - mean(dat1$Y2))^2))
  )
  
  return(data.frame(mseA = mseA,
                    mseB = mseB,
                    mseC = mseC,
                    mseintercept = mseintercept,
                    R2T = R2T,
                    R2A = R2A,
                    R2B = R2B,
                    R2C = R2C))
  
}


# Run ---------------------------------------------------------------------

outY2 <- "g2ttotalss"
outY1vec <- c("g1ttotalss","g1tmathss","gktmathss","gkpresent")
all.res <- list()
all.prop.perms <- seq(from = 0, to = 1, by = 0.05)

for (k in 1:length(outY1vec)){
  
  outY1 <- outY1vec[k]
  
  for (i in 1:length(all.prop.perms)){
    
    this.prop.perm <- all.prop.perms[i]
    n.sims <- 1000
    hold.mses <- list()
    for (j in 1:n.sims){
      
      this.dat <- prep_data(dat = dat, 
                            outY2 = outY2, 
                            outY1 = outY1, 
                            prop.perm = this.prop.perm,
                            this.seed = j)
      these.mses <- run_analysis(dat = this.dat$dat, 
                                 dat1 = this.dat$dat1, 
                                 dat2 = this.dat$dat2, 
                                 dat3 = this.dat$dat3,
                                 this.seed = j)
      hold.mses[[j]] <- these.mses
      rm(this.dat,these.mses)
      
    }
    
    all.mses <- do.call("rbind",hold.mses)
    
    all.res[[length(all.res) + 1]] <-     
      cbind(data.frame(outY2 = outY2,
                       outY1 = outY1),
            data.frame(t(c(propperm = this.prop.perm, 
                           colMeans(all.mses),
                           setNames(apply(all.mses,2,var),
                                    paste0(colnames(all.mses),"_","var")),
                           n_run = n.sims))))
    
    rm(this.prop.perm,hold.mses,all.mses)
    print(i)
    
  }
  
  rm(outY1)
  print(k)
  
}

all.res <- as.data.frame(do.call("rbind",all.res))
save(all.res, file = "bpr_project_star_results.Rdata")


# Evaluate ----------------------------------------------------------------

load("bpr_project_star_results.Rdata")

#MSE Results
all.res.melt.mse <- as.data.frame(
  all.res %>% select(outY2,outY1,propperm,n_run,starts_with("mse") & !ends_with("var")) %>%
    pivot_longer(cols = starts_with("mse"),
                 names_to = "Method", 
                 names_prefix = "mse",
                 values_to = c("MSE"))
)
add.res.melt.mse.var <- as.data.frame(
  all.res %>% select(outY2,outY1,propperm,starts_with("mse") & ends_with("var")) %>%
    pivot_longer(cols = starts_with("mse") & ends_with("var"),
                 names_to = "Method", 
                 names_prefix = "mse",
                 values_to = c("MSE_var"))
)
add.res.melt.mse.var$Method <- gsub("_var","",add.res.melt.mse.var$Method)

all.res.melt.mse <- merge(all.res.melt.mse,add.res.melt.mse.var,
                          by = c("outY2","outY1","propperm","Method"))

all.res.melt.mse$low <- all.res.melt.mse$MSE + 
  qnorm(0.005)*sqrt(all.res.melt.mse$MSE_var/all.res.melt.mse$n_run)
all.res.melt.mse$high <- all.res.melt.mse$MSE + 
  qnorm(0.995)*sqrt(all.res.melt.mse$MSE_var/all.res.melt.mse$n_run)

all.res.melt.mse$Method <- factor(all.res.melt.mse$Method,
                                  levels = c("A","B","C","intercept"))

all.res.melt.mse$pcor <- NA
all.res.melt.mse$pcor[all.res.melt.mse$outY1 == "g1ttotalss"] <- 
  paste0("G1 SAT Test Proxy\ncor(Y1,Y2) = ",round(cor(dat$g2ttotalss,dat$g1ttotalss),2))
all.res.melt.mse$pcor[all.res.melt.mse$outY1 == "g1tmathss"] <- 
  paste0("G1 Math Test Proxy\ncor(Y1,Y2) = ",round(cor(dat$g2ttotalss,dat$g1tmathss),2))
all.res.melt.mse$pcor[all.res.melt.mse$outY1 == "gktmathss"] <- 
  paste0(" KG Math Test Proxy \ncor(Y1,Y2) = ",round(cor(dat$g2ttotalss,dat$gktmathss),2))
all.res.melt.mse$pcor[all.res.melt.mse$outY1 == "gkpresent"] <- 
  paste0("  KG Attendance Proxy  \ncor(Y1,Y2) = ",round(cor(dat$g2ttotalss,dat$gkpresent),2))

pdf(file = paste0("ps_mse.pdf"),width=8.5,height=3)
ggplot(all.res.melt.mse, aes(x = propperm, y = MSE, color = Method, shape = Method)) +
  geom_point() + geom_line() +
  theme_bw() + facet_wrap(~ pcor, nrow = 1) +
  scale_color_manual(values = c("dodgerblue3","magenta4","seagreen4","grey40"),
                     labels = c('A' = expression(tau^A),
                                'B' = expression(tilde(tau)^B),
                                'C' = expression(tau^C),
                                'intercept' = "constant")) +
  scale_shape_manual(values = c(0,1,17,3),
                     labels = c('A' = expression(tau^A),
                                'B' = expression(tilde(tau)^B),
                                'C' = expression(tau^C),
                                'intercept' = "constant")) +
  xlab("Degree of Distribution Shift") +
  theme(panel.spacing = unit(1, "lines")) +
  theme(legend.position = "top") + 
  ylab(expression("MSE (w.r.t. "*tau*")"))
dev.off()

pdf(file = paste0("ps_mse_CI.pdf"),width=8.5,height=3)
ggplot(all.res.melt.mse, aes(x = propperm, y = MSE, color = Method, shape = Method)) +
  geom_line(linewidth = 0.2) + geom_linerange(aes(ymin = low, ymax = high)) +
  theme_bw() + facet_wrap(~ pcor, nrow = 1) +
  scale_color_manual(values = c("dodgerblue3","magenta4","seagreen4","grey40"),
                     labels = c('A' = expression(tau^A),
                                'B' = expression(tilde(tau)^B),
                                'C' = expression(tau^C),
                                'intercept' = "constant")) +
  scale_shape_manual(values = c(0,1,17,3),
                     labels = c('A' = expression(tau^A),
                                'B' = expression(tilde(tau)^B),
                                'C' = expression(tau^C),
                                'intercept' = "constant")) +
  xlab("Degree of Distribution Shift") +
  theme(panel.spacing = unit(1, "lines")) +
  theme(legend.position = "top") + 
  ylab(expression("MSE (w.r.t. "*tau*")"))
dev.off()

#R2 Results
all.res.melt.R2 <- as.data.frame(
  all.res %>% select(outY2,outY1,propperm,n_run,starts_with("R2") & !ends_with("var")) %>%
    pivot_longer(cols = starts_with("R2"),
                 names_to = "Method", 
                 names_prefix = "R2",
                 values_to = c("R2"))
)
add.res.melt.R2.var <- as.data.frame(
  all.res %>% select(outY2,outY1,propperm,starts_with("R2") & ends_with("var")) %>%
    pivot_longer(cols = starts_with("R2") & ends_with("var"),
                 names_to = "Method", 
                 names_prefix = "R2",
                 values_to = c("R2_var"))
)
add.res.melt.R2.var$Method <- gsub("_var","",add.res.melt.R2.var$Method)

all.res.melt.R2 <- merge(all.res.melt.R2,add.res.melt.R2.var,
                          by = c("outY2","outY1","propperm","Method"))

all.res.melt.R2$low <- all.res.melt.R2$R2 + 
  qnorm(0.005)*sqrt(all.res.melt.R2$R2_var/all.res.melt.R2$n_run)
all.res.melt.R2$high <- all.res.melt.R2$R2 + 
  qnorm(0.995)*sqrt(all.res.melt.R2$R2_var/all.res.melt.R2$n_run)

all.res.melt.R2$Method <- factor(all.res.melt.R2$Method,
                                 levels = c("T","A","B","C"))

all.res.melt.R2$pcor <- NA
all.res.melt.R2$pcor[all.res.melt.R2$outY1 == "g1ttotalss"] <- 
  paste0("G1 SAT Test Proxy\ncor(Y1,Y2) = ",round(cor(dat$g2ttotalss,dat$g1ttotalss),2))
all.res.melt.R2$pcor[all.res.melt.R2$outY1 == "g1tmathss"] <- 
  paste0("G1 Math Test Proxy\ncor(Y1,Y2) = ",round(cor(dat$g2ttotalss,dat$g1tmathss),2))
all.res.melt.R2$pcor[all.res.melt.R2$outY1 == "gktmathss"] <- 
  paste0(" KG Math Test Proxy \ncor(Y1,Y2) = ",round(cor(dat$g2ttotalss,dat$gktmathss),2))
all.res.melt.R2$pcor[all.res.melt.R2$outY1 == "gkpresent"] <- 
  paste0("  KG Attendance Proxy  \ncor(Y1,Y2) = ",round(cor(dat$g2ttotalss,dat$gkpresent),2))

pdf(file = paste0("ps_r2.pdf"),width=8.5,height=3.25)
ggplot(all.res.melt.R2, aes(x = propperm, y = R2, color = Method, shape = Method)) +
  geom_point() + geom_line() +
  theme_bw() + facet_wrap(~ pcor, nrow = 1) +
  scale_color_manual(values = c("black","dodgerblue3","magenta4","seagreen4"),
                     labels = c('T' = expression(tau),
                                'A' = expression(tau^A),
                                'B' = expression(tilde(tau)^B),
                                'C' = expression(tau^C))) +
  scale_shape_manual(values = c(16,0,1,17,3),
                     labels = c('T' = expression(tau),
                                'A' = expression(tau^A),
                                'B' = expression(tilde(tau)^B),
                                'C' = expression(tau^C))) +
  xlab("Degree of Distribution Shift") +
  theme(panel.spacing = unit(1, "lines")) +
  theme(legend.position = "top") + 
  ylab("R-Squared\n(on Target Data)") +
  geom_hline(yintercept = 0, color = "grey35", linetype = 2)
dev.off()

pdf(file = paste0("ps_r2_CI.pdf"),width=8.5,height=3.25)
ggplot(all.res.melt.R2, aes(x = propperm, y = R2, color = Method, shape = Method)) +
  geom_line(linewidth = 0.2) + geom_linerange(aes(ymin = low, ymax = high)) +
  theme_bw() + facet_wrap(~ pcor, nrow = 1) +
  scale_color_manual(values = c("black","dodgerblue3","magenta4","seagreen4"),
                     labels = c('T' = expression(tau),
                                'A' = expression(tau^A),
                                'B' = expression(tilde(tau)^B),
                                'C' = expression(tau^C))) +
  scale_shape_manual(values = c(16,0,1,17,3),
                     labels = c('T' = expression(tau),
                                'A' = expression(tau^A),
                                'B' = expression(tilde(tau)^B),
                                'C' = expression(tau^C))) +
  xlab("Degree of Distribution Shift") +
  theme(panel.spacing = unit(1, "lines")) +
  theme(legend.position = "top") + 
  ylab("R-Squared\n(on Target Data)") +
  geom_hline(yintercept = 0, color = "grey35", linetype = 2)
dev.off()

