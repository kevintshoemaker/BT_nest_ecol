### Macey et al. nest protection -------------

# Questions --------------

#  Were nests re-protected after a predation event?

# Clear workspace -----------------

rm(list=ls())

# load packages --------------

library(Hmisc)
library(brms)
library(glmmTMB)
library(DHARMa)
library(tidyverse)
library(lubridate)
library(rstan)
library(jagsUI)

# remove.packages(c("StanHeaders", "rstan"))
# install.packages(
#   "StanHeaders",
#   repos = c(
#     "https://mc-stan.org/r-packages/",
#     getOption("repos")
#   )
# )
# install.packages(
#   "rstan",
#   repos = c(
#     "https://mc-stan.org/r-packages/",
#     getOption("repos")
#   )
# )


# read data --------------

eggdf <- read.csv("egg_success_150624.csv")
temp <- as.Date(eggdf$Day_Nested, origin = "1899-12-30")   # issues need to be fixed

eggdf$DATE <- ymd(sprintf("%s-01-01", eggdf$Year)) + yday(temp)

names(eggdf)

eggdf <- eggdf[,!grepl("X",names(eggdf))]
names(eggdf)

eggdf$NESTTEMP <- (eggdf$Inc_Mean_nest_daily_min_temp_fifty_five_days + eggdf$Inc_Mean_nest_daily_max_temp_fifty_five_days)/2

eggdf$MINTEMP <- eggdf$Inc_Mean_nest_daily_min_temp_fifty_five_days

keep <- c(
  "Nest","Merged_ID","Site","Year","DATE",
  "NESTTEMP","MINTEMP",
  "Excluders","Mammal_Predated_Egg","Hatched_Egg","Unsucc_Egg",
  "Estimated_Mass","Carapace_Length","BCI"
)
lu <- data.frame(
  full=keep,
  short=c("NEST","femID","SITE","YEAR","DATE","NESTTEMP","MINTEMP","EXCLUD","PRED","HATCH","UNSUCC","MASS","LENGTH","BCI")
)
eggdf <- eggdf[,lu$full]
names(eggdf) <- lu$short

eggdf <- eggdf[order(eggdf$NEST),]

table(eggdf$SITE,eggdf$YEAR)
table(eggdf$SITE,eggdf$YEAR,eggdf$EXCLUD)
eggdf$SITEYEAR <- paste0(eggdf$SITE,eggdf$YEAR)

# make nest df ------------------

nestdf <- eggdf %>% 
  group_by(NEST) %>% 
  summarise(SITE=first(SITE),
            YEAR=first(YEAR),
            DATE=first(DATE),
            femID=first(femID),            
            femCL=first(LENGTH),           # fem body size
            EXCLUD=first(EXCLUD),
            NESTTEMP=first(NESTTEMP),
            MINTEMP=first(MINTEMP),
            EGGS=n(),                      # clutch size
            HATCH=sum(HATCH),              # total hatched successfully
            PREDtot=sum(PRED),                  
            PRED=ifelse(any(PRED==1),1,0)    # pred or no
            # HATCH=ifelse(any(HATCH==1),1,0),
            # UNSUCCrate=sum(UNSUCC)/n(),
            # UNSUCC=ifelse(any(UNSUCC==1),1,0)
            )

mean(nestdf$EGGS)   # mean clutch size of 3.34
mean(nestdf$HATCH)    # mean of 1.18 successfully hatched per nest

table(nestdf$SITE,nestdf$YEAR)
table(nestdf$SITE,nestdf$YEAR,nestdf$EXCLUD)
nestdf$SITEYEAR <- paste0(nestdf$SITE,nestdf$YEAR)

table(nestdf$femID)

nestdf$nestID <- paste0(nestdf$femID,"_",nestdf$YEAR)

# incubation df ------------------

incdf <- read.csv("incubation_no_temp_excluder_150702a.csv")
incdf$Day_Nested
temp <- as.Date(incdf$Day_Nested, origin = "1899-12-30")   # issues need to be fixed

incdf$DATE <- ymd(sprintf("%s-01-01", incdf$Year)) + yday(temp)
incdf$YEAR <- year(incdf$DATE)

incdf$nestID <- paste0(incdf$Merged_ID,"_",incdf$YEAR)
names(incdf)

keep <- c(
  "nestID",
  "DATE",
  "YEAR",
  "Site",
  "Merged_ID",
  "Mean_Incubation",
  "Excluders",
  "Eggs_Laid",
  "Predated_Y_N",
  "Mammal_Predated_Eggs",
  "Hatched_Eggs"
)

lu2 <- data.frame(
  full=keep,
  short=c("nestID","DATE","YEAR","SITE","femID","INCDAYS",
          "EXCLUD","EGGS","PRED","PREDtot","HATCH")
)
incdf <- incdf[,lu2$full]
names(incdf) <- lu2$short

# add temperature data to incubation dataframe
incdf$NESTTEMP <- nestdf$NESTTEMP[match(incdf$nestID,nestdf$nestID)]
incdf$MINTEMP <- nestdf$MINTEMP[match(incdf$nestID,nestdf$nestID)]

# add incubation info to nest dataframe

nestdf$INCDAYS <- incdf$INCDAYS[match(nestdf$nestID,incdf$nestID)]

# prelim analyses --------------------

t1 <- table(eggdf$EXCLUD,eggdf$PRED)     # egg level
t1_row <- rowSums(t1)
t1_col <- colSums(t1)
t1[,2]/t1_row   # 29% egg pred with excluders, 60% egg pred without excluders
t1    # matches table 2
chisq.test(t1)   # highly significant

with(nestdf,tapply(PREDtot,EXCLUD,mean,na.rm=T))


t1 <- table(nestdf$EXCLUD,nestdf$PRED)    # nest level
t1_row <- rowSums(t1)
t1_col <- colSums(t1)
t1[,2]/t1_row   # 37.7% nest pred with excluders, 62.5% nest pred without excluders
t1    # matches table 2
chisq.test(t1)   # not significant (due to low sample size)


t2 <- table(eggdf$EXCLUD,eggdf$HATCH)
t2_row <- rowSums(t2)
t2_col <- colSums(t2)
t2[,2]/t2_row   # 37.8% hatching success with excluders, 30% hatching success without excluders
t2    # matches table 2
chisq.test(t2)   # not significant


t2 <- table(nestdf$EXCLUD,ifelse(nestdf$HATCH>0,1,0))
t2_row <- rowSums(t2)
t2_col <- colSums(t2)
t2[,2]/t2_row   # 54.7% nest 'success' with excluders, 37.5% nest success without excluders (meaning at least one hatched)
t2    # matches table 2
chisq.test(t2)   # not significant

rm(keep,t1,t1_col,t1_row,t2,t2_col,t2_row,temp)

 # only non-predated eggs
eggdf_nopred <- subset(eggdf,PRED==0)

t3 <- table(eggdf_nopred$EXCLUD,eggdf_nopred$UNSUCC)
t3_row <- rowSums(t3)
t3_col <- colSums(t3)
t3[,2]/t3_row   # 47% unsuccessful with excluders, 25% unsuccessful without excluders
t3    # matches table 2
chisq.test(t3)   # significant (barely!)
# NOTE: there seems to be a negative effect of excluders on egg hatching rate... 

tapply(nestdf$NESTTEMP,nestdf$EXCLUD,mean,na.rm=T)
tapply(nestdf$NESTTEMP,nestdf$EXCLUD,sd,na.rm=T)

summary(lm(NESTTEMP~EXCLUD,data=nestdf))  # strong effect of nest excluders on nest temperature??

rm(eggdf_nopred,t3,t3_col,t3_row)

summary(lm(EGGS~femCL,nestdf))
summary(lm(EGGS~femCL,nestdf))

# Incubation prelim analyses ----------------

plot(incdf$INCDAYS~incdf$NESTTEMP)
plot(incdf$INCDAYS~incdf$MINTEMP)

ggplot(incdf,aes(MINTEMP,INCDAYS))+
  geom_point(aes(col=EXCLUD))

incmod1 <- lm(INCDAYS~NESTTEMP,incdf)
summary(incmod1)

incmod2 <- lm(INCDAYS~MINTEMP,incdf)
summary(incmod2)     # I can't replicate suzanne's result...

tapply(incdf$INCDAYS,incdf$EXCLUD,mean)

tt1 <- t.test(INCDAYS~EXCLUD,incdf)
tt1  # not significant- slightly longer incubation for nests with excluder

tapply(incdf$NESTTEMP,incdf$EXCLUD,mean,na.rm=T)  # warmer nest temps for excluder- opposite of previous result!

tapply(nestdf$NESTTEMP,nestdf$EXCLUD,mean,na.rm=T)  # opposite effect in nest dataframe- more info

t.test(NESTTEMP~EXCLUD,nestdf)  # significant result for nests

summary(lm(NESTTEMP~EXCLUD+YEAR,nestdf))

tmeans <- with(nestdf,tapply(NESTTEMP,YEAR,mean,na.rm=T)  )
# tsds <- with(nestdf,tapply(NESTTEMP,YEAR,sd,na.rm=T)  )

with(nestdf,tapply(NESTTEMP,EXCLUD,mean,na.rm=T)  )

incdf2 <- incdf[!is.na(incdf$NESTTEMP),]  # only one nest without excluder that has both temp and incubation timing...

summary(lm(INCDAYS~SITE+EXCLUD,incdf))

table(incdf$SITE)


# Try JAGS ------------------------

cat("
model{

# carapace length (interpolation only)

cl.mean ~ dnorm(0,1)  # standardized...
cl.sd ~ dunif(0,5)
cl.prec <- pow(cl.sd,-2)
for(n in 1:nnests){
  cl[n] ~ dnorm(cl.mean,cl.prec)  # data node
}

# nest temperature process (interpolation)

for(y in 1:nyears){
  nt.base[y] ~ dnorm(0,1)
}

nt.siteeff.sd ~ dgamma(3,3)
nt.siteeff.prec <- pow(nt.siteeff.sd,-2)
for(s in 1:nsites){
  nt.siteeff[s] ~ dnorm(0,nt.siteeff.prec)
}

# nt.base ~ dnorm(0,1)   
nt.excl.eff ~ dnorm(0,0.1)  # effect of excluder on temperature
nt.sd ~ dunif(0,5)
nt.prec <- pow(nt.sd,-2)
for(n in 1:nnests){
  nt.exp[n] <- nt.base[yr_f[n]] + nt.excl.eff * is_excl[n] + nt.siteeff[site_f[n]]
  nt[n] ~ dnorm(nt.exp[n],nt.prec)    #   data node
}

# nest predation process

   # random effect
pred.syeff.sd ~ dgamma(3,3)
pred.syeff.prec <- pow(pred.syeff.sd,-2)
for(t in 1:nsy){
  pred.syeff[t] ~ dnorm(0,pred.syeff.prec)
}

pred.base ~ dunif(0,1)    # probability of nest being predated
pred.base.l <- log(pred.base/(1-pred.base))
pred.excl.eff ~ dnorm(0,.1)
for(n in 1:nnests){  # loop through nests
  logit(pred.exp[n]) <-  pred.base.l + pred.excl.eff*is_excl[n] + pred.syeff[sy[n]]   # predation process- at the nest level
}

# number of eggs predated...

eggpred ~ dunif(0,1)   # probability of egg being predated if nest is predated

# clutch size  (total eggs laid)

neggs.base ~ dunif(0,10)
neggs.base.l <- log(neggs.base)
neggs.cl.eff ~ dnorm(0,.1)
for(n in 1:nnests){  # loop through nests
  log(neggs.exp[n]) <- neggs.base.l + neggs.cl.eff*cl[n]
}

# probability of hatching successfully if not predated

hatch.ebv.sd ~ dunif(0.25,5)
hatch.ebv.prec <- pow(hatch.ebv.sd,-2)
hatch.base ~ dunif(0,1)
hatch.base.l <- log(hatch.base/(1-hatch.base))
hatch.nt.eff ~ dnorm(0,1)   # effect of nest temp on the prob of hatching successfully
hatch.excl.eff ~ dnorm(0,.1)    # effect of pred excl on prob of hatching successfully
for(n in 1:nnests){  # loop through nests
  hatch.ebv[n] ~ dnorm(0,hatch.ebv.prec)
  logit(hatch.exp[n]) <- hatch.base.l + hatch.nt.eff*nt[n] + 
                          hatch.excl.eff*is_excl[n] + 
                          hatch.ebv[n]
}

# Incubation
inc.base ~ dunif(50,100)
inc.sd ~ dunif(0,10)
inc.prec <- pow(inc.sd,-2)
inc.excl.eff ~ dnorm(0,.1)
inc.nt.eff ~ dnorm(0,1)   # effect of nest temp on incubation duration

   # random effect
# inc.syeff.sd ~ dunif(0,5)
# inc.syeff.prec <- pow(inc.syeff.sd,-2)
# for(t in 1:nsy){
#   inc.syeff[t] ~ dnorm(0,inc.syeff.prec)
# }

for(n in 1:nnests){
  inc.exp[n] <- inc.base + inc.nt.eff*nt[n] +
                 inc.excl.eff*is_excl[n]    #+ 
                 # inc.syeff[sy[n]]
}

# likelihood

for(n in 1:nnests){  # loop through nests
  pred[n] ~ dbern(pred.exp[n])  # data node: is the nest predated or no?
  neggs[n] ~ dpois(neggs.exp[n])  # data node: clutch size
  npred[n] ~ dbin(eggpred,neggs[n])   # data node: number of eggs predated
  inc[n] ~ dnorm(inc.exp[n],inc.prec)   # data node: incubation timing
}

for(n2 in 1:nnests2){  # loop through nests with at least one non-predated egg
  hatched[ndx[n2]] ~ dbin(hatch.exp[ndx[n2]],neggs.rem[ndx[n2]])   # data.node: number of eggs successfully hatched
}  

# posterior predictive checks

for(n in 1:nnests){  # loop through nests
  hatched.sim[n] ~ dbin(hatch.exp[n],neggs.rem[n])
  sq.err.obs[n] <- pow(hatched[n]-(hatch.exp[n]*neggs.rem[n]),2)
  sq.err.sim[n] <- pow(hatched.sim[n]-(hatch.exp[n]*neggs.rem[n]),2)
}

sse.obs <- sum(sq.err.obs[1:nnests])
sse.sim <- sum(sq.err.sim[1:nnests])

# derived parameters: excl vs no

  # nesttemp
nt.excl <- nt.base + nt.excl.eff
nt.noexcl <- nt.base

  # predation, nest level
logit(pred.excl) <-  pred.base.l + pred.excl.eff
logit(pred.noexcl) <-  pred.base.l

  # clutch size
neggs.this <- neggs.base

  # eggs remaining after predation
neggsrem.excl <- neggs.this * (pred.excl*eggpred + (1-pred.excl))
neggsrem.noexcl <- neggs.this * (pred.noexcl*eggpred + (1-pred.noexcl))

logit(phatch.excl) <- hatch.base.l + hatch.excl.eff # + hatch.nt.eff*nt.excl
logit(phatch.noexcl) <- hatch.base.l #+ hatch.nt.eff*nt.noexcl

hatch.excl <- neggsrem.excl * phatch.excl
hatch.noexcl <- neggsrem.noexcl * phatch.noexcl

}
    
",file="jagsmod5.txt")


# Prepare for JAGS ---------------

nt2 <- (nestdf$NESTTEMP-mean(nestdf$NESTTEMP,na.rm=T))/sd(nestdf$NESTTEMP,na.rm=2)
# hist(nt2)
# nt2[is.na(nt2)] = 0
# nt2[1] <- -0.5
yr_f <- as.numeric(as.factor(nestdf$YEAR))
nyears <- length(unique(yr_f))

site_f <- as.numeric(as.factor(nestdf$SITE))
nsites <- max(site_f)

# temp <- subset(nestdf,is.na(YEAR))  

dfj <- list(
  nsy = length(sort(unique(nestdf$SITEYEAR))),    # 21 siteyears
  sy = as.numeric(as.factor(nestdf$SITEYEAR)),
  cl = (nestdf$femCL-mean(nestdf$femCL,na.rm=T))/sd(nestdf$femCL,na.rm=2),
  is_excl = ifelse(nestdf$EXCLUD=="Y",1,0),
  nnests = nrow(nestdf),
  nnests2 = length(which(nestdf$EGGS-nestdf$PREDtot>0)),
  ndx = which(nestdf$EGGS-nestdf$PREDtot>0),
  nt = nt2,
  pred = nestdf$PRED,
  neggs = nestdf$EGGS,
  npred = nestdf$PREDtot,
  inc = nestdf$INCDAYS,
  neggs.rem = nestdf$EGGS-nestdf$PREDtot,
  hatched = nestdf$HATCH,
  nyears = nyears,
  yr_f = yr_f,
  nsites = nsites,
  site_f = site_f
)
dfj
# any(dfj$npred>dfj$neggs)
# any(dfj$hatched>dfj$neggs.rem)
# any((dfj$hatched+dfj$npred)>dfj$neggs)
# length(which(((dfj$hatched+dfj$npred)==dfj$neggs)))

ntinit <- dfj$nt
ntinit[is.na(ntinit)] <- 0
ntinit[!is.na(dfj$nt)] <- NA
ifj <- function(){
  list(
    cl.mean=runif(1,-.01,0.01),
    cl.sd=runif(1,0.1,0.2),
    nt.base=runif(nyears,-.01,0.01),
    nt.sd=runif(1,0.1,0.2),
    nt.excl.eff=runif(1,-.1,0.1),
    nt.siteeff.sd=runif(1,0.1,0.2),
    pred.base=runif(1,0.3,0.4),
    pred.excl.eff=runif(1,-.1,0.1),
    pred.syeff.sd=runif(1,0.05,0.1),
    eggpred=runif(1,0.4,0.5),
    neggs.base=runif(1,2.5,3.5),
    neggs.cl.eff=runif(1,0.1,0.2),
    hatch.base=runif(1,0.4,0.5),
    hatch.ebv.sd=runif(1,0.3,0.4),
    hatch.ebv=runif(dfj$nnests,-0.01,0.01),
    hatch.nt.eff=runif(1,-.1,0.1),
    hatch.excl.eff=runif(1,-.1,0.1),
    inc.base=rnorm(1,79,2),
    inc.sd=runif(1,0.3,0.5),
    inc.excl.eff=runif(1,-.1,0.1),
    inc.nt.eff=runif(1,-.1,0.1),
    # inc.syeff.sd=runif(1,0.05,0.1),
    nt=ntinit
  )
}
ifj()

pts <- c(
  "cl.mean",
  "cl.sd",
  "nt.base",
  "nt.sd",
  "nt.siteeff.sd",
  "nt.siteeff",
  "nt.excl.eff",
  "pred.base",
  "pred.excl.eff",
  "pred.syeff.sd",
  "pred.syeff",
  "eggpred",
  "neggs.base",
  "neggs.cl.eff",
  "hatch.base",
  "hatch.ebv.sd",
  "hatch.nt.eff",
  "hatch.excl.eff",
  "sse.obs",
  "sse.sim",
  "nt.excl",
  "pred.excl",
  "neggsrem.excl",
  "phatch.excl",
  "hatch.excl",
  "nt.noexcl",
  "pred.noexcl",
  "inc.base",
  "inc.sd",
  "inc.excl.eff",
  "inc.nt.eff",
  "neggsrem.noexcl",
  "phatch.noexcl",
  "hatch.noexcl"
)
pts

# ?jags
ni=500000    # takes 15 mins or so to run
nb=250000
nt=100
nc=3
mod=jags(dfj, ifj, pts, "jagsmod5.txt",
     n.chains=nc, n.adapt=100, n.iter=ni, n.burnin=nb, n.thin=nt,
     parallel=TRUE, n.cores=nc)
saveRDS(mod,"output/JAGSout5.rds")

# mod=readRDS("output/JAGSout5.rds")

sims <- mod$sims.list

# make parameter table -------------------------

pts

t=pts[1]
tabfunc <- function(t){
  x=mod$sims.list[[t]]
  ret <- NULL
  if(length(dim(x))>1){
    for(i in 1:dim(x)[2]){
      this <-t(as.data.frame(quantile(x[,i],c(0.5,0.025,0.975))))
      rownames(this) <- paste0(t,"_",i)
      ret <- rbind(ret,this)
    }
  }else{
    this <- t(as.data.frame(quantile(x,c(0.5,0.025,0.975))))
    rownames(this) <- t
    ret <- rbind(ret,this)
    
  }
  return(ret)
}
# tabfunc(t)

temp <- lapply(pts,function(t) tabfunc(t)  )
table1 <- do.call(rbind,temp)

write.csv(table1,"output/table1_5.csv",row.names=T)

# visualize posterior ---------------------

# hist(sims$nt.base)
traceplot(mod,"nt.base")
densityplot(mod,"nt.base")

# summary(lm(NESTTEMP~as.factor(YEAR)+EXCLUD,nestdf))  # confirm result

hist(sims$nt.sd)
traceplot(mod,"nt.sd")
densityplot(mod,"nt.sd")

hist(sims$nt.siteeff.sd)
densityplot(mod,"nt.siteeff.sd")

hist(sims$nt.excl.eff)
traceplot(mod,"nt.excl.eff")
densityplot(mod,"nt.excl.eff")   # strong effect of excluder on temperature (disappears after accounting for among-year differences)
quantile(sims$nt.excl.eff,c(0.5,0.025,0.975))   # -1.16 -1.83 -0.49    (no overlap with zero...)

hist(sims$cl.mean)
# traceplot(mod,"cl.mean")
densityplot(mod,"cl.mean")

hist(sims$cl.sd)
# traceplot(mod,"cl.sd")
densityplot(mod,"cl.sd")

hist(sims$pred.base)
# traceplot(mod,"pred.base")
densityplot(mod,"pred.base")

hist(sims$pred.excl.eff)
# traceplot(mod,"pred.excl.eff")
densityplot(mod,"pred.excl.eff")    # pred excluders prevent predation!
quantile(sims$pred.excl.eff,c(0.5,0.025,0.975))  # -1.14 -2.38  0.029     

hist(sims$pred.syeff.sd)
traceplot(mod,"pred.syeff.sd")
densityplot(mod,"pred.syeff.sd")  

# hist(sims$pred.syeff[,1])

hist(sims$eggpred)
# traceplot(mod,"eggpred")
densityplot(mod,"eggpred")   # about 40% of eggs are destroyed in a predated nest

hist(sims$neggs.base)
# traceplot(mod,"neggs.base")
densityplot(mod,"neggs.base")
quantile(sims$neggs.base,c(0.5,0.025,0.975))

hist(sims$neggs.cl.eff)
# traceplot(mod,"neggs.cl.eff")
densityplot(mod,"neggs.cl.eff")  # no effect of CL on egg production (maybe remove this?)
quantile(sims$neggs.cl.eff,c(0.5,0.025,0.975))

hist(sims$hatch.base)
# traceplot(mod,"hatch.base")
densityplot(mod,"hatch.base")   # close to 100%

hist(sims$hatch.ebv.sd)
# traceplot(mod,"hatch.ebv.sd")
densityplot(mod,"hatch.ebv.sd")   # converged

hist(sims$hatch.nt.eff)
# traceplot(mod,"hatch.nt.eff")
densityplot(mod,"hatch.nt.eff")   # no effect of nest temperature on hatch success? [not in model currently]

hist(sims$hatch.excl.eff)
# traceplot(mod,"hatch.excl.eff")
densityplot(mod,"hatch.excl.eff")  # neg effect of excluder on nest success
quantile(sims$hatch.excl.eff,c(0.5,0.025,0.975))   # -1.32 -3.48  0.74    (overlaps zero but...)

hist(sims$sse.obs)
# traceplot(mod,"sse.obs")
densityplot(mod,"sse.obs")

# hist(sims$sse.sim)
# # traceplot(mod,"sse.sim")
# densityplot(mod,"sse.sim")

hist(sims$nt.excl)
# traceplot(mod,"nt.excl")
densityplot(mod,"nt.excl[3]")

hist(sims$nt.noexcl)
# traceplot(mod,"nt.noexcl")
densityplot(mod,"nt.noexcl[3]")

hist(sims$pred.excl)
# traceplot(mod,"pred.excl")
densityplot(mod,"pred.excl")

hist(sims$pred.noexcl)
# traceplot(mod,"pred.noexcl")
densityplot(mod,"pred.noexcl")

hist(sims$neggsrem.excl)
# traceplot(mod,"neggsrem.excl")
densityplot(mod,"neggsrem.excl")

hist(sims$neggsrem.noexcl)
# traceplot(mod,"neggsrem.noexcl")
densityplot(mod,"neggsrem.noexcl")

hist(sims$phatch.excl)
# traceplot(mod,"phatch.excl")
densityplot(mod,"phatch.excl")

hist(sims$phatch.noexcl)
# traceplot(mod,"phatch.noexcl")
densityplot(mod,"phatch.noexcl")

hist(sims$hatch.excl)
# traceplot(mod,"hatch.excl")
densityplot(mod,"hatch.excl")

hist(sims$hatch.noexcl)
# traceplot(mod,"hatch.noexcl")
densityplot(mod,"hatch.noexcl")

png("figs/pvalfig.png",4,4,units="in",res=600)
plot(sims$sse.obs,sims$sse.sim,xlab="SSE, observed", ylab="SSE, simulated")    # good model fit...
abline(1,1,lwd=3,col=gray(0.7))
dev.off()

pval = sum(sims$sse.sim>sims$sse.obs)/nmcmc


hist(sims$inc.base)
# traceplot(mod,"inc.base")
densityplot(mod,"inc.base")
round(quantile(sims$inc.base,c(0.5,0.025,0.975)),1)

hist(sims$inc.sd)
# traceplot(mod,"inc.sd")
densityplot(mod,"inc.sd")

hist(sims$inc.nt.eff)
densityplot(mod,"inc.nt.eff")

hist(sims$inc.excl.eff)
# traceplot(mod,"inc.excl.eff")
densityplot(mod,"inc.excl.eff")
round(quantile(sims$inc.excl.eff,c(0.5,0.025,0.975)),1)

# visualize effects ---------------------- 

library(ggplot2)
nmcmc <- length(sims$cl.mean)

# nesttemp
nt.excl <- sims$nt.base[,3] + sims$nt.excl.eff
nt.noexcl <- sims$nt.base[,3]

# predation, nest level
pred.excl <-  plogis(qlogis(sims$pred.base) + sims$pred.excl.eff)
pred.noexcl <-  sims$pred.base

# clutch size
neggs.this <- sims$neggs.base

# eggs remaining after predation
neggsrem.excl <- neggs.this * (pred.excl*sims$eggpred + (1-pred.excl))
neggsrem.noexcl <- neggs.this * (pred.noexcl*sims$eggpred + (1-pred.noexcl))

# hatching probability
phatch.excl <- plogis(qlogis(sims$hatch.base) + sims$hatch.excl.eff + sims$hatch.nt.eff*nt.excl) 
phatch.noexcl <- plogis(qlogis(sims$hatch.base) + sims$hatch.nt.eff*nt.noexcl)

hatch.excl <- neggsrem.excl * phatch.excl
hatch.noexcl <- neggsrem.noexcl * phatch.noexcl

inc.excl <- sims$inc.base+sims$inc.excl.eff+sims$inc.nt.eff*nt.excl
inc.noexcl <- sims$inc.base + sims$inc.nt.eff*nt.noexcl

# incubation duration

df <- data.frame(
  value = c(inc.noexcl,inc.excl),
  excl = rep(c("N","Y"),each=nmcmc)
)

gginc <- ggplot(df,aes(excl,value))+
  # geom_boxplot(col=gray(0.2),fill="white",size=0.5,width=0.75) +
  geom_violin(col=alpha("black",0.5),fill=alpha("white",0.5),lwd=1.4,draw_quantiles = TRUE) +
  labs(x="Predator excluder?",y="Incubation duration",title = "")+
  theme_bw()
gginc

   # inc duration vs temperature

temps <- seq(min(nestdf$NESTTEMP,na.rm=T),max(nestdf$NESTTEMP,na.rm=T),length=20)
temps2 <- (temps-mean(nestdf$NESTTEMP,na.rm=T))/sd(nestdf$NESTTEMP,na.rm=T)

preds1 <- lapply(temps2,function(t) sims$inc.base + sims$inc.nt.eff*t  )
preds2 <- lapply(temps2,function(t) sims$inc.base + sims$inc.excl.eff + sims$inc.nt.eff*t  )

df1 <- data.frame(
  med = sapply(preds1,function(t) quantile(t,0.5)),
  lb = sapply(preds1,function(t) quantile(t,0.025)),
  ub = sapply(preds1,function(t) quantile(t,0.975)),
  temp = temps
)
df2 <- data.frame(
  med = sapply(preds2,function(t) quantile(t,0.5)),
  lb = sapply(preds2,function(t) quantile(t,0.025)),
  ub = sapply(preds2,function(t) quantile(t,0.975)),
  temp = temps
)
df1$excl <- "N"
df2$excl <- "Y"
df <- rbind(df1,df2)

gginc_temp <- ggplot(df1,aes(x=temp,y=med))+
  geom_ribbon(aes(x=temp,ymin=lb,ymax=ub),alpha=0.5) +
  geom_line(lwd=2) +
  # geom_point(data=incdf, aes(NESTTEMP,INCDAYS)) + 
  labs(x="Nest temperature (C)",y="Incubation duration",title = "") +
  theme_bw()
gginc_temp

png("figs/incubfig.png",5.5,3.5,units="in",res=600)
cowplot::plot_grid(gginc,gginc_temp,labels="auto")
dev.off()


# nest temperature

df <- data.frame(
  value = c(nt.noexcl,nt.excl)*sd(nestdf$NESTTEMP,na.rm=T)+mean(nestdf$NESTTEMP,na.rm=T),
  excl = rep(c("N","Y"),each=nmcmc)
)

ggtemp <- ggplot(df,aes(excl,value))+
  # geom_boxplot(col=gray(0.2),fill="white",size=0.5,width=0.75) +
  geom_violin(col=alpha("black",0.5),fill=alpha("white",0.5),lwd=1.4,draw_quantiles = TRUE) +
  labs(x="Predator excluder?",y="Nest temperature",title = "Nest temperature") +
  theme_bw()
ggtemp

# as a function of year

years <- 2009:2012
years2 <- 1:4

preds1 <- lapply(years2,function(t) sims$nt.base[,t]  )
preds2 <- lapply(years2,function(t) sims$nt.base + sims$nt.excl.eff )

df1 <- data.frame(
  value = do.call("c",preds1),
  year = rep(years,each=nmcmc)
)
df2 <- data.frame(
  value = do.call("c",preds2),
  year = rep(years,each=nmcmc)
)
df1$excl <- "N"
df2$excl <- "Y"
df <- rbind(df1,df2)
df$year <- as.factor(df$year)

ggtemp_year <- ggplot(df,aes(x=year,y=value))+
  geom_violin(aes(fill=excl)) +
  labs(x="Year",y="Nest temperature (C)",title = "") +
  scale_fill_grey(start = 0, end = .9) +
  theme_bw()
ggtemp_year

png("figs/ntfig_year.png",5,3.5,units="in",res=600)
ggtemp_year
dev.off()



# nest predation

df <- data.frame(
  value = c(pred.noexcl,pred.excl),
  excl = rep(c("N","Y"),each=nmcmc)
)

ggpred <- ggplot(df,aes(excl,value))+
  # geom_boxplot(col=gray(0.2),fill="white",size=0.5,width=0.75) +
  geom_violin(col=alpha("black",0.5),fill=alpha("white",0.5),lwd=1.4) +
  labs(x="Predator excluder?",y="Prob. nest pred.",title = "Prob. nest predation")+
  theme_bw()


# eggs remaining after predation

df <- data.frame(
  value = c(neggsrem.noexcl,neggsrem.excl),
  excl = rep(c("N","Y"),each=nmcmc)
)

ggrem <- ggplot(df,aes(excl,value))+
  # geom_boxplot(col=gray(0.2),fill="white",size=0.5,width=0.75) +
  geom_violin(col=alpha("black",0.5),fill=alpha("white",0.5),lwd=1.4) +
  labs(x="Predator excluder?",y="Eggs remaining",title = "Eggs remaining after pred") +
  theme_bw()

# prob of hatching successfully

df <- data.frame(
  value = c(phatch.noexcl,phatch.excl),
  excl = rep(c("N","Y"),each=nmcmc)
)

ggphatch <- ggplot(df,aes(excl,value))+
  # geom_boxplot(col=gray(0.2),fill="white",size=0.5,width=0.75) +
  geom_violin(col=alpha("black",0.5),fill=alpha("white",0.5),lwd=1.4) +
  labs(x="Predator excluder?",y="Prob. hatch",title = "Hatching success") +
  theme_bw()

# number of eggs hatching successfully

df <- data.frame(
  value = c(hatch.noexcl,hatch.excl),
  excl = rep(c("N","Y"),each=nmcmc)
)

gghatch = ggplot(df,aes(excl,value))+
  # geom_boxplot(col=gray(0.2),fill="white",size=0.5,width=0.75) +
  geom_violin(col=alpha("black",0.5),fill=alpha("white",0.5),lwd=1.4) +
  labs(x="Predator excluder?",y="# Eggs",title = "# Eggs Hatching Successfully") +
  theme_bw()

library(cowplot)

png("figs/pred_excl2.png",6,5,units="in",res=600)
plot_grid(ggpred,ggrem,ggphatch,gghatch,byrow=T,labels="auto")
dev.off()

















# models in brms ---------------------------

# # structural equation models?
# names(eggdf)
# pred_mod <- bf(PRED ~ EXCLUD + (1|YEAR) + (1|SITE) )
# succ_mod <- bf(cover ~ firesev)
# 
# k_fit_brms <- brm(rich_mod +
#                     cover_mod +
#                     set_rescor(FALSE),
#                   data=keeley,
#                   cores=4, chains = 2)

# ?bf

# 
# ## model predation as a function of excluded or not, accounting for random effects of year and site
# pred_mod <- bf(PRED ~ EXCLUD + (1|YEAR) + (1|SITE) )
# 
# mod1 <- brm(pred_mod,data=eggdf,family=bernoulli() )
# 
# mod1
# brms::conditional_effects(mod1)
# 
# ## model hatch success as a function of excluded or not
# pred_mod <- bf(HATCH ~ EXCLUD + (1|YEAR) + (1|SITE) )
# 
# mod1 <- brm(pred_mod,data=eggdf,family=bernoulli() )
# 
# mod1
# brms::conditional_effects(mod1)
# 
# 
# ## model unsuccess as a function of excluded or not
# 
# eggdf_nopred$NESTTEMPs <- scale(eggdf_nopred$NESTTEMP)[,1]
# eggdf_nopred$LENGTHs <- scale(eggdf_nopred$LENGTH)[,1]
# 
# 
# # form <- bf(UNSUCC ~ EXCLUD + mi(NESTTEMPs) + mi(LENGTHs) ) +
# #   bf(NESTTEMPs | mi() ~ 1 + gaussian()) + bf(LENGTHs | mi() ~ 1 + gaussian()) + set_rescor(FALSE)
# # mod1 <- brm(form, data = eggdf_nopred,chains = 2, cores = 2)
# 
# pred_mod <- bf(UNSUCC ~ EXCLUD  )  #+ (1|YEAR) + (1|SITE) )  + (1|SITEYEAR)
# mod1 <- brm(pred_mod,data=eggdf_nopred,family=bernoulli(),chains = 2, cores = 2)
# mod1
# brms::conditional_effects(mod1)
# 
# 
# form <- bf(HATCH ~ (1-PRED) * (1-UNSUCC),
#              PRED ~ EXCLUD + (1|YEAR) + (1|SITE), 
#              UNSUCC ~ EXCLUD + (1|YEAR) + (1|SITE), nl=TRUE)
# 
# fit_loss1 <- brm(formula = form, data = eggdf, family = bernoulli(),
#                  chains = 2, cores = 2)
# 
# fit_loss1
# conditional_effects(fit_loss1)

# random brms code -----------------------

fit <-
  brm(data = my_data, 
      family = binomial,
      y | trials(1) ~ 1 + x1 + x2,
      prior = c(prior(normal(0, 2), class = Intercept),
                prior(normal(0, 2), class = b)))

b2 <- brm (Species ~ Petal.Length + Petal.Width + Sepal.Length + Sepal.Width, data=ir,family="categorical", 
    n.chains=3, n.iter=3000, n.warmup=600, prior=c(set_prior ("normal (0, 8)")))


model_chocolate_brms <- brm(
  bf(choice ~ dark + soft + nuts + (1 | subj)),
  data = chocolate,
  family = poisson(),
  prior = c(
    prior(normal(0, 3), class = Intercept),
    prior(normal(0, 3), class = b),
    prior(exponential(1), class = sd)
  ),
  chains = 4, cores = 4, iter = 2000, seed = 1234,
  backend = "cmdstanr", threads = threading(2), refresh = 0,
  file = "models/model_chocolate_brms"
)

pp_check(fit1, resp = "tarsus")

bf_tarsus <- bf(tarsus ~ sex + (1|p|fosternest) + (1|q|dam))
bf_back <- bf(back ~ hatchdate + (1|p|fosternest) + (1|q|dam))
fit2 <- brm(bf_tarsus + bf_back + set_rescor(TRUE), 
            data = BTdata, chains = 2, cores = 2)



conditional_effects(fit_zinb1)

# 2 chains on 2 cores...
fit_rent1 <- brm(rentsqm ~ t2(area, yearc) + (1|district), data = rent99,
                 chains = 2, cores = 2)


  # visualize the conditional smooths
conditional_smooths(fit_rent2)


## combined model formula. Nonlinear model specification followed by the regression models
    # nl = TRUE sets it to treat the non-linear formula literally

nlform <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
             ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1, nl = TRUE)
nlprior <- c(prior(normal(5000, 1000), nlpar = "ult"),
             prior(normal(1, 2), nlpar = "omega"),
             prior(normal(45, 10), nlpar = "theta"))
fit_loss1 <- brm(formula = nlform, data = loss, family = gaussian(),
                 prior = nlprior, control = list(adapt_delta = 0.9))






# Map figure ----------------------

# use Qgis...












