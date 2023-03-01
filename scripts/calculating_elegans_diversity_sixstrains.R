#theta and pi from full population 
ce.fullpop.df <- read.csv("ce_full_population/sub_region_diversity.csv", header = TRUE, sep = ",")

#theta and pi from six strains in anthelmintic DR 
ce.anthDR.df <- read.csv("ce_dose_response/sub_region_diversity.csv", header = TRUE, sep = ",")

#subset to retain "full" subregions 
ce.fullpop.df.subset <- subset(ce.fullpop.df, ce.fullpop.df$sub_region == "full")


ce.anthDR.df.subset <- subset(ce.anthDR.df, ce.anthDR.df$sub_region == "full")

#sum theta in 6 strains for DR and full elegans strains
sum.anthDR <- sum(ce.anthDR.df.subset$theta)
sum.full <- sum(ce.fullpop.df.subset$theta)

#theta calculation
theta.calc <- (sum.anthDR/sum.full)
theta.calc

#mean pi 
ce.fullpop.pi.mean <- mean(ce.fullpop.df$pi)
ce.anthDR.pi.mean <- mean(ce.anthDR.df$pi)

pi.calc <- (ce.anthDR.pi.mean/ce.fullpop.pi.mean)