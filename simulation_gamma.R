#############################
# Start of simulation
#############################


library(fitdistrplus)
library(ggplot2)
set.seed(10)

###############################################
# Lets make a nice graph gamma 0.5, 0.75, 1.0, 1.25, 1.5

theme_andy <- function(){
  theme_gray() %+replace% theme(
     panel.grid.major= element_line(linetype = "longdash"),
     panel.grid.minor= element_blank(),
     axis.title.y=element_text(margin=margin(0,30,0,0))
) }

x <- seq(0.01,2,length.out=500)
s <- seq(0.8,1.6,length.out=5)
df <- data.frame(expand.grid(x=x,s=s))
df$d <- dgamma(df$x,df$s)
sl <- paste0(s)
sl[2] <- "1.0"
df$s <- factor(df$s,levels=s, labels=sl)

p <- ggplot(data=df, aes(x=x,y=d,color=s)) + 
     geom_line(size=1.2) + scale_color_discrete(name = "Shape") +
     labs(x=NULL,y=NULL,title='Density of Gamma Distribution') +
     theme_andy()

png(file = "GammaDist.png", bg = "transparent", height=5, width=9, units="in", res=1000)
p
dev.off()

###############################################

###############################################
# Lets start with a simple scenario
# sample sizes of 5, 10, 20
# estimate distribution

shape <- 0.5
rate <- 2

d5 <- rgamma(5,shape,rate)
summary(fitdist(d5,"gamma"))

d10 <- rgamma(10,shape,rate)
summary(fitdist(d10,"gamma"))

d20 <- rgamma(20,shape,rate)
summary(fitdist(d20,"gamma"))

###############################################

# Function to create graph

# Function to estimate bias/variance
# given sample size for gamma

## try/catch if error
#out <- tryCatch(
#   { # part with the potential error
#     #r <- ???? #whatever code steps you want
#    r
#   }, error=function(cond){ 
#      print("Function Failed, Error message is \n\n")
#      print(Cond)
#      return(-1)
#      } )
#return(out)

simg <- function(n,shape,rate,method="mle"){
    # if failed, return NULL
    out <- tryCatch({
        disg <- rgamma(n,shape,rate)
        fitg <- fitdist(disg, "gamma", method=method)
        if (method=="mle") {
        se_shape <- fitg$sd[1]
        } else {
        se_shape <- NA
        }
        bias_shape <- fitg$estimate[1] - shape
        bias_rate <- fitg$estimate[2] - rate
        r <- c(bias_shape,bias_rate,se_shape)
    }, error=function(cond){
        return(c(NA,NA,NA))
    })
    return(out)
}

simr <- function(rep,n,shape,rate,method="mle"){
    res <- replicate(rep,simg(n,shape,rate,method))
    shape_res <- res[1,]
    rate_res <- res[2,]
    shapese_res <- res[3,]
    mean_shapebias <- mean(shape_res,na.rm = TRUE)
    std_shapebias <- sd(shape_res,na.rm = TRUE)
    mean_ratebias <- mean(rate_res,na.rm = TRUE)
    std_ratebias <- sd(rate_res,na.rm = TRUE)
    mean_shapese <- mean(shapese_res,na.rm = TRUE)
    tot_sim <- sum(!is.na(shape_res))
    return(c(mean_shapebias,std_shapebias,mean_ratebias,std_ratebias,mean_shapese,tot_sim))
}

# Now do the actual simulation
reps <- 10000  # when I do the real analysis, will do something like 10k
samp_sizes <- c(10,20,30,40,50,60,70,80,90,100) #5:100
shape <- c(0.6,0.8,1.0,1.2,1.4) #seq(0.5,1.5,0.1)
rate <- 1 # rate doesn't really matter, will change depending on unit of distance

sim_full <- expand.grid(Reps=reps,Rate=rate,SampSize=samp_sizes,Shape=shape)
#names(sim_full) <- c('Reps','SampSize','Shape','Rate')

b <- Sys.time()
print(b)
sim_res <- mapply(simr,sim_full$Reps,sim_full$SampSize,sim_full$Shape,sim_full$Rate)
e <- Sys.time()
print(e)
print(e - b) # 15 minutes

sim_res <- data.frame(t(sim_res))
names(sim_res) <- c('ShapeBias','ShapeStd','RateBias','StdRate','ShapeSE','TotSim')
sim_full$ShapeSE <- sim_res$ShapeSE
sim_full$ShapeBias <- sim_res$ShapeBias
sim_full$TotSim <- sim_res$TotSim

write.csv(sim_full,"SimulationTable.csv",row.names=FALSE)

# Same table, but using method of moments fit

b <- Sys.time()
print(b)
sim_res2 <- mapply(simr,sim_full$Reps,sim_full$SampSize,sim_full$Shape,sim_full$Rate,"mme")
e <- Sys.time()
print(e)
print(e - b) # 15 minutes

sim_res2 <- data.frame(t(sim_res2))
names(sim_res2) <- c('ShapeBias','ShapeStd','RateBias','StdRate','ShapeSE','TotSim')
#sim_full$ShapeSE <- sim_res2$ShapeSE
sim_full$ShapeBias <- sim_res2$ShapeBias
sim_full$TotSim <- sim_res2$TotSim

write.csv(sim_full,"SimulationTableMoM.csv",row.names=FALSE)


#############################

# Graph for exact SE

seG <- function(samp,shape){
    vr <- shape/(samp*(trigamma(shape)*shape - 1))
    return(sqrt(vr))
}

samp <- 10:100
shapes <- c(0.5,1.0,1.5,2.0)
se_df <- data.frame(expand.grid(Ssize=samp,Shape=shapes))
se_df$SE <- seG(se_df$Ssize,se_df$Shape)
sl <- paste0(shapes)
sl[2] <- "1.0"
sl[4] <- "2.0"
se_df$Shape <- factor(se_df$Shape,labels=sl)

p <- ggplot(data=se_df, aes(x=Ssize,y=SE,color=Shape)) + 
     geom_line(size=1.2) + scale_color_discrete(name = "Shape") +
     labs(x='Sample Size',y=NULL,title='Standard Error of Gamma Shape') +
     scale_y_continuous(breaks=seq(0.1,0.9,0.1)) + 
     scale_x_continuous(breaks=seq(10,100,10)) +
     theme_andy()

png(file = "GammaDist_SE.png", bg = "transparent", height=5, width=9, units="in", res=1000)
p
dev.off()


#############################
# Notes

# For non-parametric KDE, just use gamma kernel
#https://www.ism.ac.jp/editsec/aism/pdf/052_3_0471.pdf
#https://cran.r-project.org/web/packages/DELTD/DELTD.pdf

# Alt ways to transform KDE estimates to obey bounds
# Log would be easy
#https://thirdorderscientist.org/homoclinic-orbit/2013/10/24/kernel-density-estimation-for-random-variables-with-bounded-support-mdash-the-transformation-trick
#https://andrewpwheeler.com/2015/07/20/transforming-kde-estimates-from-logistic-to-probability-scale-in-r/

# To get error intervals for KDE, can use bootstrap
#https://stats.stackexchange.com/a/207455/1036

# Maybe check out
# https://www.tandfonline.com/doi/full/10.1080/00949655.2023.2173194
# for bandwidth selection

#############################

library(fitdistrplus)
set.seed(10)

x <- c(1.5,0.4,0.8,1.9,0.3,2.5)
gamma_fit <- fitdist(x,"gamma")
summary(gamma_fit)