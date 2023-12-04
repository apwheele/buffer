# Analysis of 

library(fitdistrplus)
library(ggplot2)
set.seed(10)

################
# Functions

# axis.title.y=element_text(margin=margin(0,30,0,0))

theme_andy <- function(){
  theme_bw() %+replace% theme(
     text = element_text(size = 16),
     panel.grid.major= element_line(linetype = "longdash"),
     panel.grid.minor= element_blank()
) }


# idea via https://stats.stackexchange.com/questions/6588/how-can-i-estimate-the-density-of-a-zero-inflated-parameter-in-r/6596#6596
# https://www.stata-journal.com/sjpdf.html?articlenum=gr0003


# Note BW here is on "meter" scale
bootdens <- function(x,samples=99,bw=1000,n=2000,from=0.00001){
    # calc original
    y <- log(x)
    lbw <- max(abs(log(bw)),0.1)
    frv <- log(max(0.05,min(x) - 2*bw))
    tov <- log(max(x) + 2*bw)
    g <- density(y,bw=lbw,n=n,from=frv,to=tov)
    xgrid <- exp(g$x)
    g$y <- g$y/xgrid
    g$x <- xgrid
    # calculate bootstrap resamples
    bm <- matrix(0,n,samples)
    for (i in 1:samples){
        by <- sample(y,replace=TRUE)
        bd <- density(by,bw=lbw,n=n,from=frv,to=tov)
        bm[,i] <- bd$y/xgrid
    }
    return(list(data.frame(x=g$x,y=g$y),bm,x))
}

plotdens <- function(bd,title=NULL,quants=c(0.025,0.975),x='Distance',y=NULL){
    df <- bd[[1]]
    df$low <- apply(bd[[2]],1,quantile,probs=quants[1])
    df$hig <- apply(bd[[2]],1,quantile,probs=quants[2])
    # add in gamma estimates
    gamma_fit <- fitdist(bd[[3]],"gamma")
    df$gamd <- dgamma(df$x,gamma_fit$estimate[1],gamma_fit$estimate[2])
    pint <- 100*(quants[2] - quants[1])
    cap <- paste0("Confidence intervals are pointwise ",pint,"% via bootstrap")
    dfr = data.frame(x=bd[[3]])
    p <- ggplot(data=df, aes(x=x,y=y,ymin=low,ymax=hig)) + 
      geom_ribbon(fill = "blue", alpha=0.1) +
      geom_line(size=1.5, color="blue") + 
      geom_line(aes(x=x,y=gamd),color="black",size=1.5) +
      geom_rug(data=dfr, aes(x=x), sides="b", alpha=0.6, inherit.aes = FALSE) +
      labs(x='Distance',y='Density',title=title) +
      theme_andy()
    return(p)
}
################

# First individual, paperbag

d <- c(24.15,
       12.56,
       10.95,
       11.27,
       11.91,
       12.56,
       9.66 ,
       9.98 ,
       8.05 ,
       9.34 ,
       5.8  ,
       7.73 ,
       6.12 ,
       15.78,
       14.81,
       11.91,
       3.22 ,
       2.58 ,
       12.56,
       1.93 ,
       6.44 ,
       1.29 ,
       2.58 ,
       8.05 ,
       8.69 ,
       0.97 ,
       8.69 ,
       1.93 ,
       12.88,
       27.05,
       5.15 ,
       29.3 ,
       27.05,
       28.34,
       4.19 ,
       11.27,
       17.07,
       4.83 ,
       8.69 ,
       8.05 ,
       16.42,
       8.05 ,
       13.52,
       9.34 ,
       10.3 ,
       10.95,
       13.52,
       13.2 ,
       13.85,
       21.57,
       14.17,
       14.17,
       24.47,
       15.78,
       25.12,
       17.07,
       16.42,
       28.01,
       17.07,
       16.74,
       20.29,
       22.86,
       39.61,
       29.95,
       28.98,
       28.66,
       21.25,
       26.73,
       24.47,
       27.37,
       27.05,
       25.76,
       29.3 ,
       26.4 ,
       30.59,
       31.88,
       42.5 ,
       35.1 ,
       51.52)


gamma_fit <- fitdist(d,"gamma")

# > summary(gamma_fit)
# Fitting of the distribution ' gamma ' by maximum likelihood
# Parameters :
#        estimate Std. Error
# shape 2.1031211 0.31169449
# rate  0.1262443 0.02111655

res <- bootdens(d,bw=3)
p <- plotdens(res,quant=c(0.005,0.995),title='Gamma Fit')


# Save graph
png(file = "JTC_01V2.png", bg = "transparent", height=5, width=9, units="in", res=1000)
p
dev.off()


# second individual, serial exposer
d2 <- c(0.37,
        0.92,
        1.19,
        1.44,
        1.49,
        1.52,
        1.61,
        1.65,
        1.65,
        1.70,
        1.81,
        1.91,
        2.05,
        2.20,
        2.20,
        2.20,
        2.27,
        2.27,
        2.35,
        2.36,
        2.56,
        2.57,
        2.64,
        2.65,
        2.93,
        3.13,
        3.16,
        3.16,
        3.25,
        3.41,
        3.41,
        3.42,
        3.43,
        3.47,
        3.50,
        3.51,
        3.76,
        3.76,
        3.82,
        3.83,
        3.91,
        4.16,
        4.16,
        4.29,
        4.36,
        4.41,
        4.54,
        4.93,
        4.94,
        4.94,
        5.15,
        5.38,
        5.95,
        6.10,
        6.24,
        6.24,
        6.29,
        6.33,
        6.55,
        6.67,
        6.80,
        6.91,
        7.03,
        7.04,
        7.10,
        7.13,
        7.17,
        7.17,
        7.22,
        7.23,
        7.30,
        7.32,
        7.42,
        7.55,
        7.57,
        7.95,
        8.01,
        8.03,
        8.06,
        8.20,
        8.34,
        8.52,
        8.67,
        8.76,
        8.82,
        8.82,
        8.92,
        9.21,
        9.95,
        10.52,
        10.87,
        10.89,
        11.18,
        11.36,
        11.46,
        11.47,
        11.47,
        11.54,
        11.57,
        12.93,
        13.04,
        13.16,
        13.17,
        13.34,
        13.82,
        13.90,
        13.91,
        13.91,
        14.80,
        14.98,
        14.98,
        15.07,
        15.18,
        15.18,
        15.31,
        15.46,
        15.55,
        15.71,
        15.98,
        16.54,
        16.56,
        16.99,
        17.08,
        17.46,
        18.10,
        19.96,
        20.31,
        20.52,
        22.65)

res2 <- bootdens(d2,bw=2)
p <- plotdens(res2, quant=c(0.005,0.995)) # 99% CIs

# Save graph
png(file = "JTC_02.png", bg = "transparent", height=5, width=9, units="in", res=1000)
p
dev.off()

gamma_fit2 <- fitdist(d2,"gamma")

# > summary(gamma_fit2)
#Fitting of the distribution ' gamma ' by maximum likelihood
#Parameters :
#       estimate Std. Error
#shape 2.0415169  0.2363680
#rate  0.2576694  0.0337945 

length(d) #79
length(d2) #129


# Figuring out max according to bootstrap

mc <- apply(res2[[2]],2,max)
dr <- res2[[2]]
x <- 0

for (i in 1:length(mc)) {
    max_row <- res2[[1]]$x[dr[,i] == mc[i]]
    x <- x + max_row
}

x <- x/length(mc) # [1] 1.813224


mc <- apply(res[[2]],2,max)
dr <- res[[2]]
x <- 0

for (i in 1:length(mc)) {
    max_row <- res[[1]]$x[dr[,i] == mc[i]]
    x <- x + max_row
}

x <- x/length(mc) # [1] 1.739261


# test my function with real gamma data

set.seed(10)
n <- 50
d3 <- rgamma(n,0.7)

res3 <- bootdens(d3,bw=2)
p1 <- plotdens(res3,title='Simulated Gamma with Shape=0.7')

d4 <- rgamma(n,1.5)
res4 <- bootdens(d4,bw=2)
p2 <- plotdens(res4,title='Simulated Gamma with Shape=1.5')

png(file = "JTC_SimDensity.png", bg = "transparent", height=5, width=9, units="in", res=1000)
p1 + p2
dev.off()
