setwd("/Users/Niladri/Documents/Research/Permutation/metrics-paper")

library(nullabor)
library(plyr)
library(ggplot2)

####Modified Binned Distance: No indexing

bdist_mod <- function(X,PX, nbin.X = 5, nbin.Y = 5) {
	if(!is.numeric(X[,1])){
	X[,1] <- as.numeric(X[,1])
	nij <- as.numeric(table(cut(X[,1], breaks=seq(min(X[,1]), max(X[,1]),length.out = length(unique(X[,1])) + 1), include.lowest = TRUE),cut(X[,2], breaks=seq(min(lineup.dat[,2]), max(lineup.dat[,2]),length.out = nbin.Y + 1), include.lowest = TRUE)))
	}else
		nij <- as.numeric(table(cut(X[,1], breaks=seq(min(lineup.dat[,1]), max(lineup.dat[,1]),length.out = nbin.X + 1), include.lowest = TRUE),cut(X[,2], breaks=seq(min(lineup.dat[,2]), max(lineup.dat[,2]),length.out = nbin.Y + 1), include.lowest = TRUE)))
	if(!is.numeric(PX[,1])){
	PX[,1] <- as.numeric(PX[,1])
	mij <- as.numeric(table(cut(PX[,1], breaks=seq(min(X[,1]), max(X[,1]),length.out = length(unique(X[,1])) + 1), include.lowest = TRUE),cut(PX[,2], breaks=seq(min(lineup.dat[,2]), max(lineup.dat[,2]),length.out = nbin.Y + 1), include.lowest = TRUE)))
	}else
	mij <- as.numeric(table(cut(PX[,1], breaks=seq(min(lineup.dat[,1]), max(lineup.dat[,1]),length.out = nbin.X + 1), include.lowest = TRUE),cut(PX[,2], breaks=seq(min(lineup.dat[,2]), max(lineup.dat[,2]),length.out = nbin.Y + 1), include.lowest = TRUE)))
	sqrt(sum((nij-mij)^2))
}


diff_bin <- function(lineup.dat, m, n){
	dat.bin <- NULL
	for(i in 1:20){
			bin.j <- NULL
		for(j in 1:20){
			X <- lineup.dat[lineup.dat$.sample == i, ]
			PX <- lineup.dat[lineup.dat$.sample == j, ]
			bin <- bdist_mod(X, PX, nbin.X = m, nbin.Y = n)
			bin.dat1 <- data.frame(p = m, q = n, pos.1 = i, pos.2 = j, bin = bin)
			bin.j <- rbind(bin.j, bin.dat1)
						}
				dat.bin <- rbind(dat.bin, bin.j)
				}
	dat.bin <- subset(dat.bin, pos.1 != pos.2 & pos.2 != pos)
	dat.bin.mean <- ddply(dat.bin, .(pos.1, p, q), summarize, bin.m = mean(bin), len = length(bin))
	diff <- with(dat.bin.mean, bin.m[len == 19] - max(bin.m[len != 19]))
return(diff)
dat.bin
}

calc_diff <- function(lineup.dat, m, n, pos){
	d <- sapply(1:20, function(i){
			X <- lineup.dat[lineup.dat$.sample == i, ]
		sapply(1:20, function(j){
			PX <- lineup.dat[lineup.dat$.sample == j, ]
			dis <- bdist_mod(X, PX, nbin.X = m, nbin.Y = n)
		})
	})
	require(reshape)
	d.m <- melt(d)
	names(d.m) <- c("pos.2", "plotno", "bin")
	dat.bin <- subset(d.m, plotno != pos.2 & pos.2 != pos)
	require(plyr)
	dat.bin.mean <- ddply(dat.bin, .(plotno), summarize, bin.m = mean(bin), len = length(bin))
 with(dat.bin.mean, bin.m[len == 19] - max(bin.m[len != 19]))
}


bin_diff <- function(lineup.dat, xlow, xhigh, ylow, yhigh, pos, plot = FALSE){
		d <- sapply(xlow:xhigh, function(m){
		sapply(ylow:yhigh, function(n){
		calc_diff(lineup.dat, m, n, pos)
		})
	})
	d.m <- melt(d)
	names(d.m) <- c("q", "p", "Diff")
	d.m$p <- d.m$p + xlow - 1
	d.m$q <- d.m$q + ylow - 1
	d.m <- data.frame(p = d.m$p, q = d.m$q, Diff = d.m$Diff)
	if(plot){
		require(ggplot2)
		p <- ggplot(d.m, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill
		 	=   Diff)) + scale_fill_gradient(high ="blue", low ="white") + 
		 	xlab("p") + ylab("q")
		return(list(dat = d.m, p = p))
	}else{
	return(dat = d.m)
	}
}

pt <- proc.time()
bin_diff(lineup.dat, 2, 10, 2, 10, pos, plot = TRUE)
proc.time() - pt
	
### Data 1 -- Cluster Data

set.seed(2000)
X1 <- data.frame(x = rnorm(100), y = rnorm(100))
X1[1:33, 1] <- X1[1:33, 1] - 5
X1[34:66, 1] <- X1[34:66, 1] + 5
X1[67:100, 2] <- X1[67:100, 2] + sqrt(75)

dat <- X1

qplot(x,y, data = dat, geom = "point", size = I(5), xlab = "X", ylab = "Y")

ggsave("data1.pdf", height = 4, width = 4)

pos <- sample(20, 1)
	
lineup.dat <- lineup(null_permute(names(dat)[1]), true = dat, pos = pos)

qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)
	
ggsave("bin-select-lineup1.pdf", height = 4, width = 4)
	
null <- subset(lineup.dat, .sample != pos & .sample == sample(.sample, 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
ggsave("null1.pdf", height = 4, width = 4)
	

ptm <- proc.time()
p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}
proc.time() - ptm

ptm1 <- proc.time()
D <- bin_diff(lineup.dat, 8, 11, 7, 12, pos)
proc.time() - ptm1


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("bin-select-plot1.pdf", height = 4, width = 5)

max.min <- Diff[c(which.min(Diff$diff), which.max(Diff$diff) ), ]
max.min

### Data 2 -- Categorical Variables

set.seed(2000)
y <- rep(1:10, each = 10)
mu <- sample(10, 10, replace = T)
mu <- sort(mu)
x <- NULL
for(i in mu){
xx <- rpois(10, i)
x <- c(x, xx)
}

X3 <- data.frame(x, y)

dat <- X3

qplot(x,y, data = dat, geom = "point", size = I(5), xlab = "X", ylab = "Y")

#ggsave("data2.pdf", height = 4, width = 4)

	pos <- sample(20, 1)
	
	sample.dat <- NULL
	for(i in 1:19){
		x <- rpois(100, 10)
		y <- NULL
		for(j in 1:100){
		y[j] <- 1.6643 + 0.5803*x[j] + rnorm(1)
		}
		new.dat <- data.frame(x = x, y = y, .n = i)
		sample.dat <- rbind(sample.dat, new.dat)
	}
	lineup.dat <- lineup(true = dat, samples = sample.dat, pos = pos)
	
	qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)
	
#	ggsave("bin-select-lineup2.pdf", height = 4, width = 4)
	
null <- subset(lineup.dat, .sample != pos & .sample == sample(.sample, 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
#ggsave("null2.pdf", height = 4, width = 4)
	
	
	
p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

#ggsave("bin-select-plot2.pdf", height = 4, width = 5)

max.min <- Diff[c(which.min(Diff$diff), which.max(Diff$diff) ), ]
max.min

### Data 3 -- Nonlinear relation with outliers

set.seed(2000)
x <- runif(100, 1, 10)
x <- sort(x, decreasing = F)
y <- rexp(100, 10)/(x^2)

xx <- rep(max(x), 5)
yy <- rep(max(y), 5)

x <- c(x, xx + rnorm(5, 0, 0.2))
y <- c(y, yy + rnorm(5, 0, 0.002))

dat3 <- data.frame(x = x, y = y)

qplot(x, y, data = dat3, size = I(5), xlab = "X", ylab = "Y")

ggsave("data3.pdf", height = 4, width = 4)


	pos <- sample(20, 1)
	
	sample.dat <- NULL
	for(i in 1:19){
		x <- rpois(100, 10)
		y <- NULL
		for(j in 1:100){
		y[j] <- 1.6643 + 0.5803*x[j] + rnorm(1)
		}
		new.dat <- data.frame(x = x, y = y, .n = i)
		sample.dat <- rbind(sample.dat, new.dat)
	}
	lineup.dat <- lineup(true = dat, samples = sample.dat, pos = pos)
	lineup.dat <- lineup(null_permute(names(dat3)[1]), true = dat3, pos = pos)
	qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)
	
	ggsave("bin-select-lineup3.pdf", height = 4, width = 4)
	
null <- subset(lineup.dat, .sample != pos & .sample == sample(.sample, 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
ggsave("null3.pdf", height = 4, width = 4)
	
	
	
p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("bin-select-plot3.pdf", height = 4, width = 5)

max.min <- Diff[c(which.min(Diff$diff), which.max(Diff$diff) ), ]
max.min


#### Data 4 -- Linear Relationship with Outliers

set.seed(2011)
x <- rexp(100, 2)
y <- rexp(100, 2)

xx <- 3
yy <- 3

x <- c(x, xx)
y <- c(y, yy)

dat <- data.frame(x = x, y = y)

qplot(x, y, data = dat, size = I(5), xlab = "X", ylab = "Y")

ggsave("data5.pdf", height = 4, width = 4) 

mod <- lm(y ~ x)

pos <- sample(20, 1)
lineup.dat <- lineup(null_permute(names(dat)[1]), true = dat, pos = pos)
qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)

ggsave("bin-select-lineup5.pdf", height = 4, width = 4)

null <- subset(lineup.dat, .sample == 14)

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
ggsave("null5.pdf", height = 4, width = 4)


p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("bin-select-plot5.pdf", height = 4, width = 5)

max.min <- Diff[c(which.min(Diff$diff), which.max(Diff$diff) ), ]
max.min

### Data 5 -- Curved Pattern in the residual plot

set.seed(2000)
x <- rnorm(100, 5, 2)
y <- NULL
for(i in 1:100){
y[i] <- 34 + 80*x[i] + 100*(x[i])^2 + rnorm(1, 0, 150)
}
dat <- data.frame(x = x, y = y)

#qplot(x, y , data = dat, geom = "point", size = I(3), xlab = "X", ylab = "Y")

mod <- lm(y ~ x, data = dat)

dat4 <- data.frame(x = x, y = mod$resid)

qplot(x, y, data = dat4, geom = "point", size = I(5), xlab = "X", ylab = "Y")

ggsave("data4.pdf", height = 4, width = 4)

sample.dat <- NULL
for(k in 1:19){
	y.samp <- NULL
	for(i in 1:100){
		y.samp[i] <- rnorm(1, mod$coef[1] + mod$coef[2]*x[i], sd(mod$resid))
	}
	mod.samp <- lm(y.samp ~ x)
	dat <- data.frame(x = x, y = mod.samp$resid, .n = k)
	sample.dat <- rbind(sample.dat, dat)
}
pos <- sample(20, 1)
lineup.dat <- lineup(true = dat4, samples = sample.dat, pos = pos)
qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)
	
ggsave("bin-select-lineup4.pdf", height = 4, width = 4)

null <- subset(lineup.dat, .sample != pos & .sample == sample(.sample, 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
ggsave("null4.pdf", height = 4, width = 4)
	
p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("bin-select-plot4.pdf", height = 4, width = 5)

max.min <- Diff[c(which.min(Diff$diff), which.max(Diff$diff) ), ]
max.min

### Data 6 -- Residual plot for Nonconstant variance

set.seed(2005)
x<-rnorm(100,5,2)
y <- NULL
for(i in 1:100){
y[i] <-rnorm(1,34+80*x[i],10 + 30*x[i])
}

dat <- data.frame(x = x, y = y)
qplot(x,y,data = dat, xlab="X",ylab="Y", size = I(3))
#lines(x,mod4$fitted)

mod4<-lm(y~x, data = dat)
summary(mod4)

dat6 <- data.frame(x = x, y = mod4$resid)
qplot(x,y,data = dat6, xlab="X",ylab="Y", size = I(5))

ggsave("data6.pdf", height = 4, width = 4) 

sample.dat <- NULL
for(k in 1:19){
y.samp <- NULL
for(i in 1:100){
y.samp[i] <- rnorm(1, mod4$coef[1] + mod4$coef[2]*x[i], sd(mod4$resid))
}
mod <- lm(y.samp ~ x)
samp.dat <- data.frame(x = x, y = mod$resid, .n = k)
sample.dat <- rbind(sample.dat, samp.dat)
}

pos <-  sample(20, 1)
lineup.dat <- lineup(true = dat6, samples = sample.dat, pos = pos)
qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)

ggsave("bin-select-lineup6.pdf", height = 4, width = 4)

null <- subset(lineup.dat, .sample == sample(.sample[.sample != pos], 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
ggsave("null6.pdf", height = 4, width = 4)


p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("bin-select-plot6.pdf", height = 4, width = 5)

max.min <- Diff[c(which.min(Diff$diff), which.max(Diff$diff) ), ]
max.min


### Data 7 - Spiral data

set.seed(2000)
f.spiral <- function(n=100, mx=8*pi, a=0, b=1) {
  theta <- runif(n, min=0, max=mx)
  radius <- a + b * theta
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  return(df = data.frame(x, y))
}

s <- f.spiral(mx=10*pi, n=500)

qplot(x, y, data = s, geom = "point", size = I(5), xlab = "X", ylab = "Y")

ggsave("data7.pdf",  height = 4, width = 4)

pos <- sample(20, 1)
lineup.dat <- lineup(null_permute(names(s)[1]), true = s, pos = pos)
qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)

ggsave("bin-select-lineup7.pdf", height = 4, width = 4)

null <- subset(lineup.dat, .sample == sample(.sample[.sample != pos], 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")

ggsave("null7.pdf", height = 4, width = 4)

p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("bin-select-plot7.pdf", height = 4, width = 5)

max.min <- Diff[c(which.min(Diff$diff), which.max(Diff$diff) ), ]
max.min

### Data 8 - Data with contamination

set.seed(2000)
x1 <-  rnorm(100)
y1 <- NULL
for(i in 1:100){
	y1[i] <- 5 + 5*x1[i] + rnorm(1, 0, 5)
}
x2 <- rnorm(15, -1.75, 1/3)
y2 <- NULL
for(i in 1:15){
	y2[i] <- 15 + rnorm(1, 0, 5/3)
}

x <- c(x1, x2)
y <- c(y1, y2)

dat <- data.frame(x = x, y = y)

qplot(x, y, data = dat, geom = "point", size = I(5), xlab = "X", ylab = "Y")

ggsave("data8.pdf", height = 4, width = 4)

mod <- lm(y ~ 1)

sample.dat <- NULL
for(k in 1:19){
	y.samp <- rnorm(115, mod$coef, sd(mod$resid))
	samp <- data.frame(x = x, y = y.samp, .n = k)
	sample.dat <- rbind(sample.dat, samp)
}

pos <- sample(20, 1)

lineup.dat <- lineup(samples = sample.dat, true = dat, pos = pos)
qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)

ggsave("bin-select-lineup8.pdf", height = 4, width = 4)

null <- subset(lineup.dat, .sample == sample(.sample[.sample != pos], 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")

ggsave("null8.pdf", height = 4, width = 4)

p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("bin-select-plot8.pdf", height = 4, width = 5)

max.min <- Diff[c(which.min(Diff$diff), which.max(Diff$diff) ), ]
max.min

#=========================================================================================
# Francis Anscombe's quartet dataset
#=========================================================================================

###Dataset 1

x <- c(10.0, 8.0, 13.0, 9.0, 11.0, 14.0, 6.0, 4.0, 12.0, 7.0, 5.0)
y <- c(8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68)

dat1 <- data.frame(x = x, y = y)

qplot(x, y, data = dat1, size = I(3), xlim = c(2, 19.5), ylim = c(2, 13)) + geom_smooth(method = "lm", se = FALSE)

## Simulated Dataset 1

set.seed(2010)
#set.seed(2000)
sim.x <- rnorm(50, 5)
#sim.x <- rpois(20, 9)
sim.y <- NULL
for(i in 1:50){
sim.y[i] <- 3.0 + 0.7*sim.x[i] + rnorm(1)
}

sim.dat1 <- data.frame(x = sim.x, y = sim.y)

qplot(x, y, data = sim.dat1, size = I(5), xlab = "X", ylab = "Y") 
ggsave("anscombe-1.pdf", height = 4, width = 4)

pos <- sample(20, 1)
lineup.dat <- lineup(null_permute(names(sim.dat1)[1]), true = sim.dat1, pos = pos)
qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)

ggsave("anscombe-lineup-1.pdf", height = 4, width = 4)

null <- subset(lineup.dat, .sample == sample(.sample[.sample != pos], 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
ggsave("anscombe-null-1.pdf", height = 4, width = 4)


p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("anscombe-nbin-1.pdf", height = 4, width = 5)

max.min <- rbind(Diff[which.min(Diff$diff),], Diff[which.max(Diff$diff),])
max.min

###Dataset 2

x <- c(10.0, 8.0, 13.0, 9.0, 11.0, 14.0, 6.0, 4.0, 12.0, 7.0, 5.0)
y <- c(9.14, 8.14, 8.74, 8.77, 9.26, 8.10, 6.13, 3.10, 9.13, 7.26, 4.74)

dat2 <- data.frame(x = x, y = y)
qplot(x, y, data = dat2, size = I(3), xlim = c(2, 19.5), ylim = c(2, 13)) + geom_smooth(method = "lm", se = FALSE)

## Simulated Dataset 2

set.seed(2010)
#sim.x <- rpois(20, 9)
sim.x <- rnorm(50, 5)
sim.y <- NULL
for(i in 1:50){
sim.y[i] <- -5.9957 + 1.4808*sim.x[i] - 0.1267*sim.x[i]*sim.x[i]
}

sim.dat2 <- data.frame(x = sim.x, y = sim.y)

qplot(x, y, data = sim.dat2, size = I(5), xlab = "X", ylab = "Y")

ggsave("anscombe-2.pdf", height = 4, width = 4)

pos <- sample(20, 1)
lineup.dat <- lineup(null_permute(names(sim.dat2)[1]), true = sim.dat2, pos = pos)
qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)

ggsave("anscombe-lineup-2.pdf", height = 4, width = 4)

null <- subset(lineup.dat, .sample == sample(.sample[.sample != pos], 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
ggsave("anscombe-null-2.pdf", height = 4, width = 4)


p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("anscombe-nbin-2.pdf", height = 4, width = 5)

max.min <- rbind(Diff[which.min(Diff$diff),], Diff[which.max(Diff$diff),])
max.min

### Dataset3

x <- c(10.0, 8.0, 13.0, 9.0, 11.0, 14.0, 6.0, 4.0, 12.0, 7.0, 5.0)
y <- c(7.46, 6.77, 12.74, 7.11, 7.81, 8.84, 6.08, 5.39, 8.15, 6.42, 5.73)

dat3 <- data.frame(x = x, y = y)
qplot(x, y, data = dat3, size = I(3), xlim = c(2, 19.5), ylim = c(2, 13)) + geom_smooth(method = "lm", se = FALSE)

### Simulated Dataset 3

set.seed(2010)
#sim.x <- rpois(20, 9)
sim.x <- rnorm(50, 5)
sim.xx <- sim.x[-which.max(sim.x)]

sim.yy <- NULL
for(i in 1:49){
sim.yy[i] <- 4.0056 + 0.5454*sim.xx[i] 
}

sim.x <- c(sim.xx, sim.x[which.max(sim.x)])
sim.y <- c(sim.yy, 12.74)
sim.dat3 <- data.frame(x = sim.x, y = sim.y)
qplot(x, y, data = sim.dat3, size = I(5), xlab = "X", ylab = "Y") 
ggsave("anscombe-3.pdf", height = 4, width = 4)

pos <- sample(20, 1)
lineup.dat <- lineup(null_permute(names(sim.dat3)[1]), true = sim.dat3, pos = pos)
qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)

ggsave("anscombe-lineup-3.pdf", height = 4, width = 4)

null <- subset(lineup.dat, .sample == sample(.sample[.sample != pos], 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
ggsave("anscombe-null-3.pdf", height = 4, width = 4)


p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("anscombe-nbin-3.pdf", height = 4, width = 5)

max.min <- rbind(Diff[which.min(Diff$diff),], Diff[which.max(Diff$diff),])
max.min

###Dataset 4

x <- c(8,8,8,8,8,8,8,19,8,8,8)
y <- c(6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.50, 5.56, 7.91, 6.89)

dat4 <- data.frame(x = x, y = y)
qplot(x, y, data = dat4, size = I(3), xlim = c(2, 19.5), ylim = c(2, 13)) + geom_smooth(method = "lm", se = FALSE)

### Simulated Dataset 4

set.seed(2010)
#sim.x <- c(rep(8,40), rep(19,10))
sim.x <- c(rep(8,49), 19)
sim.y <- NULL
for(i in 1:49){
sim.y[i] <- 3 + 0.5*sim.x[i] + rnorm(1, 0, 1.4)
}
sim.y[50] <- 3 + 0.5*sim.x[50]
#for(i in 41:50){
#sim.y[i] <- rnorm(1, 12.5, 0.25) 
#}

sim.dat4 <- data.frame(x = sim.x, y = sim.y)
qplot(x, y, data = sim.dat4, size = I(5), xlab = "X", ylab = "Y")

ggsave("anscombe-4.pdf", height = 4, width = 4)


sim.x <- rpois(50, 9)
sample.dat <- NULL
for(k in 1:19){
	sim.y <- NULL
for(i in 1:50){
	sim.y[i] <- 3.0 + 0.5*sim.x[i] + rnorm(1)
}
dat <- data.frame(x = sim.x, y = sim.y, .n = k)
sample.dat <- rbind(sample.dat, dat)
}

pos <- sample(20, 1)
lineup.dat <- lineup(true = sim.dat4, samples = sample.dat, pos = pos)
qplot(x,y, data = lineup.dat, size = I(2.5), alpha = I(0.6), xlab = "X", ylab = "Y") + facet_wrap(~.sample)

ggsave("anscombe-lineup-4.pdf", height = 4, width = 4)

null <- subset(lineup.dat, .sample == sample(.sample[.sample != pos], 1))

qplot(x, y, data = null, geom = "point", size = I(5), xlab = "X", ylab = "Y")
	
ggsave("anscombe-null-4.pdf", height = 4, width = 4)


p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("anscombe-nbin-4.pdf", height = 4, width = 5)

max.min <- rbind(Diff[which.min(Diff$diff),], Diff[which.max(Diff$diff),])
max.min


#===========================############===============================================

##====================== Not in the paper==============================

## Dataset 1

X1 <- data.frame(x = rnorm(100), y = rnorm(100))
X1[1:33, 1] <- X1[1:33, 1] - 5
X1[34:66, 1] <- X1[34:66, 1] + 5
X1[67:100, 2] <- X1[67:100, 2] + sqrt(75)

library(ggplot2)
qplot(x,y, data = X1, geom = "point", size = I(5), xlab = "X", ylab = "Y")

ggsave("data1.pdf", height = 4, width = 4)

## Dataset 2

rho <- 0.6
z <- rnorm(50)
zz <- NULL
for(i in 1:length(z)){
zz[i] <- rnorm(1, 3 + rho*2*z[i], sqrt(1 - rho^2)*2)
}
X2 <- data.frame(x = z, y = zz)

qplot(x,y, data = X2, geom = "point", size = I(3), xlab = "X", ylab = "Y")

ggsave("data2.pdf", height = 4, width = 4)

## Dataset 3

#x <- rep(1:20, sample(5, 20, replace = TRUE))

x <- rep(1:20, each = 10)
y <- rnorm(length(x))

X3 <- data.frame(x, y)

qplot(x,y, data = X3, geom = "point", size = I(2.5), xlab = "X", ylab = "Y")

ggsave("data3.pdf", height = 4, width = 4)


bdist <- function(X,PX, nbins=10) {
	nij <- as.numeric(table(cut(X[,1], breaks=nbins),cut(X[,2], breaks=nbins)))
	mij <- as.numeric(table(cut(PX[,1], breaks=nbins),cut(PX[,2], breaks=nbins)))
	sqrt(sum((nij-mij)^2))
}

#qplot(x,z,data = X)

bin.d <- NULL
bin.dat1 <- NULL

for(m in 1:500){

#rho <- runif(1)	
#p <- rnorm(length(X2$x), 5 + rho*2*X2$x, sqrt((1 - rho^2)*4))

Y <- data.frame(X3$x, sample(X3$y))	

for(l in 2:20){

bin.dist <- bdist(X3, Y, nbins = l)

bin.dat <- data.frame(nbin = l, b = bin.dist) 

bin.dat1 <- rbind(bin.dat1, bin.dat)

}

bin.d <- rbind(bin.d, bin.dat1)

	
}

qplot(factor(nbin), b, data = bin.d, geom = "boxplot", xlab = "number of bins", ylab = "Binned distance")

ggsave("data3-nbin.pdf", height = 4, width = 7)

#ggsave("nbin-100.pdf", height = 3.5, width = 8.5)
#qplot(b, geom = "histogram", binwidth = 0.25, data = bin.d, facets = nbin ~ .)

summary.b <- ddply(bin.d, .(nbin), summarize, med = quantile(b, prob = 0.5), q1 = quantile(b, prob = 0.25), q3 = quantile(b, prob = 0.75),  IQR = quantile(b, prob = 0.75) - quantile(b, prob = 0.25), se = sd(b))

qplot(nbin, med, data = summary.b) + geom_errorbar(aes(ymin = q1, ymax = q3)) + geom_smooth(se = FALSE)

qplot(nbin, IQR, data = summary.b) + geom_smooth(se = FALSE)

##############################################################################

files.txt <- dir("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/large-p-exp","*.txt")

pos <- 5
i <- 61
lineup.dat <- read.table(paste("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/large-p-exp/",files.txt[i], sep = ""), header = T)

qplot(x, 0, data = lineup.dat, col = factor(cl), geom = "jitter", ylim = c(-2,2), ylab = "") + facet_wrap(~ .sample)



plot_large_p_small_n_30_100_0_1_3.png  13 pos = 15

plot_large_p_small_n_30_60_0_1_1.png   47 pos = 18

plot_large_p_small_n_30_80_0_1_3.png   61 pos = 5

#### Tips Data

library(reshape2)

qplot(total_bill, tip, data = tips)

samples <- NULL
for(i in 1:19){
	x <- rexp(dim(tips)[1], 1/8)
    y <- 0.1*x + rnorm(dim(tips)[1], 3)
#    x <- 60*rbeta(dim(tips)[1], 5, 1)
#    y <- 0.08*x + rnorm(dim(tips)[1])
    dat <- data.frame(total_bill = x, tip = y, .n = i)
    samples <- rbind(samples, dat)   
}

dat <- data.frame(total_bill = tips$total_bill, tip = tips$tip)

pos <- sample(20, 1)
lineup.dat <- lineup(true = dat, samples = samples, pos = pos)

qplot(total_bill, tip, data = lineup.dat) + facet_wrap(~ .sample)

p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

### Dataset with outliers

x <- rnorm(40)
y <- NULL
for(i in 1:40){
y[i] <- rnorm(1, 2 + 0.1*2*x[i], sqrt(4*(1 - 0.1^2)))
}

xx <- rnorm(10, 10, 0.1)
yy <- rnorm(10, 12, 0.3)

x <- c(x, xx)
y <- c(y, yy)

dat <- data.frame(x = x,y = y)

qplot(x,y, data = dat, geom = "point", size = I(3), xlab = "X", ylab = "Y")

ggsave("bin-select-data6.pdf", height = 4, width = 4)

pos <- sample(20, 1)
lineup.dat <- lineup(null_permute("x"), true = dat, pos = pos)

qplot(x, y, data = lineup.dat) + facet_wrap(~ .sample)

p <- 2:10
q <- 2:10

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}


ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

samples <- NULL
for(j in 1:19){
	x <- rnorm(50, 5, 2.5)
	y <- NULL
	for(i in 1:50){
	y[i] <- rnorm(1, 2 + 0.1*2*x[i], sqrt(4*(1 - 0.1^2)))
	}
	dat1 <- data.frame(x = x, y = y, .n = j)
	samples <- rbind(samples, dat1)
}

pos <- sample(20, 1)
lineup.dat <- lineup(true = dat, samples = samples, pos = pos)

qplot(x, y, data = lineup.dat) + facet_wrap(~ .sample)

ggsave("bin-select-lineup6.pdf", height = 4, width = 4)

p <- 2:15
q <- 2:15

Diff <- NULL
for(r in p){
	dif1 <- NULL
	for(s in q){
		diff <- diff_bin(lineup.dat, r, s)
		dif <- data.frame(diff = diff, p = r, q = s)
		dif1 <- rbind(dif1, dif)
	}
	Diff <- rbind(Diff, dif1)
}

ggplot(Diff, aes(x = factor(p), y = factor(q))) + geom_tile(aes(fill = diff)) + scale_fill_gradient(high ="blue", low ="white") + xlab("p") + ylab("q")

ggsave("bin-select-plot5.pdf", height = 4, width = 4)

