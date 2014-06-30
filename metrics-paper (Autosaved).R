library(nullabor)
library(ggplot2)
library(plyr)
library(reshape)
library(fpc)
library(tourr)

##====================================Distance Metrics===========================================

## Distance based on Boxplots with indexing

box_dist_indx <- function(i, j){
	X <- lineup.dat[lineup.dat$.sample == i, ]
	PX <- lineup.dat[lineup.dat$.sample == j, ]
	X.sum <- ddply(X, .(group), summarize, sum.stat = quantile(val, c(0.25, 0.5, 0.75)))
	PX.sum <- ddply(PX, .(group), summarize, sum.stat = quantile(val, c(0.25, 0.5, 0.75)))
	abs.diff.X <- abs(X.sum$sum.stat[X.sum$group == levels(X.sum$group)[1]] - X.sum$sum.stat[X.sum$group == levels(X.sum$group)[2]])
	abs.diff.PX <- abs(PX.sum$sum.stat[PX.sum$group == levels(PX.sum$group)[1]] - PX.sum$sum.stat[PX.sum$group == levels(PX.sum$group)[2]])
	sqrt(sum((abs.diff.X - abs.diff.PX)^2))
}

## Distance based on Boxplots: No indexing

box_dist <- function(X, PX){
	X.sum <- ddply(X, .(group), summarize, sum.stat = quantile(val, c(0.25, 0.5, 0.75)))
	PX.sum <- ddply(PX, .(group), summarize, sum.stat = quantile(val, c(0.25, 0.5, 0.75)))
	abs.diff.X <- abs(X.sum$sum.stat[X.sum$group == levels(X.sum$group)[1]] - X.sum$sum.stat[X.sum$group == levels(X.sum$group)[2]])
	abs.diff.PX <- abs(PX.sum$sum.stat[PX.sum$group == levels(PX.sum$group)[1]] - PX.sum$sum.stat[PX.sum$group == levels(PX.sum$group)[2]])
	sqrt(sum((abs.diff.X - abs.diff.PX)^2))
}

### Regression based distance with indexing

reg_bin_indx <- function(i, j, nbins = 1){
	X <- lineup.dat[lineup.dat$.sample == i, ]
	PX <- lineup.dat[lineup.dat$.sample == j, ]
	ss <- seq(min(X[,1]), max(X[,1]), length = nbins + 1)
	beta.X <- NULL ; beta.PX <- NULL
	for(k in 1:nbins){
		X.sub <- subset(X, X[,1] >= ss[k] & X[,1] <= ss[k + 1])
		PX.sub <- subset(PX, X[,1] >= ss[k] & X[,1] <= ss[k + 1])
		b.X <- as.numeric(coef(lm(X.sub[,2] ~ X.sub[,1])))
		b.PX <- as.numeric(coef(lm(PX.sub[,2] ~ PX.sub[,1])))
		beta.X <- rbind(beta.X, b.X)
		beta.PX <- rbind(beta.PX, b.PX)
	}
	beta.X <- subset(beta.X, !is.na(beta.X[,2]))
	beta.PX <- subset(beta.PX, !is.na(beta.PX[,2]))
	sum((beta.X[,1] - beta.PX[,1])^2 + (beta.X[,2] - beta.PX[,2])^2)
}


### Regression based distance: No indexing

reg_bin <- function(X, PX, nbins = 1){
#	X <- lineup.dat[lineup.dat$.sample == i, ]
#	PX <- lineup.dat[lineup.dat$.sample == j, ]
	ss <- seq(min(X[,1]), max(X[,1]), length = nbins + 1)
	beta.X <- NULL ; beta.PX <- NULL
	for(k in 1:nbins){
		X.sub <- subset(X, X[,1] >= ss[k] & X[,1] <= ss[k + 1])
		PX.sub <- subset(PX, X[,1] >= ss[k] & X[,1] <= ss[k + 1])
		b.X <- as.numeric(coef(lm(X.sub[,2] ~ X.sub[,1])))
		b.PX <- as.numeric(coef(lm(PX.sub[,2] ~ PX.sub[,1])))
		beta.X <- rbind(beta.X, b.X)
		beta.PX <- rbind(beta.PX, b.PX)
	}
	beta.X <- subset(beta.X, !is.na(beta.X[,2]))
	beta.PX <- subset(beta.PX, !is.na(beta.PX[,2]))
	sum((beta.X[,1] - beta.PX[,1])^2 + (beta.X[,2] - beta.PX[,2])^2)
}

### Distance for Univariate Data

dist_uni_indx <- function(i, j){
	xx <- lineup.dat[lineup.dat$.sample == i, 1]
	yy <- lineup.dat[lineup.dat$.sample == j, 1]
	stat.xx <- c(mean(xx), sd(xx), moments::skewness(xx), moments::kurtosis(xx))
	stat.yy <- c(mean(yy), sd(yy), moments::skewness(yy), moments::kurtosis(yy))
	sqrt(sum((stat.xx - stat.yy)^2))
}

####Modified Binned Distance with indexing

bdist_mod_indx <- function(i,j, nbin.X = 5, nbin.Y = 5) {
	X <- lineup.dat[lineup.dat$.sample == i,]
	PX <- lineup.dat[lineup.dat$.sample == j,]
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

####Modified Binned Distance: No indexing

bin_dist <- function(X,PX, nbin.X = 5, nbin.Y = 5) {
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



##Weighted Bin Distance with indexing

wbdist_indx <- function(i,j, nbins=10) {
	X <- lineup.dat[lineup.dat$.sample == i,]
	PX <- lineup.dat[lineup.dat$.sample == j,]
	d1 <- MASS::kde2d(X[,1],X[,2],n=nbins,lims=c(range(X[,1]), range(X[,2])))
	d2 <- MASS::kde2d(PX[,1],PX[,2],n=nbins,lims=c(range(PX[,1]), range(PX[,2])))
	
	sqrt(sum((d1$z-d2$z)^2)/(sum(d1$z^2) * sum(d2$z^2)))	
}

## Distances based on separation: No Indexing

min_sep_dist <- function(X, PX, clustering = FALSE, nclust = 3){
	require(fpc)
	dX <- dist(X[,1:2])
	dPX <- dist(PX[,1:2])
	if(clustering){
			X$cl <- X[,3]
			PX$cl <- PX[,3]
			X.clus <- sort(cluster.stats(dX, clustering = X$cl)$separation)
			PX.clus <- sort(cluster.stats(dPX, clustering = X$cl)$separation)
	}
	else{
	complete.X <- cutree(hclust(dX), nclust)
	complete.PX <- cutree(hclust(dPX), nclust)
	X.clus <- sort(cluster.stats(dX, complete.X)$separation)
	PX.clus <- sort(cluster.stats(dPX, complete.PX)$separation)
	}
	sqrt(sum((X.clus - PX.clus)^2))
}

## Distances based on separation with indexing
	
min_sep_dist_indx <- function(i, j, clustering = FALSE, nclust = 3){
	X <- lineup.dat[lineup.dat$.sample == i,]
	PX <- lineup.dat[lineup.dat$.sample == j,]
	require(fpc)
	dX <- dist(X[,1:2])
	dPX <- dist(PX[,1:2])
	if(clustering){
			X$cl <- X[,3]
			PX$cl <- PX[,3]
	X.clus <- sort(cluster.stats(dX, clustering = X$cl)$separation)
	PX.clus <- sort(cluster.stats(dPX, clustering = X$cl)$separation)
	}
	else{
	complete.X <- cutree(hclust(dX), nclust)
	complete.PX <- cutree(hclust(dPX), nclust)
	X.clus <- sort(cluster.stats(dX, complete.X)$separation)
	PX.clus <- sort(cluster.stats(dPX, complete.PX)$separation)
	}
	sqrt(sum((X.clus - PX.clus)^2))
}

### Distances based on average separation 

ave_sep_dist <- function(X, PX, clustering = FALSE, nclust = 3) {
    dX <- dist(X[, 1:2])
    dPX <- dist(PX[, 1:2])
    if (clustering) {
        X$cl <- X[, 3]
        PX$cl <- PX[, 3]
        X.clus <- sort(cluster.stats(dX, clustering = X$cl)$average.toother)
        PX.clus <- sort(cluster.stats(dPX, clustering = PX$cl)$average.toother)
    } else {
        complete.X <- cutree(hclust(dX), nclust)
        complete.PX <- cutree(hclust(dPX), nclust)
       X.clus <- sort(cluster.stats(dX, complete.X)$average.toother)
        PX.clus <- sort(cluster.stats(dPX, complete.PX)$average.toother)
    }
    sqrt(sum((X.clus - PX.clus)^2))
} 

### Distances based on average separation with indexing

ave_sep_dist_indx <- function(i, j, clustering = FALSE, nclust = 3){
	X <- lineup.dat[lineup.dat$.sample == i,]
	PX <- lineup.dat[lineup.dat$.sample == j,]
	dX <- dist(X[,1:2])
	dPX <- dist(PX[,1:2])
	if(clustering){
			X$cl <- X[,3]
			PX$cl <- PX[,3]
        X.clus <- sort(cluster.stats(dX, clustering = X$cl)$average.toother)
        PX.clus <- sort(cluster.stats(dPX, clustering = PX$cl)$average.toother)
        	}
	else{
	complete.X <- cutree(hclust(dX), nclust)
	complete.PX <- cutree(hclust(dPX), nclust)
	X.clus <- sort(cluster.stats(dX, complete.X)$average.toother)
	PX.clus <- sort(cluster.stats(dPX, complete.PX)$average.toother)
	}
	sqrt(sum((X.clus - PX.clus)^2))
}	


### hausdorff Distance

##====================================Application to the Turk Experiment=====================================================

pos.1 <- 1:20
pos.2 <- 1:20

dat.pos <- expand.grid(pos.1 = pos.1, pos.2 = pos.2)

###====================================Turk 1 Experiment===============================================================


files <- dir("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/exp1")

results <- NULL

for(k in files){	
### Reading the lineup data

dat <- read.table(paste("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/exp1/",k, sep = ""), header = T)
#dat <- read.table(paste("U:/Documents/Research/Permutation/exp1/",i, sep = ""), header = T)

### Melting the data

dat.m <- melt(dat, id = c("age", "grp", "weight"))


### Changing the categorical variable to a numerical variable

dat.m$x <- as.numeric(dat.m$grp) 

### Breaking the variable name to get the position and type of plot (null or obs)

dat.m$plot <- substring(dat.m$variable, 1, 3)

dat.m$position <- substring(dat.m$variable, 4, 5)

### Finding the observed data

obs <- dat.m[dat.m$plot == "obs", c("x", "value") ]

### Storing the pic name

file.name <- k

split.file.name <- unlist(strsplit(file.name,"\\."))

split.2 <- unlist(strsplit(split.file.name[1],"\\_"))

pic_name <- paste("plot_", split.2[2],"_", split.2[3], "_", split.2[4], "_", split.2[5], "_", split.2[6], ".png", sep = "" )

### Calculating the distance metrics on the lineups. Only binned, weighted bin, canberra distance 
### and hausdorff are calculated.

dat.m <- dat.m[, c("x", "value", "position")]

names(dat.m) <- c("group", "val", ".sample")

lineup.dat <- dat.m

lineup.dat$group <- as.factor(lineup.dat$group)

stat <- ddply(dat.pos,.(pos.1, pos.2),summarize,box.dist = box_dist_indx(pos.1, pos.2), bin.dist = bdist_mod_indx(pos.1, pos.2, nbin.X = 8, nbin.Y = 8) )

res <- data.frame(pic_name, stat)

results <- rbind(results, res)

}

lineup.dat$.sample <- as.numeric(lineup.dat$.sample)

#write.table(results, "turk1-metrics.txt", row.names = F)

#results <- read.table("turk1-metrics.txt", header = T)
qplot(factor(group), val, data = lineup.dat, geom = "boxplot", col = group, ylab = "", xlab = "group") + facet_wrap(~.sample)

#ggsave("turk1-diff-box-prop.pdf", height = 5, width = 5.5)

res.exp1 <- read.csv("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/raw_data_turk1.csv")

res.exp1 <- subset(res.exp1, select = c(pic_name, response, plot_location, time_taken))

res.dat <- ddply(res.exp1, .(pic_name), summarize, prop = sum(response)/length(response), pos = mean(plot_location), m.time = median(time_taken))

metrics.sub <- subset(results, pos.1 != pos.2)

dat.merge <- merge(metrics.sub, res.dat, by = "pic_name")

dat.merge <- subset(dat.merge, pos.2 != pos)

dd <- ddply(dat.merge, .(pic_name, pos.1), summarize, box.mean = mean(box.dist), bin.mean = mean(bin.dist), len = length(box.dist), prop = mean(prop), m.time = mean(m.time))

prop.dist <- ddply(dd, .(pic_name), summarize, diff.box = box.mean[len == 19] - max(box.mean[len == 18]), grtr.box = sum(box.mean[len == 18] > box.mean[len == 19]) , diff.bin = bin.mean[len == 19] - max(bin.mean[len == 18]), grtr.bin = sum(bin.mean[len == 18] > bin.mean[len == 19]) , prop = mean(prop), m.time = mean(m.time))

prop.dist$special <- ifelse(prop.dist$pic_name == "plot_turk1_100_8_12_2.png", 1, 0)

### Faceting

prop.diff <- subset(prop.dist, select = c(pic_name, diff.box, diff.bin, prop, m.time, special))
prop.diff.m <- melt(prop.diff, id = c("pic_name", "prop", "special", "m.time"))

levels(prop.diff.m$variable) <- c("Boxplot Based Distance", "Binned Distance")

qplot(value, prop, data = prop.diff.m, geom = "point", size = I(3), ylim = c(0, 1), xlab = "Difference", ylab = "Detection Rate", shape = factor(special)) + geom_smooth(se = FALSE) + facet_wrap( ~ variable, scales = "free_x") + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("turk1-prop-box-bin.pdf", height = 4, width = 8.5)

qplot(value, m.time, data = prop.diff.m, geom = "point", size = I(3), xlab = "Difference", ylab = "Mean Time to Respond", shape = factor(special)) + facet_wrap( ~ variable, scales = "free_x") + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("turk1-mtime-box-bin.pdf", height = 4, width = 8.5)

grtr.diff <- subset(prop.dist, select = c(pic_name, grtr.box, grtr.bin, prop, special))
grtr.diff.m <- melt(grtr.diff, id = c("pic_name", "prop", "special"))

levels(grtr.diff.m$variable) <- c("Boxplot Based Distance", "Binned Distance")

qplot(value, prop, data = grtr.diff.m, geom = "point", size = I(3), ylim = c(0, 1), xlab = "Greater than Observed Plot", ylab = "Detection Rate", shape = factor(special)) + facet_wrap( ~ variable) + theme(legend.position = "none")

ggsave("turk1-grtr-box-bin.pdf", height = 4, width = 8.5)


### Individual Plots

qplot(diff.box, prop, data = prop.dist, ylim = c(0, 1), size = I(3), xlab = "Difference", ylab = "Detection Rate", shape = factor(special)) + geom_smooth( se = FALSE) + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("turk1-diff-box-prop.pdf", height = 5, width = 5.5)

qplot(factor(grtr.box), prop, data = prop.dist, ylim = c(0, 1), size = I(3), alpha = I(0.6), xlab = "Greater than observed plot", ylab = "Detection Rate", shape = factor(special)) + theme(legend.position = "none")

ggsave("turk1-grtr-box-prop.pdf", height = 5, width = 5.5)

qplot(diff.bin, prop, data = prop.dist, ylim = c(0, 1), size = I(3), xlab = "Difference", ylab = "Detection Rate", shape = factor(special)) + theme(legend.position = "none") + geom_smooth( se = FALSE) + geom_vline(xintercept = 0, col = "red")

ggsave("turk1-diff-bin-prop-8.pdf", height = 5, width = 5.5)

qplot(factor(grtr.bin), prop, data = prop.dist, ylim = c(0, 1), size = I(3), alpha = I(0.6), xlab = "greater than observed plot", ylab = "Detection Rate", shape = factor(special)) + theme(legend.position = "none")

ggsave("turk1-grtr-bin-prop-8.pdf", height = 5, width = 5.5)

qplot(diff.box, m.time, data = prop.dist, size = I(3), xlab = "difference", ylab = "median time to respond") + geom_smooth(method = "lm", se = FALSE) + geom_vline(xintercept = 0, col = "red")

ggsave("turk1-diff-box-mtime.pdf", height = 5, width = 5.5)

qplot(diff.bin, m.time, data = prop.dist, size = I(3), xlab = "difference", ylab = "median time to respond") + geom_smooth(method = "lm", se = FALSE) + geom_vline(xintercept = 0, col = "red")

ggsave("turk1-diff-bin-mtime.pdf", height = 5, width = 5.5)


###============================================================================
###Turk 2 Experiment
###============================================================================

files.png <- dir("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/exp2","*.png")

files.txt <- dir("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/exp2","*.txt")

metrics <- NULL
for(i in 1:length(files.txt)){
	dat <- read.table(paste("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/exp2/",files.txt[i], sep = ""), header = T)
	dat.m <- melt(dat, id = "X")
	dat.m$.sample <- substring(dat.m$variable, 2)
	lineup.dat <- data.frame(x = dat.m$X, z = dat.m$value, .sample = dat.m$.sample)
	metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, reg.bin = reg_bin_indx(pos.1, pos.2), bin.dist = bdist_mod_indx(pos.1, pos.2, nbin.X = 2, nbin.Y = 2), reg.no.int = reg_no_int_indx(pos.1, pos.2))
	metrics.dat <- data.frame(metrics.dat, pic_name = files.png[i])
	metrics <- rbind(metrics, metrics.dat)
}

#write.table(metrics, "turk2-metrics.txt", row.names = F)

#lineup.dat$.sample <- as.numeric(lineup.dat$.sample)

#qplot(x,z, data = lineup.dat, geom = "point", alpha = I(0.5), xlab = "X", ylab = "Y") + geom_smooth(method = "lm", se = FALSE) + facet_wrap(~.sample)


res.exp2 <- read.csv("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/raw_data_turk2.csv")

res.exp2 <- subset(res.exp2, select = c(pic_name, response, plot_location, time_taken))

res.dat <- ddply(res.exp2, .(pic_name), summarize, prop = sum(response)/length(response), pos = mean(plot_location), m.time = median(time_taken))

metrics.sub <- subset(metrics, pos.1 != pos.2)

dat.merge <- merge(metrics.sub, res.dat, by = "pic_name")

dat.merge <- subset(dat.merge, pos.2 != pos)

dd <- ddply(dat.merge, .(pic_name, pos.1), summarize, reg.mean = mean(reg.bin), bin.mean = mean(bin.dist), reg.no.int.mean = mean(reg.no.int), len = length(reg.bin), prop = mean(prop), m.time = mean(m.time))

prop.dist <- ddply(dd, .(pic_name), summarize, diff.reg = reg.mean[len == 19] - max(reg.mean[len == 18]), grtr.reg = sum(reg.mean[len == 18] > reg.mean[len == 19]), diff.reg.no.int = reg.no.int.mean[len == 19] - max(reg.no.int.mean[len == 18]), grtr.reg.no.int = sum(reg.no.int.mean[len == 18] > reg.no.int.mean[len == 19]), diff.bin = bin.mean[len == 19] - max(bin.mean[len == 18]), grtr.bin = sum(bin.mean[len == 18] > bin.mean[len == 19]), prop = mean(prop), m.time = mean(m.time))

prop.dist$special <- ifelse(prop.dist$pic_name == "plot_turk2_100_350_12_3.png", 1, 0)

### Facetted Plots

prop.diff <- subset(prop.dist, select = c(pic_name, diff.reg, diff.bin, prop, m.time, special))

prop.diff.m <- melt(prop.diff, id = c("pic_name", "prop", "special", "m.time"))

levels(prop.diff.m$variable) <- c("Regression Based Distance", "Binned Distance")

qplot(value, prop, data = prop.diff.m, geom = "point", size = I(3), xlab = "Difference", ylab = "Detection Rate", shape = factor(special)) + geom_smooth(se = FALSE) + facet_wrap( ~ variable, scales = "free_x") + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("turk2-prop-reg-bin.pdf", height = 4, width = 8.5)

qplot(value, m.time, data = prop.diff.m, geom = "point", size = I(3), xlab = "Difference", ylab = "Mean Time to Respond", shape = factor(special)) + facet_wrap( ~ variable, scales = "free_x") + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("turk2-mtime-reg-bin.pdf", height = 4, width = 8.5)

grtr.diff <- subset(prop.dist, select = c(pic_name, grtr.reg, grtr.bin, prop, special))
grtr.diff.m <- melt(grtr.diff, id = c("pic_name", "prop", "special"))

levels(grtr.diff.m$variable) <- c("Regression Based Distance", "Binned Distance")

qplot(value, prop, data = grtr.diff.m, geom = "point", size = I(3), ylim = c(0, 1), xlab = "Greater than Observed Plot", ylab = "Detection Rate", shape = factor(special)) + facet_wrap( ~ variable) + theme(legend.position = "none")

ggsave("turk2-grtr-reg-bin.pdf", height = 4, width = 8.5)

### Individual Plots

qplot(diff.reg, prop, data = prop.dist, size = I(3), xlab = "Difference", ylab = "Detection Rate", shape = factor(special)) + geom_smooth(se = FALSE) + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("turk2-diff-reg-prop.pdf", height = 5, width = 5.5)

qplot(factor(grtr.reg), prop, data = prop.dist, size = I(3), xlab = "Greater than observed plot", ylab = "Detection Rate", shape = factor(special))  + theme(legend.position = "none")

ggsave("turk2-grtr-reg-prop.pdf", height = 5, width = 5.5)

qplot(diff.bin, prop, data = prop.dist, size = I(3), xlab = "Difference", ylab = "Detection Rate", shape = factor(special)) + geom_smooth(se = FALSE) + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("turk2-diff-bin-prop-2.pdf", height = 5, width = 5.5)

qplot(factor(grtr.bin), prop, data = prop.dist, size = I(3), xlab = "Greater than observed plot", ylab = "Detection Rate", shape = factor(special))  + theme(legend.position = "none")

ggsave("turk2-grtr-bin-prop-2.pdf", height = 5, width = 5.5)

qplot(diff.reg, m.time , data = prop.dist, size = I(3), xlab = "difference", ylab = "median time to respond") + geom_smooth(method = "lm", se = FALSE) + geom_vline(xintercept = 0, col = "red")

ggsave("turk2-diff-reg-mtime.pdf", height = 5, width = 5.5)

qplot(diff.bin, m.time , data = prop.dist, size = I(3), xlab = "difference", ylab = "median time to respond") + geom_smooth(method = "lm", se = FALSE) + geom_vline(xintercept = 0, col = "red")

ggsave("turk2-diff-bin-mtime.pdf", height = 5, width = 5.5)


###===================================================================================================
###Large p, small n
###===================================================================================================

files.png <- dir("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/large-p-exp","*.png")

files.txt <- dir("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/large-p-exp","*.txt")

metrics1 <- NULL
metrics2 <- NULL
for(i in 1:length(files.txt)){
	dat <- read.table(paste("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/large-p-exp/",files.txt[i], sep = ""), header = T)
	if(dim(dat)[2] == 4){
	lineup.dat <- data.frame(x = dat$x, z = dat$cl, cl = dat$cl, .sample = dat$.sample)
	metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, b.mod = bdist_mod_indx(pos.1, pos.2, nbin.X = 10, nbin.Y = 10), sep.dist = min_sep_dist_indx(pos.1, pos.2, clustering = TRUE, nclust = 2), ave.dist = ave_sep_dist_indx(pos.1, pos.2, clustering = TRUE, nclust = 2))
	metrics.dat1 <- data.frame(metrics.dat, pic_name = files.png[i])
	metrics1 <- rbind(metrics1, metrics.dat1)
	}
	if(dim(dat)[2] == 6){
	lineup.dat <- data.frame(x = dat$X1, z = dat$X2, cl = dat$cl, .sample = dat$.sample)
	metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize,  b.mod = bdist_mod_indx(pos.1, pos.2, nbin.X = 5, nbin.Y = 5), sep.dist = min_sep_dist_indx(pos.1, pos.2, clustering = TRUE, nclust = 3), ave.dist = ave_sep_dist_indx(pos.1, pos.2, clustering = TRUE, nclust = 2))
	metrics.dat2 <- data.frame(metrics.dat, pic_name = files.png[i])
	metrics2 <- rbind(metrics2, metrics.dat2)
	}
	metrics <- rbind(metrics1, metrics2)
}

#write.table(metrics, "largep-metrics.txt", row.names = F)

#qplot(X1, X2, data = dat, geom = "point", alpha = I(0.7), color = factor(cl)) + facet_wrap(~.sample) + scale_colour_discrete(name = "Group")

res.exp.lp <- read.csv("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/raw_data_turk7.csv")

res.exp.lp <- subset(res.exp.lp, pic_name %in% files.png, select = c(pic_name, response, plot_location, time_taken))

res.dat <- ddply(res.exp.lp, .(pic_name), summarize, prop = sum(response)/length(response), pos = mean(plot_location), m.time = mean(time_taken))

metrics.sub <- subset(metrics, pos.1 != pos.2)

dat.merge <- merge(metrics.sub, res.dat, by = "pic_name")

dat.merge <- subset(dat.merge, pos.2 != pos)

dd <- ddply(dat.merge, .(pic_name, pos.1), summarize, bin.mean = mean(b.mod), sep.mean = mean(sep.dist), ave.mean = mean(ave.dist), len = length(b.mod), prop = mean(prop), m.time = mean(m.time))

prop.dist <- ddply(dd, .(pic_name), summarize, diff.bin = bin.mean[len == 19] - max(bin.mean[len == 18]), diff.sep = sep.mean[len == 19] - max(sep.mean[len == 18]), diff.ave = ave.mean[len == 19] - max(ave.mean[len == 18]), prop = mean(prop), m.time = mean(m.time))

#write.csv(prop.dist, "prop.dist.csv", row.names = FALSE) ## This includes min separation, 
## average separation, wb.ratio and average between matrix.

#prop.dist <- read.csv("prop.dist.csv")
#ggpairs(prop.dist[, 2:7])

library(GGally)
ggpairs(prop.dist[, -1])

nomatch.sep.ave <- subset(prop.dist, (diff.sep > 0 & diff.ave < 0) | (diff.sep < 0 & diff.ave > 0))

prop.dist <- ddply(dd, .(pic_name), summarize, diff.bin = bin.mean[len == 19] - max(bin.mean[len == 18]), grtr.bin = sum(bin.mean[len == 18] > bin.mean[len == 19]), diff.sep = sep.mean[len == 19] - max(sep.mean[len == 18]), grtr.sep = sum(sep.mean[len == 18] > sep.mean[len == 19]), diff.ave = ave.mean[len == 19] - max(ave.mean[len == 18]), grtr.ave = sum(ave.mean[len == 18] > ave.mean[len == 19]), prop = mean(prop), m.time = mean(m.time))

prop.dist$shape <- ifelse(prop.dist$pic_name == "plot_large_p_small_n_30_100_0_2_3.png", 1, 0)

### Facetted Plots

prop.diff <- subset(prop.dist, select = c(pic_name, diff.sep, diff.ave, diff.bin,  prop, m.time, shape))

prop.diff.m <- melt(prop.diff, id = c("pic_name", "prop", "m.time", "shape"))

levels(prop.diff.m$variable) <- c("Minimum Separation", "Average Separation","Binned Distance")

qplot(value, prop, data = prop.diff.m, geom = "point", size = I(3), xlab = "Difference", ylab = "Detection Rate", shape = factor(shape)) + geom_smooth(se = FALSE) + facet_wrap( ~ variable, scales = "free_x") + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("largep-prop-sep-bin.pdf", height = 4, width = 8.5)

qplot(value, m.time, data = prop.diff.m, geom = "point", size = I(3), xlab = "Difference", ylab = "Mean Time to Respond", shape = factor(shape)) + facet_wrap( ~ variable, scales = "free_x") + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("largep-mtime-sep-bin.pdf", height = 4, width = 8.5)

grtr.diff <- subset(prop.dist, select = c(pic_name, grtr.sep, grtr.ave, grtr.bin, prop, shape))
grtr.diff.m <- melt(grtr.diff, id = c("pic_name", "prop", "shape"))

levels(grtr.diff.m$variable) <- c("Minimum Separation", "Average Separation","Binned Distance")

qplot(value, prop, data = grtr.diff.m, geom = "point", size = I(3), ylim = c(0, 1), xlab = "Greater than Observed Plot", ylab = "Detection Rate") + facet_wrap( ~ variable) + theme(legend.position = "none")

ggsave("largep-grtr-sep-bin.pdf", height = 4, width = 8.5)

### Individual Plots

qplot(diff.bin, prop, data = prop.dist, size = I(3), xlab = "Difference", ylab = "Detection Rate", shape = factor(special)) + geom_smooth(se = FALSE) + geom_vline(xintercept = 0, col = "red") + theme(legend.position = "none")

ggsave("largep-diff-bin-prop-10-5.pdf", height = 5, width = 5.5)

qplot(diff.bin, m.time, data = prop.dist, size = I(3), xlab = "Difference", ylab = "Median time to respond") + geom_smooth(method = "lm",se = FALSE) + geom_vline(xintercept = 0, col = "red")

ggsave("largep-diff-bin-mtime.pdf", height = 5, width = 5.5)

qplot(diff.sep, prop, data = prop.dist, size = I(3), xlab = "Difference", ylab = "Detection Rate")+ geom_smooth(se = FALSE) + geom_vline(xintercept = 0, col = "red")+ theme(legend.position = "none")

ggsave("largep-diff-clus-prop.pdf", height = 5, width = 5.5)

qplot(diff.sep, m.time, data = prop.dist, size = I(3), xlab = "Difference", ylab = "Median time to respond") + geom_smooth(method = "lm",se = FALSE) + geom_vline(xintercept = 0, col = "red")

ggsave("largep-diff-clus-mtime.pdf", height = 5, width = 5.5)

qplot(factor(grtr.bin), prop, data = prop.dist, size = I(3), xlab = "Greater than observed plot", ylab = "Detection Rate", alpha = I(0.6), shape = factor(special)) + theme(legend.position = "none")

ggsave("largep-grtr-bin-prop-10-5.pdf", height = 5, width = 5.5)

qplot(factor(grtr.sep), prop, data = prop.dist, size = I(3), xlab = "greater than observed plot", ylab = "Detection Rate", alpha = I(0.6), shape = factor(special)) + theme(legend.position = "none")

ggsave("largep-grtr-clus-prop.pdf", height = 5, width = 5.5)

qplot(factor(grtr.sep), m.time, data = prop.dist, size = I(3), xlab = "greater than observed plot", ylab = "prop correct", alpha = I(0.6))

time.dist <- merge(prop.dist, res.exp.lp, by = "pic_name")

qplot(diff.bin, time_taken, data = subset(time.dist, time_taken < 250), size = I(3), alpha = I(0.3)) + geom_smooth(method = "lm", se = FALSE)

qplot(diff.sep, time_taken, data = subset(time.dist, time_taken < 250), size = I(3), alpha = I(0.3)) + geom_smooth(method = "lm", se = FALSE)


####============================================================================
## Lendie's Data - Turk 10 Experiment
####============================================================================

files.csv <- dir("/Users/Niladri/Documents/Research/Permutation/Adam's Data/Adam's and Lendie's data/turk10/lineups/data/", "*.csv")

pic_details <- read.csv("/Users/Niladri/Documents/Research/Permutation/Adam's Data/Adam's and Lendie's data/turk10/lineups/picture-details.csv")

files.svg.all <- dir("/Users/Niladri/Documents/Research/Permutation/Adam's Data/Adam's and Lendie's data/turk10/lineups/images/", "*.svg")

files.svg <- matrix(unlist(strsplit(as.character(pic_details$pic_name), "/")), ncol = 2, byrow = TRUE)[,2]

res.lendie <- read.csv("/Users/Niladri/Documents/Research/Permutation/Adam's Data/Adam's and Lendie's data/turk10/raw_data_turk10.csv")

res.lendie <- subset(res.lendie, pic_name %in% files.svg, select = c(pic_name, response, response_no, plot_location, time_taken))

res.dat <- ddply(res.lendie, .(pic_name), summarize, prop = sum(response)/length(response), pos = mean(plot_location), m.time = median(time_taken))

metrics <- NULL
for(i in 1:length(files.csv)){
	dat <- read.table(paste("/Users/Niladri/Documents/Research/Permutation/Adam's Data/Adam's and Lendie's data/turk10/lineups/data/",files.csv[i], sep = ""), header = T, sep = ",")
	lineup.dat <- data.frame(x = dat$naive1.qq.x, y = dat$naive1.qq.y, .sample = dat$.n)
	metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, reg.bin = reg_bin_indx(pos.1, pos.2), bin.dist = bdist_mod_indx(pos.1, pos.2, nbin.X = 5, nbin.Y = 5))
	metrics.dat <- data.frame(metrics.dat, data_name = files.csv[i])
	metrics <- rbind(metrics, metrics.dat)
}

pic_details <- subset(pic_details, select = c("pic_name","data_name"))
pic_details$pic_name <- files.svg

metrics.sub <- subset(metrics, pos.1 != pos.2)

metrics.m <- merge(metrics.sub, pic_details, by = "data_name")

res.mer <- merge(metrics.m, res.dat, by = "pic_name")

dat.merge <- subset(res.mer, pos.2 != pos)

dd <- ddply(dat.merge, .(pic_name, data_name, pos.1), summarize, reg.mean = mean(reg.bin), bin.mean = mean(bin.dist), len = length(bin.dist), prop = mean(prop), m.time = mean(m.time))

prop.dist <- ddply(dd, .(pic_name, data_name), summarize, diff.bin = bin.mean[len == 19] - max(bin.mean[len == 18]), grtr.bin = sum(bin.mean[len == 18] > bin.mean[len == 19]), diff.reg = reg.mean[len == 19] - max(reg.mean[len == 18]), grtr.reg = sum(reg.mean[len == 18] > reg.mean[len == 19]), prop = mean(prop), m.time = mean(m.time))

qplot(diff.reg, prop, data = prop.dist)

qplot(diff.bin, prop, data = prop.dist)


### Metrics Example

set.seed(1500)
X1 <- rnorm(50, 10, 2)
X2 <- NULL
for(i in 1:50){
	X2[i] <- rnorm(1, mean = 3 + 0.7*(X1[i] - 10), sd = sqrt(4*(1 - 0.7^2)))
}
true.dat <- data.frame(X1 = X1, X2 = X2)
cor(X1, X2)
qplot(X1, X2, geom = "point", size = I(3))

ggsave("dat-example-1.pdf", height = 3.5, width = 3.5)

binplot <- function(x, y, nbins= 8, plot=TRUE) {
	xd <- cut(x, breaks=nbins, labels=as.character(1:nbins))
	yd <- cut(y, breaks=nbins, labels=as.character(1:nbins))
	
	ndf <- as.data.frame(xtabs(~yd+xd))
	X <- data.frame(x=x, y=y)
	X$xnew <- (x-min(x))/(max(x)-min(x))*nbins + 0.5
	X$ynew <- (y-min(y))/(max(y)-min(y))*nbins + 0.5
	
	if (plot) {
	print(ggplot() + geom_tile(aes(xd,yd,fill=Freq), colour="grey50", data=ndf) + scale_fill_gradient2(name = "Count") + xlab("p") + ylab("q") )
	}
	invisible(ndf)
}
      
      
nij <- binplot(X1,X2) 

ggsave("bin-example-1.pdf", height = 3.5, width = 4.2) 

freqplot <- function(x, y, nbins= 8, plot=TRUE) {
	xd <- cut(x, breaks=nbins, labels=as.character(1:nbins))
	yd <- cut(y, breaks=nbins, labels=as.character(1:nbins))
	
	ndf <- as.data.frame(xtabs(~yd+xd))
	X <- data.frame(x=x, y=y)
	
	if (plot) {
	print(ggplot() + geom_tile(aes(xd,yd,fill=0.5), colour="grey50", data=ndf) + scale_fill_gradient2(name = "Count") + xlab("p") + ylab("q") + geom_text(aes(xd, yd, label = Freq), data = ndf) + theme(legend.position = "none") )
	}
	invisible(ndf)
}
      
      
mij <- freqplot(X1,X2) 

ggsave("freq-example-1.pdf", height = 3.5, width = 3.5) 

X1 <- sample(X1)
samp.dat <- data.frame(X1 = X1, X2 = X2)
qplot(X1, X2, geom = "point", size = I(3), xlab = "Permuted X1")

ggsave("dat-example-2.pdf", height = 3.5, width = 3.5)

nij <- binplot(X1,X2) 

ggsave("bin-example-2.pdf", height = 3.5, width = 4.2) 

mij <- freqplot(X1,X2) 

ggsave("freq-example-2.pdf", height = 3.5, width = 3.5) 



###================================================================
#### Exp1 : data generation
###================================================================

set.seed(1000)

b0 <- 5
b1 <- 15
b2 <- 8
sigma <- 12

x1 <- rpois(100, 30)
x2 <- rep(c(1,2), c(51, 49))
eps <- rnorm(100, 0, sigma)

y = b0 + b1*x1 + b2*x2 + eps

mod1 <- lm(y ~ x1)

summary(mod1)$sigma ##11.77

qplot(factor(x2), mod1$resid, geom = "boxplot")

qplot(group, val, data = obs.dat, geom = "boxplot")

###==================================================================
### Turk 1 Example
###==================================================================

dat <- read.table(file.choose(), header = T) ## plot_turk1_100_16_12_3

dat.m <- melt(dat, id = c("age", "grp", "weight"))

dat.m$x <- as.numeric(dat.m$grp) 

dat.m$plot <- substring(dat.m$variable, 1, 3)

dat.m$position <- substring(dat.m$variable, 4, 5)

dat.m <- dat.m[, c("x", "value", "position")]

names(dat.m) <- c("group", "val", ".sample")

lineup.dat <- dat.m

lineup.dat$group <- as.factor(lineup.dat$group)

lineup.dat$.sample <- as.numeric(as.character(lineup.dat$.sample))

levels(lineup.dat$group) <- c("A", "B")

qplot(group, val, data = lineup.dat, geom = "boxplot", col = group, ylab = "", xlab = "Group") + facet_wrap(~ .sample) + scale_color_discrete(name = "Group")

ggsave("turk1-example.pdf", height = 5, width = 5.5)


###=========================================================================
### Exp1 -- generation of distribution of distance metric
###=========================================================================

dat <- read.table(file.choose(), header = T)  ## plot_turk1_100_8_12_2

### The detection rate for the above lineup is 0.28. But the difference for both binned distance and regression based distance are large negative.

### Melting the data

dat.m <- melt(dat, id = c("age", "grp", "weight"))

### Changing the categorical variable to a numerical variable

dat.m$x <- as.numeric(dat.m$grp) 

### Breaking the variable name to get the position and type of plot (null or obs)

dat.m$plot <- substring(dat.m$variable, 1, 3)

dat.m$position <- substring(dat.m$variable, 4, 5)

### Finding the observed data

obs <- dat.m[dat.m$plot == "obs", c("x", "value") ]

dat.m <- dat.m[, c("x", "value", "position")]

names(dat.m) <- c("group", "val", ".sample")

lineup.dat <- dat.m

lineup.dat$group <- as.factor(lineup.dat$group)

lineup.dat$.sample <- as.numeric(as.character(lineup.dat$.sample))

obs.dat <- subset(lineup.dat, .sample == 20)

levels(lineup.dat$group) <- c("A", "B")

qplot(group, val, data = lineup.dat, geom = "boxplot", col = group, ylab = "", xlab = "Group") + facet_wrap(~ .sample) + scale_color_discrete(name = "Group") + theme(legend.position = "none")

ggsave("lineup-high-prop-neg-diff.pdf", height = 4, width = 4.5)


dat.bin <- dat.bin.28 <- dat.box <- NULL
for (i in 1:1000){
	samp.dat <- data.frame(group = obs.dat$group, val = rnorm(dim(obs.dat)[1], 0, 11.77))
	dat1 <- sapply(1:18, function(k){
		null.dat <- data.frame(group = obs.dat$group, val = rnorm(dim(obs.dat)[1], 0, 11.77))
 		b1 = bin_dist(samp.dat, null.dat, lineup.dat = lineup.dat, X.bin = 2, Y.bin = 2)
 		b2 = bin_dist(samp.dat, null.dat, lineup.dat = lineup.dat, X.bin = 2, Y.bin = 8)
 		s = box_dist(samp.dat, null.dat)
 		return(list(b1 = b1, b2 = b2, s = s))
 		})
 		dat1 <- matrix(unlist(dat1), nrow = 3 )
 		dat.bin <- c(dat.bin, mean(dat1[1,]))
 		dat.bin.28 <- c(dat.bin.28, mean(dat1[2,]))
 		dat.box <- c(dat.box, mean(dat1[3,]))
 }
 
df <- as.data.frame(t(dat1))
names(df) <- c("Binned-2-2", "Binned-2-8", "Boxplot Dist")

write.csv(df, "distr-turk1.csv", row.names = FALSE)
 
pos.1 <- 1:20
pos.2 <- 1:20

dat.pos <- expand.grid(pos.1 = pos.1, pos.2 = pos.2)

metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, bin.dist = bdist_mod_indx(pos.1, pos.2, nbin.X = 2, nbin.Y = 2))

pos <- 20

metrics.dat <- subset(metrics.dat, pos.1 != pos.2 & pos.2 != pos)

dd1 <- ddply(metrics.dat, .(pos.1), summarize, bin.mean = mean(bin.dist), len = length(bin.dist))

dat.bin <- as.data.frame(dat.bin)

ggplot()  + geom_density(data = dat.bin, aes(x = dat.bin), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd1, len == 19), aes(x= bin.mean, xend = bin.mean, y=0.02*max(density(dat.bin$dat.bin)$y), yend = 0.2*max(density(dat.bin$dat.bin)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd1, len != 19), aes(x = bin.mean, xend = bin.mean, y = rep(0.02*max(density(dat.bin$dat.bin)$y),19), yend = rep(0.1*max(density(dat.bin$dat.bin)$y),19)), size=1, alpha = I(0.7)) + xlab("Binned Distance (p = q = 2)") + ylab("") + geom_text(data = dd1, y = - 0.03*max(density(dat.bin$dat.bin)$y), size = 2.5, aes(x = bin.mean, label = pos.1)) + ylim(c(- 0.04*max(density(dat.bin$dat.bin)$y), max(density(dat.bin$dat.bin)$y) + 0.1*max(density(dat.bin$dat.bin)$y)))
        
ggsave("distribution-bin-dist-2-2-exp1.pdf", height = 4, width = 4.5)   

### bin_dist with x.bin = 2, y.bin = 8

metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, bin.dist = bdist_mod_indx(pos.1, pos.2, nbin.X = 2, nbin.Y = 8))

pos <- 20

metrics.dat <- subset(metrics.dat, pos.1 != pos.2 & pos.2 != pos)

dd1 <- ddply(metrics.dat, .(pos.1), summarize, bin.mean = mean(bin.dist), len = length(bin.dist))

dat.bin <- as.data.frame(dat.bin.28)

ggplot()  + geom_density(data = dat.bin, aes(x = dat.bin.28), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd1, len == 19), aes(x= bin.mean, xend = bin.mean, y=0.02*max(density(dat.bin$dat.bin.28)$y), yend = 0.2*max(density(dat.bin$dat.bin.28)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd1, len != 19), aes(x = bin.mean, xend = bin.mean, y = rep(0.02*max(density(dat.bin$dat.bin.28)$y),19), yend = rep(0.1*max(density(dat.bin$dat.bin.28)$y),19)), size=1, alpha = I(0.7)) + xlab("Binned Distance (p = 2, q = 8)") + ylab("") + geom_text(data = dd1, y = - 0.03*max(density(dat.bin$dat.bin.28)$y), size = 2.5, aes(x = bin.mean, label = pos.1)) + ylim(c(- 0.04*max(density(dat.bin$dat.bin.28)$y), max(density(dat.bin$dat.bin.28)$y) + 0.1*max(density(dat.bin$dat.bin.28)$y)))
        
ggsave("distribution-bin-dist-2-8-exp1.pdf", height = 4, width = 4.5)   

### box_dist

metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, box.dist = box_dist_indx(pos.1, pos.2))

pos <- 20

metrics.dat <- subset(metrics.dat, pos.1 != pos.2 & pos.2 != pos)

dd3 <- ddply(metrics.dat, .(pos.1), summarize, box.mean = mean(box.dist), len = length(box.dist))

dat.box <- as.data.frame(dat.box)

ggplot()  + geom_density(data = dat.box, aes(x = dat.box), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= box.mean, xend = box.mean, y=0.02*max(density(dat.box$dat.box)$y), yend = 0.2*max(density(dat.box$dat.box)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = box.mean, xend = box.mean, y = rep(0.02*max(density(dat.box$dat.box)$y),19), yend = rep(0.1*max(density(dat.box$dat.box)$y),19)), size=1, alpha = I(0.7)) + xlab("Boxplot Based Distance") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dat.box$dat.box)$y), size = 2.5, aes(x = box.mean, label = pos.1)) + ylim(c(- 0.04*max(density(dat.box$dat.box)$y), max(density(dat.box$dat.box)$y) + 0.1*max(density(dat.box$dat.box)$y)))
        
ggsave("distribution-box-dist-exp1.pdf", height = 4, width = 4.5)   

###==================================================================
### Turk 2 Example
###==================================================================

dat <- read.table(file.choose(), header = TRUE) #plot_turk2_100_600_12_2
dat.m <- melt(dat, id = "X")
dat.m$.sample <- substring(dat.m$variable, 2)
lineup.dat <- data.frame(x = dat.m$X, z = dat.m$value, .sample = dat.m$.sample)

lineup.dat$.sample <- as.numeric(as.character(lineup.dat$.sample))

qplot(x, z, data = lineup.dat, alpha = I(0.1), xlab = "X1", ylab = "X2") + geom_smooth(method = "lm", se = FALSE, size = 1) + facet_wrap(~ .sample)

ggsave("turk2-example.pdf", height = 5, width = 5.5)



###=================================================================
### Exp 2: data generation
###=================================================================

b0 <- 6
b1 <- -3.5
sigma <- 12
n <- 100

x1 <- rnorm(n, 0, 1)

y <- b0 + b1*x1 + rnorm(100, 0, sigma)

qplot(x1, y, geom = "point") + geom_smooth(method = "lm", se = FALSE) 

### Using Lineup

dat <- read.table(file.choose(), header = TRUE) #plot_turk2_100_350_12_3
dat.m <- melt(dat, id = "X")
dat.m$.sample <- substring(dat.m$variable, 2)
lineup.dat <- data.frame(x = dat.m$X, z = dat.m$value, .sample = dat.m$.sample)

lineup.dat$.sample <- as.numeric(as.character(lineup.dat$.sample))

qplot(x, z, data = lineup.dat, alpha = I(0.1), xlab = "X1", ylab = "X2") + geom_smooth(method = "lm", se = FALSE, size = 1) + facet_wrap(~ .sample)

ggsave("lineup-exp2-neg-diff-large-prop.pdf", height = 4, width = 4.5)

 
## From the lineup data

obs.dat <- lineup.dat[lineup.dat$.sample == 10, ]   

#qplot(x, z, data = obs.dat, geom = "point") + geom_smooth(method = "lm", se = FALSE) 

mod2 <- lm(z ~ 1, data = obs.dat)

mean.null <- predict(mod2)

sd.null <- summary(mod2)$sigma

###=============================================================================
### generation of distribution of distance metric
###=============================================================================


### using bin_dist with x.bin = 2, y.bin = 2

dat.bin <- dat.bin.82 <- dat.reg <- dat.reg.no.int <-  NULL
for (i in 1:1000){
	zz <- NULL
	for(i in 1:dim(obs.dat)[1]){
		zz[i] <- rnorm(1, mean.null[i], sd.null)
	}
	samp.dat <- data.frame(group = obs.dat$x, z = zz )
	dat1 <- sapply(1:18, function(k){
		yy <- NULL
	for(i in 1:dim(obs.dat)[1]){
		yy[i] <- rnorm(1, mean.null[i], sd.null)
		}
		null.dat <- data.frame(group = obs.dat$x, z = yy)
		b1 = bin_dist(samp.dat, null.dat, lineup.dat = lineup.dat, X.bin = 2, Y.bin = 2)
 		b2 = bin_dist(samp.dat, null.dat, lineup.dat = lineup.dat, X.bin = 8, Y.bin = 2)
 		r.int = reg_dist(samp.dat, null.dat)
 		r.no.int <- reg_no_int_dist(samp.dat, null.dat)
 		return(list(b1 = b1, b2 = b2, r.int = r.int, r.no.int = r.no.int))
 		})
	 	dat1 <- matrix(unlist(dat1), nrow = 4 )
 		dat.bin <- c(dat.bin, mean(dat1[1,]))
 		dat.bin.82 <- c(dat.bin.82, mean(dat1[2,]))
 		dat.reg <- c(dat.reg, mean(dat1[3,]))
 		dat.reg.no.int <- c(dat.reg.no.int, mean(dat1[4,]))
}

df <- as.data.frame(t(dat1))
names(df) <- c("Binned-2-2", "Binned-8-2", "Regression Dist", "Regression No Intercept")

write.csv(df, "distr-turk2.csv", row.names = FALSE)


#opt_diff(lineup.dat, var = c('x', 'z'), 2, 10, 2, 10, 19, plot = TRUE) 

metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, bin.dist = bdist_mod_indx(pos.1, pos.2, nbin.X = 2, nbin.Y = 2), bin.dist.82 = bdist_mod_indx(pos.1, pos.2, nbin.X = 8, nbin.Y = 2), reg.dist = reg_bin_indx(pos.1, pos.2), reg.dist.no.int = reg_no_int_indx(pos.1, pos.2))

pos <- 10

metrics.dat <- subset(metrics.dat, pos.1 != pos.2 & pos.2 != pos)

dd3 <- ddply(metrics.dat, .(pos.1), summarize, bin.mean = mean(bin.dist), bin.mean.82 = mean(bin.dist.82), reg.mean = mean(reg.dist), reg.mean.no.int = mean(reg.dist.no.int), len = length(bin.dist))

dat.bin <- as.data.frame(dat.bin)

ggplot()  + geom_density(data = dat.bin, aes(x = dat.bin), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= bin.mean, xend = bin.mean, y=0.02*max(density(dat.bin$dat.bin)$y), yend = 0.2*max(density(dat.bin$dat.bin)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = bin.mean, xend = bin.mean, y = rep(0.02*max(density(dat.bin$dat.bin)$y),19), yend = rep(0.1*max(density(dat.bin$dat.bin)$y),19)), size=1, alpha = I(0.7)) + xlab("Binned Distance (p = q = 2)") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dat.bin$dat.bin)$y), size = 2.5, aes(x = bin.mean, label = pos.1)) + ylim(c(- 0.04*max(density(dat.bin$dat.bin)$y), max(density(dat.bin$dat.bin)$y) + 0.1*max(density(dat.bin$dat.bin)$y)))


ggsave("distribution-bin-dist-2-2-exp2.pdf", height = 4, width = 4.5)   

dat.bin.82 <- as.data.frame(dat.bin.82)

ggplot()  + geom_density(data = dat.bin.82, aes(x = dat.bin.82), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= bin.mean.82, xend = bin.mean.82, y=0.02*max(density(dat.bin.82$dat.bin.82)$y), yend = 0.2*max(density(dat.bin.82$dat.bin.82)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = bin.mean.82, xend = bin.mean.82, y = rep(0.02*max(density(dat.bin.82$dat.bin.82)$y),19), yend = rep(0.1*max(density(dat.bin.82$dat.bin.82)$y),19)), size=1, alpha = I(0.7)) + xlab("Binned Distance (p = 8, q = 2)") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dat.bin.82$dat.bin.82)$y), size = 2.5, aes(x = bin.mean.82, label = pos.1)) + ylim(c(- 0.04*max(density(dat.bin.82$dat.bin.82)$y), max(density(dat.bin.82$dat.bin.82)$y) + 0.1*max(density(dat.bin.82$dat.bin.82)$y)))
        
ggsave("distribution-bin-dist-8-2-exp2.pdf", height = 4, width = 4.5)         

dat.reg <- as.data.frame(dat.reg)

ggplot()  + geom_density(data = dat.reg, aes(x = dat.reg), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= reg.mean, xend = reg.mean, y=0.02*max(density(dat.reg$dat.reg)$y), yend = 0.2*max(density(dat.reg$dat.reg)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = reg.mean, xend = reg.mean, y = rep(0.02*max(density(dat.reg$dat.reg)$y),19), yend = rep(0.1*max(density(dat.reg$dat.reg)$y),19)), size=1, alpha = I(0.7)) + xlab("Regression Based Distance") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dat.reg$dat.reg)$y), size = 2.5, aes(x = reg.mean, label = pos.1)) + ylim(c(- 0.04*max(density(dat.reg$dat.reg)$y), max(density(dat.reg$dat.reg)$y) + 0.1*max(density(dat.reg$dat.reg)$y)))

ggsave("distribution-reg-dist-exp2.pdf", height = 4, width = 4.5)

dat.reg.no.int <- as.data.frame(dat.reg.no.int)

ggplot()  + geom_density(data = dat.reg.no.int, aes(x = dat.reg.no.int), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= reg.mean.no.int, xend = reg.mean.no.int, y=0.02*max(density(dat.reg.no.int$dat.reg.no.int)$y), yend = 0.2*max(density(dat.reg.no.int$dat.reg.no.int)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = reg.mean.no.int, xend = reg.mean.no.int, y = rep(0.02*max(density(dat.reg.no.int$dat.reg.no.int)$y),19), yend = rep(0.1*max(density(dat.reg.no.int$dat.reg.no.int)$y),19)), size=1, alpha = I(0.7)) + xlab("Regression Based Distance \n (only slope)") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dat.reg.no.int$dat.reg.no.int)$y), size = 2.5, aes(x = reg.mean.no.int, label = pos.1)) + ylim(c(- 0.04*max(density(dat.reg.no.int$dat.reg.no.int)$y), max(density(dat.reg.no.int$dat.reg.no.int)$y) + 0.1*max(density(dat.reg.no.int$dat.reg.no.int)$y)))


ggsave("distribution-reg-no-int-dist-exp2.pdf", height = 4, width = 4.5) 

### P-value distance

pval_dist <- function(X) {
        as.numeric(summary(lm(X[,2] ~ X[,1], data = X))$coefficients[,4][2])
}


dat <- NULL
for (i in 1:1000){
	zz <- NULL
	for(k in 1:dim(obs.dat)[1]){
		zz[k] <- rnorm(1, mean = mean.null[k], sd = sd.null)
	}
	samp.dat <- data.frame(group = obs.dat$x, z = zz )
	qplot(group, z, data = samp.dat) + geom_smooth(method = "lm", se = FALSE)
	dat <- c(dat, pval_dist(samp.dat))
}


pval <- ddply(lineup.dat, .(.sample), summarize, p = pval_dist(data.frame(lineup.dat[lineup.dat$.sample == .sample,1], lineup.dat[lineup.dat$.sample == .sample,2])))

pos <- 10
m <- 20

qplot(dat, geom = "density", fill = I("grey80"), colour = I("grey80"), 
        xlab = "p-value", ylab = "") + geom_segment(aes(x = pval$p[pval$.sample != 
        10], xend = pval$p[pval$.sample != pos], y = rep(0.01 * min(density(dat)$y), 
        (m - 1)), yend = rep(0.1 * max(density(dat)$y), (m - 1))), size = 1, alpha = I(0.7)) + 
        geom_segment(aes(x = pval$p[pval$.sample == 
        10], xend = pval$p[pval$.sample == 10], y = 0.01 * min(density(dat)$y), yend = 0.2 * max(density(dat)$y)), 
            colour = "darkorange", size = 1) + geom_text(data = pval, y = -0.03 * max(density(dat)$y), 
        size = 2.5, aes(x = p, label = .sample)) + ylim(c(-0.04 * max(density(dat)$y), 
        max(density(dat)$y) + 0.1))

ggsave("distribution-pval-exp2.pdf", height = 4, width = 4.5) 

### Large p, Small n Data

dist <- read.csv("dat-dist-100-1.csv")

lineup.dat <- read.table(file.choose(), header = TRUE)  # plot_large_p_small_n_30_100_0_2_3

qplot(X1, X2, data = lineup.dat, geom = "point", alpha = I(0.8), size = I(1.5), col = factor(cl), xlab = "PD1", ylab = "PD2") + facet_wrap(~.sample) + scale_color_discrete(name = "Group")


ggsave("lineup-large-p-small-n.pdf", height = 4, width = 4.5)

metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, bin.dist = bdist_mod_indx(pos.1, pos.2, nbin.X = 6, nbin.Y = 4), min.sep.dist = min_sep_dist_indx(pos.1, pos.2), ave.sep.dist = ave_sep_dist_indx(pos.1, pos.2))

pos <- 20

metrics.dat <- subset(metrics.dat, pos.1 != pos.2 & pos.2 != pos)

dd3 <- ddply(metrics.dat, .(pos.1), summarize, bin.mean = mean(bin.dist),  min.sep.mean = mean(min.sep.dist), ave.sep.mean = mean(ave.sep.dist), len = length(bin.dist))

ggplot()  + geom_density(data = dist, aes(x = dat.b), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= bin.mean, xend = bin.mean, y=0.02*max(density(dist$dat.b)$y), yend = 0.2*max(density(dist$dat.b)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = bin.mean, xend = bin.mean, y = rep(0.02*max(density(dist$dat.b)$y),19), yend = rep(0.1*max(density(dist$dat.b)$y),19)), size=1, alpha = I(0.7)) + xlab("Binned Distance (p = 6, q = 4)") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dist$dat.b)$y), size = 2.5, aes(x = bin.mean, label = pos.1)) + ylim(c(- 0.04*max(density(dist$dat.b)$y), max(density(dist$dat.b)$y) + 0.1*max(density(dist$dat.b)$y)))

ggsave("distribution-bin-dist-6-4-lpexp.pdf", height = 4, width = 4.5) 

ggplot()  + geom_density(data = dist, aes(x = dat.s), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= min.sep.mean, xend = min.sep.mean, y=0.02*max(density(dist$dat.s)$y), yend = 0.2*max(density(dist$dat.s)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = min.sep.mean, xend = min.sep.mean, y = rep(0.02*max(density(dist$dat.s)$y),19), yend = rep(0.1*max(density(dist$dat.s)$y),19)), size=1, alpha = I(0.7)) + xlab("Minimum Separation Based Distance") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dist$dat.s)$y), size = 2.5, aes(x = min.sep.mean, label = pos.1)) + ylim(c(- 0.04*max(density(dist$dat.s)$y), max(density(dist$dat.s)$y) + 0.1*max(density(dist$dat.s)$y)))

ggsave("distribution-min-sep-dist-lpexp.pdf", height = 4, width = 4.5) 

ggplot()  + geom_density(data = dist, aes(x = dat.ms), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= ave.sep.mean, xend = ave.sep.mean, y=0.02*max(density(dist$dat.ms)$y), yend = 0.2*max(density(dist$dat.ms)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = ave.sep.mean, xend = ave.sep.mean, y = rep(0.02*max(density(dist$dat.ms)$y),19), yend = rep(0.1*max(density(dist$dat.ms)$y),19)), size=1, alpha = I(0.7)) + xlab("Average Separation Based Distance") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dist$dat.ms)$y), size = 2.5, aes(x = ave.sep.mean, label = pos.1)) + ylim(c(- 0.04*max(density(dist$dat.ms)$y), max(density(dist$dat.ms)$y) + 0.1*max(density(dist$dat.ms)$y)))

ggsave("distribution-ave-sep-dist-lpexp.pdf", height = 4, width = 4.5) 

opt_bin_diff(lineup.dat, var = c('X1', 'X2'), 2, 10, 2, 10, pos = 16,plot = TRUE, m = 20)

prop.dist <- read.csv("prop.dist.csv")

####========================================================================================================     
        

### using bin_dist with x.bin = 8, y.bin = 2

opt_diff(lineup.dat, var = c('x', 'z'), 2, 10, 2, 10, 10, plot = TRUE) 

dat <- NULL
for (i in 1:1000){
	zz <- NULL
	for(i in 1:dim(obs.dat)[1]){
		zz[i] <- rnorm(1, mean.null[i], sd.null)
	}
	samp.dat <- data.frame(group = obs.dat$x, z = zz )
	dat1 <- replicate(18, {
		yy <- NULL
	for(i in 1:dim(obs.dat)[1]){
		yy[i] <- rnorm(1, mean.null[i], sd.null)
		}
		null.dat <- data.frame(group = obs.dat$x, z = yy)
		bin_dist(samp.dat, null.dat, lineup.dat = lineup.dat, X.bin = 8, Y.bin = 2)
	})
	dat <- c(dat, mean(dat1))
}



ddd <- distmet(lineup.dat, var = c("x", "z"), 'bin_dist', null_permute("x"), pos = 10, dist.arg = list(X.bin = 8, Y.bin = 2))

pos <- 10
m <- 20        

### Using reg_dist

dat <- NULL
for (i in 1:1000){
	zz <- NULL
	for(i in 1:dim(obs.dat)[1]){
		zz[i] <- rnorm(1, mean.null[i], sd.null)
	}
	samp.dat <- data.frame(group = obs.dat$x, z = zz )
	dat1 <- replicate(18, {
		yy <- NULL
	for(i in 1:dim(obs.dat)[1]){
		yy[i] <- rnorm(1, mean.null[i], sd.null)
		}
		null.dat <- data.frame(group = obs.dat$x, z = yy)
		reg_dist(samp.dat, null.dat)
	})
	dat <- c(dat, mean(dat1))
}


ddd <- distmet(lineup.dat, var = c("x", "z"), 'reg_dist', null_permute("x"), pos = 10)

pos <- 10
m <- 20

dat.reg <- as.data.frame(dat.reg)

ggplot()  + geom_density(data = dat.reg, aes(x = dat.reg), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= reg.mean, xend = reg.mean, y=0.02*max(density(dat.reg$dat.reg)$y), yend = 0.2*max(density(dat.reg$dat.reg)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = reg.mean, xend = reg.mean, y = rep(0.02*max(density(dat.reg$dat.reg)$y),19), yend = rep(0.1*max(density(dat.reg$dat.reg)$y),19)), size=1, alpha = I(0.7)) + xlab("Distance") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dat.reg$dat.reg)$y), size = 2.5, aes(x = reg.mean, label = pos.1)) + ylim(c(- 0.04*max(density(dat.reg$dat.reg)$y), max(density(dat.reg$dat.reg)$y) + 0.1*max(density(dat.reg$dat.reg)$y)))

ggsave("distribution-reg-dist-exp2.pdf", height = 4, width = 4.5) 

### Reg distance: no intercept

reg_no_int_dist <- function(X, PX, nbins = 1) {
    ss <- seq(min(X[, 1]), max(X[, 1]), length = nbins + 1)
    beta.X <- NULL
    beta.PX <- NULL
    for (k in 1:nbins) {
        X.sub <- subset(X, X[, 1] >= ss[k] & X[, 1] <= ss[k + 1])
        PX.sub <- subset(PX, X[, 1] >= ss[k] & X[, 1] <= ss[k + 1])
        b.X <- as.numeric(coef(lm(X.sub[, 2] ~ X.sub[, 1])))
        b.PX <- as.numeric(coef(lm(PX.sub[, 2] ~ PX.sub[, 1])))
        beta.X <- rbind(beta.X, b.X)
        beta.PX <- rbind(beta.PX, b.PX)
    }
    beta.X <- subset(beta.X, !is.na(beta.X[, 2]))
    beta.PX <- subset(beta.PX, !is.na(beta.PX[, 2]))
    sum((beta.X[, 2] - beta.PX[, 2])^2)
}

### Reg distance: no intercept

reg_no_int_indx <- function(i, j, nbins = 1){
	X <- lineup.dat[lineup.dat$.sample == i, ]
	PX <- lineup.dat[lineup.dat$.sample == j, ]
	ss <- seq(min(X[,1]), max(X[,1]), length = nbins + 1)
	beta.X <- NULL ; beta.PX <- NULL
	for(k in 1:nbins){
		X.sub <- subset(X, X[,1] >= ss[k] & X[,1] <= ss[k + 1])
		PX.sub <- subset(PX, X[,1] >= ss[k] & X[,1] <= ss[k + 1])
		b.X <- as.numeric(coef(lm(X.sub[,2] ~ X.sub[,1])))
		b.PX <- as.numeric(coef(lm(PX.sub[,2] ~ PX.sub[,1])))
		beta.X <- rbind(beta.X, b.X)
		beta.PX <- rbind(beta.PX, b.PX)
	}
	beta.X <- subset(beta.X, !is.na(beta.X[,2]))
	beta.PX <- subset(beta.PX, !is.na(beta.PX[,2]))
	sum((beta.X[,2] - beta.PX[,2])^2)
}



dat <- NULL
for (i in 1:1000){
	zz <- NULL
	for(i in 1:dim(obs.dat)[1]){
		zz[i] <- rnorm(1, mean.null[i], sd.null)
	}
	samp.dat <- data.frame(group = obs.dat$x, z = zz )
	dat1 <- replicate(18, {
		yy <- NULL
	for(i in 1:dim(obs.dat)[1]){
		yy[i] <- rnorm(1, mean.null[i], sd.null)
		}
		null.dat <- data.frame(group = obs.dat$x, z = yy)
		reg_no_int_dist(samp.dat, null.dat)
	})
	dat <- c(dat, mean(dat1))
}


ddd <- distmet(lineup.dat, var = c("x", "z"), 'reg_no_int_dist', null_permute("x"), pos = 10)

dat.reg.no.int <- as.data.frame(dat.reg.no.int)

ggplot()  + geom_density(data = dat.reg.no.int, aes(x = dat.reg.no.int), fill = "grey80", col = "grey80" ) + geom_segment(data = subset(dd3, len == 19), aes(x= reg.mean.no.int, xend = reg.mean.no.int, y=0.02*max(density(dat.reg.no.int$dat.reg.no.int)$y), yend = 0.2*max(density(dat.reg.no.int$dat.reg.no.int)$y)), colour="darkorange", size=1)  + geom_segment(data = subset(dd3, len != 19), aes(x = reg.mean.no.int, xend = reg.mean.no.int, y = rep(0.02*max(density(dat.reg.no.int$dat.reg.no.int)$y),19), yend = rep(0.1*max(density(dat.reg.no.int$dat.reg.no.int)$y),19)), size=1, alpha = I(0.7)) + xlab("Distance") + ylab("") + geom_text(data = dd3, y = - 0.03*max(density(dat.reg.no.int$dat.reg.no.int)$y), size = 2.5, aes(x = reg.mean.no.int, label = pos.1)) + ylim(c(- 0.04*max(density(dat.reg.no.int$dat.reg.no.int)$y), max(density(dat.reg.no.int$dat.reg.no.int)$y) + 0.1*max(density(dat.reg.no.int$dat.reg.no.int)$y)))


ggsave("distribution-reg-no-int-dist-exp2.pdf", height = 4, width = 4.5) 


###=======================================================================================================================

res.dat <- read.csv(file.choose()) 
is.character(res.dat$pic_name)
res.lineup <- subset(res.dat, pic_name == "plot_turk2_100_350_12_3.png") 

pval <- function(t){
	summary(lm(z ~ x, data = subset(lineup.dat, lineup.dat$.sample == t)))$coefficients[,4][2]
}

dat <- read.table("dat_turk2_100_350_12_3.txt", header = T)
dat.m <- melt(dat, id = "X")
dat.m$.sample <- substring(dat.m$variable, 2)
lineup.dat <- data.frame(x = dat.m$X, z = dat.m$value, .sample = dat.m$.sample)
metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, reg.bin = reg_bin_indx(pos.1, pos.2), bin.dist = bdist_mod_indx(pos.1, pos.2, nbin.X = 8, nbin.Y = 2))

metrics.sub <- subset(metrics.dat, pos.1 != pos.2 & pos.2 != 10)

dd <- ddply(metrics.sub, .(pos.1), summarize, reg.mean = mean(reg.bin), bin.mean = mean(bin.dist))

dd$pic_name = "plot_turk2_100_350_12_3.png"

dat.merge <- merge(dd, res.lineup, by = "pic_name")

dat.merge <- subset(dat.merge, select = c(reg.mean, bin.mean, response_no, pos.1))

rel <- ddply(dat.merge, .(pos.1), summarize, rel.freq = sum(response_no == pos.1)/dim(res.lineup)[1], reg.mean = mean(reg.mean), bin.mean = mean(bin.mean))


pvalue.dat <- ddply(lineup.dat, .(.sample), summarize, p = pval(.sample))

pval.merge <- merge(rel, pvalue.dat, by.x = "pos.1", by.y = ".sample")

pval.merge$plot_loc <- ifelse(pos.1 == 10, 1, 0)

qplot(log(p), rel.freq, data = pval.merge, geom = "point", col = plot_loc) + geom_linerange(aes(x = log(p), ymin = 0, ymax = rel.freq )) + theme(legend.position="none") + scale_colour_continuous(high = "red", low = "black") + scale_x_continuous("p-value" )  + scale_y_continuous("Relative Frequency", limits = c(0,1), breaks = c(0, 0.5, 1))

ggsave("rel-pvalue.pdf", height = 2.5, width = 4)

rel.m <- melt(rel, id = c("pos.1", "rel.freq"))

rel.m$plot_loc <- ifelse(pos.1 == 10, 1, 0)

levels(rel.m$variable) <- c("Regression Based", "Binned")

qplot(value, rel.freq, data = rel.m, geom = "point", col = plot_loc) + geom_linerange(aes(x = value, ymin = 0, ymax = rel.freq )) + theme(legend.position="none") + scale_colour_continuous(high = "red", low = "black") + facet_grid(variable ~ ., scales = "free_x") + scale_x_continuous("Mean Distances" )  + scale_y_continuous("Relative Frequency", limits = c(0,1), breaks = c(0, 0.5, 1))

# + theme(legend.position="none", axis.ticks = element_blank(), axis.text.x = element_blank())
      
 
 ### Turk Experiment : Large p, Small n
 

lineup.dat <- read.table(file.choose(), header = TRUE)  # plot_large_p_small_n_30_80_0_2_3

qplot(X1, X2, data = lineup.dat, geom = "point", alpha = I(0.6), size = I(3), col = factor(cl)) + facet_wrap(~.sample) + scale_color_discrete(name = "Group")

ggsave("lineup-large-p-small-n.pdf", height = 5, width = 5.5)

#### New Lineup

generate_plot_2d<-function(n=30,p, noise=1, m=20){
	x<-matrix(rnorm(p*n),ncol=p)
	if(noise==0){
x[1:10,(p-1)]<-x[1:10,(p-1)]+3
x[11:20,(p-1)]<-x[11:20,(p-1)]-3
x[21:30,p]<-x[21:30,p]+sqrt(27)
}
colnames(x)<-paste("X",1:(p),sep="")
x<-scale(x)
x<-data.frame(x, cl=factor(c(rep(1,n/3),rep(2,n/3),rep(3,n/3))))
d=2

optima <- save_history(x[,-(p+1)], tour_path=guided_tour(index_f=pda_pp(cl=x[,(p+1)], lambda=0.2), max.tries=1000), max_bases=100, rescale=F)
nbases<-dim(optima)[3]
optima.global<-unclass(optima)[,,nbases]

projdata.true<-data.frame(as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)], nbases=rep(nbases,n))

projdata.samples<-NULL
flag <- 0
 while(flag < 19) {
 x[,(p+1)]<-sample(x[,(p+1)])
 optima <- save_history(x[,-(p+1)], guided_tour(index_f=pda_pp(cl=x[,(p+1)], lambda=0.2),max.tries=100), max_bases=100, rescale=F)
 nbases<-dim(optima)[3]
 optima.global<-unclass(optima)[,,nbases]
 projdata<-data.frame(as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)], nbases=rep(nbases,30))
 lamb<-summary(manova(cbind(X1, X2)~cl, data=projdata), test="Wilks")[[4]][3]
 if(lamb < 0.001){
 projdata.samples<-rbind(projdata.samples, projdata)
 flag <- flag + 1
 }else
 {  
 	projdata <- NULL
 	projdata.samples <- rbind(projdata.samples, projdata)
 	flag = flag
 	}
 cat(flag, lamb, "\n") 
}
projdata.samples$.n <- rep(1:19, each = 30)
#pos<-sample(m,1)
lineup.data<-lineup(true=projdata.true, samples=projdata.samples,pos=1)
return(lineup.data)
}

dat <- generate_plot_2d(p = 80, noise = 0)
qplot(X1, X2, data=dat, colour=cl) + facet_wrap(~ .sample) + scale_colour_discrete(name="Group") + scale_y_continuous("X2") + scale_x_continuous("X1")

ggsave("largep-lineup-new-1.pdf", height = 5, width = 5.5)

n <- 30; p <- 80
dat.b <- dat.s <- NULL
for (i in 1:15){
		flag <- 0
		while(flag < 1){
		x<-matrix(rnorm(p*n),ncol=p)
		x[1:10,(p-1)]<-x[1:10,(p-1)]+3
		x[11:20,(p-1)]<-x[11:20,(p-1)]-3
		x[21:30,p]<-x[21:30,p]+sqrt(27)
		colnames(x)<-paste("X",1:(p),sep="")
		x<-scale(x)
		x<-data.frame(x, cl=factor(c(rep(1,n/3),rep(2,n/3),rep(3,n/3))))
		x[,(p+1)]<-sample(x[,(p+1)])
		optima <- save_history(x[,-(p+1)], tour_path=guided_tour(index_f=pda_pp(cl=x[,(p+1)], 		lambda=0.2), max.tries=100), max_bases=100, rescale=F)
		nbases<-dim(optima)[3]
		optima.global<-unclass(optima)[,,nbases]
		projdata.null<-data.frame(as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)], 				nbases=rep(nbases,n))
		lamb<-summary(manova(cbind(X1, X2)~cl, data=projdata.null), test="Wilks")[[4]][3]
		 if(lamb < 0.001){
 		flag <- flag + 1
 		}
 		cat(flag, lamb, "\n")
 		}
		
		dat1 <- sapply(1:5, function(k){
		flag1 <- 0
		while(flag1 < 1){	
		 x[,(p+1)]<-sample(x[,(p+1)])
 		optima <- save_history(x[,-(p+1)], guided_tour(index_f=pda_pp(cl=x[,(p+1)], 				lambda=0.2),max.tries=100), max_bases=100, rescale=F)
 		nbases<-dim(optima)[3]
 		optima.global<-unclass(optima)[,,nbases]
 		projdata<-data.frame(as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)], 					nbases=rep(nbases,30))
		lamb1 <- summary(manova(cbind(X1, X2)~cl, data=projdata), test="Wilks")[[4]][3]
		 if(lamb1 < 0.001){
 		flag1 <- flag1 + 1
 		}
 		cat("flag1 =", flag1, "i = ", i, lamb1, "\n")
 		}
 		projdata.null$cl <- as.numeric(projdata.null$cl)
 		projdata$cl <- as.numeric(projdata$cl)
 		b = bin_dist(projdata.null, projdata, lineup.dat = dat, X.bin = 6, Y.bin = 4)
 		s = sep_dist(projdata.null, projdata, clustering = TRUE, nclust = 3)
 		return(list(b = b, s = s))
 		})
 		dat1 <- matrix(unlist(dat1), nrow = 2 )
 		dat.b <- c(dat.b, mean(dat1[1,]))
 		dat.s <- c(dat.s, mean(dat1[2,]))
} 

#opt_diff(lineup.dat, var = c('X1', 'X2'), 2, 10, 2, 10, 20, plot = TRUE)    

#detach(package:plyr)
library(dplyr)
ddd <- distmet(dat, var = c("X1", "X2"), 'bin_dist', null_permute("cl"), pos = 1, dist.arg = list(X.bin = 6, Y.bin = 4))

pos <- 1
m <- 20

qplot(dat.b, geom = "density", fill = I("grey80"), colour = I("grey80"), 
        xlab = "Permutation distribution", ylab = "") + geom_segment(aes(x = ddd$lineup$mean.dist[ddd$lineup$plotno != 
        pos], xend = ddd$lineup$mean.dist[ddd$lineup$plotno != pos], y = rep(0.01 * min(density(dat.b)$y), 
        (m - 1)), yend = rep(0.05 * max(density(dat.b)$y), (m - 1))), size = 1, alpha = I(0.7)) + 
        geom_segment(aes(x = ddd$lineup$mean.dist[ddd$lineup$plotno == pos], xend = ddd$lineup$mean.dist[ddd$lineup$plotno == 
            pos], y = 0.01 * min(density(dat.b)$y), yend = 0.1 * max(density(dat.b)$y)), 
            colour = "darkorange", size = 1) + geom_text(data = ddd$lineup, y = -0.03 * max(density(dat.b)$y), 
        size = 2.5, aes(x = mean.dist, label = plotno)) + ylim(c(-0.04 * max(density(dat.b)$y), 
        max(density(dat.b)$y) + 0.1))
        
ggsave("bin-dist-largep-6-4-new-1.pdf", height = 5, width = 5.5)        
        
### Sep dist

dat$cl <- as.numeric(dat$cl)
ddd <- distmet(dat, var = c("X1", "X2", "cl"), 'sep_dist', null_permute("cl"), pos = 1, dist.arg = list(clustering = TRUE, nclust = 3))  


pos <- 1; m = 20

qplot(dat.s, geom = "density", fill = I("grey80"), colour = I("grey80"), 
        xlab = "Permutation distribution", ylab = "") + geom_segment(aes(x = ddd$lineup$mean.dist[ddd$lineup$plotno != 
        pos], xend = ddd$lineup$mean.dist[ddd$lineup$plotno != pos], y = rep(0.01 * min(density(dat.s)$y), 
        (m - 1)), yend = rep(0.05 * max(density(dat.s)$y), (m - 1))), size = 1, alpha = I(0.7)) + 
        geom_segment(aes(x = ddd$lineup$mean.dist[ddd$lineup$plotno == pos], xend = ddd$lineup$mean.dist[ddd$lineup$plotno == 
            pos], y = 0.01 * min(density(dat.s)$y), yend = 0.1 * max(density(dat.s)$y)), 
            colour = "darkorange", size = 1) + geom_text(data = ddd$lineup, y = -0.03 * max(density(dat.s)$y), 
        size = 2.5, aes(x = mean.dist, label = plotno)) + ylim(c(-0.04 * max(density(dat.s)$y), 
        max(density(dat.s)$y) + 0.1))

      
ggsave("sep-dist-largep-new-1.pdf", height = 5, width = 5.5)  


### Case Study

## Separation Distance

X <- lineup.dat[lineup.dat$.sample == 20, ]
PX <- lineup.dat[lineup.dat$.sample == 7, ]

	

M2sep_dist_indx <- function(i, j, clustering = FALSE, nclust = 3){
	X <- lineup.dat[lineup.dat$.sample == i,]
	PX <- lineup.dat[lineup.dat$.sample == j,]
	require(fpc)
	dX <- dist(X[,1:2])
	dPX <- dist(PX[,1:2])
	if(clustering){
			X$cl <- X[,3]
			PX$cl <- PX[,3]
        X.clus <- as.vector(cluster.stats(dX, clustering = X$cl)$ave.between.matrix)
        PX.clus <- as.vector(cluster.stats(dPX, clustering = PX$cl)$ave.between.matrix)
        	}else{
	complete.X <- cutree(hclust(dX), nclust)
	complete.PX <- cutree(hclust(dPX), nclust)
	X.clus <- as.vector(cluster.stats(dX, complete.X)$ave.between.matrix)
	PX.clus <- as.vector(cluster.stats(dPX, complete.PX)$ave.between.matrix)
	}
	sqrt(sum((X.clus - PX.clus)^2))
}	

M3sep_dist_indx <- function(i, j, clustering = FALSE, nclust = 3){
	X <- lineup.dat[lineup.dat$.sample == i,]
	PX <- lineup.dat[lineup.dat$.sample == j,]
	require(fpc)
	dX <- dist(X[,1:2])
	dPX <- dist(PX[,1:2])
	if(clustering){
			X$cl <- X[,3]
			PX$cl <- PX[,3]
        X.clus <- cluster.stats(dX, clustering = X$cl)$wb.ratio
        PX.clus <- cluster.stats(dPX, clustering = PX$cl)$wb.ratio
        	}else{
	complete.X <- cutree(hclust(dX), nclust)
	complete.PX <- cutree(hclust(dPX), nclust)
	X.clus <- cluster.stats(dX, complete.X)$wb.ratio
	PX.clus <- cluster.stats(dPX, complete.PX)$wb.ratio
	}
	sqrt(sum((X.clus - PX.clus)^2))
}


lineup.dat <- read.table(file.choose(), header = TRUE)  # plot_large_p_small_n_30_80_0_2_3

qplot(X1, X2, data = lineup.dat, geom = "point", alpha = I(0.6), size = I(2), col = factor(cl)) + facet_wrap(~.sample) + scale_color_discrete(name = "Group")

part.largep <- distmet(lineup.dat, var = c("X1", "X2", "cl"), 'M3sep_dist', null_permute("cl"), pos = 20, repl = 1, dist.arg = list(clustering = TRUE, nclust = 3))

part.largep$lineup$plotno[order(part.largep$lineup$mean.dist, decreasing = TRUE)]   

files.png <- dir("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/large-p-exp","*.png")

files.txt <- dir("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/large-p-exp","*.txt")

metrics1 <- NULL
metrics2 <- NULL
for(i in 1:length(files.txt)){
	dat <- read.table(paste("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/large-p-exp/",files.txt[i], sep = ""), header = T)
	if(dim(dat)[2] == 4){
	lineup.dat <- data.frame(x = dat$x, z = dat$cl, cl = dat$cl, .sample = dat$.sample)
	metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, sep.dist = M3sep_dist_indx(pos.1, pos.2, clustering = TRUE, nclust = 2))
	metrics.dat1 <- data.frame(metrics.dat, pic_name = files.png[i])
	metrics1 <- rbind(metrics1, metrics.dat1)
	}
	if(dim(dat)[2] == 6){
	lineup.dat <- data.frame(x = dat$X1, z = dat$X2, cl = dat$cl, .sample = dat$.sample)
	metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, sep.dist = M3sep_dist_indx(pos.1, pos.2, clustering = TRUE, nclust = 3))
	metrics.dat2 <- data.frame(metrics.dat, pic_name = files.png[i])
	metrics2 <- rbind(metrics2, metrics.dat2)
	}
	metrics <- rbind(metrics1, metrics2)
}

#write.table(metrics, "largep-metrics.txt", row.names = F)

#qplot(X1, X2, data = dat, geom = "point", alpha = I(0.7), color = factor(cl)) + facet_wrap(~.sample) + scale_colour_discrete(name = "Group")

res.exp.lp <- read.csv("/Users/Niladri/Documents/Research/Permutation/paper-metrics-data-code/Mahbub's data/raw_data_turk7.csv")

res.exp.lp <- subset(res.exp.lp, pic_name %in% files.png, select = c(pic_name, response, response_no, plot_location, time_taken))

res.dat <- ddply(res.exp.lp, .(pic_name), summarize, prop = sum(response)/length(response), pos = mean(plot_location), m.time = median(time_taken))

metrics.sub <- subset(metrics, pos.1 != pos.2)

dat.merge <- merge(metrics.sub, res.dat, by = "pic_name")

dat.merge <- subset(dat.merge, pos.2 != pos)

dd <- ddply(dat.merge, .(pic_name, pos.1), summarize, sep.mean = mean(sep.dist), len = length(sep.dist), prop = mean(prop), m.time = mean(m.time))

prop.dist <- ddply(dd, .(pic_name), summarize, diff.sep = sep.mean[len == 19] - max(sep.mean[len == 18]), grtr.sep = sum(sep.mean[len == 18] > sep.mean[len == 19]), prop = mean(prop), m.time = mean(m.time))

qplot(diff.sep, prop, data = prop.dist, size = I(3), xlab = "Difference", ylab = "Detection Rate")+ geom_smooth(se = FALSE) + geom_vline(xintercept = 0, col = "red")+ theme(legend.position = "none")

######===================================================================================          
wilks <- function(x){
	summary(manova(cbind(X1, X2)~cl, data = subset(lineup.dat, lineup.dat$.sample == x)), test="Wilks")[[4]][3]
}          
dddd <- ddply(lineup.dat, .(.sample), summarize,  w = wilks(.sample))  
sum(dddd$w - dddd$w[dddd$.sample == 19])/19

wilks_dist <- function(X, PX, clustering = FALSE, nclust = 3) {
    dX <- dist(X[, 1:2])
    dPX <- dist(PX[, 1:2])
    if (clustering) {
        X$cl <- X[, 3]
        PX$cl <- PX[, 3]
        X.clus <- summary(manova(cbind(X[,1], X[,2])~cl, data = X), test="Wilks")[[4]][3]
        PX.clus <- summary(manova(cbind(PX[,1], PX[,2])~cl, data = PX), test="Wilks")[[4]][3]
    } else {
        complete.X <- cutree(hclust(dX), nclust)
        complete.PX <- cutree(hclust(dPX), nclust)
       X.clus <- summary(manova(cbind(X[,1], X[,2])~cl, data = X), test="Wilks")[[4]][3]
        PX.clus <- summary(manova(cbind(PX[,1], PX[,2])~cl, data = PX), test="Wilks")[[4]][3]
    }
    sn.cl <- sign(X.clus - PX.clus)
    ifelse(sn.cl == -1, sqrt(sum((X.clus - PX.clus)^2)),0)
} 

eden <- ddply(dat.pos, .(pos.1, pos.2), summarize, d = wilks_dist(lineup.dat[lineup.dat$.sample == pos.1, ], lineup.dat[lineup.dat$.sample == pos.2, ], clustering = TRUE)) 

eden <- subset(eden, pos.1 != pos.2 & pos.2 != 20)  

ddply(eden, .(pos.1), summarize, mean(d))   

wilks_dist_indx <- function(i, j, clustering = FALSE, nclust = 3){
	X <- lineup.dat[lineup.dat$.sample == i,]
	PX <- lineup.dat[lineup.dat$.sample == j,]
	require(fpc)
	dX <- dist(X[,1:2])
	dPX <- dist(PX[,1:2])
        X$cl <- X[, 3]
        PX$cl <- PX[, 3]
        X.clus <- summary(manova(cbind(X[,1], X[,2])~cl, data = X), test="Wilks")[[4]][3]
        PX.clus <- summary(manova(cbind(PX[,1], PX[,2])~cl, data = PX), test="Wilks")[[4]][3]
    sn.cl <- sign(X.clus - PX.clus)
    ifelse(sn.cl == -1, sqrt(sum((X.clus - PX.clus)^2)),0)
    }	
  