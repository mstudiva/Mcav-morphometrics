setwd("/path/to/working/directory")
getwd()

# install.packages("Rmisc")
# install.packages("stringr")
# install.packages("gridExtra")

library(ggplot2)
library(Rmisc)
library(stringr)
library(gridExtra)

#-------------------
avg=read.csv("template_AVG.csv", head=T)
str(avg)
avg$Site=str_wrap(avg$Site, width=6) # sets site names to maximum width so they wrap on the graph
# importing average dataset 
avg$Site=factor(avg$Site, levels=unique(avg$Site)) # maintains site order as imported
avg$Depth=factor(avg$Depth, levels=unique(avg$Depth)) # maintains depth order as imported


sem=read.csv("template_SEM.csv", head=T)
str(sem)
# importing SEM dataset
sem$Site=factor(sem$Site, levels=unique(sem$Site)) # maintains site order as imported
sem$Depth =factor(sem$Depth, levels=unique(sem$Depth)) # maintains depth order as imported

# creates a manual color palette for graphs
depth_colors= c(shallow="dodgerblue", mesophotic="darkorange")
# depth_colors= c(shallow="grey", mesophotic="grey45")


# plots figure of CD means +/- SEM
# scale_fill_manual forces colors from above color palette
# theme() removes x axis numbers and label
CD=ggplot(avg, aes(x=Site, y=avg$CD, fill=Depth)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=avg$CD-sem$CD, ymax=avg$CD+sem$CD), width=.2, position=position_dodge(.9)) + labs(x="Bank", y="Size (mm)", title="Corallite Diameter (CD)") + theme_bw() + scale_fill_manual(values=depth_colors) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5))

#takes the legend and saves it as a separate object (grob)
get_legend<-function(CD){
  tmp <- ggplot_gtable(ggplot_build(CD))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(CD)

#replots CD, without legend using guides()
CD=ggplot(avg, aes(x=Site, y=avg$CD, fill=Depth)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=avg$CD-sem$CD, ymax=avg$CD+sem$CD), width=.2, position=position_dodge(.9)) + labs(x="Bank", y="Size (mm)", title="Corallite Diameter (CD)") + theme_bw() + scale_fill_manual(values=depth_colors) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + guides(fill=FALSE)

#has both axes labels
TH=ggplot(avg, aes(x=Site, y=avg$TH, fill=Depth)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=avg$TH-sem$TH, ymax=avg$TH +sem$TH), width=.2, position=position_dodge(.9)) + labs(x="Bank", y="Size (mm)", title="Theca Height (TH)") + theme_bw() + scale_fill_manual(values=depth_colors) + theme(plot.title = element_text(hjust=0.5)) + guides(fill=FALSE)

#removes y axis labels
CH=ggplot(avg, aes(x=Site, y=avg$CH, fill=Depth)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=avg$CH-sem$CH, ymax=avg$CH+sem$CH), width=.2, position=position_dodge(.9)) + labs(x="Bank", y="Size (mm)", title="Corallite Height (CH)") + theme_bw() + scale_fill_manual(values=depth_colors) + theme(axis.title.y=element_blank(), plot.title = element_text(hjust=0.5)) + guides(fill=FALSE)

#has no axes labels
CS=ggplot(avg, aes(x=Site, y=avg$CS, fill=Depth)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=avg$CS-sem$CS, ymax=avg$CS +sem$CS), width=.2, position=position_dodge(.9)) + labs(x="Bank", y="Size (mm)", title="Corallite Spacing (CS)") + theme_bw() + scale_fill_manual(values=depth_colors) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(hjust=0.5)) + guides(fill=FALSE)

#creates plot object with a 2x2 grid of plots and legend
#forces layout to be CD, CS, TH, CH, then legend taking up all of column 3 in center
plot=grid.arrange(CD, TH, CS, CH, legend, ncol=3, nrow=2, layout_matrix=cbind(c(1,2),c(3,4),c(5,5)), widths=c(3,2.8,1), heights=c(2.6,3))

#saves plot as PDF and EPS
ggsave("CD_CS_TH_CH_meansfig.pdf", plot= plot, width=10, height=7, units="in", dpi=300)
ggsave("CD_CS_TH_CH_meansfig.eps", plot= plot, width=10, height=7, units="in", dpi=300)

# CW=ggplot(avg, aes(x=Site, y=avg$CW, fill=Depth)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=avg$CW-sem$CW, ymax=avg$CW +sem$CW), width=.2, position=position_dodge(.9)) + labs(x="Bank", y="Size (mm)", title="Columella Width (CW)") + theme_bw() + scale_fill_manual(values=depth_colors) + guides(fill=FALSE)

#-------------------
# normality tests

data = read.csv("template_means.csv", head=TRUE)
str(data)

shapiro.test(data$cd)
shapiro.test(data$cw)
shapiro.test(data$l1s)
shapiro.test(data$l4s)
shapiro.test(data$t1c)
shapiro.test(data$th)
shapiro.test(data$ch)
shapiro.test(data$cs)

# only t1c is normal

# transformations
data$cd.log = log(data$cd)
data$cw.log = log(data$cw)
data$l1s.log = log(data$l1s)
data$l4s.log = log(data$l4s)
data$t1c.log = log(data$t1c)
data$th.log = log(data$th)
data$ch.log = log(data$ch)
data$cs.log = log(data$cs)

shapiro.test(data$cd.log)
shapiro.test(data$cw.log)
shapiro.test(data$l1s.log)
shapiro.test(data$l4s.log)
shapiro.test(data$t1c.log)
shapiro.test(data$th.log)
shapiro.test(data$ch.log)
shapiro.test(data$cs.log)

# l1s, t1c, th, ch, and cs still not normal

#------------------
# histograms for size distribution

# exporting EPS
setEPS()
postscript("metrics_hist.eps", height=6, width=8)
par(mfrow=c(2,3))
hist(data$cd, breaks=50, xlab= NA, main="CD")
hist(data$cw, breaks= 50, xlab=NA, ylab=NA, main ="CW")
hist(data$cs, breaks= 50, xlab= NA,  ylab=NA, main ="CS")
abline(v = 8.0, col="red", lwd=2)
hist(data$th, breaks= 50, xlab="Size in mm",  ylab=NA, main ="TH")
hist(data$ch, breaks= 50, xlab= NA,  ylab=NA, main ="CH")
hist(data$l1s, breaks= 50, xlab="size in mm", main ="L1S")
abline(v = 1.15, col="red", lwd=2)
dev.off()

# exporting PDF
pdf("metrics_hist.pdf", height=6, width=8)
par(mfrow=c(2,3))
hist(data$cd, breaks=50, xlab= NA, main="CD")
hist(data$cw, breaks= 50, xlab=NA, ylab=NA, main ="CW")
hist(data$cs, breaks= 50, xlab= NA,  ylab=NA, main ="CS")
abline(v = 8.0, col="red", lwd=2)
hist(data$th, breaks= 50, xlab="Size in mm",  ylab=NA, main ="TH")
hist(data$ch, breaks= 50, xlab= NA,  ylab=NA, main ="CH")
hist(data$l1s, breaks= 50, xlab="size in mm", main ="L1S")
abline(v = 1.15, col="red", lwd=2)
dev.off()

# hist(data$l4s, breaks= 100, xlab="size in mm", main ="L4S")
# hist(data$t1c, breaks= 100, xlab="size in mm", main ="T1C")

