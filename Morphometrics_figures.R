setwd("/path/to/working/directory")
getwd()

# install.packages("Rmisc")
# install.packages("stringr")
# install.packages("gridExtra")
# install.packages("ggpubr")
# install.packages("FSA")
# install.packages("rcompanion")
# install.packages("dichromat")

library(ggplot2)
library(stringr)
library(gridExtra)
library(ggpubr)
library(Rmisc)
library(FSA)
library(rcompanion)
library(RColorBrewer)

# generates citations for theses, manuscripts, etc.
# citation(package="FSA")
# citation(package="ggpubr")

#-------------------
# data import

data = read.csv("morphometrics_means.csv", head=TRUE)
str(data)

# sets site name to maximum width so they wrap on the graph
# data$site=str_wrap(data$site, width=6) 
# maintains site order as imported
# data$site=factor(data$site, levels=unique(data$site)) 
# maintains depth order as imported
data$depth=factor(data$depth, levels=unique(data$depth)) 

# sets bank name to maximum width so they wrap on the graph
# data$bank=str_wrap(data$bank, width=6) 
# maintains bank order as imported
data$bank=factor(data$bank, levels=unique(data$bank)) 

#------------------
# Kruskal-Wallis tests for each metric
cd.kw <- compare_means(cd ~ bank, data=data, method = "kruskal.test")
cs.kw <- compare_means(cs ~ bank, data=data, method = "kruskal.test")
th.kw <- compare_means(th ~ bank, data=data, method = "kruskal.test")
ch.kw <- compare_means(ch ~ bank, data=data, method = "kruskal.test")

# not used in figure, but reported in manuscript
l1s.kw <- compare_means(l1s ~ bank, data=data, method = "kruskal.test")
l4s.kw <- compare_means(l4s ~ bank, data=data, method = "kruskal.test")
t1c.kw <- compare_means(t1c ~ bank, data=data, method = "kruskal.test")
cw.kw <- compare_means(cw ~ bank, data=data, method = "kruskal.test")

# merging result dataframes and exporting as .csv
kruskal <- rbind(cd.kw,cw.kw,l1s.kw,l4s.kw,t1c.kw,th.kw,ch.kw,cs.kw)
write.csv(kruskal, file="Kruskal_outputs.csv")

# pairwise Dunn tests by bank

cd.dunn <- dunnTest(cd ~ bank, data=data, method = "bh")$res
cs.dunn <- dunnTest(cs ~ bank, data=data, method = "bh")$res
th.dunn <- dunnTest(th ~ bank, data=data, method = "bh")$res
ch.dunn <- dunnTest(ch ~ bank, data=data, method = "bh")$res

# not used in figure, but reported in manuscript
l1s.dunn <- dunnTest(l1s ~ bank, data=data, method = "bh")$res
l4s.dunn <- dunnTest(l4s ~ bank, data=data, method = "bh")$res
t1c.dunn <- dunnTest(t1c ~ bank, data=data, method = "bh")$res
cw.dunn <- dunnTest(cw ~ bank, data=data, method = "bh")$res

# merging result dataframes and exporting as .csv
dunn <- rbind(cd.dunn,cw.dunn,l1s.dunn,l4s.dunn,t1c.dunn,th.dunn,ch.dunn,cs.dunn)
write.csv(dunn, file="Dunn_outputs.csv")

#-------------------
# compact letter display for figure

cd.list <- cldList(P.adj ~ Comparison, data = cd.dunn, threshold = 0.05)
cs.list <- cldList(P.adj ~ Comparison, data = cs.dunn, threshold = 0.05)
th.list <- cldList(P.adj ~ Comparison, data = th.dunn, threshold = 0.05)
ch.list <- cldList(P.adj ~ Comparison, data = ch.dunn, threshold = 0.05)

# not used in figure, but reported in manuscript
l1s.list <- cldList(P.adj ~ Comparison, data = l1s.dunn, threshold = 0.05)
l4s.list <- cldList(P.adj ~ Comparison, data = l4s.dunn, threshold = 0.05)
t1c.list <- cldList(P.adj ~ Comparison, data = t1c.dunn, threshold = 0.05)
cw.list <- cldList(P.adj ~ Comparison, data = cw.dunn, threshold = 0.05)

#-------------------
# creates a manual color palette for graphs
depth_colors= c(shallow=rgb(0.7,0.95,1), mesophotic="dodgerblue")

# creates a function that forces y axis labels to have 1 decimal place
scale <- function(x) sprintf("%.1f", x)

# plots boxplot of cd with sample points
# palette forces colors from above color palette
# theme() removes both axis labels
cd <- ggboxplot(data, x="bank", y="cd", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Corallite Diameter (CD)", fill='Depth') + theme_bw() + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=cd.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
cd

#takes the legend and saves it as a separate object (grob)
get_legend<-function(cd){
  tmp <- ggplot_gtable(ggplot_build(cd))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(cd)

# replots cd, without legend using guides()
cd <- ggboxplot(data, x="bank", y="cd", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Corallite Diameter (CD)") + theme_bw() + theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=cd.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
cd

# has x axis labels
th <- ggboxplot(data, x="bank", y="th", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Theca Height (TH)") + theme_bw() + theme(axis.title.y=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=th.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
th

# has y axis labels
cs <- ggboxplot(data, x="bank", y="cs", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Corallite Spacing (CS)") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=cs.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
cs

# has y axes labels
ch <- ggboxplot(data, x="bank", y="ch", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Corallite Height (CH)") + theme_bw() + theme(plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=ch.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
ch

#creates plot object with a 2x2 grid of plots and legend
#forces layout to be CD, CS, TH, CH, then legend taking up all of column 3 in center
plot=grid.arrange(cs, ch, cd, th, legend, ncol=3, nrow=2, layout_matrix=cbind(c(1,2),c(3,4),c(5,5)), widths=c(3,2.8,1), heights=c(2.6,3))

#saves plot as PDF and EPS
ggsave("CD_CS_TH_CH_meansfig_boxplot.pdf", plot= plot, width=10, height=7, units="in", dpi=300)
ggsave("CD_CS_TH_CH_meansfig_boxplot.eps", plot= plot, width=10, height=7, units="in", dpi=300)

#-------------------
# other metrics for supplemental

# plots boxplot of cw with sample points
# palette forces colors from above color palette
# theme() removes x axis labels
cw <- ggboxplot(data, x="bank", y="cw", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Columella Width (CW)", fill='Depth') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=cw.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
cw

#takes the legend and saves it as a separate object (grob)
get_legend<-function(cw){
  tmp <- ggplot_gtable(ggplot_build(cw))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend=get_legend(cw)

# replots cw, without legend using guides()
cw <- ggboxplot(data, x="bank", y="cw", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Columella Width (CW)", fill='Depth') + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=cw.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
cw

# has no axes labels
t1c <- ggboxplot(data, x="bank", y="t1c", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Thickness 1st Cycle Costa (T1C)") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=t1c.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
t1c

# has both axes labels
l1s <- ggboxplot(data, x="bank", y="l1s", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Length 1st Cycle Septa (L1S)") + theme_bw() + theme(plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=l1s.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
l1s

# has x axis labels
l4s <- ggboxplot(data, x="bank", y="l4s", color="grey30", fill= "depth", palette= depth_colors, add="jitter", width=0.7, size=0.75) + labs(x="Site", y="Size (mm)", title="Length 4th Cycle Septa (L4S)") + theme_bw() + theme(axis.title.y=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), label.x=4.5, hjust=0.5, label.y=0, vjust=-1) +geom_text(data=l4s.list, aes(x = Group, y=0, vjust=-2.5, label=Letter)) + scale_x_discrete(breaks=c("west.mesophotic","west.shallow","east.mesophotic","east.shallow","bright.mesophotic","mcgrail.mesophotic","pulley.mesophotic","tortugas.shallow"),labels=c("WFGB","","EFGB","","BRT","MCG","PRG","DRT")) + scale_y_continuous(labels=scale)
l4s

#creates plot object with a 2x2 grid of plots and legend
#forces layout to be CD, CS, TH, CH, then legend taking up all of column 3 in center
plot2=grid.arrange(cw, l1s, t1c, l4s, legend, ncol=3, nrow=2, layout_matrix=cbind(c(1,2),c(3,4),c(5,5)), widths=c(3,2.8,1), heights=c(2.6,3))

#saves plot as PDF and EPS
ggsave("CW_T1C_L1S_L4S_supp_boxplot.pdf", plot= plot2, width=10, height=7, units="in", dpi=300)
ggsave("CW_T1C_L1S_L4S_supp_boxplot.eps", plot= plot2, width=10, height=7, units="in", dpi=300)

#-------------------
# normality tests

data = read.csv("normality_test.csv", head=TRUE)
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
postscript("metrics_hist.eps", height=6, width=12)
par(mfrow=c(2,3))
hist(data$cd, breaks=50, xlab= NA, main="CD", col="grey50")
hist(data$cw, breaks= 50, xlab=NA, ylab=NA, main ="CW", col="grey50")
hist(data$cs, breaks= 50, xlab= NA,  ylab=NA, main ="CS", col="grey50")
abline(v = 8.4, col="darkorange", lwd=4)
hist(data$th, breaks= 50, xlab="Size (mm)", main ="TH", col="grey50")
hist(data$ch, breaks= 50, xlab= "Size (mm)",  ylab=NA, main ="CH", col="grey50")
hist(data$l1s, breaks= 50, xlab="Size (mm)", main ="L1S", col="grey50")
abline(v = 1.15, col="darkorange", lwd=4)
dev.off()

# exporting PDF
pdf("metrics_hist.pdf", height=6, width=12)
par(mfrow=c(2,3))
hist(data$cd, breaks=50, xlab= NA, main="CD", col="grey50")
hist(data$cw, breaks= 50, xlab=NA, ylab=NA, main ="CW", col="grey50")
hist(data$cs, breaks= 50, xlab= NA,  ylab=NA, main ="CS", col="grey50")
abline(v = 8.4, col="darkorange", lwd=4)
hist(data$th, breaks= 50, xlab="Size (mm)", main ="TH", col="grey50")
hist(data$ch, breaks= 50, xlab= "Size (mm)",  ylab=NA, main ="CH", col="grey50")
hist(data$l1s, breaks= 50, xlab="Size (mm)", main ="L1S", col="grey50")
abline(v = 1.15, col="darkorange", lwd=4)
dev.off()

# hist(data$l4s, breaks= 100, xlab="size in mm", main ="L4S")
# hist(data$t1c, breaks= 100, xlab="size in mm", main ="T1C")

h<-hist(data$cs, breaks= 200, xlab= NA,  ylab=NA, main ="CS", col="grey50")
h
abline(v = 8.25, col="darkorange", lwd=5)
par(mfrow=c(1,1))

d<- density(data$cs)
plot(d)
abline(v = 8.45, col="darkorange", lwd=5)
