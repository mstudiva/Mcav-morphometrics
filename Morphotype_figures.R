setwd("/path/to/working/directory")
getwd()

# install.packages("Rmisc")
# install.packages("stringr")
# install.packages("gridExtra")
# install.packages("ggpubr")

library(ggplot2)
library(stringr)
library(gridExtra)
library(ggpubr)

#-------------------
# data import

data = read.csv("morphotype_means.csv", head=TRUE)
str(data)

# converting zoox_cm to millions of cells
data$zoox_cm <- (data$zoox_cm/1e6)

# maintains morphotype order as imported
#data$morphotype =factor(data$morphotype, levels(data$morphotype)[c(2,1)]) 
#data$morphotype =factor(data$morphotype, levels(data$morphotype)[c(2,1)]) 

#------------------
# Wilcoxon test for p values

# morphometrics
cs.p <- compare_means(cs ~ morphotype, data=data)
cd.p <- compare_means(cd ~ morphotype, data=data)
ch.p <- compare_means(ch ~ morphotype, data=data)
l1s.p <- compare_means(l1s ~ morphotype, data=data)

# zoox metrics
zoox_cm.p <- compare_means(zoox_cm ~ morphotype, data=data)
chla_cm.p <- compare_means(chla_cm ~ morphotype, data=data)
chlc_cm.p <- compare_means(chlc_cm ~ morphotype, data=data)
chla_cell.p <- compare_means(chla_cell ~ morphotype, data=data)
chlc_cell.p <- compare_means(chlc_cell ~ morphotype, data=data)
chla_c.p <- compare_means(chla_c ~ morphotype, data=data)

# not used in figs
th.p <- compare_means(th ~ morphotype, data=data)
cw.p <- compare_means(cw ~ morphotype, data=data)
l4s.p <- compare_means(l4s ~ morphotype, data=data)
t1c.p <- compare_means(t1c ~ morphotype, data=data)

#------------------
#creates a manual color palette for graphs
morphotype_colors = c("depth-generalist"="coral", "shallow-only"="cyan3")

# creates a function that forces y axis labels to have 1 decimal place
scale <- function(x) sprintf("%.1f", x)

# plots boxplot of cs with sample points
# palette forces colors from above color palette
# theme() removes x axis label
cs <- ggboxplot(data, x="morphotype", y="cs", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.75, size=0.75) + labs(x="Morphotype", y="Size (mm)", title="Corallite Spacing (CS)") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
cs

# theme() removes y axis label
l1s <- ggboxplot(data, x="morphotype", y="l1s", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.75, size=0.75) + labs(x="Morphotype", y="Size (mm)", title="Length 1st Cycle Septa (L1S)") + theme_bw() + theme(plot.title = element_text(hjust=0.5), legend.position = "none", axis.title.y=element_blank()) + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
l1s

# theme() removes x axis label and y axis label
cd <- ggboxplot(data, x="morphotype", y="cd", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.75, size=0.75) + labs(x="Morphotype", y="Size (mm)", title="Corallite Diameter (CD)") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
cd

# theme() removes x axis label and y axis label
ch <- ggboxplot(data, x="morphotype", y="ch", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.75, size=0.75) + labs(x="Morphotype", y="Size (mm)", title="Corallite Height (CH)") + theme_bw() + theme(plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
ch

# creates plot object with a 2x2 grid of plots
# forces layout to be cs, cd, ch, l1s
plot1=grid.arrange(cs, ch, cd, l1s, ncol=2, nrow=2, layout_matrix=cbind(c(1,2),c(3,4)), widths=c(2.2,2), heights=c(2.6,3))

#saves plot as PDF and EPS
ggsave("morphotype_meansfig_box.pdf", plot= plot1, width=6, height=6, units="in", dpi=300)
ggsave("morphotype_meansfig_box.eps", plot= plot1, width=6, height=6, units="in", dpi=300)

#-------------------------
# not used in fig

# both axis labels
th <- ggboxplot(data, x="morphotype", y="th", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.75, size=0.75) + labs(x="Morphotype", y="Size (mm)", title="Theca Height (TH)") + theme_bw() + theme(plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
th

# theme() removes x axis label
cw <- ggboxplot(data, x="morphotype", y="cw", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.75, size=0.75) + labs(x="Morphotype", y="Size (mm)", title="Columella Width (CW)") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
cw

# theme() removes x axis label and y axis label
t1c <- ggboxplot(data, x="morphotype", y="t1c", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.75, size=0.75) + labs(x="Morphotype", y="Size (mm)", title="Thickness 1st Cycle Costa (T1C)") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
t1c

# theme() removes y axis label
l4s <- ggboxplot(data, x="morphotype", y="l4s", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.75, size=0.75) + labs(x="Morphotype", y="Size (mm)", title="Length 4th Cycle Septa (L4S)") + theme_bw() + theme(axis.title.y=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1)+ scale_y_continuous(labels=scale)
l4s

# creates plot object with a 2x2 grid of plots
# forces layout to be cw, th, l4s, t1c
plot2=grid.arrange(cw, th, t1c, l4s, ncol=2, nrow=2, layout_matrix=cbind(c(1,2),c(3,4)), widths=c(2.2,2), heights=c(2.6,3))

#saves plot as PDF and EPS
ggsave("morphotype_supp_box.pdf", plot= plot2, width=6, height=6, units="in", dpi=300)
ggsave("morphotype_supp_box.eps", plot= plot2, width=6, height=6, units="in", dpi=300)

#------------------------
# zoox metrics

# theme() removes x axis label
zoox_cm <- ggboxplot(data, x="morphotype", y="zoox_cm", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.7, size=0.75) + labs(x="Morphotype", y="Symbiont Density \n(millions cells / sq. cm)", title="Symbiont Density") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
zoox_cm

# theme() removes x axis label
chla_cm <- ggboxplot(data, x="morphotype", y="chla_cm", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.7, size=0.75) + labs(x="Morphotype", y="Areal Chl \n(ug / sq. cm)", title="Areal Chl a") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale) 
chla_cm

# theme() removes x axis label and y axis label
chlc_cm <- ggboxplot(data, x="morphotype", y="chlc_cm", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.7, size=0.75) + labs(x="Morphotype", y="Areal Chl \n(ug / sq. cm)", title="Areal Chl c2") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
chlc_cm

# theme() removes x axis label and y axis label
chla_c <- ggboxplot(data, x="morphotype", y="chla_c", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.7, size=0.75) + labs(x="Morphotype", y="Chl a:c2", title="Chl a:c2") + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
chla_c

# both axis labels
chla_cell <- ggboxplot(data, x="morphotype", y="chla_cell", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.7, size=0.75) + labs(x="Morphotype", y="Cellular Chl \n(pg / cell)", title="Cellular Chl a") + theme_bw() + theme(plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale)
chla_cell

# both axis labels
chlc_cell <- ggboxplot(data, x="morphotype", y="chlc_cell", color="grey30", fill= "morphotype", palette= morphotype_colors, add="jitter", width=0.7, size=0.75) + labs(x="Morphotype", y="Cellular Chl \n(pg / cell)", title="Cellular Chl c2") + theme_bw() + theme(axis.title.y=element_blank(), plot.title = element_text(hjust=0.5), legend.position = "none") + stat_compare_means(aes(label = paste0("p=", ..p.format..)), paired=FALSE, label.x=1.5, hjust=0.5, label.y=0, vjust=-1) + scale_y_continuous(labels=scale) 
chlc_cell

# creates plot object with a 2x3 grid of plots
# forces layout to be zoox_cm, chla_c, chla_cm, chlc_cm, chla_cell, chlc_cell
plot3=grid.arrange(zoox_cm, chla_cm, chla_cell, chla_c, chlc_cm, chlc_cell, ncol=2, nrow=3, layout_matrix=cbind(c(1,2,3),c(4,5,6)), widths=c(2.2,2), heights=c(2.6,2.6,3))

#saves plot as PDF and EPS
ggsave("morphotype_zoox_box.pdf", plot= plot3, width=6, height=8, units="in", dpi=300)
ggsave("morphotype_zoox_box.eps", plot= plot3, width=6, height=8, units="in", dpi=300)

