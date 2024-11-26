library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(dplyr)
library(viridis)
library(scales)
library(reshape2)
library(stringr)


setwd(" tarabren/singlecell/sc/paper_finalfigures")

int_obj = readRDS(" tarabren/singlecell/weekfeb1/integrated.object_split_feb1.RDS")
int_obj <- RenameIdents(object = int_obj,
                        "0" = "1",
                        "1" = "2",
                        "2" = "3",
                        "3" = "4",
                        "4" = "5",
                        "5" = "6",
                        "6" = "7",
                        "7" = "8",
                        "8" = "9",
                        "9_0" = "10",
                        "9_1" = "11",
                        "10_0" = "12",
                        "10_1" = "13",
                        "10_2_0" = "14",
                        "10_2_1" = "15",
                        "10_3" = "16",
                        "11_0" = "17",
                        "11_1" = "18",
                        "11_2" = "19",
                        "12" = "20",
                        "13" = "21",
                        "14_0" = "22",
                        "14_1" = "23",
                        "15_0" = "24",
                        "15_1" = "25",
                        "16" = "26",
                        "17" = "27")

DefaultAssay(int_obj) <- "RNA"

## cell types ##
int_obj <- RenameIdents(object = int_obj,
                        "1" = "Epithelial",
                        "2" = "Epithelial",
                        "3" = "Epithelial",
                        "4" = "Epithelial",
                        "5" = "Epithelial",
                        "6" = "Fibroblast",
                        "7" = "Epithelial",
                        "8" = "Epithelial",
                        "9" = "Epithelial",
                        "10" = "Epithelial",
                        "11" = "Epithelial",
                        "12" = "Epithelial",
                        "13" = "Epithelial",
                        "14" = "Epithelial",
                        "15" = "Epithelial",
                        "16" = "Epithelial",
                        "17" = "Epithelial",
                        "18" = "Epithelial",
                        "19" = "Epithelial",
                        "20" = "Immune",
                        "21" = "Endothelial",
                        "22" = "Epithelial",
                        "23" = "Epithelial",
                        "24" = "Muscle",
                        "25" = "Muscle",
                        "26" = "Epithelial",
                        "27" = "Muscle")

int_obj$seurat_clusters = Idents(int_obj)

############################# Cell Type Heat Map ############################# 
epimarkers = read.csv("singlecell/sc/celltype_finalmarkers_sep2/epithelial_vsother_markers_thresholds.csv")
epimarkers = sort(epimarkers$X)
endomarkers = read.csv("singlecell/sc/celltype_finalmarkers_sep2/endothelial_vsother_markers_thresholds.csv")
endomarkers = sort(endomarkers$X)
immunemarkers = read.csv("singlecell/sc/celltype_finalmarkers_sep2/immune_vsother_markers_thresholds.csv")
immunemarkers = sort(immunemarkers$X)
fibromarkers = read.csv("singlecell/sc/celltype_finalmarkers_sep2/fibroblast_vsother_markers_thresholds.csv")
fibromarkers = sort(fibromarkers$X)
muscmarkers = read.csv("singlecell/sc/celltype_finalmarkers_sep2/muscle_vsother_markers_thresholds.csv")
muscmarkers = sort(muscmarkers$X)

allmarkers <- c(epimarkers, fibromarkers, immunemarkers, endomarkers, muscmarkers)

# only do d0 cells #
d0subset <- subset(int_obj, subset = grouped == "D0")

heatmap_counts = DoHeatmap(subset(d0subset, downsample = 250), features = allmarkers, group.by="seurat_clusters", angle = 0, slot = 'counts', raster = FALSE, size = 3, disp.max = 6) + scale_fill_gradientn(colors = c("blue", "yellow"))
heatmap_scaledata = DoHeatmap(subset(d0subset, downsample = 250), features = allmarkers, group.by="seurat_clusters", angle = 0, slot = 'scale.data', raster = FALSE, size = 3) #+ scale_fill_gradientn(colors = c("blue","red", "yellow"))

ggsave(heatmap_counts, filename = "celltype_markers_heatmap_counts_D0_6max_250cells_sep2.pdf", width = 10, height = 10, units = "in")
ggsave(heatmap_scaledata, filename = "celltype_markers_heatmap_scaledata_D0_250cells_sep2.pdf", width = 10, height = 10, units = "in")

#d0subsetmuscle = subset(d0subset, subset = seurat_clusters == "Muscle") #ncol to get number of muscle cells



############################# Cell Types Module Score UMAPS ############################# 
setwd("tarabren/singlecell/sc/paper_finalfigures/modulescore_celltype")

epithelial_markers = read.csv("tarabren/singlecell/sc/celltype_finalmarkers_sep2/epithelial_vsother_markers_thresholds.csv")
X_values_list <- list(as.character(epithelial_markers$X))
int_obj <- AddModuleScore(int_obj,features = X_values_list,name = "Epithelial")

endothelial_markers = read.csv("tarabren/singlecell/sc/celltype_finalmarkers_sep2/endothelial_vsother_markers_thresholds.csv")
X_values_list <- list(as.character(endothelial_markers$X))
int_obj <- AddModuleScore(int_obj,features = X_values_list,name = "Endothelial")

immune_markers = read.csv("tarabren/singlecell/sc/celltype_finalmarkers_sep2/immune_vsother_markers_thresholds.csv")
X_values_list <- list(as.character(immune_markers$X))
int_obj <- AddModuleScore(int_obj,features = X_values_list,name = "Immune")

fibroblast_markers = read.csv("/tarabren/singlecell/sc/celltype_finalmarkers_sep2/fibroblast_vsother_markers_thresholds.csv")
X_values_list <- list(as.character(fibroblast_markers$X))
int_obj <- AddModuleScore(int_obj,features = X_values_list,name = "Fibroblast")

muscle_markers = read.csv("tarabren/singlecell/sc/celltype_finalmarkers_sep2/muscle_vsother_markers_thresholds.csv")
X_values_list <- list(as.character(muscle_markers$X))
int_obj <- AddModuleScore(int_obj,features = X_values_list,name = "Muscle")

int_obj@meta.data$pass_thresh <- NA
int_obj@meta.data$pass_thresh <- ifelse(int_obj@meta.data$Epithelial1 > 0, "Epithelial", int_obj@meta.data$pass_thresh)
int_obj@meta.data$pass_thresh <- ifelse(int_obj@meta.data$Endothelial1 > 1, "Endothelial", int_obj@meta.data$pass_thresh)
int_obj@meta.data$pass_thresh <- ifelse(int_obj@meta.data$Immune1 > 1.5, "Immune", int_obj@meta.data$pass_thresh)
int_obj@meta.data$pass_thresh <- ifelse(int_obj@meta.data$Fibroblast1 > .7, "Fibroblast", int_obj@meta.data$pass_thresh)
int_obj@meta.data$pass_thresh <- ifelse(int_obj@meta.data$Muscle1 > 2, "Muscle", int_obj@meta.data$pass_thresh)

d0 <- subset(x = int_obj, subset = grouped == 'D0')
d2 <- subset(x = int_obj, subset = grouped == 'D2')
d4 <- subset(x = int_obj, subset = grouped == 'D4')
d8 <- subset(x = int_obj, subset = grouped == 'D8')

#FeaturePlot(object = d0, features = 'Epithelial1')

d0@meta.data$pass_thresh[is.na(d0@meta.data$pass_thresh)] <- "Unclassified"

bigplot = DimPlot(d0, label=F, group.by="pass_thresh", pt.size =  0.5, cols = c("darkorange", "deeppink", "turquoise4", "#B95FBB", "chartreuse3", "darkgrey"), combine = TRUE)
bigplot = bigplot + ggtitle('Cell Types by Module Score at PBS')
bigplot <- bigplot + theme(
  legend.title = element_text(size = 25),
  legend.position = "none",
  legend.text = element_text(size = 22),
  legend.key.size = unit(1.5, 'cm'),
  plot.title = element_text(size = 39, hjust = 0.5),
  axis.title.x = element_text(size = 45),
  axis.title.y = element_text(size = 45),
  axis.text.x = element_text(size = 45),
  axis.text.y = element_text(size = 45),
)

ggsave(bigplot, filename = "celltypesbymodscore_umap_PBS_intobj.pdf", width = 11, height = 10, units = "in")


############Over Time################
d0threshplot = DimPlot(d0, label=F, group.by="pass_thresh", pt.size =  0.5, cols = c("darkorange", "deeppink", "turquoise4", "#B95FBB", "chartreuse3"), combine = TRUE)
d0threshplot = d0threshplot + ggtitle('PBS')

d2threshplot = DimPlot(d2, label=F, group.by="pass_thresh", pt.size =  0.5, cols = c("darkorange", "deeppink", "turquoise4", "#B95FBB", "chartreuse3"), combine = TRUE)
d2threshplot = d2threshplot + ggtitle('Day2')

d4threshplot = DimPlot(d4, label=F, group.by="pass_thresh", pt.size =  0.5, cols = c("darkorange", "deeppink", "turquoise4", "#B95FBB", "chartreuse3"), combine = TRUE)
d4threshplot = d4threshplot + ggtitle('Day4')

d8threshplot = DimPlot(d8, label=F, group.by="pass_thresh", pt.size =  0.5, cols = c("darkorange", "deeppink", "turquoise4", "#B95FBB", "chartreuse3"), combine = TRUE)
d8threshplot = d8threshplot + ggtitle('Day8')

myTheme <- theme(axis.text = element_text(size = 16),
                 legend.title = element_text(size = 25),
                 legend.position = "none",
                 legend.text = element_text(size = 22),
                 legend.key.size = unit(1.5, 'cm'),
                 plot.title = element_text(size = 35, hjust = 0.5),
                 axis.title.x = element_text(size = 28),
                 axis.title.y = element_text(size = 28),
                 axis.text.x = element_text(size = 27),
                 axis.text.y = element_text(size = 27),
)

d0threshplot1 = d0threshplot + myTheme
d2threshplot1 = d2threshplot +myTheme
d4threshplot1 = d4threshplot + myTheme
d8threshplot1 = d8threshplot + myTheme

finalplot = ggarrange(d0threshplot1, d2threshplot1, d4threshplot1, d8threshplot1, ncol=2, nrow=2, common.legend = T, legend="right")
ggsave(finalplot, filename = "celltypesbymodscore_alldays_upmaps_intobj.pdf", width = 16, height = 13)

############################# Stacked Bar Plots of Celltype ############################# 
pt <- table(int_obj$orig.ident, int_obj$pass_thresh)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

library(stringr)
pt <- pt %>%
  mutate(Var1 = case_when(
    str_starts(Var1, "Day0") ~ "PBS",
    str_starts(Var1, "Day2") ~ "Day2",
    str_starts(Var1, "Day8") ~ "Day8",
    TRUE ~ Var1
  ))

summarized_data <- pt %>%
  group_by(Var1, Var2) %>%
  summarise(Freq = sum(Freq))

group_totals <- summarized_data %>%
  group_by(Var1) %>%
  summarise(total = sum(Freq))

overall_total <- sum(group_totals$total)
library(RColorBrewer)

custom_colors <- c("darkorange", "deeppink", "turquoise4", "#B95FBB", "chartreuse3")
summarized_data$Var1 <- factor(summarized_data$Var1, levels = c("Day8", "Day4", "Day2", "PBS"))
summarized_data <- summarized_data[order(summarized_data$Var1), ]

plotd0 = ggplot((summarized_data), aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = custom_colors) + 
  theme(legend.title = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(size = 47),
        axis.title.x = element_text(size = 45),
        axis.text.y = element_text(size = 47),
        axis.title.y = element_text(size = 45),
        legend.text = element_text(size = 39),
        legend.key.size = unit(1.2, 'cm')) + coord_flip()

#finalplot = plotd0 + geom_text(aes(Var1, total, label = paste("n=", total), fill = NULL), data = group_totals, position = position_fill(vjust = 1.02), size = 7)

ggsave(plotd0, filename = "stackedbarplot_celltype_nolabels.pdf", width = 14, height = 8, units = "in")

############################# Epithelial PBS UMAP ############################# 
epiobj = subset(int_obj, idents = c("1","2","3","4","5","7","8","9","10","11","12","13","14","15","16","17","18","19","22","23","26"))
my_cols <- c('1'='#2FF18B','2'='#FF61C5','3'='#1FA195','4'='#ff9a36','5'='#D4D915',
             '6'='seagreen4','7'='mistyrose4','8'='#F68282','9'='#B95FBB','10'='#31C53F',
             '11'='#CCB1F9','12'='skyblue4','13'='brown4','14'='magenta1','15'='royalblue4',
             '16'='#E6C122', '17'='#4B4BF7','18'='#28CECA','19'='bisque4','20'='#6A9AFF','21'='purple4',
             '22'='yellow4','23'='turquoise4','24'='tan2','25'='deeppink2','26'='purple' ,'27'='firebrick3')

my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
d0 <- subset(x = epiobj, subset = grouped == 'D0')
d0plot = DimPlot(d0,
                 cols = my_cols2, label=F, pt.size = .5, combine = T)
d0plot <- d0plot + ylim(-12, 11) + xlim(-10,10)
d0plot <- d0plot + ggplot2::ggtitle("PBS")
#d0plot = LabelClusters(d0plot, id = "ident",  fontface = "bold", color = "black",size = )

myTheme <- theme(legend.title = element_text(size = 25),
                 legend.position = "none",
                 legend.text = element_text(size = 22),
                 legend.key.size = unit(1.5, 'cm'),
                 plot.title = element_text(size = 39, hjust = 0.5),
                 axis.title.x = element_text(size = 45),
                 axis.title.y = element_text(size = 45),
                 axis.text.x = element_text(size = 45),
                 axis.text.y = element_text(size = 45),
)

d0plot = d0plot + myTheme
ggsave(d0plot, filename = "epithelial_umap_nolabels.pdf", width = 11, height = 11, units = "in")

############################# Module Score Basal Suprabasal ############################# 
setwd("tarabren/singlecell/sc/nov2/supra_vs_all_analysis/")
int_obj = subset(int_obj, idents = c("1","2","3","4","5","7","8","9","10","11","12","13","14","15","16","17","18","19","22","23","26"))

basal_markers = list(c("Itga6","Itgb4","Prxl2a","Stmn1","Lamb3","Wfdc2"))
suprabasal_markers = list(c("Krt4","Tgm3","Dsg1a","Dsc2","Lypd3","Sbsn","Rab11a","Gltp","Chit1","Degs1"))

#basal_markers = list(c("Krt14","Krt5","Itgb4"))
#suprabasal_markers = list(c("Krt13","Krt4"))
int_obj <- AddModuleScore(int_obj,features = basal_markers,name = "BasalMarkers")
int_obj <- AddModuleScore(int_obj,features = suprabasal_markers,name = "SuprabasalMarkers")

basalplot = FeaturePlot(object = int_obj, features = "BasalMarkers1") 
suprabasalplot = FeaturePlot(object = int_obj, features = "SuprabasalMarkers1") 

int_obj@meta.data$pass_thresh = NA
int_obj@meta.data$pass_thresh <- ifelse(int_obj@meta.data$BasalMarkers1 > -.6, "BasalMarkers", int_obj@meta.data$pass_thresh)
int_obj@meta.data$pass_thresh <- ifelse(int_obj@meta.data$SuprabasalMarkers1 > -.05, "SuprabasalMarkers", int_obj@meta.data$pass_thresh)

int_obj@meta.data$pass_thresh[is.na(int_obj@meta.data$pass_thresh)] <- "Unclassified"

d0 <- subset(x = int_obj, subset = grouped == 'D0')
d2 <- subset(x = int_obj, subset = grouped == 'D2')
d4 <- subset(x = int_obj, subset = grouped == 'D4')
d8 <- subset(x = int_obj, subset = grouped == 'D8')

threshplot = DimPlot(d0, label=F, group.by="pass_thresh", pt.size =  0.7, cols = c("#0072B2", "maroon", "darkgrey"), combine = TRUE)

myTheme <- theme(legend.title = element_text(size = 17),
                 legend.position = "none",
                 legend.text = element_text(size = 22),
                 legend.key.size = unit(1.5, 'cm'),
                 plot.title = element_text(size = 31, hjust = 0.5),
                 axis.title.x = element_text(size = 45),
                 axis.title.y = element_text(size = 45),
                 axis.text.x = element_text(size = 45),
                 axis.text.y = element_text(size = 45),
)

threshplot = threshplot + ggtitle("Basal and Suprabasal by Modscore on PBS") + myTheme + ylim(-12, 11) + xlim(-10,10)
ggsave(threshplot, filename = "basalsuprabasal_bymodscore_PBS.pdf", width = 11, height = 11, units = "in")

############Over Time################
d0threshplot = DimPlot(d0, label=F, group.by="pass_thresh", pt.size =  0.5, cols = c("#0072B2", "maroon", "darkgrey"), combine = TRUE)
d0threshplot = d0threshplot + ggtitle('PBS') + ylim(-12, 11) + xlim(-10,10)

d2threshplot = DimPlot(d2, label=F, group.by="pass_thresh", pt.size =  0.5, cols = c("#0072B2", "maroon", "darkgrey"), combine = TRUE)
d2threshplot = d2threshplot + ggtitle('Day 2') + ylim(-12, 11) + xlim(-10,10)

d4threshplot = DimPlot(d4, label=F, group.by="pass_thresh", pt.size =  0.5, cols = c("#0072B2", "maroon", "darkgrey"), combine = TRUE)
d4threshplot = d4threshplot + ggtitle('Day 4') + ylim(-12, 11) + xlim(-10,10)

d8threshplot = DimPlot(d8, label=F, group.by="pass_thresh", pt.size =  0.5, cols = c("#0072B2", "maroon", "darkgrey"), combine = TRUE)
d8threshplot = d8threshplot + ggtitle('Day 8') + ylim(-12, 11) + xlim(-10,10)

myTheme <- theme(axis.text = element_text(size = 16),
                 legend.title = element_text(size = 25),
                 legend.position = "none",
                 legend.text = element_text(size = 22),
                 legend.key.size = unit(1.5, 'cm'),
                 plot.title = element_text(size = 35, hjust = 0.5),
                 axis.title.x = element_text(size = 28),
                 axis.title.y = element_text(size = 28),
                 axis.text.x = element_text(size = 27),
                 axis.text.y = element_text(size = 27),
)

d0threshplot1 = d0threshplot + myTheme
d2threshplot1 = d2threshplot +myTheme
d4threshplot1 = d4threshplot + myTheme
d8threshplot1 = d8threshplot + myTheme

finalplot = ggarrange(d0threshplot1, d2threshplot1, d4threshplot1, d8threshplot1, ncol=2, nrow=2, common.legend = T, legend="right")
ggsave(finalplot, filename = "basalsupra_bymodscore_alldays_upmaps.pdf", width = 16, height = 13)

############################# Stacked bar plot for Module Score Basal Suprabasal ############################# 
pt <- table(int_obj$orig.ident, int_obj$pass_thresh)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

library(stringr)
pt <- pt %>%
  mutate(Var1 = case_when(
    str_starts(Var1, "Day0") ~ "PBS",
    str_starts(Var1, "Day2") ~ "Day2",
    str_starts(Var1, "Day8") ~ "Day8",
    TRUE ~ Var1
  ))

summarized_data <- pt %>%
  group_by(Var1, Var2) %>%
  summarise(Freq = sum(Freq))

group_totals <- summarized_data %>%
  group_by(Var1) %>%
  summarise(total = sum(Freq))

overall_total <- sum(group_totals$total)

custom_colors <- c("#0072B2", "maroon")
summarized_data$Var1 <- factor(summarized_data$Var1, levels = c("Day8", "Day4", "Day2", "PBS"))
summarized_data <- summarized_data[order(summarized_data$Var1), ]

plotd0 = ggplot((summarized_data), aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = custom_colors) + 
  theme(legend.title = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(size = 47),
        axis.title.x = element_text(size = 45),
        axis.text.y = element_text(size = 47),
        axis.title.y = element_text(size = 45),
        legend.text = element_text(size = 39),
        legend.key.size = unit(1.2, 'cm')) + coord_flip()

ggsave(plotd0, filename = "stackedbarplot_basalsupra_modscore.pdf", width = 14, height = 8, units = "in")

############################# Bubble Plots of basal/suprabasal markers ############################# 
setwd("tarabren/singlecell/sc/nov2/supra_vs_all_analysis/")

d2 <- subset(x = int_obj, subset = grouped == 'D2')
d4 <- subset(x = int_obj, subset = grouped == 'D4')
d8 <- subset(x = int_obj, subset = grouped == 'D8')

basal_markers = unlist(basal_markers)
suprabasal_markers = unlist(suprabasal_markers)

allmarkers = c(basal_markers, suprabasal_markers)

allbubble = DotPlot(d0, features = allmarkers) + RotatedAxis()
ggsave(allbubble, file = "allmarkes_supravsbasal_bubbleplot.pdf", width = 8, height = 6)

# dot plots over times
d0plot = DotPlot(d0, features = allmarkers) + RotatedAxis() +ggtitle("PBS")
d2plot = DotPlot(d2, features = allmarkers) + RotatedAxis() +ggtitle("Day 2")
d4plot = DotPlot(d4, features = allmarkers) + RotatedAxis() + ggtitle("Day 4")
d8plot = DotPlot(d8, features = allmarkers) + RotatedAxis() + ggtitle("Day 8")

finalplot = ggarrange(d0plot, d2plot, d4plot, d8plot, ncol=2, nrow=2, common.legend = TRUE, legend="right")
finalplot
ggsave(finalplot, file = "allmarkes_throughdays_bubbleplot.pdf", width = 13, height = 11)

############################# Heatmaps of basal/suprabasal markers ############################# 
basal_markers = c("Krt14","Krt5","Itgb4")
suprabasal_markers = c("Krt13","Krt4")

allmarkers <- c(basal_markers, suprabasal_markers)

d0subset <- subset(int_obj, subset = grouped == "D0")
d0subset = subset(d0subset, idents = c("1","2","3","4","5","7","8","9","10","11","12","13","14","15","16","17","18","19","22","23","26"))
d0subset$seurat_clusters = Idents(d0subset)

heatmap_counts = DoHeatmap(subset(d0subset, downsample = 200), features = allmarkers, group.by="seurat_clusters", angle = 0, slot = 'counts', raster = FALSE, size = 3, disp.max = 10) + scale_fill_gradientn(colors = c("blue", "yellow"))
heatmap_counts = heatmap_counts + theme(legend.position = "none",
                                        axis.text.y = element_text(size = 29),
                                        axis.title.y = element_text(size = 29),
                                        ) 

ggsave(heatmap_counts, filename = "basalsupra_markers_heatmap_rawcounts_10max.pdf", width = 11, height = 4, units = "in")

############################# Heatmap of Cell count proportions by day ############################# 
int_obj = subset(int_obj, idents = c("1","2","3","4","5","7","8","9","10","11","12","13","14","15","16","17","18","19","22","23","26"))
pt <- table(Idents(int_obj), int_obj$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt <- pt %>%
  mutate(Var2 = case_when(
    str_starts(Var2, "Day0") ~ "Day0",
    str_starts(Var2, "Day2") ~ "Day2",
    str_starts(Var2, "Day4") ~ "Day4",
    str_starts(Var2, "Day8") ~ "Day8",
    TRUE ~ Var1
  ))

summarized_data <- pt %>%
  group_by(Var1, Var2) %>%
  summarise(Freq = sum(Freq))

group_totals <- summarized_data %>%
  group_by(Var2) %>%
  summarise(total = sum(Freq))

overall_total <- sum(group_totals$total)

summarized_data$Var1 <- factor(summarized_data$Var1, levels = as.character(1:27))

summarized_data = as.data.frame(summarized_data)
colnames(summarized_data) = c("Cluster","Day", "N_Cells")

summarized_data$N_Cells <- summarized_data$N_Cells + 1

group_totals <- group_totals %>% rename(Day = Var2)
summarized_data <- summarized_data %>%
  left_join(group_totals, by = "Day")

summarized_data <- summarized_data %>%
  mutate(proportion = N_Cells / total)
summarized_data$proportion <- summarized_data$proportion + 0.0001

summarized_data <- summarized_data %>%
  mutate(log_proportion = log(N_Cells / total))

summarized_data <- summarized_data %>%
  group_by(Cluster) %>%
  mutate(normalized_proportion = proportion / sum(proportion)) %>%
  ungroup()

summarized_data <- summarized_data %>%
  group_by(Cluster) %>%
  mutate(log2_fold_change = log2(proportion / proportion[Day == "Day0"])) %>%
  ungroup()

summarized_data <- summarized_data %>%
  mutate(Day = ifelse(Day == "Day0", "PBS", Day))
summarized_data = as.data.frame(summarized_data)

summarized_data <- summarized_data %>%
  mutate(Day = factor(Day, levels = c("PBS", "Day2", "Day4", "Day8")))

finalplot = ggplot(summarized_data, aes(x = Day, y = factor(Cluster), fill = proportion)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Proportions by Cluster and Day",
       x = "Day",
       y = "Cluster",
       fill = "Proportion of Total Cells") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 36),
        legend.title = element_text(size = 29), 
        axis.text.x = element_text(size = 34),
        axis.text.y = element_text(size = 31),
        axis.title.y = element_text(size = 34),
        legend.text = element_text(size = 32),
        legend.key.size = unit(1.6, 'cm'))

ggsave(finalplot, filename = "epiproportions_log_proportions.pdf", width = 11, height = 10, units = "in")

########## split cluster 18 based on markers 

markers18 = list(c("Aqp5", "Kcnn4", "Bhlha15"))

int_obj <- AddModuleScore(int_obj,features = markers18, name = "Markers18_")

threshold <- 0.5

cluster_18_cells <- WhichCells(int_obj, idents = "18")
markers18_scores <- FetchData(int_obj, vars = "Markers18_1")

FeaturePlot(int_obj, features = "Markers18_1")

int_obj@meta.data$pass_thresh = NA
int_obj@meta.data$pass_thresh <- ifelse(int_obj@meta.data$Markers18 > 1, "18_1", int_obj@meta.data$pass_thresh)
int_obj@meta.data$pass_thresh <- ifelse(int_obj@meta.data$Markers18 < 1, "18_2", int_obj@meta.data$pass_thresh)

int_obj@meta.data$subcluster <- NA
int_obj@meta.data$subcluster[cluster_18_cells] <- ifelse(markers18_scores[cluster_18_cells, "Markers18"] > threshold, "18_1", "18_2")

Idents(int_obj) <- "subcluster"






