setwd("~/Documents/Research/Geukensia/salmon_quants/")

rm(list = ls())
library(edgeR)
library(data.table)        
library(ggplot2)
library(cowplot)

##FULL DATA
files = c("GA1-30_quant.sf","GA3-30_quant.sf","GA4-30_quant.sf","GA5-30_quant.sf",
          "GA6-20_quant.sf","GA7-20_quant.sf","GA8-20_quant.sf","GA9-20_quant.sf","GA10-20_quant.sf",
          "GAC1_quant.sf","GAC2_quant.sf","GAC3_quant.sf","GAC4_quant.sf","GAC5_quant.sf",
          "MA12-30_quant.sf","MA13-30_quant.sf","MA14-30_quant.sf", "MA15-30_quant.sf",
          "MA16-20_quant.sf","MA17-20_quant.sf","MA18-20_quant.sf","MA20-20_quant.sf",
          "MAC1_quant.sf","MAC2_quant.sf","MAC3_quant.sf","MAC5_quant.sf")
geukData <- readDGE(files = files, path = "~/Documents/Research/Geukensia/salmon_quants",
                    columns = c(1,5), group = c("GA30","GA30","GA30","GA30",
                                                "GA20","GA20","GA20","GA20","GA20",
                                                "GAC","GAC","GAC","GAC","GAC",
                                                "MA30","MA30","MA30", "MA30",
                                                "MA20","MA20","MA20","MA20",
                                                "MAC","MAC","MAC","MAC"))

keep <- filterByExpr(geukData, group = geukData$samples$group)
keep <- rowSums(cpm(geukData) > 1) >= 2
#GLM approach
geukData <- geukData[keep, , keep.lib.sizes=FALSE] 
geukData <- calcNormFactors(geukData, method="TMM")
#generate design matrix
design <- model.matrix(~0 + geukData$samples$group)
colnames(design) <- levels(geukData$samples$group)
geukData <- estimateDisp(geukData,design)
geukData$common.dispersion
fit <-glmQLFit(geukData, design)

##plot showing relationship of samples to each other
pch <- rep(c(0,1,2,15,16,17))
colors <- rep(c("darkgreen","red","blue"),2)
par(mar=c(5,4,4,9),xpd=TRUE)
recPCA <- plotMDS(geukData, col=colors[geukData$samples$group], pch=pch[geukData$samples$group], 
                  top = nrow(geukData), gene.selection = 'common')
legend("topright", inset=c(-0.6,0), legend=levels(geukData$samples$group), pch=pch, col=colors, ncol=2)

#make comparisons across the different treatments
my.contrasts <- makeContrasts(
  GACvsMAC = GAC-MAC,
  GA20vsMA20 = GA20-MA20,
  GA30vsMA30 = GA30-MA30,
  GACvsGA20 = GAC-GA20,
  GACvsGA30 = GAC-GA30,
  MACvsMA20 = MAC-MA20,
  MACvsMA30 = MAC-MA30,
  GA20vsGA30 = GA20-GA30,
  MA20vsMA30 = MA20-MA30, 
  levels=design
)

thresh_low <- 0.05
#find genes that are different between the control groups

tr_GACvsMAC <- glmTreat(fit, contrast=my.contrasts[,"GACvsMAC"], lfc=log2(1.5)) #looks at genes above a certain log fold change threshold
tt_GACvsMAC <- topTags(tr_GACvsMAC, n=Inf)
#topids_GACvsMAC <- tt_GACvsMAC$table[tt_GACvsMAC$table$FDR<thresh_low]
#topids_GACvsMAC.df <- data.frame(topids_GACvsMAC)
tt_GACvsMAC.df <- data.frame(tt_GACvsMAC)
#tr_GACvsMAC.df <- data.frame(tr_GACvsMAC)
GACvsMAC_logFDR <- -log10(tt_GACvsMAC.df$FDR)
GACvsMAC_logFDR.df <- data.frame(GACvsMAC_logFDR)
tt_GACvsMAC.df <- cbind(tt_GACvsMAC.df,GACvsMAC_logFDR.df)
write.table(topTags(tr_GACvsMAC,n=30),
            file = "~/Documents/Research/Geukensia/GACvsMAC_tt_summary.txt",
            quote=FALSE)
is.de_GACvsMAC <- decideTestsDGE(tr_GACvsMAC, adjust.method="fdr",p.value=0.05)
summary(is.de_GACvsMAC)
GACvsMACplot <- ggplot(tt_GACvsMAC.df) +aes(x = logFC,y = GACvsMAC_logFDR,colour = GACvsMAC_logFDR) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_gradient(low = "#000000", high = "#000000") +
  labs(x = "log2 Fold Change",y = "-log10(FDR)") +
  theme_bw() + geom_hline(yintercept=1.3) +
  theme_minimal(base_size = 22) + theme(legend.position = "none")


#find genes that are responding differently in the 20 degree temp for GA and MA

tr_GA20vsMA20 <- glmTreat(fit, contrast=my.contrasts[,"GA20vsMA20"], lfc=log2(1.5)) #looks at genes above a certain log fold change threshold
tt_GA20vsMA20 <- topTags(tr_GA20vsMA20, n=Inf)
#topids_GA20vsMA20 <- tt_GA20vsMA20$table[tt_GA20vsMA20$table$FDR<thresh_low]
#topids_GA20vsMA20.df <- data.frame(topids_GA20vsMA20)
tt_GA20vsMA20.df <- data.frame(tt_GA20vsMA20)
tr_GA20vsMA20.df <- data.frame(tr_GA20vsMA20)
GA20vsMA20_logFDR <- -log10(tt_GA20vsMA20.df$FDR)
GA20vsMA20_logFDR.df <- data.frame(GA20vsMA20_logFDR)
tt_GA20vsMA20.df <- cbind(tt_GA20vsMA20.df,GA20vsMA20_logFDR.df)
write.table(topTags(tr_GA20vsMA20,n=30),
            file = "~/Documents/Research/Geukensia/GA20vsMA20_tt_summary.txt",
            quote=FALSE)
summary(decideTestsDGE(tr_GA20vsMA20, adjust.method="fdr",p.value=0.05))
GA20vsMA20plot <- ggplot(tt_GA20vsMA20.df) +aes(x = logFC,y = GA20vsMA20_logFDR,colour = GA20vsMA20_logFDR) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_gradient(low = "#000000", high = "#000000") +
  labs(x = "log2 Fold Change",y = "-log10(FDR)") +
  theme_bw() + geom_hline(yintercept=1.3) + 
  theme_minimal(base_size = 22) + theme(legend.position = "none")


#find genes that are responding differently in 30 degree temp for GA and MA

tr_GA30vsMA30 <- glmTreat(fit, contrast=my.contrasts[,"GA30vsMA30"], lfc=log2(1.5)) #looks at genes above a certain log fold change threshold
tt_GA30vsMA30 <- topTags(tr_GA30vsMA30, n=Inf)
#topids_GA30vsMA30 <- tt_GA30vsMA30$table[tt_GA30vsMA30$table$FDR<thresh_low]
#topids_GA30vsMA30.df <- data.frame(topids_GA30vsMA30)
tt_GA30vsMA30.df <- data.frame(tt_GA30vsMA30)
#tr_GA30vsMA30.df <- data.frame(tr_GA30vsMA30)
GA30vsMA30_logFDR <- -log10(tt_GA30vsMA30.df$FDR)
GA30vsMA30_logFDR.df <- data.frame(GA30vsMA30_logFDR)
tt_GA30vsMA30.df <- cbind(tt_GA30vsMA30.df,GA30vsMA30_logFDR.df)
write.table(topTags(tr_GA30vsMA30,n=30),
            file = "~/Documents/Research/Geukensia/GA20vsMA20_tt_summary.txt",
            quote=FALSE)
summary(decideTestsDGE(tr_GA30vsMA30, adjust.method="fdr",p.value=0.05))
GA30vsMA30plot <- ggplot(tt_GA30vsMA30.df) +aes(x = logFC,y = GA30vsMA30_logFDR,colour = GA30vsMA30_logFDR) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_gradient(low = "#000000", high = "#000000") +
  labs(x = "log2 Fold Change",y = "-log10(FDR)") +
  theme_bw() + geom_hline(yintercept=1.3) +
  theme_minimal(base_size = 22) + theme(legend.position = "none")

#save plots for constitutive differences figure
cons_plots <- plot_grid(GACvsMACplot,GA20vsMA20plot,GA30vsMA30plot,labels="AUTO")
save_plot("Cons_Differences_withMA13-30.png",cons_plots,ncol=2,nrow=2)

#find genes responding differently between GA temps

tr_GA20vsGA30 <- glmTreat(fit, contrast=my.contrasts[,"GA20vsGA30"], lfc=log2(1.5)) #looks at genes above a certain log fold change threshold
tt_GA20vsGA30 <- topTags(tr_GA20vsGA30, n=Inf)
#topids_GA20vsGA30 <- tt_GA20vsGA30$table[tt_GA20vsGA30$table$FDR<thresh_low]
#topids_GA20vsGA30.df <- data.frame(topids_GA20vsGA30)
tt_GA20vsGA30.df <- data.frame(tt_GA20vsGA30)
#tr_GA20vsGA30.df <- data.frame(tr_GA20vsGA30)
GA20vsGA30_logFDR <- -log10(tt_GA20vsGA30.df$FDR)
GA20vsGA30_logFDR.df <- data.frame(GA20vsGA30_logFDR)
tt_GA20vsGA30.df <- cbind(tt_GA20vsGA30.df,GA20vsGA30_logFDR.df)
write.table(topTags(tr_GA20vsGA30,n=30),
            file = "~/Documents/Research/Geukensia/GA20vsGA30_tt_summary.txt",
            quote=FALSE)
summary(decideTestsDGE(tr_GA20vsGA30, adjust.method="fdr",p.value=0.05))
GA20vsGA30plot <- ggplot(tt_GA20vsGA30.df) +aes(x = logFC,y = GA20vsGA30_logFDR,colour = GA20vsGA30_logFDR) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_gradient(low = "#000000", high = "#000000") +
  labs(x = "log2 Fold Change",y = "-log10(FDR)") +
  theme_bw() + geom_hline(yintercept=1.3) +
  theme_minimal(base_size = 22) + theme(legend.position = "none")

#find genes responding differently between MA temps

tr_MA20vsMA30 <- glmTreat(fit, contrast=my.contrasts[,"MA20vsMA30"], lfc=log2(1.5)) #looks at genes above a certain log fold change threshold
tt_MA20vsMA30 <- topTags(tr_MA20vsMA30, n=Inf)
#topids_MA20vsMA30 <- tt_MA20vsMA30$table[tt_MA20vsMA30$table$FDR<thresh_low]
#topids_MA20vsMA30.df <- data.frame(topids_MA20vsMA30)
tt_MA20vsMA30.df <- data.frame(tt_MA20vsMA30)
#tr_MA20vsMA30.df <- data.frame(tr_MA20vsMA30)
MA20vsMA30_logFDR <- -log10(tt_MA20vsMA30.df$FDR)
MA20vsMA30_logFDR.df <- data.frame(MA20vsMA30_logFDR)
tt_MA20vsMA30.df <- cbind(tt_MA20vsMA30.df,MA20vsMA30_logFDR.df)
write.table(topTags(tr_MA20vsMA30,n=30),
            file = "~/Documents/Research/Geukensia/MA20vsMA30_tt_summary.txt",
            quote=FALSE)
summary(decideTestsDGE(tr_MA20vsMA30, adjust.method="fdr",p.value=0.05))
MA20vsMA30plot <- ggplot(tt_MA20vsMA30.df) +aes(x = logFC,y = MA20vsMA30_logFDR,colour = MA20vsMA30_logFDR) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_gradient(low = "#000000", high = "#000000") +
  labs(x = "log2 Fold Change",y = "-log10(FDR)") +
  theme_bw() + geom_hline(yintercept=1.3) +
  theme_minimal(base_size = 22) + theme(legend.position = "none")

temp_plots <- plot_grid(GA20vsGA30plot,MA20vsMA30plot, labels="AUTO")
save_plot("Temp_Differences_withMA13-30_all.png",temp_plots,device="png",ncol=2)

