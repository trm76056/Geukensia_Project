setwd("~/Documents/Research/Geukensia/")

library(PopGenome)
library(tidyverse)
library(ggplot2)
library(coin)
library(cowplot)

#read in the fasta alignments
GENOME.class <- readData("OrthoFinder_mafft")

#number of sites analyzed in each alignment
GENOME.class@n.sites

#get summary information from the alignments
get.sum.data(GENOME.class)

#set the populations
GENOME.class <- set.populations(GENOME.class,list(
  c("GA1-30","GA10-20","GA3-30","GA4-30","GA5-30","GA6-20","GA7-20","GA8-20","GA9-20","GAC1","GAC2","GAC3","GAC4","GAC5"),
  c("MA12-30","MA13-30","MA15-30","MA16-20","MA17-20","MA18-20","MA20-20","MAC1","MAC2","MAC3","MAC5")))

#neutrality statistics
GENOME.class <- neutrality.stats(GENOME.class)
get.neutrality(GENOME.class)
get.neutrality(GENOME.class)[[1]]
tajd <- GENOME.class@Tajima.D
tajD.df <- data.frame(tajd)
tajD_GAplot <- ggplot(tajD.df) +aes(x = pop.1) +geom_histogram(bins = 30L, fill = "#000000") +
  labs(x = "Tajima's D",y = "Number")+
  theme_minimal(base_size = 22)
tajD_MAplot <- ggplot(tajD.df) +aes(x = pop.2) +geom_histogram(bins = 30L, fill = "#000000") +
  labs(x = "Tajima's D",y = "Number") +
  theme_minimal(base_size = 22)

#fu and li's D instead of Tajima's D
#fuliD <- GENOME.class@Fu.Li.D
#fuliD.df <- data.frame(fuliD)
#fuliD_GAplot <- ggplot(fuliD.df) +aes(x = pop.1) +geom_histogram(bins = 30L, fill = "#000000") +
#  labs(x = "Fu and Li's D*",y = "Number")+
#  theme_minimal(base_size = 22)
#fuliD_MAplot <- ggplot(fuliD.df) +aes(x = pop.2) +geom_histogram(bins = 30L, fill = "#000000") +
#  labs(x = "Fu and Li's D*",y = "Number") +
#  theme_minimal(base_size = 22)

#calculate pi dxy and fst
GENOME.class <- diversity.stats(GENOME.class, pi=TRUE) #diversity within populations
#plot(GENOME.class@Pi)
pi_persite <- GENOME.class@Pi / GENOME.class@n.valid.sites
pi_persite.df <- data.frame(pi_persite)
pi_persite.df.2 <- data.frame(pi_persite = c(pi_persite.df$pop.1, pi_persite.df$pop.2), 
                              population = factor(rep(c("Georgia","Massachusetts"),c(length(pi_persite.df$pop.1),length(pi_persite.df$pop.2)))))

mean(pi_persite.df$pop.1, na.rm = TRUE)
max(pi_persite.df$pop.1, na.rm = TRUE)
mean(pi_persite.df$pop.2, na.rm = TRUE)
max(pi_persite.df$pop.2, na.rm = TRUE)
GENOME.class <- F_ST.stats(GENOME.class, mode="nucleotide")
#plot(GENOME.class@nucleotide.F_ST)
pi <- GENOME.class@Pi
pi.df <- data.frame(pi)
fst <- t(GENOME.class@nuc.F_ST.pairwise)
fst[fst<0]<-0
GENOME.class <- F_ST.stats.2(GENOME.class,snn=TRUE,Phi_ST=TRUE, new.populations=GENOME.class@populations)
snn <- GENOME.class@Hudson.Snn
snn[snn<0]<-0
snn.df <- data.frame(snn)
mean(snn.df$Snn, na.rm = TRUE)
sapply(snn.df, function(x) head(row.names(snn.df)[order(x, decreasing = TRUE)], 10)) #get top two orthologs for snn
snnplot <- ggplot(as.data.frame(snn)) + aes(x = `Snn`) + geom_histogram(bins = 30L, fill = "#000000") +
  labs(x = expression(paste(italic("S")[nn])),y = "Number") +
  theme_minimal(base_size = 22)
phi_st <- GENOME.class@Phi_ST
phi_st[phi_st<0]<-0
phi_st.df <- data.frame(phi_st)
mean(phi_st.df$Phi_ST, na.rm = TRUE)
max(phi_st.df$Phi_ST, na.rm = TRUE)
phi_stPlot <- ggplot(as.data.frame(phi_st)) + aes(x = `Phi_ST`) + geom_histogram(bins = 30L, fill = "#000000") +
  labs(x = expression(paste(Phi)[ST]),y="Number") +
  theme_minimal(base_size = 22)
pi_persiteGA <- ggplot(pi_persite.df) +aes(x = pop.1) +geom_histogram(bins = 30L, fill = "#000000") +
  labs(x = expression(paste(Pi," per site")),y = "Number") +
  theme_minimal() + theme(axis.text=element_text(size=10),axis.title=element_text(size=22)) + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
pi_persiteMA <- ggplot(pi_persite.df) +aes(x = pop.2) +geom_histogram(bins = 30L, fill = "#000000") +
  labs(x = expression(paste(Pi," per site")),y = "Number") +
  theme_minimal() + theme(axis.text=element_text(size=10),axis.title=element_text(size=22)) + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))

sapply(pi_persite.df, function(x) head(row.names(pi_persite.df)[order(x, decreasing = TRUE)], 10)) #get top two orthologs for pi
sapply(fuliD.df, function(x) head(row.names(fuliD.df)[order(x, decreasing = TRUE)], 10))
sapply(fuliD.df, function(x) head(row.names(fuliD.df)[order(x, decreasing = FALSE)], 10))
sapply(phi_st.df, function(x) head(row.names(phi_st.df)[order(x, decreasing = TRUE)], 10))

#from john to get net nucleotide divergence per base pair
GENOME.class <- diversity.stats.between(GENOME.class)
dxy<-GENOME.class@nuc.diversity.between
dxy<-dxy/GENOME.class@n.sites
dxdy<-GENOME.class@nuc.diversity.within
dxdy<-dxdy/GENOME.class@n.sites
dxdy2<-((dxdy[,1])+(dxdy[,2]))/2
d_int<-((dxdy[,1])+(dxdy[,2]))
net<-dxy-dxdy2
net[net<0]<-0
#net<-t(net)
netdf <- data.frame(net)
net_pop.df <- data.frame(net)
mean(net_pop.df$pop1.pop2)
max(net_pop.df$pop1.pop2)
netPlot <- ggplot(net_pop.df) + aes(x = pop1.pop2) +geom_histogram(bins = 30L, fill = "#000000") +
  labs(x = expression(paste(italic("d")[A])),y = "Number") +
  theme_minimal(base_size = 22) + geom_vline(xintercept = mean(netdf$pop1.pop2)) + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
sapply(netdf, function(x) head(row.names(netdf)[order(x, decreasing = TRUE)], 10)) #get top two orthologs for dxy

#create figure for within population diversity stats
within_plots <- plot_grid(pi_persiteGA,pi_persiteMA,tajD_GAplot,tajD_MAplot,labels="AUTO")
save_plot("WithinPopDiversityPlots_withMA13-30.png",within_plots,ncol=2,nrow=2)

#create figure for between population diversity stats
between_plots <- plot_grid(netPlot,snnplot,labels="AUTO")
save_plot("BetweenPopDiversityPlots_withMA13-30.png",between_plots,ncol=2,nrow=1)

#create data frame of all values
all.df <- data.frame(snn,net,pi_persite,tajd,fst)
colnames(all.df) <- c("Snn","Da","Pi_GA","Pi_MA","TajimaD_GA","TajimaD_MA","FST")

library(dplyr)

#from John
randiv<-readData("~/Documents/Research/Geukensia/OrthoFinder_mafft/")
#ranSnnout = NULL
#testSnnout = NULL

all_samples <- c("GA1-30","GA10-20","GA3-30","GA4-30","GA5-30","GA6-20","GA7-20","GA8-20","GA9-20","GAC1","GAC2","GAC3","GAC4","GAC5",
                 "MA12-30","MA15-30","MA16-20","MA17-20","MA18-20","MA20-20","MAC1","MAC2","MAC3","MAC5")
bigJ <- 1000
ga_length <- 14
ma_length <- 10


#permutations for snn
for (i in 1:bigJ) {
  
  testrand1<-c(as.character(unlist(sample(all_samples,(ga_length),replace=F))))
  #  testrand1<-c(as.character(unlist(sample(randiv$Tip,(40-i),replace=F))))
  #rand1 grabs sizeN sequences to set populations at random, keeping size correct. now how to do more than 2 regions?
  #randRem should be whatever is not picked from data into rand1, otherwise biases the values down, cannot sample w replace, repeat for all locations
  
  randRem2<-c(as.character(unlist(setdiff(all_samples,testrand1))))
  #this is putting all of the remainder into 'randRem' and works
  
  randiv<-set.populations(randiv,list(testrand1,randRem2))
  randiv@region.data@populations2
  randivSnn<-F_ST.stats.2(randiv,snn=TRUE)
  #testSnnout[i]<-rantestSnn@Hudson.Snn
  randivSnn.df <- data.frame(randivSnn@Hudson.Snn)
  snn.df <- snn.df %>% add_column(i = randivSnn.df$Snn)
}

snn.df <- snn.df %>% add_column(p_value = rowSums(snn.df[,2:ncol(snn.df)]>snn.df[,1])/bigJ)
rownames(snn.df)[which(snn.df$p_value <0.05)] #get transcripts that had a p-value less than 0.05

#permutations for dA
for (i in 1:bigJ) {
  testrand1<-c(as.character(unlist(sample(all_samples,(ga_length),replace=F))))
  #  testrand1<-c(as.character(unlist(sample(randiv$Tip,(40-i),replace=F))))
  #rand1 grabs sizeN sequences to set populations at random, keeping size correct. now how to do more than 2 regions?
  #randRem should be whatever is not picked from data into rand1, otherwise biases the values down, cannot sample w replace, repeat for all locations
  
  randRem2<-c(as.character(unlist(setdiff(all_samples,testrand1))))
  #this is putting all of the remainder into 'randRem' and works
  
  #calculate dA for the permutation
  randiv<-set.populations(randiv,list(testrand1,randRem2))
  randiv <- neutrality.stats(randiv)
  randiv <- diversity.stats(randiv)
  randiv <- diversity.stats.between(randiv)
  randxy<-randiv@nuc.diversity.between
  randxy<-randxy/randiv@n.sites
  randxdy<-randiv@nuc.diversity.within
  randxdy<-randxdy/randiv@n.sites
  randxdy2<-((randxdy[,1])+(randxdy[,2]))/2
  rand_int<-((randxdy[,1])+(randxdy[,2]))
  rannet<-randxy-randxdy2
  rannet[rannet<0]<-0
  rannet.df <- data.frame(rannet)
  netdf <- netdf %>% add_column(i = rannet.df$pop1.pop2)
}

netdf <- netdf %>% add_column(p_value = rowSums(netdf[,2:ncol(netdf)]>netdf[,1])/bigJ)
rownames(netdf)[which(netdf$p_value <0.05)] #get transcripts that had a p-value less than 0.05

#permutations for phi st
for (i in 1:bigJ) {
  
  testrand1<-c(as.character(unlist(sample(all_samples,(ga_length),replace=F))))
  #  testrand1<-c(as.character(unlist(sample(randiv$Tip,(40-i),replace=F))))
  #rand1 grabs sizeN sequences to set populations at random, keeping size correct. now how to do more than 2 regions?
  #randRem should be whatever is not picked from data into rand1, otherwise biases the values down, cannot sample w replace, repeat for all locations
  
  randRem2<-c(as.character(unlist(setdiff(all_samples,testrand1))))
  #this is putting all of the remainder into 'randRem' and works
  
  randiv<-set.populations(randiv,list(testrand1,randRem2))
  randiv <- F_ST.stats.2(randiv,Phi_ST=TRUE)
  #rantestdiv@region.data@populations2
  randivphi_st <- randiv@Phi_ST
  randivphi_st[randivphi_st<0]<-0
  randivphi_st.df <- data.frame(randivphi_st)
  phi_st.df <- phi_st.df %>% add_column(i = randivphi_st.df$Phi_ST)
}

phi_st.df <- phi_st.df %>% add_column(p_value = rowSums(phi_st.df[,2:ncol(phi_st.df)]>phi_st.df[,1])/bigJ)
rownames(phi_st.df)[which(phi_st.df$p_value <0.05)] #get transcripts that had a p-value less than 0.05


#permutations for pi per site in GA and MA
pi_persite_ga.df <- data.frame(pi_persite.df$pop.1)
pi_persite_ma.df <- data.frame(pi_persite.df$pop.2)
for (i in 1:bigJ) {
  
  testrand1<-c(as.character(unlist(sample(all_samples,(ga_length),replace=F))))
  #  testrand1<-c(as.character(unlist(sample(randiv$Tip,(40-i),replace=F))))
  #rand1 grabs sizeN sequences to set populations at random, keeping size correct. now how to do more than 2 regions?
  #randRem should be whatever is not picked from data into rand1, otherwise biases the values down, cannot sample w replace, repeat for all locations
  
  randRem2<-c(as.character(unlist(setdiff(all_samples,testrand1))))
  #this is putting all of the remainder into 'randRem' and works
  
  
  randiv<-set.populations(randiv,list(testrand1,randRem2))
  randiv@region.data@populations2
  randiv <- diversity.stats(randiv, pi=TRUE) #diversity within populations
  rantestpi_persite <- randiv@Pi / randiv@n.valid.sites
  rantestpi_persite.df <- data.frame(rantestpi_persite)
  pi_persite_ga.df <- pi_persite_ga.df %>% add_column(i = rantestpi_persite.df$pop.1)
  pi_persite_ma.df <- pi_persite_ma.df %>% add_column(i = rantestpi_persite.df$pop.2)
}

pi_persite_ga.df <- pi_persite_ga.df %>% add_column(p_value = rowSums(pi_persite_ga.df[,3:ncol(pi_persite_ga.df)]>pi_persite_ga.df[,1])/bigJ)
rownames(pi_persite_ga.df)[which(pi_persite_ga.df$p_value <0.05)] #get transcripts that had a p-value less than 0.05
pi_persite_ma.df <- pi_persite_ma.df %>% add_column(p_value = rowSums(pi_persite_ma.df[,3:ncol(pi_persite_ma.df)]>pi_persite_ma.df[,2])/bigJ)
rownames(pi_persite_ma.df)[which(pi_persite_ma.df$p_value <0.05)] #get transcripts that had a p-value less than 0.05


finalall.df <- cbind(all.df,snn.df$p_value)


#Calculate ka/ks / dn/ds
library(ape)
library(adegenet)
try <- fasta2DNAbin("OrthoFinder_mafft/new_OG0017815.fa")
dnds(try)

