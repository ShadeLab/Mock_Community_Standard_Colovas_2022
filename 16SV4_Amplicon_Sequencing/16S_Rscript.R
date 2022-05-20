#################################### Mock_Community_Standar_16S #####################################
##
# Date: May 19th 2022
# By : Colovas et al

library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(viridis)
library(microbiome)
library(reshape2)
library(ggpubr)
library(broom)
library(ggfortify)
library(dplyr)
library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(metagenomeSeq)
library(devtools)
library(metagMisc)

### Set working directory
#setwd("/Volumes/ShadeLab/WorkingSpace/MarcoMechan_WorkingSpace/Mucilage_TX08001_2020/Controls/phyloseq")

#### Panel A: Expected Mock community
Controls.metadata=read.table("Mock_comm_metadata.txt", sep="\t", header = T)

Control.expected.ord <- transform(Controls.metadata, 
                                  Control.expected.ord  = factor(
                                    Genus,
                                    levels=c("Bacillus","Escherichia-Shigella","Salmonella","Pseudomonas","Ochrobactrum","Rhizobium","Staphylococcus","Streptomyces"),
                                    ordered =TRUE))
col_vector <- c("Bacillus"="#ffeda0", "Escherichia-Shigella"="#feb24c","Salmonella"="#fee6ce","Pseudomonas"="#a1d99b","Ochrobactrum"="#756bb1","Rhizobium"="#bcbddc","Staphylococcus"="#fa9fb5","Streptomyces"="#7095c5")

Control.expected.plot <- ggplot(Control.expected.ord, aes(y=Abundance, x=Source,order=Control.expected.ord)) + 
  geom_bar(stat="identity", aes(fill=Control.expected.ord)) + 
  scale_fill_manual(values =col_vector) + 
  guides(fill=guide_legend(keywidth = 1,keyheight = 1)) + 
  labs(fill="Bacterial Genus")+
  ylab("Relative Abundance (Genus) \n") +
  xlab("Expected Mock Community") +
  theme(axis.text.x=element_text(size=0,angle=0,hjust =0.5,colour = "black", face = "bold"),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=14,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=14, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=10),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#ggtitle("Expected Mock Community diversity")
Control.expected.plot

#### Panel B: Normalized Mock community
#### Import Qiime2 files

otu=read.table("otu_table.txt", header = T, row.names = 1)
tax=read.delim("taxonomy.tsv",row.names = 1)
map=read.table("metadata.csv", header = T, row.names = 1, sep=",")

tax=tax[-2]
tax_df <- colsplit(tax$Taxon, '; ', names =  c("Kingdom", "Phylum", "Class", 
                                               "Order", "Family", "Genus", "Species"))
tax_df[1:7] <- lapply(tax_df[1:7], function(x) gsub(".*__", "", x))
rownames(tax_df) <- rownames(tax)

OTU=otu_table(as.matrix(otu), taxa_are_rows = T)
TAX=tax_table(as.matrix(tax_df))
MAP=sample_data(map)

otuPhyloseq=phyloseq(OTU,TAX,MAP)
sample_sums(otuPhyloseq)

#Filtering mitochondria, Chloroplast and Unclassified taxa
otuPhyloseq_filt <- otuPhyloseq %>%
  subset_taxa(Kingdom != "Unassigned")
otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Family != "Mitochondria") | is.na(Family))
otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Order != "Chloroplast") | is.na(Class))
sample_sums(otuPhyloseq_filt) 
filtered_otus <- phyloseq_to_df(otuPhyloseq_filt, addtax=T, addtot=T)

####### Plotting all Control samples ################
###### OTU table was generated with 4 samples: ShadeLab Mock Community, Zymo Mock Community, Shadelab Negative Control, RTSF(Core) Negative Control
cont.otu <- prune_taxa(taxa_sums(otuPhyloseq_filt) > 1, otuPhyloseq_filt)
print(cont.otu)
sample_sums(cont.otu)
#### ShadeLab Mock: 147356 reads, Zymo Mock:207140 reads, Shadelab Neg.Control: 60 reads, RTSF Neg.Control:69 reads

all_contr_genus <- cont.otu %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.001) %>%
  arrange(Genus)

n <- dim(all_contr_genus )[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vect = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(all_contr_genus ,aes(x=Sample,y=Abundance,fill=Genus)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vect) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus) \n") +
  xlab("Positive controls") +
  scale_x_discrete(expand = c(.3, 0),
                   limit = c("Mock.ShadeLab","Mock.Zymo", "Muc.NC1", "RTSF.NC"),
                   labels = c("ShadeLab Mock","Zymo Mock", "SL NC", "RTSF NC"))+
  theme(axis.text.x=element_text(size=14,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Genus Composition")

#### Removing Negative controls
Control.samples <- prune_samples(sample_sums(cont.otu) >= 140000, cont.otu)
Control.samples <- prune_taxa(taxa_sums(Control.samples) > 0, Control.samples)
print(Control.samples)
sample_sums(Control.samples)

#### Rarefying Positive controls 
set.seed(13)
Control.samples.rare = rarefy_even_depth(Control.samples, rngseed=1, sample.size=min(sample_sums(Control.samples)), replace=F)
print (Control.samples.rare) ### Samples were rarefied to 147356 reads
sample_sums(Control.samples.rare)

### Subsetting ShadeLab Mock Community
Shadelab.Comm <- subset_samples(Control.samples.rare, Sample.name%in%c("Mock.ShadeLab"))
Shadelab.Comm <- prune_taxa(taxa_sums(Shadelab.Comm) > 0, Shadelab.Comm)

Shadelab.Community.genus <- Shadelab.Comm %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  arrange(Genus)

Shadelab.Community.genus.ord <- transform(Shadelab.Community.genus, 
                                          Shadelab.Community.genus.ord  = factor(
                                            Genus,
                                            levels=c("Bacillus","Escherichia-Shigella","Salmonella","Pseudomonas","Ochrobactrum","Staphylococcus","Streptomyces"),
                                            ordered =TRUE))
col.vector <- c("Bacillus"="#ffeda0", "Escherichia-Shigella"="#feb24c","Salmonella"="#fee6ce","Pseudomonas"="#a1d99b","Ochrobactrum"="#756bb1","Staphylococcus"="#fa9fb5","Streptomyces"="#7095c5")

Shadelab.Community.genus.plot <- ggplot(Shadelab.Community.genus.ord, aes(y=Abundance, x=Sample,order=Shadelab.Community.genus.ord)) + 
  geom_bar(stat="identity", aes(fill=Shadelab.Community.genus.ord)) + 
  scale_fill_manual(values =col.vector) + 
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  labs(fill="Bacterial Genus")+
  ylab("Relative Abundance (Genus) \n") +
  xlab("Rarefied Mock Community") +
  theme(axis.text.x=element_text(size=0,angle=0,hjust =0.5, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 0),
        strip.text.x = element_text(size=14,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=14, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title.y =element_text(size=0,face="bold", vjust = 10),
        axis.title.x =element_text(size=14,face="bold", vjust = 0),
        legend.text = element_text(size=10),
        legend.title = element_text(size =12, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#ggtitle("Genus Composition")
Shadelab.Community.genus.plot

#### Panel C: Correcting Abundance with 16S copy number
Shadelab.Community.otu <- Shadelab.Comm %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt() %>%
  filter(Abundance > 22) %>%
  arrange(Genus)

Merged.Df <- merge (x=Shadelab.Community.otu,y=Controls.metadata[,c(5:7)], by='Genus')

Merged.Df2 <- Merged.Df %>%
  mutate(Abundance = (Abundance/Copy.16S),
         Relativ.Abundance = sum(as.numeric(Abundance)),
         Relativ.Abundance = ((Abundance*100)/Relativ.Abundance)/100)


Norm.Comm.genus.ord <- transform(Merged.Df2,
                                 Norm.Comm.genus.ord  = factor(
                                   Genus,
                                   levels=c("Bacillus","Escherichia-Shigella","Salmonella","Pseudomonas","Ochrobactrum","Staphylococcus","Streptomyces"),
                                   ordered =TRUE))
col.vector <- c("Bacillus"="#ffeda0", "Escherichia-Shigella"="#feb24c","Salmonella"="#fee6ce","Pseudomonas"="#a1d99b","Ochrobactrum"="#756bb1","Staphylococcus"="#fa9fb5","Streptomyces"="#7095c5")

Correct.Comm.plot <- ggplot(Norm.Comm.genus.ord, aes(y=Relativ.Abundance, x=Sample,order=Norm.Comm.genus.ord)) + 
  geom_bar(stat="identity", aes(fill=Norm.Comm.genus.ord)) + 
  scale_fill_manual(values =col.vector) + 
  guides(fill=guide_legend(reverse=F,keywidth = 1,keyheight = 1)) + 
  labs(fill="Bacterial Genus")+
  ylab("Relative Abundance (Genus) \n") +
  xlab("Corrected Mock Community") +
  theme(axis.text.x=element_text(size=0,angle=0,hjust =0.5),
        axis.text.y = element_text(size = 0),
        strip.text.x = element_text(size=14,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=0, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title.y =element_text(size=0,face="bold", vjust = 10),
        axis.title.x =element_text(size=14,face="bold", vjust = 0),
        legend.text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#ggtitle("Genus Composition")
Correct.Comm.plot

### Merging panels
legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

Shadelab_legend = legend(Control.expected.plot)
Figure.Mock.Comm <- cowplot::plot_grid(Control.expected.plot + theme(legend.position = "none"),
                                       Shadelab.Community.genus.plot + theme(legend.position = "none"), 
                                       Correct.Comm.plot + theme(legend.position = "none"),
                                       ncol=3, align = "v",
                                       axis = "lr") 
SL.legend = cowplot::plot_grid(Shadelab_legend, ncol=1)
Figure.Mock.Comm = cowplot::plot_grid(Figure.Mock.Comm, SL.legend, nrow=1, ncol=2, rel_widths = c(0.8,0.2))

Figure.Mock.Comm
                      
