# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# March 20, 2023
# Rscript: I4512ivac_alpha.R
# R v4.2.2/R v2022.12.0

# PART 3: Using the phyloseq R package to calculate alpha diversity
## anything ending with f is for Fecal samples
## anything ending with f2 is for Fecal samples without MVMC treatment
## anything ending with cc is for Cecal Content samples


# ALPHA DIVERSITY: RICHNESS & EVENNESS___________________________________________________
library(phyloseq) #v1.42.0
library(tidyverse) #v1.3.2
library(vegan) #v2.6.4
library(ggplot2) #3.4.0
library(dplyr)
library(readr)
library(reshape2)

set.seed(19760620)

## Open I4512ivac_phylo_all.rds and create phylo subset objects
PHYLO <- read_rds("./outputs/I4512ivac_phylo_all.rds")
PHYLOf <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Fecal", x = PHYLO)
PHYLOf2 <- prune_samples(samples = PHYLOf@sam_data$Challenge=="I4512i", x = PHYLOf)
PHYLOcc <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Cecal_contents", x = PHYLO)

#Create Taxonomy R matrix
TAXf <-
  PHYLOf@tax_table %>%
  as(., Class = 'matrix') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV')
TAXcc <-
  PHYLOcc@tax_table %>%
  as(., Class = 'matrix') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV')

#Subset and prune_taxa
EFGf <-
  PHYLOf %>%
  subset_samples(PHYLOf@sam_data$Vaccine %in% c('Mock', 'Enterisol') &
                   PHYLOf@sam_data$Challenge %in% c('Mock', 'I4512i')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %>%
  prune_taxa(taxa=colSums(.@otu_table > 0) > 5) # taxa detected in more than 5 samples
EFGcc <-
  PHYLOcc %>%
  subset_samples(PHYLOcc@sam_data$Vaccine %in% c('Mock', 'Enterisol') &
                   PHYLOcc@sam_data$Challenge %in% c('Mock', 'I4512i')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %>%
  prune_taxa(taxa=colSums(.@otu_table > 0) > 5) # taxa detected in more than 5 samples

#Create dpc, set, and Group columns
EFGf@sam_data$dpc <- EFGf@sam_data$Days_post_challenge
EFGf@sam_data$set <- paste(EFGf@sam_data$Days_post_challenge, EFGf@sam_data$Vaccine,EFGf@sam_data$Challenge, sep = '_')
EFGf@sam_data$Group <- paste(EFGf@sam_data$Vaccine,EFGf@sam_data$Challenge, sep = '_')

EFGcc@sam_data$dpc <- EFGcc@sam_data$Days_post_challenge
EFGcc@sam_data$set <- paste(EFGcc@sam_data$Days_post_challenge, EFGcc@sam_data$Vaccine,EFGcc@sam_data$Challenge, sep = '_')
EFGcc@sam_data$Group <- paste(EFGcc@sam_data$Vaccine,EFGcc@sam_data$Challenge, sep = '_')

# Plot_richness (phyloseq) - total # OTUs/taxa in sample
## See how spread out your data actually is overall (all points not constrained by 25%/75% percentile)
plot_richness(PHYLOf, x="Treatment", color="Treatment", measures=c("Shannon", "Simpson", "Chao1"))
plot_richness(PHYLOcc, x="Treatment", color="Treatment", measures=c("Shannon", "Simpson", "Chao1"))
    ### Shannon
tiff(filename="./graphics/I4512ivac_fecal_plotrichness_shannon_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
richshanf <- plot_richness(PHYLOf, x="Treatment", color="Treatment", measures=c("Shannon")) + 
  geom_point(aes(color=Treatment)) +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
richf + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
tiff(filename="./graphics/I4512ivac_cecalcontents_plotrichness_shannon_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
richshancc <- plot_richness(PHYLOcc, x="Treatment", color="Treatment", measures=c("Shannon")) + 
  geom_point(aes(color=Treatment)) +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
richcc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
    ### Simpson
tiff(filename="./graphics/I4512ivac_fecal_plotrichness_simpson_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
richsimpf <- plot_richness(PHYLOf, x="Treatment", color="Treatment", measures=c("Simpson")) + 
  geom_point(aes(color=Treatment)) +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
richsimpf + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
tiff(filename="./graphics/I4512ivac_cecalcontents_plotrichness_simpson_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
richsimpcc <- plot_richness(PHYLOcc, x="Treatment", color="Treatment", measures=c("Simpson")) + 
  geom_point(aes(color=Treatment)) +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
richsimpcc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
    ### Chao1
tiff(filename="./graphics/I4512ivac_fecal_plotrichness_chao1_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
richchaof <- plot_richness(PHYLOf, x="Treatment", color="Treatment", measures=c("Chao1")) + 
  geom_point(aes(color=Treatment)) +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
richchaof + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
tiff(filename="./graphics/I4512ivac_cecalcontents_plotrichness_chao1_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
richchaocc <- plot_richness(PHYLOcc, x="Treatment", color="Treatment", measures=c("Chao1")) + 
  geom_point(aes(color=Treatment)) +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
richchaocc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

##Rarefy_even_depth
EFGf_rare <- rarefy_even_depth(EFGf)
EFGcc_rare <- rarefy_even_depth(EFGcc)

# Creating boxplots for each alpha diversity measure
    ### Shannon
EFGf_rare@sam_data$Shannon <- diversity(EFGf_rare@otu_table)
tiff(filename="./graphics/I4512ivac_fecal_shannon_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
shanf <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Group, y=Shannon, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
shanf + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
tiff(filename="./graphics/I4512ivac_fecal_shannon_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
shan1f <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Shannon, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
shan1f + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

EFGcc_rare@sam_data$Shannon <- diversity(EFGcc_rare@otu_table)
tiff(filename="./graphics/I4512ivac_cecalcontents_shannon_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
shancc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Group, y=Shannon, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
shancc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
tiff(filename="./graphics/I4512ivac_cecalcontents_shannon_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
shan1cc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Shannon, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
shan1cc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

    ### Simpson
EFGf_rare@sam_data$Simpson <- diversity(EFGf_rare@otu_table)
tiff(filename="./graphics/I4512ivac_fecal_simpson_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
simpf <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Group, y=Simpson, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
simpf + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
tiff(filename="./graphics/I4512ivac_fecal_simpson_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
simp1f <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Simpson, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
simp1f + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

EFGcc_rare@sam_data$Simpson <- diversity(EFGcc_rare@otu_table)
tiff(filename="./graphics/I4512ivac_cecalcontents_simpson_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
simpcc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Group, y=Simpson, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
simpcc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
tiff(filename="./graphics/I4512ivac_cecalcontents_simpson_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
simp1cc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Simpson, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
simp1cc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

    ### Chao1
EFGf_rare@sam_data$Chao1 <- diversity(EFGf_rare@otu_table)
tiff(filename="./graphics/I4512ivac_fecal_chao1_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
chaof <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Group, y=Chao1, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
chaof + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
tiff(filename="./graphics/I4512ivac_fecal_chao1_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
chao1f <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Chao1, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
chao1f + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

EFGcc_rare@sam_data$Chao1 <- diversity(EFGcc_rare@otu_table)
tiff(filename="./graphics/I4512ivac_cecalcontents_chao1_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
chaocc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Group, y=Chao1, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
chaocc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
tiff(filename="./graphics/I4512ivac_cecalcontents_chao1_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
chao1cc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Chao1, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
chao1cc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()









# Linear model
mod_datf <-
  EFGf_rare@sam_data %>%
  data.frame() %>%
  mutate(Group=factor(Group, levels=c('Enterisol_I4512i', 'Mock_I4512i', 'Mock_Mock')),
         Days_post_challenge = factor(Days_post_challenge, levels = c('-28', '-1', '2', '3', '7','10', '14')))
colnames(EFGf_rare@sam_data)

shannon_modf <- lm(data = mod_datf, formula= Shannon ~ Days_post_challenge * Group)
shannon_modf %>% summary()

mod_datcc <-
  EFGcc_rare@sam_data %>%
  data.frame() %>%
  mutate(Group=factor(Group, levels=c('Enterisol_I4512i', 'Mock_I4512i', 'Mock_Mock')),
         Days_post_challenge = factor(Days_post_challenge, levels = c('14')))
colnames(EFGcc_rare@sam_data)

shannon_modcc <- lm(data = mod_datcc, formula= Shannon ~ Group)
shannon_modcc %>% summary()
simpson_modf <- lm(data = mod_datf, formula= Simpson ~ Group)


# EMMEANS - estimated marginal means (aka least squares means) for summarizing a fitted linear model that include factors
#Fecal
library(emmeans)
EMMEANSf <- emmeans(shannon_modf, ~ Group | Days_post_challenge)
shan_contrastsf <- emmeans::contrast(EMMEANSf, method='pairwise', adjust='fdr')
contrasts_tablef <- shan_contrastsf %>% as.data.frame()
plot_datf <- EMMEANSf %>% as.data.frame()

tiff(filename="./graphics/I4512ivac_fecal_shannon_alpha_div_over_time.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
plot_datf %>%
  ggplot(aes(x=Days_post_challenge,
             y=emmean,ymin=lower.CL,
             ymax=upper.CL,
             group=Group, color=Group)) +
  geom_pointrange(position = position_dodge(width = .2)) +
  geom_line() + ggtitle('Alpha diversity over time') +
  ylab('Shannon') +
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_I4512i = 'skyblue',
                              Enterisol_I4512i = 'orange')) +
  theme_bw()
dev.off()

#Cecal_contents
EMMEANScc <- emmeans(shannon_modcc, ~ Group)
shan_contrastscc <- emmeans::contrast(EMMEANScc, method='pairwise', adjust='fdr')
contrasts_tablecc <- shan_contrastscc %>% as.data.frame()
plot_datcc <- EMMEANScc %>% as.data.frame()





# Alpha contrasts table
#Fecal
alpha_contrasts_tablef <-
  contrasts_tablef %>%
  mutate(across(where(is.numeric), signif, digits=2),
         sig=ifelse(p.value < 0.05, TRUE, FALSE))
write_tsv(alpha_contrasts_tablef, "./outputs/I4512ivac_fecal_alpha_tests.tsv")

#Cecal_contents
alpha_contrasts_tablecc <-
  contrasts_tablecc %>%
  mutate(across(where(is.numeric), signif, digits=2),
         sig=ifelse(p.value < 0.05, TRUE, FALSE))
write_tsv(alpha_contrasts_tablecc, "./outputs/I4512ivac_cecalcontents_alpha_tests.tsv")


# Beta dispersion
BETA_DISPERf <- betadisper(DISTf, group = sample_data(EFGf)$set)
EFGf_rare@sam_data$Distance_from_centroid <- BETA_DISPERf$distances
EFGf_rare@sam_data %>% ggplot(aes(x=Group, y=Distance_from_centroid)) + geom_boxplot() +
  facet_wrap(~Days_post_challenge, ncol=1)

BETA_DISPERcc <- betadisper(DISTcc, group = sample_data(EFGcc)$set)
EFGcc_rare@sam_data$Distance_from_centroid <- BETA_DISPERcc$distances
EFGcc_rare@sam_data %>% ggplot(aes(x=Group, y=Distance_from_centroid)) + geom_boxplot() +
  facet_wrap(~Days_post_challenge, ncol=1)


EMMEANSf <- emmeans(shannon_modf, ~ Group | Days_post_challenge)
shan_contrastsf <- emmeans::contrast(EMMEANSf, method='pairwise', adjust='fdr')
contrasts_tablef <- shan_contrastsf %>% as.data.frame()
plot_datf <- EMMEANSf %>% as.data.frame()



tiff(filename="./graphics/I4512ivac_fecal_alpha_div_over_time.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
EFGf_rare@sam_data %>%
  ggplot(aes(x=Days_post_challenge,
             y=Distance_from_centroid,
             group=Group, color=Group)) +
  geom_line() + ggtitle('Beta dispersion over time') +
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_I4512i = 'skyblue',
                              Enterisol_I4512i = 'orange')) +
  theme_bw()
dev.off()