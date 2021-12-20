# BiocManager::install("phyloseq")
library(phyloseq)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)

# setwd("C:/Users/kevin/Dropbox/SNU_medical/Microbiome")
setwd("D:/Dropbox/SNU_medical/Microbiome")



## Phylum ANOSIM, PERNOVA, 0511 ------------------

load("Microbiome_EDC/Pediatrics/data/phyloseq_data_twin_omit_0602.RData")

tax_dat <- tax_table(data_twin_omit) %>% as.data.frame
# tax_dat <- tax_table(data)
otu_dat <- otu_table(data_twin_omit) %>% as.data.frame
# otu_dat_rel <- apply(otu_dat, 2, function(x) x/sum(x)) %>% as.data.frame
# otu_dat_prev <- apply(otu_dat, 1, function(x) sum(x>0)) %>% as.data.frame
# sum(otu_dat_rel[,1])

tax_dat1<- tax_dat %>% select(-Kingdom) %>% mutate_all(substring, first = 3)
str(tax_dat1)
tax_dat2 <- tax_dat1 %>% transmute_all(function(x) case_when(x %>% is.na ~ "Uncharacterized",
                                                             x == "" ~ "Uncharacterized",
                                                             TRUE ~ x))
tax_dat2$OTU <- rownames(tax_dat)
otu_dat$OTU <- rownames(tax_dat)

tax_otu_dat_phylum <- tax_dat2 %>% left_join(otu_dat, by = "OTU") %>% select(-Class, -Order, -Family, -Genus, -Species, -OTU) %>%
  group_by(Phylum) %>% summarise_all(sum)

rownames(tax_otu_dat_phylum) <- tax_otu_dat_phylum$Phylum
sam <- data.frame(sample_data(data_twin_omit))

library(vegan)

## Breast_feed (Must remove NA values)
grouping <- data.frame(grouping = sam$breast_feed_cate,
                       id = rownames(sam))
na_id <- grouping$id[which(is.na(grouping$grouping))] %>% as.character()
tax_otu_dat_phylum2 <- tax_otu_dat_phylum %>% select(-Phylum, -na_id) %>% transpose()

## Anosim
bray <- vegdist(tax_otu_dat_phylum2, "bray")
test <- anosim(tax_otu_dat_phylum2, grouping$grouping %>% na.omit)

## Permanova
test %>% summary
test <- adonis(bray ~ grouping$grouping %>% na.omit, permutations = 1000)

## RDA 0511 -----------------
library(vegan)
load("Microbiome_EDC/Pediatrics/data/phyloseq_data_twin_omit_0602.RData")

tax_dat <- tax_table(data_twin_omit) %>% as.data.frame
otu_dat <- otu_table(data_twin_omit) %>% as.data.frame
tax_dat1<- tax_dat %>% select(-Kingdom) %>% mutate_all(substring, first = 3)
str(tax_dat1)
tax_dat2 <- tax_dat1 %>% transmute_all(function(x) case_when(x %>% is.na ~ "Uncharacterized",
                                                             x == "" ~ "Uncharacterized",
                                                             TRUE ~ x))
tax_dat2$OTU <- rownames(tax_dat)
otu_dat$OTU <- rownames(tax_dat)

tax_otu_dat_phylum <- tax_dat2 %>% left_join(otu_dat, by = "OTU") %>% select(-Class, -Order, -Family, -Genus, -Species, -OTU) %>%
  group_by(Phylum) %>% summarise_all(sum)

rownames(tax_otu_dat_phylum) <- tax_otu_dat_phylum$Phylum
tax_otu_dat_phylum2 <- tax_otu_dat_phylum %>% select(-Phylum) %>% t()
# tax_otu_dat_phylum2 <- tax_otu_dat_phylum %>% select(-Phylum)

sam <- data.frame(sample_data(data_twin_omit))

phylum_dat_hel <- decostand(tax_otu_dat_phylum2, "hellinger")

otu_dat1 <- otu_dat %>% select(-OTU) %>% t()
otu_dat_hel <- decostand(otu_dat1, "hellinger")

rda_fit1 <- rda(phylum_dat_hel ~ sex_cate+preterm_cate+birth_w_cate+d_method_cate+
                  breast_feed_cate + income_cate + BPA_quart_cate + DEHP_quart_cate +ets_cate+
                  BMI_Zscore_num+energy_num+exercise_num+screen_time_num, data = sam, na.action = na.omit)

rda_fit1 <- rda(phylum_dat_hel ~ d_method_cate+ BPA_85_cate + DEHP_85_cate + 
                  breast_feed_cate, data = sam, na.action = na.omit)
# plot(rda_fit1, scaling = 2)

rda_fit2 <- rda(otu_dat_hel ~ sex_cate+preterm_cate+birth_w_cate+d_method_cate+
                  breast_feed_cate + income_cate + BPA_quart_cate + DEHP_quart_cate +ets_cate+
                  BMI_Zscore_num+energy_num+exercise_num+screen_time_num, data = sam, na.action = na.omit)

rda_fit3 <- rda(otu_dat_hel ~ BMI_Zscore_num+energy_num+exercise_num+screen_time_num, data = sam, na.action = na.omit)


# smry <- summary(rda_fit3)
# smry <- summary(rda_fit2)
smry <- summary(rda_fit1)
df1  <- data.frame(smry$sites[,1:2])       # PC1 and PC2
df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
df3  <- data.frame(smry$biplot[,1:2])     # loadings for PC1 and PC2
# df3  <- data.frame(smry$centroids[,1:2])     # loadings for PC1 and PC2
rda.plot <- ggplot(df1, aes(x=RDA1, y=RDA2)) + geom_point()+
  # geom_text(aes(label=rownames(df1)),size=4) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted")
rda.plot


rda.biplot <- rda.plot +
  geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=RDA1,y=RDA2,label=rownames(df2),
                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), 
            color="red", size=3)
rda.biplot

rda.triplot <- rda.biplot +
# rda.triplot <- rda.plot +
  geom_segment(data=df3, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="blue", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df3, 
            aes(x=RDA1,y=RDA2,label=rownames(df3),
                hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), 
            color="blue", size=3)
rda.triplot
