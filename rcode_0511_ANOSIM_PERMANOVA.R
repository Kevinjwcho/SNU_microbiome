# BiocManager::install("phyloseq")
library(phyloseq)
library(dplyr)
library(data.table)
library(stringr)

library(ggplot2)
library(gridExtra)
library(grid)

setwd("C:/Users/kevin/Dropbox/SNU_medical/Microbiome")

## PERNOVA test by phylogenetic 0405 --------------------

tax_dat <- tax_table(data_twin_omit) %>% as.data.frame
# tax_dat <- tax_table(data)
otu_dat <- otu_table(data_twin_omit) %>% as.data.frame
otu_dat_rel <- apply(otu_dat, 2, function(x) x/sum(x)) %>% as.data.frame

# sum(otu_dat_rel[,1])

tax_dat1<- tax_dat %>% select(-Kingdom) %>% mutate_all(substring, first = 3)
str(tax_dat1)
tax_dat2 <- tax_dat1 %>% transmute_all(function(x) case_when(x %>% is.na ~ "Uncharacterized",
                                                             x == "" ~ "Uncharacterized",
                                                             TRUE ~ x))
tax_dat2$OTU <- rownames(tax_dat)
otu_dat_rel$OTU <- rownames(otu_dat_rel)

tax_otu_dat <- tax_dat2 %>% left_join(otu_dat_rel)

phylum_otu_dat <- tax_otu_dat %>% select(-Class, -Order, -Species, -OTU, -Family, - Genus)
family_otu_dat <- tax_otu_dat %>% select(-Class, -Order, -Species, -OTU, -Phylum, - Genus)
genus_otu_dat <- tax_otu_dat %>% select(-Class, -Order, -Species, -OTU, -Phylum, - Family)



# preval_fun <- function(phylum_otu_dat, rank, data){
#   phylum_sum_dat <- phylum_otu_dat %>% group_by_at(vars(rank)) %>% summarise_all(sum) %>% as.matrix
#   rownames(phylum_sum_dat) <- phylum_sum_dat[,1]
#   phylum_sum_dat <- t(phylum_sum_dat[,-1]) %>% as.data.frame
#   phylum_sum_dat$new_id <- c(rownames(phylum_sum_dat))
#   phylum_sum_dat_pres <- phylum_sum_dat %>% select(-new_id) %>% apply(2, function(x)sum(x>0))
#   phylum_sum_dat_pres_rel <- phylum_sum_dat_pres / nrow(sample_data(data))
#   delete_candi <- which(phylum_sum_dat_pres_rel < 0.05) %>% names()
#   return(list(Prevalence_tab = data.frame(Prevalence = phylum_sum_dat_pres,
#                                           Relative = phylum_sum_dat_pres_rel),
#               Delete_candi = delete_candi))
# }
# 
# phylum_prev <- preval_fun(phylum_otu_dat, rank = "Phylum", data)
# family_prev <- preval_fun(family_otu_dat, rank = "Family", data)
# genus_prev <- preval_fun(genus_otu_dat, rank = "Genus", data)

phylum_edc_join <- function(phylum_otu_dat, rank, data, edc_6yr_sel){
  phylum_sum_dat <- phylum_otu_dat %>% group_by_at(vars(rank)) %>% summarise_all(sum) %>% as.matrix
  rownames(phylum_sum_dat) <- phylum_sum_dat[,1]
  phylum_sum_dat <- t(phylum_sum_dat[,-1]) %>% as.data.frame
  phylum_sum_dat$new_id <- c(rownames(phylum_sum_dat))
  phylum_sum_dat_melt <- phylum_sum_dat %>% melt(id.vars = "new_id", variable.name = rank, value.name = "Frequency")
  phylum_dat <- phylum_sum_dat_melt %>% left_join(edc_6yr_sel)
  phylum_dat <- phylum_dat %>% mutate(Frequency = Frequency %>% as.numeric)
  return(phylum_dat)
}

edc_dat$new_id <- edc_dat %>% rownames()
phylum_total <- phylum_edc_join(phylum_otu_dat, rank = "Phylum", data, edc_dat)
family_total <- phylum_edc_join(family_otu_dat, rank = "Family", data, edc_dat)
genus_total <- phylum_edc_join(genus_otu_dat, rank = "Genus", data, edc_dat)

# phylum_filter <- phylum_total %>% filter(!Phylum %in% phylum_prev$Delete_candi)
# family_filter <- family_total %>% filter(!Family %in% family_prev$Delete_candi)
# genus_filter <- genus_total %>% filter(!Genus %in% genus_prev$Delete_candi)
# 
# library(ggplot2)
# ggplot(phylum_filter, aes(y = Frequency, x = Phylum, fill = d_method_cate))+geom_boxplot()
# ggplot(family_filter, aes(y = Frequency, x = Family, fill = d_method_cate))+geom_boxplot()
# ggplot(genus_filter, aes(y = Frequency, x = Genus, fill = d_method_cate))+geom_boxplot()
# 
# filter_dataset <- list(phylum_filter, family_filter, genus_filter)


phylum_cate <- phylum_total %>% select(new_id, Phylum, Frequency, ends_with("_cate")) %>% filter(Phylum != "Uncharacterized")
phylum_candi <- phylum_cate$Phylum %>% unique %>% as.character
family_cate <- family_total %>% select(new_id, Family, Frequency, ends_with("_cate")) %>% filter(Family != "Uncharacterized")
family_candi <- family_cate$Family %>% unique %>% as.character
genus_cate <- genus_total %>% select(new_id, Genus, Frequency, ends_with("_cate")) %>% filter(Genus != "Uncharacterized")
genus_candi <- genus_cate$Genus %>% unique %>% as.character

cate_candi <- colnames(phylum_cate)[-c(1:3)]

permutation_test <- function(data = phy_sub, M = 10000){
  y_var <- data[,1]
  x_var <- data[,2]
  result <- lm(y_var ~ x_var) %>% summary
  F0 <- result$fstatistic[1]
  N = nrow(data)
  F.zp = numeric(M)
  for (i in 1:M){
    xvar_z = sample(x_var,N)
    result = lm(y_var ~ xvar_z ) %>% summary
    F.zp[i] = result$fstatistic[1]
  }
  
  
  ASL = (1+sum(F.zp >= F0)) / (M+1)
  effect_n <- rownames(result$coefficients)[2]
  p_result <- data.frame(p_v = ASL, 
                         Effect = effect_n %>% substr(7, nchar(effect_n)),
                         Nrow = N)
  return(p_result)
}
cc <- c("phylum", "family", "genus")
for(k in 1:length(cc)){
  r_n <- paste0(cc[k], "_", "result")
  rr <- data.frame()
  candi_n <- paste0(cc[k], "_candi")
  cate_n <- paste0(cc[k], "_cate")
  candi <- candi_n %>% get
  cate <- cate_n %>% get
  colnames(cate)[2] <- "kingdom"
  for(j in 1:length(cate_candi)){
    r_sub <- data.frame()
    for(i in 1:length(candi)){
      kingking <- candi[i]
      catecate <- cate_candi[j]
      cat(r_n, ":", "\n",
          "Kingdom:", kingking, "(", i, "out of", length(candi), ")\n",
          "Categorical Variable:", catecate, "(", j, "out of", length(cate_candi), ")\n")
      
      king_sub <- cate %>% filter(kingdom == kingking) %>% select(Frequency, catecate) %>% na.omit
      result <- permutation_test(data = king_sub, M = 10^4)
      colnames(result) <- paste0(substr(catecate, 1, (nchar(catecate)-5)),"_", colnames(result))
      r_sub <- rbind(r_sub, result)
      
    }
    if(j == 1){rr <- r_sub
    }else{
      rr <- cbind(rr, r_sub) 
    }
  }
  rownames(rr) <- candi
  f_dir <- paste0("Microbiome_EDC/Pediatrics/result/", r_n,"_rel.csv")
  fwrite(rr, f_dir, row.names = T)
}

## Prevalance measure ------------------------------
data_twin_omit
tax_dat <- tax_table(data_twin_omit) %>% as.data.frame
# tax_dat <- tax_table(data)
otu_dat <- otu_table(data_twin_omit) %>% as.data.frame
# otu_dat_rel <- apply(otu_dat, 2, function(x) x/sum(x)) %>% as.data.frame
otu_dat_prev <- apply(otu_dat, 1, function(x) sum(x>0)) %>% as.data.frame
# sum(otu_dat_rel[,1])

tax_dat1<- tax_dat %>% select(-Kingdom) %>% mutate_all(substring, first = 3)
str(tax_dat1)
tax_dat2 <- tax_dat1 %>% transmute_all(function(x) case_when(x %>% is.na ~ "Uncharacterized",
                                                             x == "" ~ "Uncharacterized",
                                                             TRUE ~ x))
tax_dat2$OTU <- rownames(tax_dat)
otu_dat_prev$OTU <- rownames(otu_dat_prev)
tax_otu_dat <- tax_dat2 %>% left_join(otu_dat_prev)
# tax_otu_dat <- tax_otu_dat %>% filter_all(all_vars(. !="Uncharacterized"))

phylum_otu_dat <- tax_otu_dat %>% select(-Class, -Order, -Species, -OTU, -Family, - Genus)
family_otu_dat <- tax_otu_dat %>% select(-Class, -Order, -Species, -OTU, -Phylum, - Genus)
genus_otu_dat <- tax_otu_dat %>% select(-Class, -Order, -Species, -OTU, -Phylum, - Family)

phylum_prev<- phylum_otu_dat %>% group_by(Phylum) %>% summarise_all(sum)
family_prev<- family_otu_dat %>% group_by(Family) %>% summarise_all(sum)
genus_prev<- genus_otu_dat %>% group_by(Genus) %>% summarise_all(sum)

family_prev <- family_prev %>% left_join(tax_otu_dat %>% select(Phylum, Class, Order, Family) %>% unique(), by = "Family")
genus_prev <- genus_prev %>% left_join(tax_otu_dat %>% select(Phylum, Class, Order, Family, Genus) %>% unique(), by = "Genus")

phylum_rel <- fread("Microbiome_EDC/Pediatrics/Permutation_ANOVA/twin/phylum_result_rel.csv", header = T)
family_rel <- fread("Microbiome_EDC/Pediatrics/Permutation_ANOVA/twin/family_result_rel.csv", header = T)
genus_rel <- fread("Microbiome_EDC/Pediatrics/Permutation_ANOVA/twin/genus_result_rel.csv", header = T)

colnames(phylum_rel)[1] <- "Phylum"
phylum_rel_new <- phylum_rel %>% select(!(ends_with("Effect"))) %>% select(!(ends_with("Nrow")))
phylum_rel_new <- phylum_rel_new %>% left_join(phylum_prev, by = "Phylum")
phylum_rel_new <- phylum_rel_new[,c(1, ncol(phylum_rel_new), 3:(ncol(phylum_rel_new)-1))]
colnames(phylum_rel_new)[2] <- "Prevalance"


colnames(family_rel)[1] <- "Family"
family_rel_new <- family_rel %>% select(!(ends_with("Effect"))) %>% select(!(ends_with("Nrow")))
family_rel_new <- family_rel_new %>% left_join(family_prev, by = "Family")
family_rel_new <- family_rel_new[,c((ncol(family_rel_new)-2):(ncol(family_rel_new)), 1,
                                    (ncol(family_rel_new)-3), 2:(ncol(family_rel_new)-4))]
colnames(family_rel_new)[5] <- "Prevalance"


colnames(genus_rel)[1] <- "Genus"
genus_rel_new <- genus_rel %>% select(!(ends_with("Effect"))) %>% select(!(ends_with("Nrow")))
genus_rel_new <- genus_rel_new %>% left_join(genus_prev, by = "Genus")
genus_rel_new <- genus_rel_new[,c((ncol(genus_rel_new)-3):(ncol(genus_rel_new)), 1,
                                    (ncol(genus_rel_new)-4), 2:(ncol(genus_rel_new)-5))]
colnames(genus_rel_new)[6] <- "Prevalance"

fwrite(phylum_rel_new, file = "Microbiome_EDC/Pediatrics/Permutation_ANOVA/twin_with_taxa_prev/phylum_perm_result.csv")
fwrite(family_rel_new, file = "Microbiome_EDC/Pediatrics/Permutation_ANOVA/twin_with_taxa_prev/family_perm_result.csv")
fwrite(genus_rel_new, file = "Microbiome_EDC/Pediatrics/Permutation_ANOVA/twin_with_taxa_prev/genus_perm_result.csv")


## Phylum ANOSIM, PERNOVA, 0511 ------------------

load("Microbiome_EDC/Pediatrics/data/phyloseq_data_twin_omit_0510.RData")

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

bray <- vegdist(tax_otu_dat_phylum2, "bray")
test <- anosim(tax_otu_dat_phylum2, grouping$grouping %>% na.omit)

test %>% summary
test <- adonis(bray ~ grouping$grouping %>% na.omit, permutations = 1000)

## RDA 0511 -----------------
library(vegan)
load("Microbiome_EDC/Pediatrics/data/phyloseq_data_twin_omit_0510.RData")

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
