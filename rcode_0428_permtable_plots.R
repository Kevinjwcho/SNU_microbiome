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
# data<- import_biom(BIOMfilename="Microbiome_EDC/HN00117382_Reports/HN00117382_OTU/OTU_Results/biom/otu_table.blast_NCBI_16S.biom")
# data
# edc_raw <- fread("Microbiome_EDC/edc_csv_0224.csv", header = T, fill = T)
# 
# edc_1 <- edc_raw[-c(1:4), ]
# 
# for(i in seq(2,10,2)){
#   edc_n <- paste0("edc_", i, "yr")
#   if(i == 10){
#     assign(edc_n, edc_1 %>% filter(age_C == paste0("C", i*12)))
#   }else{
#     assign(edc_n, edc_1 %>% filter(age_C == paste0("C0", i*12)))
#   }
# }
# 
# 
# ## sample 데이터의 id가 rownames로 되어 있기 때문에 이를 new_id라는 변수로 새로 만들어서 저장합니다.
# sam <- sample_data(data) %>% as.matrix %>% as.data.frame
# sam$new_id<- rownames(sam)
# rownames(sam) <- NULL
# 
# 
# ## Breast feed recoding -----------------
# 
# edc_bf <- rbind(edc_2yr, edc_4yr) %>% select(bf2, total_bf_dur, exclu_bf_dur, new_id) 
# bfbf <- edc_bf[!duplicated(edc_bf[c("new_id")]),]
# 
# which(is.na(edc_bf$bf2))
# which(is.na(edc_bf$total_bf_dur))
# which(is.na(edc_bf$exclu_bf_dur))
# 
# # edc_6yr_sel2 <- edc_6yr_sel %>% left_join(edc_2yr_bf) %>% filter(bf2 != "")
# # edc_6yr_sel2 <- edc_6yr_sel %>% left_join(bfbf)
# 
# edc_6yr <- edc_6yr %>% select(-bf2, -total_bf_dur, -exclu_bf_dur) %>% left_join(bfbf)
# 
# ## edc variable processing ----------------
# 
# str(edc_6yr %>% select(re_twin))
# edc_6yr_sel <- edc_6yr %>% transmute(t_week_num = t_week %>% as.numeric,
#                                      TSH_num = TSH %>% as.numeric,
#                                      T4Free_num = T4Free %>% as.numeric,
#                                      T3_num = T3 %>% as.numeric,
#                                      birth_w_num = BJ %>% as.numeric,
#                                      
#                                      stroop_word_num = b_stroop_1 %>% as.numeric,
#                                      stroop_col_num = b_stroop_2 %>% as.numeric,
#                                      stroop_col_word_num = b_stroop_3 %>% as.numeric,
#                                      stroop_interv_num = b_stroop_4 %>% as.numeric,
#                                      stroop_wordT_num = b_stroop_5 %>% as.numeric,
#                                      stroop_colT_num = b_stroop_6 %>% as.numeric,
#                                      stroop_intervT_num = b_stroop_7 %>% as.numeric,
#                                      stroop_wordT_num = b_stroop_8 %>% as.numeric,
#                                      
#                                      # CAT_sight_1_num = b_cat_1 %>% as.numeric,
#                                      # CAT_sight_2_num = b_cat_2 %>% as.numeric,
#                                      # CAT_sight_3_num = b_cat_3 %>% as.numeric,
#                                      # CAT_sight_4_num = b_cat_4 %>% as.numeric,
#                                      CAT_sight_5_num = b_cat_5 %>% as.numeric,
#                                      CAT_sight_6_num = b_cat_6 %>% as.numeric,
#                                      CAT_sight_7_num = b_cat_7 %>% as.numeric,
#                                      CAT_sight_8_num = b_cat_8 %>% as.numeric,
#                                      
#                                      # CAT_hear_1_num = b_cat_9 %>% as.numeric,
#                                      # CAT_hear_2_num = b_cat_10 %>% as.numeric,
#                                      # CAT_hear_3_num = b_cat_11 %>% as.numeric,
#                                      # CAT_hear_4_num = b_cat_12 %>% as.numeric,
#                                      CAT_hear_5_num = b_cat_13 %>% as.numeric,
#                                      CAT_hear_6_num = b_cat_14 %>% as.numeric,
#                                      CAT_hear_7_num = b_cat_15 %>% as.numeric,
#                                      CAT_hear_8_num = b_cat_16 %>% as.numeric,
#                                      
#                                      # CAT_resist_1_num = b_cat_17 %>% as.numeric,
#                                      # CAT_resist_2_num = b_cat_18 %>% as.numeric,
#                                      # CAT_resist_3_num = b_cat_19 %>% as.numeric,
#                                      # CAT_resist_4_num = b_cat_20 %>% as.numeric,
#                                      CAT_resist_5_num = b_cat_21 %>% as.numeric,
#                                      CAT_resist_6_num = b_cat_22 %>% as.numeric,
#                                      CAT_resist_7_num = b_cat_23 %>% as.numeric,
#                                      CAT_resist_8_num = b_cat_24 %>% as.numeric,
#                                      
#                                      # CAT_interv_1_num = b_cat_25 %>% as.numeric,
#                                      # CAT_interv_2_num = b_cat_26 %>% as.numeric,
#                                      # CAT_interv_3_num = b_cat_27 %>% as.numeric,
#                                      # CAT_interv_4_num = b_cat_28 %>% as.numeric,
#                                      CAT_interv_5_num = b_cat_29 %>% as.numeric,
#                                      CAT_interv_6_num = b_cat_30 %>% as.numeric,
#                                      CAT_interv_7_num = b_cat_31 %>% as.numeric,
#                                      CAT_interv_8_num = b_cat_32 %>% as.numeric,
#                                      
#                                      
#                                      ARS_total_num = ARS_total %>% as.numeric,
#                                      SCQ_total_num = SCQtotal %>% as.numeric,
#                                      FSIQ_num = b_wisc_1 %>% as.numeric,
#                                      
#                                      BPA_num = BPA_creatinine %>% as.numeric,
#                                      MEHHP_num = MEHHP_creatinine_2 %>% as.numeric,
#                                      MEOHP_num = MEOHP_creatinine_2 %>% as.numeric,
#                                      MnBP_num = MnBP_creatinine_2 %>% as.numeric,
#                                      
#                                      ln_ARS_total_num = ARS_total %>% as.numeric %>% log,
#                                      ln_SCQ_total_num = SCQtotal %>% as.numeric %>% log,
#                                      ln_FSIQ_num = b_wisc_1 %>% as.numeric %>% log,
#                                      
#                                      ln_BPA_num = BPA_creatinine %>% as.numeric %>% log,
#                                      ln_MEHHP_num = MEHHP_creatinine_2 %>% as.numeric %>% log,
#                                      ln_MEOHP_num = MEOHP_creatinine_2 %>% as.numeric %>% log,
#                                      ln_MnBP_num = MnBP_creatinine_2 %>% as.numeric %>% log,
#                                      
#                                      # re_twin_cate = ifelse(re_twin == "1", "twin", "single"),
#                                      breast_feed_cate = ifelse(bf2 %>% as.numeric ==1, "b_milk", "others") %>% as.factor %>% factor(levels = c("others", "b_milk")),
#                                      d_method_cate = ifelse(d_method %in% c(0, 1), "natural", "CS") %>% as.factor %>% factor(levels = c("natural", "CS")),
#                                      BMI_cate = ifelse(BMI_Zscore %>% as.numeric >= 1.04, "over", "normal") %>% as.factor,
#                                      TSH_cate = ifelse(TSH_num < 4.94, "normal", "over") %>% as.factor,
#                                      
#                                      BPA_50_cate = ifelse(BPA_num > quantile(BPA_num, 0.5, na.rm = T), "over", "less") %>% as.factor,
#                                      MEHHP_50_cate = ifelse(MEHHP_num > quantile(MEHHP_num, 0.5, na.rm = T), "over", "less") %>% as.factor,
#                                      MEOHP_50_cate = ifelse(MEOHP_num > quantile(MEOHP_num, 0.5, na.rm = T), "over", "less") %>% as.factor,
#                                      MnBP_50_cate = ifelse(MnBP_num > quantile(MnBP_num, 0.5, na.rm = T), "over", "less") %>% as.factor,
#                                      
#                                      BPA_85_cate = ifelse(BPA_num >= quantile(BPA_num, 0.85, na.rm = T), "over", "less") %>% as.factor,
#                                      MEHHP_85_cate = ifelse(MEHHP_num >= quantile(MEHHP_num, 0.85, na.rm = T), "over", "less") %>% as.factor,
#                                      MEOHP_85_cate = ifelse(MEOHP_num >= quantile(MEOHP_num, 0.85, na.rm = T), "over", "less") %>% as.factor,
#                                      MnBP_85_cate = ifelse(MnBP_num >= quantile(MnBP_num, 0.85, na.rm = T), "over", "less") %>% as.factor,
#                                      
#                                      ARS_cate = ifelse(ARS_total_num >= 19, "ADHD", "normal") %>% as.factor %>% factor(levels= c("normal", "ADHD")),
#                                      SCQ_cate = ifelse(SCQ_total_num >= 16, "ASD", "normal") %>% as.factor %>% factor(levels= c("normal", "ASD")),
#                                      # FSIQ_cate = ifelse(FSIQ_num < 70, "less", "normal") %>% as.factor %>% factor(levels= c("normal", "less")),
#                                      
#                                      CAT_sight_5_cate = ifelse(CAT_sight_5_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_sight_6_cate = ifelse(CAT_sight_6_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_sight_7_cate = ifelse(CAT_sight_7_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_sight_8_cate = ifelse(CAT_sight_8_num >= 76, "over", "less") %>% as.factor,
#                                      
#                                      CAT_hear_5_cate = ifelse(CAT_hear_5_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_hear_6_cate = ifelse(CAT_hear_6_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_hear_7_cate = ifelse(CAT_hear_7_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_hear_8_cate = ifelse(CAT_hear_8_num >= 76, "over", "less") %>% as.factor,
#                                      
#                                      CAT_resist_5_cate = ifelse(CAT_resist_5_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_resist_6_cate = ifelse(CAT_resist_6_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_resist_7_cate = ifelse(CAT_resist_7_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_resist_8_cate = ifelse(CAT_resist_8_num >= 76, "over", "less") %>% as.factor,
#                                      
#                                      CAT_interv_5_cate = ifelse(CAT_interv_5_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_interv_6_cate = ifelse(CAT_interv_6_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_interv_7_cate = ifelse(CAT_interv_7_num >= 76, "over", "less") %>% as.factor,
#                                      CAT_interv_8_cate = ifelse(CAT_interv_8_num >= 76, "over", "less") %>% as.factor,
#                                       
#                                      new_id = sapply(strsplit(new_id, "-"), function(x) paste(x, collapse = ".")))
# 
# 
# str(edc_6yr_sel)
# # table(edc_6yr_sel$re_twin_cate)
# ## 합친 다음 다시 rownames로 new_id를 넣어주고 sample_data에 저장합니다.
# sam_merge <- sam %>% select(new_id) %>% left_join(edc_6yr_sel)
# str(sam_merge)
# # table(sam_merge$re_twin_cate)
# 
# rownames(sam_merge) <- sam_merge$new_id
# sam_merge <- sam_merge %>% select(-new_id)
# sample_data(data) <- sam_merge
# 
# 
# ## Phylogenetic Tree 생성방법
# # install.packages("ape")
# # install.packages("adegenet")
# library(adegenet)
# library(ape)
# 
# # dna <-fasta2DNAbin("F:/work/Microbiome_EDC/HN00117382_Reports/otus_rep.fasta")
# # class(dna)
# # dna
# 
# ## Raw 방법은 NULL 값이 없습니다. -> nj() 함수를 쓸수 있습니다.
# # D <- dist.dna(dna, model = "raw")
# # length(D)
# # tre <- nj(D)
# 
# ## TN93 방법은 NULL 값이 많습니다. -> njs() 함수를 써야합니다. 
# ## 그리고 오래걸립니다.
# # D <- dist.dna(dna, model = "TN93")
# # length(D)
# # tre <- njs(D)
# 
# ## tree를 한번 만든다음 저장해서 필요할 때마다 불러오는 것이 좋습니다.
# # save(tre, file = "F:/work/Microbiome_EDC/RData/TN93_tre.RData")
# load("Microbiome_EDC/RData/TN93_tre.RData")
# 
# 
# ## Tree label에 id가 달려있는 것을 제외합니다.
# tre$tip.label <- sapply(strsplit(tre$tip.label, "\t"), function(x)x[1])
# phy_tree(data) <- tre
# 
# 
# otu <- otu_table(data)
# sample_dat <- sample_data(data)
# tax <- tax_table(data)
# phy_tre <- phy_tree(data)
# 
# 
# ## Microbiome processing --------------
# rank_names(data)
# tax_table(data) %>% colnames
# colnames(tax_table(data))<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# rank_names(data)
# 
# 
# save(data, edc_6yr_sel, file = "Microbiome_EDC/Rcode/markdown/report_0406/phyloseq_data.RData")

# load("Microbiome_EDC/Rcode/markdown/report_0406/phyloseq_data.RData")
load("Microbiome_EDC/Pediatrics/twin_omit.RData")



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


## Phylum stack graph 0428 ------------------

load("Microbiome_EDC/Pediatrics/twin_omit.RData")
tax_dat <- tax_table(data_twin_omit) %>% as.data.frame
# tax_dat <- tax_table(data)
otu_dat <- otu_table(data_twin_omit) %>% as.data.frame
otu_dat_rel <- apply(otu_dat, 2, function(x) x/sum(x)) %>% as.data.frame
# otu_dat_prev <- apply(otu_dat, 1, function(x) sum(x>0)) %>% as.data.frame
# sum(otu_dat_rel[,1])

tax_dat1<- tax_dat %>% select(-Kingdom) %>% mutate_all(substring, first = 3)
str(tax_dat1)
tax_dat2 <- tax_dat1 %>% transmute_all(function(x) case_when(x %>% is.na ~ "Uncharacterized",
                                                             x == "" ~ "Uncharacterized",
                                                             TRUE ~ x))
tax_dat2$OTU <- rownames(tax_dat)
otu_dat$OTU <- rownames(otu_dat)
tax_otu_dat <- tax_dat2 %>% left_join(otu_dat)
# tax_otu_dat <- tax_otu_dat %>% filter_all(all_vars(. !="Uncharacterized"))

phylum_otu_dat <- tax_otu_dat %>% select(-Class, -Order, -Family, -Genus, -Species, -OTU)
phylum_otu_dat_melt <- phylum_otu_dat %>% reshape2::melt(id = 1, variable.name = "new_id", value.name = "Abundance")

edc_dat$new_id <- rownames(edc_dat)
edc_dat_cate <- edc_dat %>% select(new_id, ends_with("_cate"))
phylum_edc_join <- phylum_otu_dat_melt %>% left_join(edc_dat_cate, by = "new_id")
str(phylum_edc_join)

library(RColorBrewer)
pdf("Microbiome_EDC/Pediatrics/Phylum_plot/absolute_value.pdf")
for(i in 4:ncol(phylum_edc_join)){
  var_n <- colnames(phylum_edc_join)[i]
  p1 <- ggplot(phylum_edc_join, aes(x = var_n %>% get(), y = Abundance, fill = Phylum))+geom_col(position = "fill",width=.5)+
    ggtitle(substr(var_n, 1, nchar(var_n)-5))+xlab("") + theme_minimal()+ylab("Absolute proportion")
  print(p1)
  cat(i, "\n")
}
dev.off()


## Propotional phylum difference plot ----------------

load("Microbiome_EDC/Pediatrics/twin_omit.RData")
tax_dat <- tax_table(data_twin_omit) %>% as.data.frame
# tax_dat <- tax_table(data)
otu_dat <- otu_table(data_twin_omit) %>% as.data.frame
otu_dat_rel <- apply(otu_dat, 2, function(x) x/sum(x)) %>% as.data.frame
# otu_dat_prev <- apply(otu_dat, 1, function(x) sum(x>0)) %>% as.data.frame
# sum(otu_dat_rel[,1])

tax_dat1<- tax_dat %>% select(-Kingdom) %>% mutate_all(substring, first = 3)
str(tax_dat1)
tax_dat2 <- tax_dat1 %>% transmute_all(function(x) case_when(x %>% is.na ~ "Uncharacterized",
                                                             x == "" ~ "Uncharacterized",
                                                             TRUE ~ x))
tax_dat2$OTU <- rownames(tax_dat)
otu_dat_rel$OTU <- rownames(otu_dat_rel)
tax_otu_dat <- tax_dat2 %>% left_join(otu_dat_rel)
# tax_otu_dat <- tax_otu_dat %>% filter_all(all_vars(. !="Uncharacterized"))

phylum_otu_dat <- tax_otu_dat %>% select(-Class, -Order, -Family, -Genus, -Species, -OTU)
phylum_otu_dat_melt <- phylum_otu_dat %>% reshape2::melt(id = 1, variable.name = "new_id", value.name = "Abundance")
edc_dat$new_id <- rownames(edc_dat)
edc_dat_cate <- edc_dat %>% select(new_id, ends_with("_cate"))
phylum_edc_join <- phylum_otu_dat_melt %>% left_join(edc_dat_cate, by = "new_id")



# library(RColorBrewer)
pdf("Microbiome_EDC/Pediatrics/Phylum_plot/mean_relative_value.pdf")
for(i in 4:ncol(phylum_edc_join)){
  var_n <- colnames(phylum_edc_join)[i]
  phylum_edc_join_sum <- phylum_edc_join %>% group_by_at(vars(new_id, Phylum, var_n)) %>% summarise(sum_ab = sum(Abundance))
  phylum_edc_join_mean <- phylum_edc_join_sum %>% group_by_at(vars(Phylum, var_n)) %>% summarise(mean_ab = mean(sum_ab))
  p1 <- ggplot(phylum_edc_join_mean, aes(x = var_n %>% get(), y = mean_ab, fill = Phylum))+geom_col(position = "fill",width=.5)+
    ggtitle(substr(var_n, 1, nchar(var_n)-5))+xlab("") + theme_minimal()+ylab("Mean proportion")
  print(p1)
  cat(i, "\n")
}
dev.off()



