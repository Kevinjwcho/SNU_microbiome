## Installing required packages -----------
# install.packages("BiocManager")
.cran_packages <- c("dplyr", "data.table", "stringr", "ggplot2", "gridExtra", "grid")
.bioc_packages <- c("phyloseq")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(.bioc_packages[!.inst])
}


## Set the work directory and load .biom data -----------------
setwd("D:/Dropbox/SNU_medical/Microbiome/Microbiome_EDC/Pediatrics") # denote the location of "sending" folder.
# setwd("C:/Users/kevin/Dropbox/SNU_medical/Microbiome/Microbiome_EDC/Pediatrics") # denote the location of "sending" folder.
data<- import_biom(BIOMfilename="sending/otu_table.blast_NCBI_16S.biom")
data # You can see the structure of data.

## Load EDC data -------------------------
edc_raw <- fread("data/EDC_age8_v1_4_csv.csv", header = T)
edc_raw1 <- edc_raw[-c(1:4), ] # Remove the first to forth row of the data.

## EDC 데이터를 나이에 따라서 나눕니다. data는 edc_2yr, edc_4yr ,... 이런식으로 저장됩니다.
for(i in seq(2,10,2)){
  edc_n <- paste0("edc_", i, "yr")
  if(i == 10){
    assign(edc_n, edc_raw1 %>% filter(age_C == paste0("C", i*12)))
  }else{
    assign(edc_n, edc_raw1 %>% filter(age_C == paste0("C0", i*12)))
  }
}


## sample 데이터의 id가 rownames로 되어 있기 때문에 이를 new_id라는 변수로 새로 만들어서 저장합니다.
sam <- sample_data(data) %>% as.matrix %>% as.data.frame
sam$new_id<- rownames(sam)
rownames(sam) <- NULL

# Breast feeding 변수 중에 bf2, total_bf_dur, exclu_bf_dur 3개의 변수가 6세에 없기 때문에 2세와 4세에서 불러와야 합니다.
# 그걸 6세에 붙이는 과정이라고 보시면 됩니다.
edc_bf <- rbind(edc_2yr, edc_4yr) %>% select(bf2, total_bf_dur, exclu_bf_dur, new_id)
bfbf <- edc_bf[!duplicated(edc_bf[c("new_id")]),]

# which(is.na(edc_bf$bf2))
# which(is.na(edc_bf$total_bf_dur))
# which(is.na(edc_bf$exclu_bf_dur))

edc_6yr <- edc_6yr %>% select(-bf2, -total_bf_dur, -exclu_bf_dur) %>% left_join(bfbf)


## edc variable processing ----------------
## 여기서 edc 변수를 설정하시면 됩니다.

edc_6yr_sel <- edc_6yr %>% transmute(age_num = age %>% as.numeric,
                                     BMI_Zscore_num = BMI_Zscore %>% as.numeric,
                                     energy_num = energy %>% as.numeric,
                                     exercise_num = severity_pa_freq %>% as.numeric * time_severity_pa %>% as.numeric +serious_pa_freq %>% as.numeric * time_serious_pa %>% as.numeric,
                                     screen_time_num = time_tv_com_ph_game %>% as.numeric,
                                     ln_BPA_num = BPA_creatinine %>% as.numeric %>% log,
                                     ln_DEHP_num = log(MEHHP_creatinine_2 %>% as.numeric + MEOHP_creatinine_2 %>% as.numeric),
                                     
                                     pat_BMI_num = fa_wei %>% as.numeric*10^4/(fa_hei %>% as.numeric)^2,
                                     mat_BMI_num = mo_wei %>% as.numeric*10^4/(mo_hei %>% as.numeric)^2,
                                     
                                     TSH_num = TSH %>% as.numeric,
                                     T4Free_num = T4Free %>% as.numeric,
                                     T3_num = T3 %>% as.numeric,
                                     T_CHO_num = T_CHO %>% as.numeric,
                                     TG_num = TG %>% as.numeric,
                                     HDL_C_num = HDL_C %>% as.numeric,
                                     LDL_C_num = LDL_C %>% as.numeric,
                                     Schwartz_eGFR_1_num = Schwartz_eGFR_1 %>% as.numeric,
                                     GOT_num = GOT %>% as.numeric,
                                     GPT_num = GPT %>% as.numeric,
                                     r_GTP_num = r_GTP %>% as.numeric,
                                     BUN_num = BUN %>% as.numeric,
                                     b_CREATININE_num = b_CREATININE %>% as.numeric,
                                     Uric_acid_num = Uric_acid %>% as.numeric,
                                     Glucose_num = Glucose %>% as.numeric,
                                     Hb_num = Hb %>% as.numeric,
                                     
                                     Adiponectin_num = Adiponectin %>% as.numeric,
                                     IL_6_num = IL_6 %>% as.numeric,
                                     leptin_num = leptin %>% as.numeric,
                                     insulin_num = insulin %>% as.numeric,
                                     hsCRP_num = hsCRP %>% as.numeric,
                                     u_Creatinine_num = u_Creatinine %>% as.numeric,
                                     
                                     homair_num = Glucose_num*insulin_num/405,
                                     
                                     sex_cate = baby_gen %>% as.factor,
                                     preterm_cate = ifelse(t_week %>% as.numeric <37, "less", "normal") %>% as.factor,
                                     d_method_cate = ifelse(re_d_method == 1, "natural", "CS") %>% as.factor,
                                     birth_w_cate = ifelse(BJ %>% as.numeric <2.5, "less", "normal") %>% as.factor,
                                     breast_feed_cate = ifelse(bf2 %>% as.numeric ==1, "b_milk", "others") %>% as.factor %>% factor(levels = c("others", "b_milk")),
                                     
                                     dm_hx_cate = ifelse(par_diab %>% as.numeric > 0, "yes", "no") %>% as.factor,
                                     thyroid_hx_cate = ifelse(par_thyroid %>% as.numeric >0 | par_thy_tumor >0, "yes", "no") %>% as.factor,
                                     
                                     income_cate = ifelse(age6_parent_dis9 %>% as.numeric >= 3, "over", "less") %>% as.factor,
                                     pat_edu_cate = ifelse(edu_fa %>% as.numeric >= 3, "over", "less") %>% as.factor,
                                     mat_edu_cate = ifelse(edu_mo %>% as.numeric >= 3, "over", "less") %>% as.factor,
                                     
                                     ets_cate = ifelse(a_ets %>% as.numeric ==1, "yes", "no") %>% as.factor,
                                     
                                     BPA_50_cate = ifelse(ln_BPA_num > quantile(ln_BPA_num, 0.5, na.rm = T), "over", "less") %>% as.factor,
                                     BPA_85_cate = ifelse(ln_BPA_num >= quantile(ln_BPA_num, 0.85, na.rm = T), "over", "less") %>% as.factor,
                                     BPA_quart_cate = case_when(
                                       ln_BPA_num < quantile(ln_BPA_num, 0.25, na.rm = T) ~ "1st",
                                       ln_BPA_num < quantile(ln_BPA_num, 0.5, na.rm = T) ~ "2nd",
                                       ln_BPA_num < quantile(ln_BPA_num, 0.75, na.rm = T) ~ "3rd",
                                       TRUE ~ "4th") %>% as.factor,
                                     
                                     DEHP_50_cate = ifelse(ln_DEHP_num > quantile(ln_DEHP_num, 0.5, na.rm = T), "over", "less") %>% as.factor,
                                     DEHP_85_cate = ifelse(ln_DEHP_num >= quantile(ln_DEHP_num, 0.85, na.rm = T), "over", "less") %>% as.factor,
                                     DEHP_quart_cate = case_when(
                                       ln_DEHP_num < quantile(ln_DEHP_num, 0.25, na.rm = T) ~ "1st",
                                       ln_DEHP_num < quantile(ln_DEHP_num, 0.5, na.rm = T) ~ "2nd",
                                       ln_DEHP_num < quantile(ln_DEHP_num, 0.75, na.rm = T) ~ "3rd",
                                       TRUE ~ "4th") %>% as.factor,
                                     
                                     BMF_6m_cate = ifelse(total_bf_dur >= 6, "over", "less") %>% as.factor,
                                     BMI_cate = ifelse(BMI_Zscore_num >= 1.036, "over", "less") %>% as.factor,
                                     TSH_cate = ifelse(TSH_num > 4.94, "over", "less") %>% as.factor,
                                     
                                     # MEHHP_50_cate = ifelse(MEHHP_num > quantile(MEHHP_num, 0.5, na.rm = T), "over", "less") %>% as.factor,
                                     # MEOHP_50_cate = ifelse(MEOHP_num > quantile(MEOHP_num, 0.5, na.rm = T), "over", "less") %>% as.factor,
                                     # MnBP_50_cate = ifelse(MnBP_num > quantile(MnBP_num, 0.5, na.rm = T), "over", "less") %>% as.factor,
                                     # 
                                     # 
                                     # MEHHP_85_cate = ifelse(MEHHP_num >= quantile(MEHHP_num, 0.85, na.rm = T), "over", "less") %>% as.factor,
                                     # MEOHP_85_cate = ifelse(MEOHP_num >= quantile(MEOHP_num, 0.85, na.rm = T), "over", "less") %>% as.factor,
                                     # MnBP_85_cate = ifelse(MnBP_num >= quantile(MnBP_num, 0.85, na.rm = T), "over", "less") %>% as.factor,
                                     # 
                                     # 
                                     # BPA_num = BPA_creatinine %>% as.numeric,
                                     # MEHHP_num = MEHHP_creatinine_2 %>% as.numeric,
                                     # MEOHP_num = MEOHP_creatinine_2 %>% as.numeric,
                                     # MnBP_num = MnBP_creatinine_2 %>% as.numeric,
                                     # 
                                     # ln_BPA_num = BPA_creatinine %>% as.numeric %>% log,
                                     # ln_MEHHP_num = MEHHP_creatinine_2 %>% as.numeric %>% log,
                                     # ln_MEOHP_num = MEOHP_creatinine_2 %>% as.numeric %>% log,
                                     # ln_MnBP_num = MnBP_creatinine_2 %>% as.numeric %>% log,
                                     # 
                                     # # re_twin_cate = ifelse(re_twin == "1", "twin", "single"),
                                     # breast_feed_cate = ifelse(bf2 %>% as.numeric ==1, "b_milk", "others") %>% as.factor %>% factor(levels = c("others", "b_milk")),
                                     # total_bf_dur = total_bf_dur %>% as.numeric,
                                     # d_method_cate = ifelse(d_method %in% c(0, 1), "natural", "CS") %>% as.factor %>% factor(levels = c("natural", "CS")),
                                     # BMI_cate = ifelse(BMI_Zscore %>% as.numeric >= 1.04, "over", "normal") %>% as.factor,
                                     # TSH_cate = ifelse(TSH_num < 4.94, "normal", "over") %>% as.factor,
                                     # 
                                     # 
                                     new_id = sapply(strsplit(new_id, "-"), function(x) paste(x, collapse = ".")),
                                     mo_id = MO_id
                                     )


str(edc_6yr_sel)
# table(edc_6yr_sel$re_twin_cate)
## 합친 다음 다시 rownames로 new_id를 넣어주고 sample_data에 저장합니다.
sam_merge <- sam %>% select(new_id) %>% left_join(edc_6yr_sel)
str(sam_merge)
# table(sam_merge$re_twin_cate)

rownames(sam_merge) <- sam_merge$new_id
sam_merge <- sam_merge %>% select(-new_id)
sample_data(data) <- sam_merge


## 제가 미리 만든 phylogentic tree 파일을 불러서 data에 넣어줍니다.
load("sending/TN93_tre.RData")


## Tree label에 id가 달려있는 것을 제외합니다.
tre$tip.label <- sapply(strsplit(tre$tip.label, "\t"), function(x)x[1])
phy_tree(data) <- tre


## 종속강목계 이름을 넣어줍니다.
rank_names(data)
tax_table(data) %>% colnames
colnames(tax_table(data))<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(data)

edc_dat <- edc_6yr_sel
save(data, edc_dat, file = "data/phyloseq_data_0510.RData")

twin_mo_id <- edc_dat$mo_id[which(duplicated(edc_dat$mo_id))]
data_twin_omit <- subset_samples(data, !(mo_id %in% twin_mo_id))

edc_dat_twin_omit <- edc_dat %>% filter(!mo_id %in% twin_mo_id)
save(data_twin_omit, edc_dat_twin_omit, file = "data/phyloseq_data_twin_omit_0602.RData")
