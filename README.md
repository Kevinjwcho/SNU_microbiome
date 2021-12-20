# Microbiome Tutorial

## 필요한 Library

Library는 _phyloseq_ 가 가장 중요합니다. 이는 주석된 명령문으로 설치할 수 있습니다. 나머지 문구는 _install.packages("")_ 로 설치하시면 됩니다.

```{r library, message=FALSE}
# BiocManager::install("phyloseq")
library(phyloseq)
library(dplyr)
library(data.table)
library(stringr)
```


## .biome 데이터 불러오기
phyloseq class 데이터는 총 4개로 이루어져 있습니다.

1. otu_table : OTU count 데이터
2. sample_data : sample에 대한 데이터
3. tax_table : Taxonomy 데이터
4. phy_tree : Phylogentic Tree 데이터
이 중에서 1번 부터 3번 까지는 otu_table.blast_NCBI_16S.biome 에 저장 되어 있습니다. 
```{r, results = 'markup'}
data<- import_biom(BIOMfilename= "F:/work/Microbiome_EDC/HN00117382_Reports/HN00117382_OTU/OTU_Results/biom/otu_table.blast_NCBI_16S.biom")
data
```

## EDC 데이터 불러오기
단, sample_data와 같은 경우 오직 성별만이 들어갔기 때문에, 위 데이터를 edc 데이터로 바꿔주어야 합니다. 다음은 EDC 데이터를 불러서 sample_data에 넣는 방법입니다.

```{r, results = 'markup'}
edc_raw <- fread("F:/work/Microbiome_EDC/edc_raw.csv", header = T, fill = T)
## 분석에 필요한 변수를 정의해줍니다. ---------------------------
## 이 작업은 값들 중에 숫자가 아닌 Not Detected와 같이 character 값이 있어서 진행하는 것입니다. 
## 만약 다 숫자(numeric) 값만 존재한다면 나이변수 만들기 전까지 skip 하셔도 됩니다.
mother_var <- "BPA_M_creatinine MEHHP_M_creatinine MEOHP_M_creatinine MnBP_M_creatinine MECCP_M_cre MBzP_M_cre_2 M3_PBA_M_creatinine Pb_M Hg_M Cd_M Mn_M"
mother_var1 <- unlist(str_split(mother_var, " "))

child_var <- "BPA_creatinine BPS_cre BPF_cre MEHHP_creatinine MEOHP_creatinine MnBP_creatinine MECCP_creatinine MBzP_creatinine_2 b3_PBA_creatinine Pb_C Hg_C Cd_C Mn_C"
child_var1 <- unlist(str_split(child_var, " "))

child_outcome <- "b_wisc_1 scqtotal ars_total1 ars_total2 ars_total"
child_outcome1 <- unlist(str_split(child_outcome, " "))

select_coln <- which(colnames(edc_raw) %in% c(mother_var1, child_var1, child_outcome1))
# select_coln %>% length()
## 위에 있는 변수들은  "검체없음" 같은 요상한 값들이 많아서 다 공백으로 처리하였습니다.
# warnings 은 무시하셔도 됩니다. 그냥 다 숫자 아닌것들은 결측처리 했다고 나오는 경고문입니다.
edc_numeric <- mutate_all(edc_raw[,..select_coln], function(x) as.numeric(as.factor(as.character(x))))
str(edc_numeric)

## 바꾼 수치형 변수를 기존의 key 변수와 합치는 작업입니다.----------
edc_final <- cbind(edc_raw[,-..select_coln], edc_numeric)
# which(colnames(edc_final) %in% c(mother_var1, child_var1, child_outcome1)) %>% length()

## 변수 잘 들어갔는지 확인하는 코드입니다.-------
dim(edc_final);dim(edc_raw)
```

데이터를 1차로 정제한 다음에 나이에 따라서 데이터를 나누어 주시는 게 연구하기 편합니다. 다음은 나이로 데이터를 4등분 하는 코드입니다.

```{r, results = 'markup'}
edc_final$old <- apply(edc_final$age %>% as.data.frame() ,1,
                       function(x)ifelse(x<30, 2, ifelse(x<57, 4, ifelse(x<80, 6, 8))))

edc_final$new_id <- sapply(strsplit(edc_final$new_id, "-"), function(x) paste(x, collapse = "."))

edc_2yr <- edc_final %>% subset(old == 2);edc_2yr %>% nrow()
edc_4yr <- edc_final %>% subset(old == 4);edc_4yr %>% nrow()
edc_6yr <- edc_final %>% subset(old == 6);edc_6yr %>% nrow()
edc_8yr <- edc_final %>% subset(old == 8);edc_8yr %>% nrow()

```
위 코드를 이용해서 다음과 같이 나누고 각각의 sample 개수를 확인할 수 있습니다.

그리고 sample data에 EDC 6세 data를 id에 맞춰서 join(merge) 한다음에 넣는 과정입니다.
```{r}
## sample 데이터의 id가 rownames로 되어 있기 때문에 이를 new_id라는 변수로 새로 만들어서 저장합니다.
sam <- sample_data(data) %>% as.matrix %>% as.data.frame
sam$new_id<- rownames(sam)

## 합친 다음 다시 rownames로 new_id를 넣어주고 sample_data에 저장합니다.
sam_merge <- merge(sam, edc_6yr, by="new_id", all.x = T)
str(sam_merge)
rownames(sam_merge) <- rownames(sam)
sam_merge <- sam_merge %>% select(-Group, -new_id)
sample_data(data) <- sam_merge
```

위까지 한 경우 data는 총 3개의 데이터가 저장된 상태입니다. 보시면 sample_data의 변수 개수가 달라진 것을 알 수 있습니다.
```{r echo=FALSE, results='markup'}
data
```


## Phylogenetic Tree 생성방법
각각의 Denovo에 대한 정보가 otus_rep_fasta에 저장되어있습니다. 이를 통해 phylogenetic tree를 만들 수 있습니다.

1. adegenet package : 데이터 불러오는 _fasta2DNAbin_ 함수가 있는 패키지
2. ape : Phylogentic Tree를 만드는 여러 함수가 있는 패키지

이 두 가지 패키지를 설치해주시면 됩니다.

```{r, message=F}
# install.packages("ape")
# install.packages("adegenet")
library(adegenet)
library(ape)
```

그리고 다음과 같이 fasta 데이터를 불러와 주시면 됩니다. 시간이 좀 걸릴 수 있습니다.
```{r eval=FALSE}
dna <-fasta2DNAbin("F:/work/Microbiome_EDC/HN00117382_Reports/otus_rep.fasta")
class(dna)
dna
```

Tree를 생성하는 방법은 여러가지가 존재합니다. 이는 _?dist.dna()_라고 치시면 여러 옵션들을 확인할 수 있습니다. 단, distance matrix가 만들어지고 나서 그 안에 NULL(NA) 즉 빈 값이 존재하면, _njs()_ 함수를 통해서 트리를 생성할 수 있고, 빈 칸이 없으면 _nj()_ 함수를 쓰셔도 됩니다. 

**시간이 오래 걸리기 때문에 한 번 만드시고 저장하는 것을 추천드립니다.**
```{r, eval = F}
## Raw 방법은 NULL 값이 없습니다. -> nj() 함수를 쓸수 있습니다.
D <- dist.dna(dna, model = "raw")
length(D)
tre <- nj(D)

## TN93 방법은 NULL 값이 많습니다. -> njs() 함수를 써야합니다. 
## 그리고 오래걸립니다.
D <- dist.dna(dna, model = "TN93")
length(D)
tre <- njs(D)

## tree를 한번 만든다음 저장해서 필요할 때마다 불러오는 것이 좋습니다.
save(tre, file = "F:/work/Microbiome_EDC/RData/TN93_tre.RData")
```

저장된 tree를 불러오는 코드입니다.
```{r}
load("F:/work/Microbiome_EDC/RData/TN93_tre.RData")
```

그 다음 **data**에 tree를 넣는 방법입니다.
```{r, results='markup'}
## Tree label에 id가 달려있는 것을 제외합니다.
tre$tip.label <- sapply(strsplit(tre$tip.label, "\t"), function(x)x[1])
phy_tree(data) <- tre
data
```

마지막으로 각각의 데이터를 따로 불러오고 싶으면, 다음과 같은 명령문으로 불러올 수 있습니다.
```{r, results = 'markup'}
otu <- otu_table(data)
sample_dat <- sample_data(data)
tax <- tax_table(data)
phy_tre <- phy_tree(data)
```

## 데이터 손질하기
*rank_names()*라는 함수를 통해 taxonomy의 ranking level을 확인할 수 있습니다. taxonomy 데이터를 보면 아직 Phylogenetic Ranking이 제대로 안 되어 있는 것을 확인할 수 있습니다. 혹은 Taxonomy data만 따로 떼서 colnames를 확인하셔도 동일합니다.
```{r, results = 'markup'}
rank_names(data)
tax_table(data) %>% colnames
```
이는 다음과 같이 label을 달아줄 수 있습니다.
```{r, results='markup'}
colnames(tax_table(data))<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(data)
```

## 그림 그리기
### Bar graph 그리기

먼저 각각 샘플의 abundance를 보기위한 bar graph를 그리기 위해서는 *otu_table*이 필요합니다. 그리고 sample별 총 count 개수가 다르기 때문에, bar graph를 통해 otu의 비율을 보려고 한다면, 상대도수로 맞춰주는 작업이 필요합니다. otu table의 read count를 각 sample 별 상대도수로 바꾸려면 다음과 같은 코드로 할 수 있습니다.
```{r, results = 'markup'}
rel_sum <- function(x) round((x/sum(x)))
rel_data = transform_sample_counts(data, rel_sum)

## 상대도수로 바꿀 경우 대부분이 0이기 때문에, 그림이 이쁘게 나오지 않습니다. 이를 Sample 별 총 count의 median 값으로 scaling 해주는 것이 일종의 표준화 역할을 하면서 그림이 이뻐집니다. 그냥 10000이나 다른 큰 숫자를 total 대신에 쓰셔도 됩니다.
total = median(sample_sums(data))
standf <- function(x, t=total) round(t*(x/sum(x)))
stand_data <- transform_sample_counts(data, standf)
## 테이블 일부만 보여줍니다.
otu_table(data)[1:10, 1:5]
otu_table(rel_data)[1:10, 1:5]
```

생각보다 그림의 용량이 커서 그리는데 오래걸립니다. 그러므로 pdf로 저장할 경로를 미리 지정한다음 pdf 파일로 저장하게하는 것이 훨씬 빠릅니다. 이 그림을 그리기 위해서는 *ggplot2*라는 패키지가 필요합니다.
```{r, eval = F}
library(ggplot2)
plot_path <- "F:/work/Microbiome_EDC/plot"
pdf(file = paste0(plot_path,"/barplot_Phylum.pdf"))
plot_bar(stand_data, fill = "Phylum")+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position = "stack")+
  theme(axis.text.x = element_blank())
dev.off()
```
이를 그리면 다음과 같습니다.
```{r, echo = F, message = F}
library(ggplot2)
plot_bar(stand_data, fill = "Phylum")+
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position = "stack")+
  theme(axis.text.x = element_blank())
```

### heat map 그리기 
히트맵은 Taxonomy가 많으면 보기 불편하기 때문에, 필터링하는 것이 좋습니다. 다음은 샘플 별 총 read 수의 중앙값으로 10%보다 위에 있는 taxonomy만 발췌하는 코드입니다.
```{r, results = 'markup'}
stand_data_abund <- filter_taxa(stand_data, function(x) sum(x>total*0.1)>0, TRUE)
stand_data_abund
```
heatmap도 pdf로 저장하는 것이 빠릅니다.
```{r, eval = F}
pdf(file = paste0(plot_path, "/heatmap.pdf"))
plot_heatmap(stand_data_abund, method = "NMDS", distance = "bray")+theme(axis.text.x = element_blank())
dev.off()
```
이를 그림으로 그리면 다음과 같습니다.
```{r, echo = F}
plot_heatmap(stand_data_abund, method = "NMDS", distance = "bray")+theme(axis.text.x = element_blank())
```

### alpha diversity 그리기

alpha diversity 값을 원한다면 *estimate_richness()* 함수를 이용하면 chao1, shannon, simpson, fisher 등 다양한 alpha diversity를 계산할 수 있습니다.
```{r, results = 'markup'}
estimate_richness(data, measures = "Chao1") %>% head
estimate_richness(data, measures = "Shannon") %>% head
estimate_richness(data, measures = "Simpson") %>% head
```

이를 그림으로 그리면 다음과 같습니다.
```{r, eval = F}
pdf(file = paste0(plot_path, "/alpha_diversity.pdf"))
plot_richness(data, measures = c("Chao1", "Shannon"))+theme(axis.text.x = element_blank())
dev.off()
```
```{r, message = F}
plot_richness(data, measures = c("Chao1", "Shannon"))+theme(axis.text.x = element_blank())
```

### Beta diversity 구하기

beta diversity 값을 원한다면 *distance()* 함수를 이용하면 chao1, shannon, simpson, fisher 등 다양한 alpha diversity를 계산할 수 있습니다.
```{r, results = 'markup'}
distance(data, method = "unifrac") %>% head
distance(data, method = "wunifrac") %>% head
distance(data, method = "bray") %>% head
distance(data, method = "dpcoa") %>% head

```

이름 그림으로 그릴 때는 *ordinate* 라는 함수를 써서 ordination을 실시해야합니다. 크게 대표적인 ordination 방법은 *NMDS*와 *PCoA*가 있습니다.
```{r, results = 'markup'}
data.ord1 <- ordinate(data, "NMDS", "bray")
data.ord2 <- ordinate(data, "NMDS", "wunifrac")
data.ord3 <- ordinate(data, "PCoA", "bray")
data.ord4 <- ordinate(data, "PCoA", "wunifrac")
```

이를 그림으로 그리면 다음과 같습니다. *"baby_gen"* 대신에 다른 edc 변수를 넣으면 그에 맞게 그려줍니다.
```{r, results = 'markeup'}
plot_ordination(data, data.ord1, type="samples", title = "Samples", color = "baby_gen")+
  geom_point(size=2)
plot_ordination(data, data.ord2, type="samples", title = "Samples", color = "baby_gen")+
  geom_point(size=2)
plot_ordination(data, data.ord3, type="samples", title = "Samples", color = "baby_gen")+
  geom_point(size=2)
plot_ordination(data, data.ord4, type="samples", title = "Samples", color = "baby_gen")+
  geom_point(size=2)
```
