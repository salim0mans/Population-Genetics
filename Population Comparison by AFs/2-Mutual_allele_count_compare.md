---
title: "Allele Frequency Comparison"
author: "Saleem Mansour"
date: "2023-11-12"
output:
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r, echo=FALSE, include=FALSE }
library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(Cairo)
library(DescTools)
library(grid)
library(gridExtra)
```

## Preparation  
```{r,echo=FALSE,include=FALSE}

ac <- read_csv("AC_3pops_all.csv")
# this information of sample count is obtained from the original vcf files by bcftools query -l <file> | wc -l 
pop1_len = 89*2
pop2_len = 88*2
pop3_len = 84*2

colnames(ac) <- c("SNP","Mean_chech","Mean_tat","Mean_yak")

df1 = ac %>% mutate(pop1_ALT = Mean_chech)
df1 = df1 %>% mutate(pop1_REF = pop1_len-Mean_chech)
df1 = df1 %>% mutate(pop2_ALT = Mean_tat)
df1 = df1 %>% mutate(pop2_REF = pop2_len-Mean_tat)
df1 = df1 %>% mutate(pop3_ALT = Mean_yak)
df1 = df1 %>% mutate(pop3_REF = pop3_len-Mean_yak)
```

### Filter the rare variants 

Exclude Any SNP that -in a particular population- satisfies:

$$ \frac{Total_{ALT}\times Total_{pop}}{Total_{alleles}} < 5 $$
```{r}
## Run this chunk only if want to filter rare variants out.
df1 <-  df1 %>%  mutate(filter1 = ifelse(  ((pop1_ALT*pop1_len/  (pop1_len+pop2_len)  )<5) | ((pop2_ALT*pop2_len/(pop2_len+pop2_len))<5) | ((pop3_ALT*pop3_len/(pop3_len+pop3_len))<5),paste0("Rare"),paste0("--")))
df1 <- df1 %>% dplyr::filter(!filter1 == "Rare")
```


```{r}
## Run this chunk only if want to filter Null alleles out.
df1 <-  df1 %>%  mutate(filter1 = ifelse(pop1_ALT==0 | pop1_REF==0 | pop2_ALT==0 | pop2_REF==0 | pop3_ALT==0 | pop3_REF==0 ,paste0("Null") ,paste0("--")))
df1 <- df1 %>% dplyr::filter(!filter1 == "Null")
```




## Statistical Tests

for multiple categorical group tests, we prefer Chi-square or G-test as an alternative test.

### Constructing the statistical test functions
```{r}
g_test <- function(x1a,x1r,x2a,x2r,x3a,x3r){
  obs <- matrix(c(x1a, x2a, x3a, x1r,x2r,x3r), nrow = 2, ncol = 3, byrow = T)
  g <- GTest(obs,correct = "none")
  return(g$p.value)
}

g_test_v <- Vectorize(g_test)


chi_sqr <- function(x1a,x1r,x2a,x2r,x3a,x3r){
  obs <- matrix(c(x1a, x2a, x3a, x1r,x2r,x3r), nrow = 2, ncol = 3, byrow = T)
  chi <- chisq.test(obs,correct = T)
  return(chi$p.value)
}

chi_sqr_v <- Vectorize(chi_sqr)

```


### Performing tests (choose one convenient test)

```{r}
## Chi-square test
df_tested <- df1 %>% mutate(p_val_chi = chi_sqr_v(pop1_ALT,pop1_REF,pop2_ALT,pop2_REF,pop3_ALT,pop3_REF))
df_tested <- df_tested  %>%  mutate(level_signif = ifelse(p_val_chi<0.001,paste0("***"),
                                                   ifelse(p_val_chi<0.01,paste0("**"),
                                                   ifelse(p_val_chi<0.05,paste0("*"),paste0("ns")))))
```


```{r}
## G-test test
df_tested <- df1 %>% mutate(p_val_g = g_test_v(pop1_ALT,pop1_REF,pop2_ALT,pop2_REF,pop3_ALT,pop3_REF))
df_tested <- df_tested  %>%  mutate(level_signif = ifelse(p_val_g<0.001,paste0("***"),
                                                   ifelse(p_val_g<0.01,paste0("**"),
                                                   ifelse(p_val_g<0.05,paste0("*"),paste0("ns")))))
```

### Adjust multiple testing (BH FDR)

FDR is a necessity for multiple testing. We choose one of the most common algorithms for it, which is Benjamini-Hochberg.
  1-Order p-values from smallest to largest
  2-Add a column of Ranking (from 1 to N)
  3-Apply the algorithm on each row as follows:
$$ adj\_p = P\_val_i\times\frac{N}{i} $$


```{r}
### Change the p_val_x (x: z or chi, fisher, b, or g) according to the chosen test (chi as default)

df_adjust <- df_tested[order(df_tested$p_val_chi, decreasing = F),]
v1 <-  seq(1,length(df_adjust$SNP))
df_adjust$rank <- v1
df_adjust <- df_adjust %>% mutate(p_val_adj = ifelse(p_val_chi*length(df_adjust$SNP)/rank>=1,1,p_val_chi*length(df_adjust$SNP)/rank))
df_adjust <- df_adjust  %>%  mutate(adj_signif = ifelse(p_val_adj<0.001,paste0("***"),
                                                   ifelse(p_val_adj<0.01,paste0("**"),
                                                   ifelse(p_val_adj<0.05,paste0("*"),paste0("--")))))

## Alternative Allele frequency

df_adjust <- df_adjust %>%  mutate(AF_Chechen = (pop1_ALT/pop1_len))
df_adjust <- df_adjust %>%  mutate(AF_Tatar = (pop2_ALT/pop2_len))
df_adjust <- df_adjust %>%  mutate(AF_Yakut = (pop3_ALT/pop3_len))

```



## plotting MAF comparison

```{r, warning=FALSE}
df_long <- df_adjust %>%
  gather(key = "case", value = "AF_pop", AF_Chechen, AF_Tatar, AF_Yakut)


g1 <- ggplot(df_long, aes(x=SNP, y=AF_pop,fill=case))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("AF_3pops.png",g1,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g2 <- ggplot(df_adjust, aes(x=SNP, y=AF_Chechen,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#aba9a9","#f2a90c", "yellow","#24820d"))
ggsave("AF_chech.png",g2,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g3 <- ggplot(df_adjust, aes(x=SNP, y=AF_Tatar,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#aba9a9","#f2a90c", "yellow","#24820d"))
ggsave("AF_tat.png",g3,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g4 <- ggplot(df_adjust, aes(x=SNP, y=AF_Yakut,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#aba9a9","#f2a90c", "yellow","#24820d"))
ggsave("AF_yak.png",g4,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

########################## only signif
df_signif <- df_adjust %>% filter(!adj_signif == "--")

df_long2 <- df_signif %>%
  gather(key = "case", value = "AF_pop", AF_Chechen, AF_Tatar, AF_Yakut)


g5 <- ggplot(df_long2, aes(x=SNP, y=AF_pop,fill=case))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("AF_3pops_signif-only.png",g5,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g6 <- ggplot(df_signif, aes(x=SNP, y=AF_Chechen,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#f2a90c", "yellow","#24820d"))
ggsave("AF_chech_signif.png",g6,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g7 <- ggplot(df_signif, aes(x=SNP, y=AF_Tatar,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#f2a90c", "yellow","#24820d"))
ggsave("AF_tat_signif.png",g7,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g8 <- ggplot(df_signif, aes(x=SNP, y=AF_Yakut,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#f2a90c", "yellow","#24820d"))
ggsave("AF_yak_signif.png",g8,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)




p_tog<- grid.arrange(g2,g3,g4, nrow = 3, ncol=1, top ="AF comparison between three Russian Populations")
p_tog_signif <-  grid.arrange(g6,g7,g8, nrow = 3, ncol=1, top ="AF comparison between three Russian Populations")
ggsave("AF_3pops_together.png",p_tog,width = 12000, height = 5000, units = "px", type = "cairo", dpi = 300)
ggsave("AF_3pops_together_signif-only.png",p_tog_signif,width = 12000, height = 5000, units = "px", type = "cairo", dpi = 300)

```

## Saving the final dataframe as a separate report
```{r}
write.csv(df_adjust,"Final_Mutual_AC_Comparison_Report.csv")
```

