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
library(Barnard)
```

## Preparaion
```{r,echo=FALSE,include=FALSE}

ac <- read_csv("AC_case-controls.csv")
# this information of sample count is obtained from the original vcf files by bcftools query -l <file> | wc -l 
ctl_len = X    # replace X with number of samples in the group * 2
t2d_len = Y    # replace Y with number of samples in the group * 2
df1 = ac %>% mutate(ctl_ALT = Mean_x)
df1 = df1 %>% mutate(ctl_REF = ctl_len-Mean_x)
df1 = df1 %>% mutate(t2d_ALT = Mean_y)
df1 = df1 %>% mutate(t2d_REF = t2d_len-Mean_y)
```

### Filter the rare variants 

Exclude Any SNP that -in a particular population- satisfies:

$$ \frac{Total_{ALT}\times Total_{pop}}{Total_{alleles}} < 5 $$
```{r}
####### large samples
## Run this chunk only if want to filter rare variants out.
df1 <-  df1 %>%  mutate(filter1 = ifelse(((ctl_ALT*ctl_len/(ctl_len+t2d_len))<5) | ((t2d_ALT*t2d_len/(t2d_len+t2d_len))<5),paste0("Rare"),paste0("--")))
df1 <- df1 %>% dplyr::filter(!filter1 == "Rare")
```


```{r}
####### small samples
## Run this chunk only if want to filter Null alleles out.
df1 <-  df1 %>%  mutate(filter1 = ifelse(((ctl_ALT*ctl_len/(ctl_len+t2d_len))==0) | ((t2d_ALT*t2d_len/(t2d_len+t2d_len))<0),paste0("Null"),paste0("--")))
df1 <- df1 %>% dplyr::filter(!filter1 == "Null")
```



## Statistical Tests  

For proportional categorical data, such as AFs, the common statistical tests for sufficient number of samples are Z-test and Chi-Square independence test. G-test is a robust alternative to Chi-square.   
Fisher Exact test is considered for small samples. However, Barnard's test is a more powerful alternative.  

### Constructing the statistical test functions
```{r}
z_test <- function(x1,x2,n1,n2){
  z <- prop.test(x = c(x1, x2), n = c(n1, n2), correct = T, alternative = "two.sided")
  return(z$p.value)
}

z_test_v <- Vectorize(z_test)


chi_sqr <- function(x1a,x1r,x2a,x2r){
  obs <- matrix(c(x1a, x1r, x2a, x2r), nrow = 2, ncol = 2, byrow = T)
  chi <- chisq.test(obs,correct = T)
  return(chi$p.value)
}

chi_sqr_v <- Vectorize(chi_sqr)


fisher_exact <- function(x1a,x1r,x2a,x2r){
  obs <- matrix(c(x1a, x1r, x2a, x2r), nrow = 2, ncol = 2, byrow = T)
  f <- fisher.test(obs)
  return(f$p.value)
}

fisher_exact_v <- Vectorize(fisher_exact)


g_test <- function(x1a,x1r,x2a,x2r){
  obs <- matrix(c(x1a, x1r, x2a, x2r), nrow = 2, ncol = 3, byrow = T)
  g <- GTest(obs,correct = "none")
  return(g$p.value)
}

g_test_v <- Vectorize(g_test)


barnard <- function(x1a,x1r,x2a,x2r){
  obs <- matrix(c(x1a, x1r, x2a, x2r), nrow = 2, ncol = 2, byrow = T)
  barn <- BarnardTest(obs, alternative = 'two.sided', method = "z-unpooled")
  return(barn[["p.value"]])
}

barnard_v <- Vectorize(barnard)

```

### Performing tests (choose one convenient test)

```{r}
## Z-test
df_tested <- df1 %>% mutate(p_val_z = z_test_v(ctl_ALT,t2d_ALT,ctl_len,t2d_len))
df_tested <- df_tested  %>%  mutate(level_signif = ifelse(p_val_z<0.001,paste0("***"),
                                                   ifelse(p_val_z<0.01,paste0("**"),
                                                   ifelse(p_val_z<0.05,paste0("*"),paste0("ns")))))
```

```{r}
## Chi-square test
df_tested <- df1 %>% mutate(p_val_chi = chi_sqr_v(ctl_ALT,ctl_REF,t2d_ALT,t2d_REF))
df_tested <- df_tested  %>%  mutate(level_signif = ifelse(p_val_chi<0.001,paste0("***"),
                                                   ifelse(p_val_chi<0.01,paste0("**"),
                                                   ifelse(p_val_chi<0.05,paste0("*"),paste0("ns")))))
```

```{r}
## Fisher's exact test
df_tested <- df1 %>% mutate(p_val_fisher = fisher_exact_v(ctl_ALT,ctl_REF,t2d_ALT,t2d_REF))
df_tested <- df_tested  %>%  mutate(level_signif = ifelse(p_val_fisher<0.001,paste0("***"),
                                                   ifelse(p_val_fisher<0.01,paste0("**"),
                                                   ifelse(p_val_fisher<0.05,paste0("*"),paste0("ns")))))
```

```{r}
## G-test test
df_tested <- df1 %>% mutate(p_val_g = g_test_v(ctl_ALT,ctl_REF,t2d_ALT,t2d_REF))
df_tested <- df_tested  %>%  mutate(level_signif = ifelse(p_val_g<0.001,paste0("***"),
                                                   ifelse(p_val_g<0.01,paste0("**"),
                                                   ifelse(p_val_g<0.05,paste0("*"),paste0("ns")))))
```


```{r}
## Barnard's  test
df_tested <- df1 %>% mutate(p_val_b = barnard_v(ctl_ALT,ctl_REF,t2d_ALT,t2d_REF))
df_tested <- df_tested  %>%  mutate(level_signif = ifelse(p_val_b<0.001,paste0("***"),
                                                   ifelse(p_val_b<0.01,paste0("**"),
                                                   ifelse(p_val_b<0.05,paste0("*"),paste0("ns")))))
```


### Adjust multiple testing (BH FDR)

FDR is a necessity for multiple testing. We choose one of the most common algorithms for it, which is Benjamini-Hochberg.
  1-Order p-values from smallest to largest
  2-Add a column of Ranking (from 1 to N)
  3-Apply the algorithm on each row as follows:
$$ adj\_p = P\_val_i\times\frac{N}{i} $$


```{r}

###### Large samples
### Change the p_val_x (x: z or chi, fisher, g, or b) according to the chosen test (Chi as default)

df_adjust <- df_tested[order(df_tested$p_val_chi, decreasing = F),]
v1 <-  seq(1,length(df_adjust$SNP))
df_adjust$rank <- v1
df_adjust <- df_adjust %>% mutate(p_val_adj = ifelse(p_val_chi*length(df_adjust$SNP)/rank>=1,1,p_val_chi*length(df_adjust$SNP)/rank))
df_adjust <- df_adjust  %>%  mutate(adj_signif = ifelse(p_val_adj<0.001,paste0("***"),
                                                   ifelse(p_val_adj<0.01,paste0("**"),
                                                   ifelse(p_val_adj<0.05,paste0("*"),paste0("--")))))

df_adjust <- df_adjust %>%  mutate(cohort_ctl = (ctl_ALT/ctl_len))
df_adjust <- df_adjust %>%  mutate(cohort_t2d = (t2d_ALT/t2d_len))
```


```{r}
###### Small samples
### Change the p_val_x (x: z or chi or fisher, g, or b) according to the chosen test (barnard as default)

df_adjust <- df_tested[order(df_tested$p_val_b, decreasing = F),]
v1 <-  seq(1,length(df_adjust$SNP))
df_adjust$rank <- v1
df_adjust <- df_adjust %>% mutate(p_val_adj = ifelse(p_val_b*length(df_adjust$SNP)/rank>=1,1,p_val_b*length(df_adjust$SNP)/rank))
df_adjust <- df_adjust  %>%  mutate(adj_signif = ifelse(p_val_adj<0.001,paste0("***"),
                                                   ifelse(p_val_adj<0.01,paste0("**"),
                                                   ifelse(p_val_adj<0.05,paste0("*"),paste0("--")))))

df_adjust <- df_adjust %>%  mutate(cohort_ctl = (ctl_ALT/ctl_len))
df_adjust <- df_adjust %>%  mutate(cohort_t2d = (t2d_ALT/t2d_len))
```


## plotting AF comparison
```{r, warning=FALSE}
df_long <- df_adjust %>%
  gather(key = "case", value = "MAF", cohort_ctl, cohort_t2d)

g1 <- ggplot(df_long, aes(x=SNP, y=MAF,fill=case))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("cohort_case-control_snps.png",g1,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g2 <- ggplot(df_adjust, aes(x=SNP, y=cohort_ctl,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#aba9a9","#fa6b7e", "#fa6bdd","#726bfa"))
ggsave("cohort_ctl_snps.png",g2,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g3 <- ggplot(df_adjust, aes(x=SNP, y=cohort_t2d,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#aba9a9","#fa6b7e", "#fa6bdd","#726bfa"))
ggsave("cohort_t2d_snps.png",g3,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)
```

## Saving the final dataframe as a separate report
```{r}
write.csv(df_final,"cohort_case-control.csv")
```


