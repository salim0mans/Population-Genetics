---
title: "Allele Frequency Comparison between Tatar and Yakut Populations according to T2D associated Variants"
author: "Saleem Mansour"
date: "2023-10-24"
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
library(readxl)
```

```{r,echo=FALSE,include=FALSE}
af <- read_excel("final_AF.xlsx")
ac <- read_csv("final_AC.csv")
```

### Filter the rare variants 

Exclude Any SNP that -in a particular population- satisfies:

$$ \frac{Total_{ALT}\times Total_{pop}}{Total_{alleles}} < 5 $$  

```{r}

df1_filter <- ac %>% dplyr::filter(!filter1 == "Rare")
```


## Statistical Tests

For proportional categorical data, such as AFs, the common statistical tests for sufficient number of samples are Z-test and Chi-Square independence test.

Starting with Z-test

The equation is:  
$$ Z =\frac{p_1 - p_2}{\sqrt{p_{avg}\times(1-p_{avg})\times(\frac{1}{n_1}+\frac{1}{n_2})}} $$  
  
For Chi-Square:
Consider the table:

                                        |       | REF | ALT |       |
                                        |-------|-----|-----|-------|
                                        | Tatar | A   | B   | A+B   |
                                        | Yakut | C   | D   | C+D   |
                                        |       | A+C | B+D | Total |

We make estimated proportions assuming that independence qualifies that any slot is a joint prob (intersection) can be resulted from just multiplying the marginal probs, and since it's a proportion, it is divided by total. e.g., A' = ((A+C)*(A+B))/total and so on.

after making the estimates, we use the Chi-square test:

$$ X^2 = \sum \frac{(observed-expected)^2}{expected}  $$  
Since we are comparing the two populations for Each SNP alone, then the df = (r-1)(c-1) = 1, the test is also two-tailed.


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
```
### Performing tests
```{r}
## The total number of alleles is 2*sample_number, 88*2 for Tatars, and 84*2 for Yakuts
df_tested <- df1_filter %>% mutate(p_val_z = z_test_v(ALT_Tatar,ALT_Yakut,176,168))
df_tested <- df_tested %>% mutate(p_val_chi = chi_sqr_v(ALT_Tatar,REF_Tatar,ALT_Yakut,REF_Yakut))
df_tested <- df_tested  %>%  mutate(level_signif = ifelse(p_val_chi<0.001,paste0("***"),
                                                   ifelse(p_val_chi<0.01,paste0("**"),
                                                   ifelse(p_val_chi<0.05,paste0("*"),paste0("ns")))))
```

### Adjust multiple testing (BH FDR)

FDR is a necessity for multiple testing. We choose one of the most common algorithms for it, which is Benjamini-Hochberg.
  1-Order p-values from smallest to largest
  2-Add a column of Ranking (from 1 to N)
  3-Apply the algorithm on each row as follows:  
  
$$ adj\_p = P\_val_i\times\frac{N}{i} $$  


```{r}
df_adjust <- df_tested[order(df_tested$p_val_chi, decreasing = F),]
v1 <-  seq(1,length(df_adjust$SNP))
df_adjust$rank <- v1
df_adjust <- df_adjust %>% mutate(p_val_adj = p_val_chi*length(df_adjust$SNP)/rank)
df_adjust <- df_adjust  %>%  mutate(adj_signif = ifelse(p_val_adj<0.001,paste0("***"),
                                                   ifelse(p_val_adj<0.01,paste0("**"),
                                                   ifelse(p_val_adj<0.05,paste0("*"),paste0("ns")))))
df_final_adj <- df_adjust[-c(6,10)]
```



## Preparing the final list

```{r}
# Considering only significant variants:
## Adding the difference of MAF as a percentage
df_final_adj <- df_final_adj %>%  mutate(MAF_Tatar = (ALT_Tatar/176))
df_final_adj <- df_final_adj %>%  mutate(MAF_Yakut = (ALT_Yakut/168))

df_final <- df_final_adj %>% mutate( Dif = ifelse(!adj_signif == "ns",
                                                  round(abs(MAF_Tatar - MAF_Yakut)*100,5), paste0("--")))
## Adding the folds of change in MAF
df_final <- df_final %>% mutate( fc = ifelse(!adj_signif == "ns",
                                                 ifelse(MAF_Tatar>=MAF_Yakut,round(MAF_Tatar/MAF_Yakut,2),round(MAF_Yakut/MAF_Tatar,2)), paste0("--")))

## Adding the name of the population with the higher MAF (for ease of following)
df_final <- df_final %>% mutate("Higher MAF" = ifelse(!adj_signif == "ns", 
                                               ifelse(MAF_Tatar>MAF_Yakut,paste0("Tatar"),paste0("Yakut")),paste0("--")))
## Adding a value that take in count both  Difference% and Fold Change and highly correlated with significant p-values
df_final <- df_final %>% mutate(GDV = ifelse(!adj_signif == "ns",as.numeric(Dif)*as.numeric(Dif)*log10(as.numeric(fc)),paste0("--")))
```




### Saving the final dataframe as a separate report
```{r}
write.csv(df_final,"Final_AF_Comparison_Report.csv")
```



## plotting original MAF comparison

```{r, warning=FALSE}
df_long <- df_final %>%
  gather(key = "pop", value = "MAF", MAF_Tatar, MAF_Yakut)

df_signif <- df_tested %>% filter()

g1 <- ggplot(df_long, aes(x=SNP, y=MAF,fill=pop))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("MAF_pop.png",g1,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g2 <- ggplot(df_final, aes(x=SNP, y=MAF_Tatar,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#f2a90c", "yellow","#24820d","#aba9a9"))
ggsave("MAF_Tatar_signif.png",g2,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)

g3 <- ggplot(df_final, aes(x=SNP, y=MAF_Yakut,fill=adj_signif))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#f2a90c", "yellow","#24820d","#aba9a9"))
ggsave("MAF_Yakut_signif.png",g3,width = 12000, height = 2000, units = "px", type = "cairo", dpi = 300)
```
