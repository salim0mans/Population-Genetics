#### QCed Summary statistic by plink (--maf 0.01) 
summary <- read.table("base_ukb_gwas_qc.tsv", header = T)
summary <- summary %>% select(SNP,BETA)

#### ASsuming the data is case controls who from multiple ethnic groups, identified by ID and annotated in Group column in both Controls and T2D

controls <- read_table("Controls_hwe_filtered_SNP.tsv")
t2d <- read_table("T2D_hwe_filtered_SNP.tsv")

SNP_list <- read.csv("Final_SNP_list_T2D.csv")

controls <- merge(SNP_list,controls,by="SNP")
controls <- controls %>% mutate(case = paste0("Control"))

t2d <- merge(SNP_list,t2d,by="SNP")
t2d <- t2d %>% mutate(case = paste0("T2D"))

all_samples <- merge(controls,t2d, by = "SNP")

p_all_samples <- all_samples %>% filter(BETA >= 0)
n_all_samples <- all_samples %>% filter(BETA < 0)


p_weighted <- p_all_samples %>% rowwise() %>% mutate(across(-c(AF,ID,Group,Genes,case,SNP), ~ .* BETA))
n_weighted <- n_all_samples %>% rowwise() %>% mutate(across(-c(AF,ID,Group,Genes,case,SNP), ~ .* BETA))

p_PRS <- colSums(p_weighted[, !colnames(p_weighted) %in% c("AF","ID","Genes","Group","case","SNP","BETA")])
n_PRS <- colSums(n_weighted[, !colnames(n_weighted) %in% c("AF","ID","Genes","Group","case","SNP","BETA")])

p_PRS <- data.frame(p_PRS)
n_PRS <- data.frame(n_PRS)

p_final_PRS <- p_PRS %>% mutate(partial = paste0("Positive"))
n_final_PRS <- n_PRS %>% mutate(partial = paste0("Negative"))

final_PRS <- merge(p_final_PRS,n_final_PRS,by="ID")
final_PRS <- final_PRS %>% mutate(total_PRS = p_PRS + n_PRS)
### from this point, each sample has annotations of case, Group, positive PRS, and negative PRS

write.csv("Final_PRS.csv",final_PRS,row.names = F)
