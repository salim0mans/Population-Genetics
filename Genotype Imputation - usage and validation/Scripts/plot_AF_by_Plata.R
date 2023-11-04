#!/bin/env Rscript --no-save

library(data.table) # For fast fread()
# Input variables
args <- commandArgs(TRUE)
infrq <- args[1] # Full path to chip data .frq file
indataset <- args[2] # Dataset name prefix
sub_dataset <- basename(indataset) # Subset the dataset name
ref_panel <- args[3] # Full path to reference panel .frq file
af_diff <- as.numeric(args[4]) # Max AF difference to panel
af_fc <- as.numeric(args[5]) # Max AF fold change to panel

# Read in the frequency files
chip <- fread(infrq, header = TRUE)
panel <- fread(ref_panel, header = TRUE)

# Take an intersection of the panel and chip data 
# based on SNP column (in format CHR_POS_REF_ALT)
isec <- merge(panel, chip, by = c("CHR", "SNP", "REF", "ALT"))
colnames(isec)[c(5,6)] <- c("AF_PANEL", "AF_CHIP")

isec$AF_diff <- abs(isec$AF_PANEL - isec$AF_CHIP)
isec$AF_ratio <- isec$AF_CHIP/isec$AF_PANEL
isec$AF_FC <- log2(isec$AF_ratio)

# Check if AFs are within the pp and fold change ranges
af_ok <- ((isec$AF_diff < af_diff) & 
         ((isec$AF_FC < af_fc) & (isec$AF_FC > -af_fc)))

exclude <- !af_ok

# Generate an exclusion list for variants not in the panel
nonpanel <- chip[!(chip$SNP) %in% (isec$SNP)]
# Generate AF list for discordant variants
af_discrepant <- isec[exclude]

# Save the plots as png
png(paste0(sub_dataset, "_AFs.png"), 
    width = 1200, height = 600)
    par(mfrow = c(1,2)) 
    par(cex.axis = 1.6, cex.lab = 1.5, cex.main = 1.6)

    # Plot first all and then excludable variants
    plot(isec$AF_PANEL, isec$AF_CHIP, col=1, pch=20, 
           main=paste0("Chip data AF vs. reference panel AF\n", 
           sub_dataset), 
           xlab="Panel AF", ylab="Chip AF")
    points(isec[exclude]$AF_PANEL, isec[exclude]$AF_CHIP, 
            col=2, pch=20)
    # Draw a legend
    legend("topleft", legend=c(
           paste0("Concordant AF, n = ", nrow(isec[!exclude])), 
           paste0("Discordant AF, n = ", nrow(isec[exclude])), 
           paste0("Non-panel, n = ", nrow(nonpanel))), 
           col=c("black", "red", "white"), pch=20, cex=1.2)

    # Add chip AF histogram to the same png
    hist(isec[!exclude]$AF_CHIP, breaks=100, 
            main=paste0("Chip AF for concordant variants\n", 
            sub_dataset), 
            xlab="Chip AF")
dev.off()

# Write out the exclusion lists (SNPID for bcftools)
# Highly discrepant AF values
write.table(isec[exclude]$SNP, 
        paste0(sub_dataset, "_exclude.txt"), 
        quote=F, row.names=F, col.names=F)
# Variants noit present in the reference panel
write.table(nonpanel$SNP, 
        paste0(sub_dataset, "_nonpanel_exclude.txt"), 
        quote=F, row.names=F, col.names=F)

# Write out the isec without exclusions 
write.table(isec[!exclude]$SNP, 
        paste0(sub_dataset, "_panel_isec.txt"), 
        quote=F, row.names=F, col.names=F)

# Write out the discrepant AFs in case needed
write.table(af_discrepant, 
        paste0(sub_dataset, "_discrepant_AF.txt"), 
        quote = F, row.names = F, sep = "\t")