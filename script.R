# Author: Ehtisham Khokhar
# Email: ehtishamshakeel@gmail.com
# Affiliation: New Mexico State University

# set working directory
setwd("C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Fruit Morphology Projects/GWAS/amova")

# install and load required packages
install.packages(c("adegenet", "poppr", "ape")) 
library(adegenet)
library(poppr)
library(ape)

# import SNP data
snp_data <- import2genind("GWAS128.str",
                          onerowperind = FALSE,
                          n.ind = 128,
                          n.loc = 40709,
                          col.lab = 1,
                          ask = FALSE)

# extract individual IDs
snp_ids <- indNames(snp_data)
print(snp_ids)
cat("Number of individuals in snp_data:", length(snp_ids), "\n")

# import hierarchy file
hierarchy <- read.csv("population.csv", header = TRUE, sep=",")
hierarchy_ids <- hierarchy[[1]]
cat("Number of individuals in hierarchy:", length(hierarchy_ids), "\n")
print(hierarchy_ids)

# check for ID consistency
snp_ids <- as.character(snp_ids)
hierarchy_ids <- as.character(hierarchy_ids)

if (length(snp_ids) != length(hierarchy_ids)) {
  stop("Error: The number of individuals in snp_data (", length(snp_ids), 
       ") does not match the number of rows in hierarchy (", length(hierarchy_ids), ").")
}

mismatch <- !all(hierarchy_ids %in% snp_ids) | !all(snp_ids %in% hierarchy_ids)
if (mismatch) {
  missing_in_snp <- hierarchy_ids[!hierarchy_ids %in% snp_ids]
  missing_in_hierarchy <- snp_ids[!snp_ids %in% hierarchy_ids]

  if (length(missing_in_snp) > 0) {
    cat("Error: The following IDs in hierarchy are not in snp_data:\n")
    print(missing_in_snp)
  }
  if (length(missing_in_hierarchy) > 0) {
    cat("Error: The following IDs in snp_data are not in hierarchy:\n")
    print(missing_in_hierarchy)
  }
  stop("ID mismatch detected. Please fix the IDs before proceeding.")
} else {
  cat("All individual IDs match between snp_data and hierarchy.\n")
}

# convert to genclone object
snp_clone <- as.genclone(snp_data)

# assign population structure
strata(snp_clone) <- hierarchy
snp_strat <- snp_clone
cat("Strata successfully set.\n")

# check hierarchy column names
print(colnames(hierarchy))

# run AMOVA
amova_res <- poppr.amova(snp_strat, ~ pop, cutoff = 0.3)
print(amova_res)

# run permutation test
amova_test <- randtest(amova_res) 
plot(amova_test)
png(filename = "amova_significance_plot.png", width = 800, height = 600)
plot(amova_test)
dev.off()

# save results
amova_df <- as.data.frame(amova_res$results)
write.csv(amova_df, file = "amova_results.csv", row.names = TRUE)

amova_covar_components <- as.data.frame(amova_res$componentsofcovariance)
write.csv(amova_covar_components, file = "amova_covar_components.csv", row.names = TRUE)

amova_statphi <- as.data.frame(amova_res$statphi)
write.csv(amova_statphi, file = "amova_statphi.csv", row.names = TRUE)

