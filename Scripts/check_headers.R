library(data.table)
ilcco_tsv <- "My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv"
finngen_gz <- "My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz"

cat("ILCCO header:\n")
print(names(fread(ilcco_tsv, nrows = 0)))

cat("FinnGen header:\n")
print(names(fread(finngen_gz, nrows = 0)))
