
my_id <- "GSE82042"
gse <- getGEO(my_id)

length(gse)
gse <- gse[[1]]
gse

pData(gse)
fdata <- fData(gse)
exprs(gse)

pData(gse)$data_processing[1]
summary(exprs(gse))
head(exprs(gse))
prob_IDs <- rownames(exprs(gse))
sum(is.na(gse@featureData@data[["Symbol"]]))
head(gse@featureData@data[["Symbol"]])
identical(rownames(exprs(gse)), rownames(gse@featureData@data))


BI_GEO_Accession <- gse@phenoData@data[["geo_accession"]]
BI_Gender <- gse@phenoData@data[["characteristics_ch1.5"]]
BI_Expression_data <- exprs(gse)
BI_Sickness_Status <- gse@phenoData@data[["characteristics_ch1.4"]]


BiFullData <- data.frame(BI_GEO_Accession,BI_Gender,BI_Sickness_Status)
BI_Positive_female <- filter(BiFullData, BI_Gender == "gender: f",
                             BI_Sickness_Status == "diagnosisbp1: 2")
BI_Positive_male <- filter(BiFullData, BI_Gender == "gender: m",
                           BI_Sickness_Status == "diagnosisbp1: 2")
BI_Negative_female <- filter(BiFullData, BI_Gender == "gender: f",
                             BI_Sickness_Status == "diagnosisbp1: 0")
BI_Negative_male <- filter(BiFullData, BI_Gender == "gender: m",
                           BI_Sickness_Status == "diagnosisbp1: 0")
BI_Positive_General <-BI_Negative_female <- filter(BiFullData,
                                                   BI_Sickness_Status == "diagnosisbp1: 2")
BI_Negaitive_General <-BI_Negative_female <- filter(BiFullData,
                                                   BI_Sickness_Status == "diagnosisbp1: 0")

Expression_BPG <- exprs(gse)[, BI_Positive_General$BI_GEO_Accession]
Expression_BNG <- exprs(gse)[, BI_Negative_female$BI_GEO_Accession]
Expression_BPFemale <- exprs(gse)[, BI_Positive_female$BI_GEO_Accession]
Expression_BNFemale <- exprs(gse)[, BI_Negative_female$BI_GEO_Accession]
Expression_BPMale <- exprs(gse)[, BI_Positive_male$BI_GEO_Accession]
Expression_BNMale <- exprs(gse)[, BI_Negative_male$BI_GEO_Accession]
Gens <- gse@featureData@data[["Symbol"]]
genXID <- data.frame(Gens, prob_IDs)

BI_meta <- data.frame(SampleID = BI_GEO_Accession,
                      sex = BI_Gender)
rownames(BI_meta) <- BI_meta$SampleID

BI_dge <- DGEList(counts = BI_Expression_data,
                  genes = genXID,
                  samples = BI_meta)
#I don't need to use voom, my data is already log2 transformed
dge_v <- voom(BI_dge, plot = TRUE)
saveRDS(dge_v, file = "dge_v.rds")


BI_dge_1 <- DGEList(counts = BI_Expression_data,
                  genes = genXID,
                  samples = BI_meta,
                 group = BI_meta$sex,
                 annotation.columbs = Gens)
dge_v_1 <- voom(BI_dge_1, plot = TRUE)
saveRDS(dge_v_1, file = "dge_v_1.rds")

Diagnosis <-  BI_Sickness_Status
BI_meta_Diag <- cbind(BI_meta,Diagnosis)


BI_dge_2 <- DGEList(counts = BI_Expression_data,
                    genes = genXID,
                    samples = BI_meta_Diag,
                    group = BI_meta_Diag$Diagnosis)
dge_v_2 <- voom(BI_dge_2, plot = TRUE)
saveRDS(dge_v_2, file = "dge_v_2.rds")


Diag_Gender_combo <- paste(BI_Gender, BI_Sickness_Status)
BI_meta_combo <- cbind(BI_meta, Diag_Gender_combo)
BI_dge_3 <- DGEList(counts = BI_Expression_data,
                    genes = genXID,
                    samples = BI_meta_Diag,
                    group = BI_meta_combo$Diag_Gender_combo)
dge_v_3 <- voom(BI_dge_2, plot = TRUE)
saveRDS(dge_v_3, file = "dge_v_3.rds")
