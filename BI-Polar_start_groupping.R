
my_id <- "GSE82042"
#gse <- getGEO(my_id, returnType = 'ExpressionSet')

#gse <- gse[[1]]
#saveRDS(gse,"C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/gse")
gse <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/gse")
fdata <- fData(gse)
prob_IDs <- rownames(exprs(gse))
#pheno <- as.data.frame(gse@phenoData@data)
#saveRDS(pheno, "C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/pheno")
pheno <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/pheno")


#BI_GEO_Accession <- gse@phenoData@data[["geo_accession"]]
#saveRDS(BI_GEO_Accession, "C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_GEO_Accession")
BI_GEO_Accession <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_GEO_Accession")
#BI_Gender <- gse@phenoData@data[["characteristics_ch1.5"]]
#saveRDS(BI_Gender , "C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_Gender")
BI_Gender <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_Gender")
#BI_Expression_data <- exprs(gse)
#saveRDS(BI_Expression_data, "C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_Expression_data")
BI_Expression_data <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_Expression_data")
#BI_Sickness_Status <- gse@phenoData@data[["characteristics_ch1.4"]]
#saveRDS(BI_Sickness_Status, "C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_Sickness_Status")
BI_Sickness_Status <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_Sickness_Status")


BiFullData <- data.frame(BI_GEO_Accession,BI_Gender,BI_Sickness_Status)
BI_Positive_female <- filter(BiFullData, BI_Gender == "gender: f",
                             BI_Sickness_Status == "diagnosisbp1: 2")
BI_Positive_male <- filter(BiFullData, BI_Gender == "gender: m",
                           BI_Sickness_Status == "diagnosisbp1: 2")
BI_Negative_female <- filter(BiFullData, BI_Gender == "gender: f",
                             BI_Sickness_Status == "diagnosisbp1: 0")
BI_Negative_male <- filter(BiFullData, BI_Gender == "gender: m",
                           BI_Sickness_Status == "diagnosisbp1: 0")
BI_Positive_General <- filter(BiFullData,
                              BI_Sickness_Status == "diagnosisbp1: 2")
BI_Negaitive_General <- filter(BiFullData,
                               BI_Sickness_Status == "diagnosisbp1: 0")

Expression_BPG <- exprs(gse)[, BI_Positive_General$BI_GEO_Accession]
Expression_BNG <- exprs(gse)[, BI_Negative_female$BI_GEO_Accession]
Expression_BPFemale <- exprs(gse)[, BI_Positive_female$BI_GEO_Accession]
Expression_BNFemale <- exprs(gse)[, BI_Negative_female$BI_GEO_Accession]
Expression_BPMale <- exprs(gse)[, BI_Positive_male$BI_GEO_Accession]
Expression_BNMale <- exprs(gse)[, BI_Negative_male$BI_GEO_Accession]


#Genes <- gse@featureData@data[["Symbol"]]
#saveRDS(Genes, "C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/Genes")
Genes <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/Genes")

#genXID <- data.frame(Genes, prob_IDs)
#saveRDS(genXID, "C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/genXID")
genXID <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/genXID")

BI_meta <- data.frame(SampleID = BI_GEO_Accession,
                      sex = BI_Gender)
rownames(BI_meta) <- BI_meta$SampleID
saveRDS(BI_meta,"C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_meta")
BI_meta <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_meta")

BI_dge <- DGEList(counts = BI_Expression_data,
                  genes = genXID,
                  samples = BI_meta)



BI_dge_1 <- DGEList(counts = BI_Expression_data,
                  genes = genXID,
                  samples = BI_meta,
                 group = BI_meta$sex,
                 annotation.columbs = Genes)

Diagnosis <-  BI_Sickness_Status
BI_meta_Diag <- cbind(BI_meta,Diagnosis)

BI_dge_2 <- DGEList(counts = BI_Expression_data,
                    genes = genXID,
                    samples = BI_meta_Diag,
                    group = BI_meta_Diag$Diagnosis)


Diag_Gender_combo <- paste(BI_Gender, BI_Sickness_Status)
BI_meta_combo <- cbind(BI_meta, Diag_Gender_combo)
BI_dge_3 <- DGEList(counts = BI_Expression_data,
                    genes = genXID,
                    samples = BI_meta_Diag,
                    group = BI_meta_combo$Diag_Gender_combo)

saveRDS(BI_dge_3 , "C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_dge_3 ")
BI_dge_3  <- readRDS("C:/Users/avive/Documents/Thesis/R code/Thesis/Bachlor/Working/Reusable Files/BI_dge_3 ")
#-------------------------------------------------------------------------------------------------------#
