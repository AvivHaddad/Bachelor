try1 <- lmFit(exprs(gse))
pheno <- data.frame(
  diagnosis = gse@phenoData@data[["diagnosisbp1:ch1"]],
  sex = gse@phenoData@data[["gender:ch1"]])
design0 <- model.matrix(~ diagnosis, data = pheno)

#lmFit fits a linear model for each gene.
#eBayes improves variance estimates across genes, 
#giving more stable p‑values (this is limma’s “empirical Bayes” step).

fit0    <- lmFit(exprs(gse), design0)
fit0    <- eBayes(fit0)
design1 <- model.matrix(~ diagnosis + sex, data = pheno)
fit1    <- lmFit(exprs(gse), design1)
fit1    <- eBayes(fit1)
Table_fit0_Diagnosis <- topTable(fit0, coef = 2)
Table_fit1_Male <- topTable(fit1, coef = "sexm")

#For each gene, get statistics for the effect 
#“diagnosis = 2 versus the baseline level (diagnosis = 0)” in the two models.

tt0 <- topTable(fit0, coef = "diagnosis2", number = Inf)
tt1 <- topTable(fit1, coef = "sexm", number = Inf)
compare0 <- cor(tt0$logFC, tt1$logFC)
compare1 <- cor(tt0$P.Value, tt1$P.Value)
plot(tt0$logFC, tt1$logFC)
PCA <- prcomp(t(exprs(gse)))
diag_color <- ifelse(pheno$diagnosis == "2", "green", "red")
plot(PCA$x[, 1], PCA$x[, 2], col = diag_color)
color <- ifelse(pheno$sex == "m", "blue", "red")
plot(PCA$x[, 1], PCA$x[, 2], col = color )



shape_column <- "sex"
color_column <-"sex"
label <- T
label_size <- 4
plot_save_name <- "PCR_Plod.pdf"

meta_table <- dge_v$targets
count_table_t <- as.data.frame(t(dge_v$E))
pca_prep <- prcomp(count_table_t, scale. = TRUE)

pca_plot <- autoplot(pca_prep, label, shape=shape_column, data = meta_table, colour = color_column,
)
ggsave(plot_save_name, device = "pdf",units = "cm", width = 16, height = 14)
pca_plot




#_1
plot_save_name_1 <- "PCR_Plod_1.pdf"
meta_table_1 <- dge_v_1$targets
count_table_t_1 <- as.data.frame(t(dge_v_1$E))
pca_prep_1 <- prcomp(count_table_t_1, scale. = TRUE)

pca_plot_1 <- autoplot(pca_prep_1, label, shape=shape_column, data = meta_table_1, colour = color_column,
)
ggsave(plot_save_name_1, device = "pdf",units = "cm", width = 16, height = 14)
pca_plot_1

#_2
shape_column_2 <- "Diagnosis"
color_column_2 <-"Diagnosis"
plot_save_name_2 <- "PCR_Plod_2.pdf"
meta_table_2 <- dge_v_2$targets
count_table_t_2 <- as.data.frame(t(dge_v_2$E))
pca_prep_2 <- prcomp(count_table_t_2, scale. = TRUE)

pca_plot_2 <- autoplot(pca_prep_2, label, shape=shape_column_2, data = BI_meta_Diag, colour = color_column_2
)
ggsave(plot_save_name_2, device = "pdf",units = "cm", width = 16, height = 14)
pca_plot_2

#_3
shape_column_3 <- "Diag_Gender_combo"
color_column_3 <-"Diag_Gender_combo"
plot_save_name_3 <- "PCR_Plod_3.pdf"
meta_table_3 <- dge_v_3$targets
count_table_t_3 <- as.data.frame(t(dge_v_3$E))
pca_prep_3 <- prcomp(count_table_t_3, scale. = TRUE)

pca_plot_3 <- autoplot(pca_prep_3, label, shape=shape_column_3, data = BI_meta_combo, colour = color_column_3
)
ggsave(plot_save_name_3, device = "pdf",units = "cm", width = 16, height = 14)
pca_plot_3

#"Each point is a sample, 
#Positioned according to its overall expression pattern 
#(compressed into PC1 and PC2), and coloured/shaped by sex/diagnosis."
#Linear: I want_: expression of gene g in sample i=baseline+effect of diagnosis+effect of sex+random 
#In matrix form:
#The design matrix holds the information about your samples (diagnosis, sex, etc.).
#The expression matrix holds the response for each gene.
#Linear modeling then gives, for each gene:
  #An estimated log fold‑change for each factor (e.g. “diagnosis 2 vs 0”).
  #A standard error, t‑statistic and p‑value testing whether that effect is zero.

#Ensures the comparison “diagnosis2” means “2 vs 0” and “sexm” means “m vs f”
pheno$diagnosis <- factor(pheno$diagnosis,
                          levels = c("0", "2")) 

pheno$sex <- factor(pheno$sex, levels = c("f", "m"))

#new 
#logFC: log2 fold change (e.g. +1 means 2‑fold higher in diagnosis 2 vs 0).
#P.Value: raw p‑value for the hypothesis “no effect of diagnosis”.
#adj.P.Val: p‑value corrected for multiple testing using Benjamini–Hochberg FDR.
res_diag <- topTable(fit1, coef = "diagnosis2", number = Inf, adjust.method = "BH")
res_sex  <- topTable(fit1, coef = "sexm", number = Inf, adjust.method = "BH")
#gives genes with FDR < 5% and at least 2‑fold change
DE_diag <- subset(res_diag, adj.P.Val < 0.05 & abs(logFC) > 1)
#How many genes are differentially expressed?
HMGADE <- sum(res_diag$adj.P.Val < 0.05)
# 0L means none passed the FDR < 0.05 threshold
#This is common when the disease effect is weak relative to noise, sample size is limited, or there are strong unmodelled sources of variation

plot(res_diag$logFC, -log10(res_diag$P.Value),
     xlab = "log2 fold change (diagnosis 2 vs 0)",
     ylab = "-log10 p-value")

#x‑axis: logFC (log2 fold change).
#Values near 0: little or no change between conditions.
#Large positive: strongly up‑regulated in diagnosis 2 vs 0.
#Large negative: strongly down‑regulated.
#y‑axis: -log10(p‑value) (often the raw p‑value).
#Larger values mean smaller p‑values, i.e. more statistically significant.
#Example: p = 0.01 → −log10(0.01) = 2; p = 10⁻⁶ → 6.
#Each point is a gene. So the volcano plot shows, at a glance:
#1-Magnitude of change (how big the effect is) along the x‑axis.
#2-Strength of evidence (how unlikely it is that the effect is just noise) along the y‑axis.
#Points far right and high up: strongly up‑regulated and very significant.
#I didn't get any points that are out of the ordinary


res_diag$negLogP <- -log10(res_diag$P.Value)

ggplot(res_diag, aes(x = logFC, y = negLogP)) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), colour = "blue", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), colour = "red", linetype = "dashed")

