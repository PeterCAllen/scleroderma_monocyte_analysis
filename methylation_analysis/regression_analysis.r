# set lib path
.libPaths("/cluster/home/pcallen/miniconda3/envs/Minfi/lib/R/library")

# set working directory
setwd("/cluster/home/pcallen/projects/Paula_SSc_Methyl_Analysis/methylation_analysis/")

library(tidyverse)
library(data.table)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
library(edgeR)
library(kableExtra)

sampleSheet <- read_csv("data/methyl_sampleSheet.csv", na = "NA")

# replacing NA in disease_subset to healthy
# sampleSheet$disease_subset <- sampleSheet$disease_subset %>% replace_na("Healthy")

# #matching rowname syntax to betas
# sampleSheet$sample_name <- gsub("-", ".", sampleSheet$sample_name)

# sampleSheet <- sampleSheet %>% column_to_rownames(var = "sample_name")

sampleSheet <- sampleSheet %>%
  mutate_at(c(12:15), ~ replace_na(., "N"))

## Imputed Betas
finalBetas_auto <- data.frame(fread("data/finalBetas_autos_Imp.txt", stringsAsFactors = FALSE))
rownames(finalBetas_auto) <- finalBetas_auto[, 1]
finalBetas_auto <- finalBetas_auto[, -1]
colnames(finalBetas_auto) <- gsub("X", "", colnames(finalBetas_auto))
finalBetas_auto <- finalBetas_auto[, order(colnames(finalBetas_auto))]

finalBetas_x <- data.frame(fread("data/finalBetas_X_Imp.txt", stringsAsFactors = FALSE))
rownames(finalBetas_x) <- finalBetas_x[, 1]
finalBetas_x <- finalBetas_x[, -1]
colnames(finalBetas_x) <- gsub("X", "", colnames(finalBetas_x))
finalBetas_x <- finalBetas_x[, order(colnames(finalBetas_x))]

finalBetasImp <- rbind(finalBetas_auto, finalBetas_x)

rm(finalBetas_auto, finalBetas_x)

## Non-imputed Betas
finalBetas_auto.nonimp <- data.frame(fread("data/finalBetas_autos_NonImp.txt", stringsAsFactors = FALSE))
rownames(finalBetas_auto.nonimp) <- finalBetas_auto.nonimp[, 1]
finalBetas_auto.nonimp <- finalBetas_auto.nonimp[, -1]
colnames(finalBetas_auto.nonimp) <- gsub("X", "", colnames(finalBetas_auto.nonimp))
finalBetas_auto.nonimp <- finalBetas_auto.nonimp[, order(colnames(finalBetas_auto.nonimp))]

finalBetas_x.nonimp <- data.frame(fread("data/finalBetas_X_NonImp.txt", stringsAsFactors = FALSE))
rownames(finalBetas_x.nonimp) <- finalBetas_x.nonimp[, 1]
finalBetas_x.nonimp <- finalBetas_x.nonimp[, -1]
colnames(finalBetas_x.nonimp) <- gsub("X", "", colnames(finalBetas_x.nonimp))
finalBetas_x.nonimp <- finalBetas_x.nonimp[, order(colnames(finalBetas_x.nonimp))]

finalBetasNonImp <- rbind(finalBetas_auto.nonimp, finalBetas_x.nonimp)
rm(finalBetas_auto.nonimp, finalBetas_x.nonimp)


# ensure both dataset columns match
all(colnames(finalBetasImp) == colnames(finalBetasNonImp))

# probes from Pidsley 2016 (EPIC)
# epic.cross1 <- read.csv('data/illumina450k_filtering/EPIC/13059_2016_1066_MOESM1_ESM.csv', head = T)
# epic.cross2 <- read.csv('data/illumina450k_filtering/EPIC/13059_2016_1066_MOESM2_ESM.csv', head = T)
# epic.cross3 <- read.csv('data/illumina450k_filtering/EPIC/13059_2016_1066_MOESM3_ESM.csv', head = T)
# epic.variants1 <- read.csv('data/illumina450k_filtering/EPIC/13059_2016_1066_MOESM4_ESM.csv', head = T)
# epic.variants2 <- read.csv('data/illumina450k_filtering/EPIC/13059_2016_1066_MOESM5_ESM.csv', head = T)
# epic.variants3 <- read.csv('data/illumina450k_filtering/EPIC/13059_2016_1066_MOESM6_ESM.csv', head = T)

# probes significantly associated with smoking
epic.smoking <- read.csv("data/illumina450k_filtering/christiansen_2021_sig_smoking_cpgs.txt", head = F) %>% pull(V1)

# additional filter probes
epic.filter.probes <- c(
  # as.character(epic.cross1$X),
  # as.character(epic.cross2$X),
  # as.character(epic.cross3$X),
  # as.character(epic.variants1$PROBE),
  # as.character(epic.variants2$PROBE),
  # as.character(epic.variants3$PROBE),
  as.character(epic.smoking)
)
# final list of unique probes
epic.filter.probes <- unique(epic.filter.probes)

# length(epic.filter.probes) #142371 probes to be filtered

# removing from dataset
finalBetasNonImp <- finalBetasNonImp[!(rownames(finalBetasNonImp) %in% epic.filter.probes), ] # 720098 probes remaining
finalBetasImp <- finalBetasImp[!(rownames(finalBetasImp) %in% epic.filter.probes), ]

detP <- data.frame(fread("data/pval.txt", stringsAsFactors = FALSE)) %>%
  column_to_rownames("V1")

colnames(detP) <- gsub("X", "", colnames(detP))

detP <- detP[rownames(finalBetasNonImp), colnames(finalBetasNonImp)]


# ChAMP Filtering ----
library(ChAMP)

## Create object for champ ----
myBetas <- list(
  beta = as.matrix(finalBetasNonImp),
  pd = sampleSheet,
  detP = as.matrix(detP)
)

## Filter failed probes ----
myLoad <- champ.filter(
  beta = as.matrix(myBetas$beta),
  pd = myBetas$pd,
  detP = myBetas$detP,
  arraytype = "EPIC",
  filterBeads = TRUE,
  filterNoCG = TRUE,
  filterSNPs = FALSE,
  filterMultiHit = TRUE,
  filterXY = TRUE,
  SampleCutoff = 0.1
)

meta <- myLoad$pd
betas <- myLoad$beta[, match(meta$methyl_id, colnames(myLoad$beta))]

pca <- prcomp(t(na.omit(betas)))$x
pca <- pca[colnames(betas), ]

# regression
ssc_model <- apply(na.omit(betas), 1, function(x) summary(lm(as.numeric(x) ~ as.factor(meta$disease_status) + as.numeric(meta$age_at_enrollment) + as.numeric(pca[, 1]) + as.numeric(pca[, 2])))$coeff[2, 4])
ssc_model_cof <- apply(na.omit(betas), 1, function(x) summary(lm(as.numeric(x) ~ as.factor(meta$disease_status) + as.numeric(meta$age_at_enrollment) + as.numeric(pca[, 1]) + as.numeric(pca[, 2])))$coeff[2,1])

# Organize & Annotate CpG Regression Results -----------------------------------------
top_cpgs <- tibble(CpG = names(ssc_model), beta_cof = ssc_model_cof, pval = ssc_model) %>% arrange(pval)
top_cpgs <- top_cpgs %>% mutate(p.adjusted = p.adjust(top_cpgs$pval, method = "BH"))


# created a function to name the top 500 genes
# gene_match_500 <- function(x){
#   EPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#   intersect <- match(x$CpG, rownames(EPIC))
#   EPIC <- EPIC[intersect,]
#   EPIC <- as_tibble(EPIC)
#
#   EPIC <- EPIC %>% dplyr::select(Name, UCSC_RefGene_Name, chr, pos, Relation_to_Island)
#
#   EPIC.granges <- data.frame(EPIC[,3:4])
#   EPIC.granges$stop <- EPIC.granges$pos+5
#   EPIC.granges$CpG <- EPIC$Name
#
#   colnames(EPIC.granges) <- c("chr", "start", "stop", "CpG")
#
#   genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#   cpgs <- makeGRangesFromDataFrame(EPIC.granges)
#
#
#   epic.genes <- matchGenes(cpgs,genes) %>% as_tibble() %>%
#     add_column(CpG=x$CpG, pval=x$pval, p.adjusted=x$p.adjusted)
#   return(epic.genes)
# }
#
# top_cpgs_annotated <- gene_match_500(top_cpgs)

data.frame(cpg = rownames(probe.features)[match(top_cpgs$CpG[1:25], rownames(probe.features))], gene = probe.features$gene[match(top_cpgs$CpG[1:25], rownames(probe.features))])

top_cpgs_annotated[1:20, ] %>% dplyr::select(CpG, name, annotation, description, pval, p.adjusted)



## Regression Model using just the cases that are on immunesuppressants

# SampleSheet
sampleSheet_IMMSUP <- sampleSheet %>%
  filter(immunosuppresive_therapy == "Y" | disease_status == "Control")

## Betas
betas_IMMSUP <- finalBetasNonImp[, sampleSheet_IMMSUP$methyl_id]
betasIMP_IMMSUP <- finalBetasImp[, sampleSheet_IMMSUP$methyl_id]

pca_immsup <- prcomp(t(na.omit(betasIMP_IMMSUP)))$x

ssc_model_immsup <- apply(na.omit(betas_IMMSUP), 1, function(x) summary(lm(as.numeric(x) ~ as.factor(sampleSheet_IMMSUP$disease_status) + as.numeric(sampleSheet_IMMSUP$age_at_enrollment) + as.numeric(pca_immsup[, 1]) + as.numeric(pca_immsup[, 2])))$coeff[2, 4])
ssc_model_immsup_cof <- apply(na.omit(betas_IMMSUP), 1, function(x) summary(lm(as.numeric(x) ~ as.factor(sampleSheet_IMMSUP$disease_status) + as.numeric(sampleSheet_IMMSUP$age_at_enrollment) + as.numeric(pca_immsup[, 1]) + as.numeric(pca_immsup[, 2])))$coeff[2, 1])

# Organize & Annotate CpG Regression Results -----------------------------------------
top_cpgs_immsup <- tibble(CpG = names(ssc_model_immsup), beta_cof = ssc_model_immsup_cof, pval = ssc_model_immsup) %>% arrange(pval)

top_cpgs_immsup <- top_cpgs_immsup %>% mutate(p.adjusted = p.adjust(top_cpgs_immsup$pval, method = "BH"))


# save.image('data/_regression_smoking_cpg_removed_wImmunosuppressants.RData')
# load("data/_regression_smoking_cpg_removed_wImmunosuppressants.RData")

# write to output
# write_delim(top_cpgs_annotated, file = "output/methylated_genes.txt")

library(ComplexHeatmap)

# Combine the meta data with methylation data -----------------------
betas <- finalBetasNonImp %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sample_id")

metadata <- sampleSheet %>%
  rownames_to_column(., var = "sample_id")

beta.table <- data.table(betas, key = "sample_id")
meta.table <- data.table(metadata, key = "sample_id")

combined_data <- merge(meta.table, beta.table)

## Creating heatmaps from top CpGs -----------------------
data_heatmap_mx <- combined_data %>%
  dplyr::select(sample_id, disease_status, disease_subset, race, c(top_cpgs_annotated$CpG[1:20]))

# Removing Samples that have too many NAs ----------------
data_heatmap_mx <- data_heatmap_mx %>%
  filter(!(rowSums(is.na(data_heatmap_mx)) / ncol(data_heatmap_mx) > 0.5))

test_meta <- data_heatmap_mx %>%
  dplyr::select(disease_subset, sample_id)

test_meta$disease_subset <- factor(test_meta$disease_subset, levels = c("Healthy", "lcSSc", "dcSSc", "ssSSc"))

test_meta <- test_meta %>%
  arrange(disease_subset)

# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
colors <- list("disease_subset" = c("Healthy" = "green", "lcSSc" = "yellow", "dcSSc" = "brown", "ssSSc" = "red"))

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  df = test_meta,
  which = "col",
  col = colors,
  show_legend = c("disease_subset" = TRUE, "sample_id" = FALSE),
  annotation_label = c("Disease Subset", ""),
  annotation_legend_param = list(disease_subset = list(title = "Disease Subset", labels = c("Control", "lcSSc", "dcSSc", "ssSSc")))
)

data_heatmap_mx <- data.frame(t(data_heatmap_mx))
colnames(data_heatmap_mx) <- data_heatmap_mx["sample_id", ]
data_heatmap_mx_1 <- data_heatmap_mx[-c(1:4), ] %>% as.matrix()

data_heatmap_mx <- matrix(as.numeric(data_heatmap_mx_1), ncol = ncol(data_heatmap_mx_1))
rownames(data_heatmap_mx) <- rownames(data_heatmap_mx_1)
colnames(data_heatmap_mx) <- colnames(data_heatmap_mx_1)

data_heatmap_mx <- data_heatmap_mx[, match(test_meta$sample_id, colnames(data_heatmap_mx))]

# rename rowname to have gene name
data_heatmap_mx_named <- data_heatmap_mx

EPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
EPIC <- as_tibble(EPIC)

EPIC <- EPIC[match(top_cpgs_annotated$CpG[1:20], EPIC$Name), ]

EPIC$gene <- sapply(strsplit(EPIC$UCSC_RefGene_Name, ";"), `[`, 1)

EPIC$cpg_name <- ifelse(!is.na(EPIC$gene), paste0(EPIC$Name, "-", EPIC$gene), paste0(EPIC$Name))

rownames(data_heatmap_mx_named) <- EPIC$cpg_name[match(EPIC$Name, rownames(data_heatmap_mx_named))]

# # Combine the heatmap and the annotation
Heatmap(data_heatmap_mx_named, name = "Beta Value", top_annotation = ha, cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = FALSE, show_row_names = FALSE)

Heatmap(data_heatmap_mx_named, name = "Beta Value", top_annotation = ha, cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = FALSE, show_row_names = FALSE)

# Comparison between Immunosuppressant vs Control | Full Case vs Control ----

top_cpgs_immsup <- top_cpgs_immsup[match(top_cpgs$CpG, top_cpgs_immsup$CpG), ]

top_cpgs_immsup$pval_negLog10 <- -log10(top_cpgs_immsup$pval)
top_cpgs$pval_negLog10 <- -log10(top_cpgs$pval)

library(ggpubr)

p1 <- ggplot(NULL, aes(x = top_cpgs$pval_negLog10, y = top_cpgs_immsup$pval_negLog10)) +
  geom_point(shape = 1, size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_x_continuous(limits = c(0, 8)) +
  scale_y_continuous(limits = c(0, 8)) +
  # coord_fixed(ratio=1) +
  ggtitle("All CpGs") +
  labs(x = "Full Model: -log10(pval)", y = "Immunosuppressant Model: -log10(pval)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20)
  )

p2 <- ggplot(NULL, aes(x = top_cpgs$beta_cof, y = top_cpgs_immsup$beta_cof)) +
  geom_point(shape = 1, size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  # scale_x_continuous(limits = c(0, 8)) +
  # scale_y_continuous(limits = c(0, 1.5)) +
  # coord_fixed(ratio=1) +
  ggtitle("All CpGs") +
  labs(x = "Full Model: Beta Coefficients", y = "Immunosuppressant Model: Beta Coefficients") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20)
  )

ggarrange(p1, p2,
  ncol = 2, nrow = 1
)
ggsave("neglog10_pval_full_and_immunosupp_model_betacof.png", width = 18, height = 10, units = "in")

# # Contour Plot
# x = top_cpgs$beta_cof; y = top_cpgs_immsup$beta_cof
# z <- sqrt(outer(x ^ 2, y ^ 2, "+"))
# 
# contour(x, y, z,
#         labcex = 1.2, labels = 1:10,
#         lwd = 2, lty = 1) 

# write_csv(full_cpgs_distinct_annotated, file = "~/Downloads/cpgs_not_in_immsup_model.csv")

# Regression with medication categories as covariate ----
ssc_model_immsup_cov <- apply(na.omit(betas), 1, function(x) {
  summary(lm(as.numeric(x) ~ as.factor(sampleSheet$disease_status) +
    as.factor(sampleSheet$immunosuppresive_therapy) +
    as.factor(sampleSheet$prednisone) +
    as.factor(sampleSheet$mmf) +
    as.factor(sampleSheet$hcq) +
    as.numeric(sampleSheet$age_at_enrollment) +
    as.numeric(pca[, 1]) + as.numeric(pca[, 2])))$coeff[2, 4]
})

# Organize & Annotate CpG Regression Results -----------------------------------------
top_cpgs_immunoCov <- tibble(CpG = names(ssc_model_immsup_cov), pval = ssc_model_immsup_cov) %>% arrange(pval)

top_cpgs_immunoCov <- top_cpgs_immunoCov %>% mutate(p.adjusted = p.adjust(top_cpgs_immunoCov$pval, method = "BH"))

# save(top_cpgs_immunoCov, top_cpgs, top_cpgs_immsup, file = "medicated_regression_pplot.RData")
load("medicated_regression_pplot.RData")

top_cpgs_immunoCov <- top_cpgs_immunoCov %>% 
  mutate(rank = 1:nrow(top_cpgs_immunoCov),
         constant = nrow(top_cpgs_immunoCov)) %>% 
  mutate(theor_dist = rank/constant) %>% 
  mutate(neglog10_pval = -log10(pval),
         neglog10_rank = -log10(theor_dist))

## Creating a pp plot ----


### Top 200 ----
p1 <- ggplot(top_cpgs_immunoCov, aes(x = neglog10_rank, y = neglog10_pval)) +
  geom_point(shape = 1, size = 3) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  # scale_x_continuous(limits = c(0, 8)) +
  # scale_y_continuous(limits = c(0, 1.5)) +
  # coord_fixed(ratio=1) +
  ggtitle("PP Plot") +
  labs(x = "Theoretical Distribution", y = "Observed") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20)
  )
ggsave(p1, "pp_plot.png", width = 10, height = 10, units = "in", device = "png")

png("pp_plot.png", width=5, height=5, units="in", res=300)
plot(x=top_cpgs_immunoCov$neglog10_rank, y=top_cpgs_immunoCov$neglog10_pval,
     main="P-P Plot", xlab="Theoretical", ylab="Empirical")
abline(a=0, b=1, col="red")
dev.off()

top_cpgs <- top_cpgs %>% 
  mutate(rank = 1:nrow(top_cpgs),
         constant = nrow(top_cpgs)) %>% 
  mutate(theor_dist = rank/constant) %>% 
  mutate(neglog10_pval = -log10(pval),
         neglog10_rank = -log10(theor_dist))

png("pp_plot_wo_medication.png", width=5, height=5, units="in", res=300)
plot(x=top_cpgs$neglog10_rank, y=top_cpgs$neglog10_pval,
     main="P-P Plot", xlab="Theoretical", ylab="Empirical")
abline(a=0, b=1, col="red")
dev.off()

## Plot Beta Values for top 5 cpgs
library(ggrepel)
library(ggpubr)
give.n <- function(x){return(c(y = 0, label = length(x)))} # function to give numbers

top_cpgs_names <- top_cpgs[1:89,] %>%
  pull(CpG)

methylation_data <- finalBetasNonImp[match(top_cpgs_names, rownames(finalBetasNonImp)),] %>%
  t() %>%
  as.data.frame()  %>%
  rownames_to_column("methyl_id")

combined_data <- left_join(sampleSheet, methylation_data, by="methyl_id")

combined_data %>%
  ggplot(., aes(x = factor(x = disease_status, levels = c("Control", "SSc")),
                y = cg22805491,
                fill=factor(disease_status, levels = c("Control", "SSc")))) +
  geom_boxplot() +
  ylim(0,1) +
  ggtitle("cg22805491 - NIN") +
  labs(x="", y="Beta Score", fill="") +
  theme_bw() +
  theme(panel.background = element_blank(),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 1), size = 7)

ggsave(filename = "topCpG_nin.png",
       width = 8, height=8, units = "in", dpi = 400)

combined_data %>%
  ggplot(., aes(x = factor(x = disease_status, levels = c("Control", "SSc")),
                y = cg18507060,
                fill=factor(disease_status, levels = c("Control", "SSc")))) +
  geom_boxplot() +
  ylim(0,1) +
  ggtitle("cg18507060 - OAS3") +
  labs(x="", y="Beta Score", fill="") +
  theme_bw() +
  theme(panel.background = element_blank(),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 1), size = 7)

ggsave(filename = "topCpG_OAS3.png",
       width = 8, height=8, units = "in", dpi = 400)

 combined_data %>%
  ggplot(., aes(x = factor(x = disease_status, levels = c("Control", "SSc")),
                y = cg08653580,
                fill=factor(disease_status, levels = c("Control", "SSc")))) +
  geom_boxplot() +
  ylim(0,1) +
  ggtitle("cg08653580 - CD5") +
  labs(x="", y="Beta Score", fill="") +
  theme_bw() +
  theme(panel.background = element_blank(),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 1), size = 7)

ggsave(filename = "topCpG_CD5.png",
       width = 8, height=8, units = "in", dpi = 400)

 combined_data %>%
  ggplot(., aes(x = factor(x = disease_status, levels = c("Control", "SSc")),
                y = cg23493751,
                fill=factor(disease_status, levels = c("Control", "SSc")))) +
  geom_boxplot() +
  ylim(0,1) +
  ggtitle("cg23493751 - CCR3") +
  labs(x="", y="Beta Score", fill="") +
  theme_bw() +
  theme(panel.background = element_blank(),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 1), size = 7)

ggsave(filename = "topCpG_CCR3.png",
       width = 8, height=8, units = "in", dpi = 400)

## Creating heatmaps from top CpGs -----------------------
library(ComplexHeatmap)

zscores <- combined_data %>% column_to_rownames("methyl_id") %>%
  dplyr::select(c(top_cpgs_names))

zscores <- scale(zscores) %>% as.data.frame()

zscores <- zscores %>%
  rownames_to_column("methyl_id")

combined_data <- left_join(sampleSheet, zscores, by="methyl_id")

data_heatmap_mx <- combined_data %>%
  dplyr::select(methyl_id, disease_status, race, c(top_cpgs_names))

meta <- data_heatmap_mx %>%
  dplyr::select(disease_status, methyl_id)

meta$disease_status <- factor(meta$disease_status, levels = c("Control", "SSc"))

meta <- meta %>%
  arrange(disease_status) %>%
  as.data.frame()

# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
colors = list('disease_status' = c("Control" = "green", "SSc" = "red"))

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  df = data.frame(disease_status=meta[,1]),
  which = 'col',
  col = colors,
  show_legend = c("disease_status" = TRUE),
  annotation_label = c(" "),
  annotation_legend_param = list(disease_status = list(title = "Disease Status", labels = c("Control", "SSc"))))

data_heatmap_mx <- data.frame(t(data_heatmap_mx))
colnames(data_heatmap_mx) <- data_heatmap_mx["methyl_id",]
data_heatmap_mx_1 <- data_heatmap_mx[-c(1:3),] %>% as.matrix()

data_heatmap_mx <- matrix(as.numeric(data_heatmap_mx_1), ncol = ncol(data_heatmap_mx_1))
rownames(data_heatmap_mx) <- rownames(data_heatmap_mx_1); colnames(data_heatmap_mx) <- colnames(data_heatmap_mx_1)

data_heatmap_mx <- data_heatmap_mx[,match(meta$methyl_id, colnames(data_heatmap_mx))]

# rename rowname to have gene name
data_heatmap_mx_named <- data_heatmap_mx

# # Combine the heatmap and the annotation
library(circlize)
beta_colors <- colorRampPalette(c("darkblue", "black", "yellow3"))(5)

p1 <- Heatmap(data_heatmap_mx_named,
              name = "Z-Score",
              col = beta_colors,
              top_annotation = ha,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              column_names_max_height = unit(4, "cm"),
              show_column_names = FALSE,
              show_row_names = FALSE,
              row_gap = c(1),
              heatmap_legend_param = list(
                legend_direction = "horizontal",
                legend_width = unit(6, "cm"),
                title_gp = gpar(fontsize = 15,
                            fontface = "bold"),
                labels_gp = gpar(fontsize = 15)
              ))

pdf("top89cpgs_heatmap.pdf", width = 8, height = 8)
  draw(p1, heatmap_legend_side="bottom", annotation_legend_side="right")



