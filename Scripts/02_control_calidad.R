# Libraries
library(ggplot2)
library(R.utils)
library(stringr)
library(pheatmap)
library(limma)

# Preparación escritorio 
rm(list=ls())
setwd("../Estudio_transcriptómico_endometrio/")
results_folder = "Results/02_control_calidad/"
dir.create(results_folder)

# Carga de datos
load("Results/01_descarga_datos/datos_crudos.RData", verbose=T)


# Filtrado de datos -------------------------------------------------------

# Nos quedamos con las columnas del identificador y de los conteos de las muestras 
join_humanWG = join[, c(1, (grep("SAMPLE", colnames(join))))]
colnames(join_humanWG) = paste(gsub("[\\. ]", "_", colnames(join_humanWG)), 
                               "WG", sep = "_")
# Y guardamos por separado los p-valores de detección 
pval_humanWG = join[, c(1, grep("Pval", colnames(join)))]
colnames(pval_humanWG) = colnames(join_humanWG)

# Combinamos la info de HumanHT y HumanWG 
id_ht = datos_HumanHT_filt$ID_REF
id_wg = join_humanWG$join_ids____WG
ids = intersect(id_ht, id_wg)
length(ids)
# 25706

# Separamos la info de HumanHT
join_humanHT = datos_HumanHT_filt[ids, c(1, 10, 13, grep("SAMPLE", colnames(datos_HumanHT_filt)))]
colnames(join_humanHT) = paste(gsub("[\\. ]", "_", colnames(join_humanHT)), 
                               "HT", sep = "_")
ilmn_to_gene = datos_HumanHT_filt[, 1:29]
all(rownames(join_humanHT) == rownames(join_humanWG))    
# TRUE 

df_final = cbind(join_humanHT, join_humanWG[-1])
counts = df_final[, 4:ncol(df_final)]

# Generamos df metadata 
metadata = data.frame("Sample_ID" = colnames(counts))
metadata$chip = ifelse(endsWith(metadata$Sample_ID, "WG"), "HumanWG-6", "HumanHT-12")
metadata$enfermedad = str_split_fixed(metadata$Sample_ID, "_", 4)[, 3]
metadata$condicion = ifelse(metadata$enfermedad %in% c("CP", "CE"), "Control", "Patient")
metadata$g_enfermedad = ifelse(metadata$enfermedad == "CE", "Control_endometrium", 
                               ifelse(metadata$enfermedad == "PE", "Patient_endometrium", 
                                      ifelse(metadata$enfermedad == "OMA", "Ovarian_endometriosis", 
                                             ifelse(metadata$enfermedad %in% c("PeLR", "PeLB", "PeLW"), "Peritoneal_endometriosis", 
                                                    ifelse(metadata$enfermedad %in% c("DiEB", "DiEIn", "RE", "REV", "SuL"), "Deep_endoemtriosis", 
                                                           ifelse(metadata$enfermedad == "CP", "Control_peritoneum", "Patient_peritoneum"))))))
metadata$sheet = str_split_fixed(metadata$Sample_ID, "_", 5)[, 4]
rownames(metadata) = metadata$Sample_ID


# Análisis exploratorio ---------------------------------------------------

# Boxplot
tidy.counts = log(counts + 1)
tidy.counts$genes = rownames(tidy.counts)
tidy.counts = reshape2::melt(tidy.counts)
colnames(tidy.counts) = c("genes", "Sample_ID", "value")
tidy.counts = merge(tidy.counts, metadata, by.x = "Sample_ID", by.y = "Sample_ID")

ggplot(tidy.counts, aes(x=Sample_ID, y=value, fill=chip)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#DE369D", "#06d6a0")) + 
  theme(axis.text.x = 
          element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Boxplot") + 
  xlab("Samples") + 
  ylab("Log10 Read counts")


# PCA
pca.res <- PCAtools::pca(mat=counts, metadata=metadata,  scale=TRUE)
# Por plataforma de secuenciación 
plot_chip = PCAtools::biplot(pca.res, colby='chip', 
                             colkey = c("#DE369D", "#06d6a0"),  
                             legendPosition = "right", 
                             lab = NULL)
plot_chip 
ggsave(filename = paste0(results_folder, "Figura_1A_PCA_Chip.png"), plot = plot_chip, 
       width = 10, height = 10, dpi = 300)


#### Solo HumanWG-6 ####
# Como hay mucha diferencia entre chips nos quedamos solo con HUMANWG-6
table(metadata$chip, metadata$condicion)
#            Control Patient
# HumanHT-12       0      72
# HumanWG-6       67     269

human_wg = rownames(metadata)[metadata$chip == "HumanWG-6"]
counts_wg = counts[, human_wg]
dim(counts_wg)
# 25706   336
metadata_wg = metadata[human_wg, ]
dim(metadata_wg)
# 336   6  

pca_res_wg <- PCAtools::pca(mat=counts_wg, metadata=metadata_wg, scale=TRUE)
# Por control y paciente
PCAtools::biplot(pca_res_wg, colby='condicion', 
                 colkey = c("#ffd166", "#6AD8F6"),  
                 legendPosition = "right", lab = NULL)

plot_sheet = PCAtools::biplot(pca_res_wg, colby='sheet', 
                             legendPosition = "right", lab = NULL)
plot_sheet
ggsave(filename = paste0(results_folder, "Figura1B_PCA_Sheet.png"), 
       plot = plot_sheet, width = 6, height = 5, dpi = 300)

# Eliminamos sheet 7, Réplicas y SAMPLE_137_SuL_sheet_3_WG
metadata_wg_filt = metadata_wg[metadata_wg$sheet != "7", ]
metadata_wg_filt = metadata_wg_filt[metadata_wg_filt$sheet != "Replicate", ]
metadata_wg_filt = metadata_wg_filt[metadata_wg_filt$Sample_ID != "SAMPLE_137_SuL_3_WG", ]
counts_wg_filt = counts_wg[, metadata_wg_filt$Sample_ID]
pval_humanWG_filt = pval_humanWG[, metadata_wg_filt$Sample_ID]


# Volvemos a hacer el PCA 
pca_res_filt = PCAtools::pca(mat = counts_wg_filt, metadata = metadata_wg_filt, scale=TRUE)
PCAtools::biplot(pca_res_filt, colby='condicion', 
                 colkey = c("#ffd166", "#6AD8F6"),  
                 legendPosition = "right", lab = NULL)


# Normalización -----------------------------------------------------------

# Transformación logarítmica de los conteos
counts_log = log2(counts_wg_filt + 1)

# Boxplot pre-normalización
png(filename = paste0(results_folder, "Figura_sup_1A_Boxplot_PreNorm.png"), 
    width = 2400, height = 1800, res = 300)
boxplot(counts_log, main="Datos Sin Normalizar (Log2)", las=2, outline=FALSE,
        names = F)
dev.off()

# Normalización por cuantiles 
counts_norm = normalizeBetweenArrays(counts_log, method = "quantile")
dim(counts_norm)
# 25706   305

# Boxplot post-normalización 
png("Results/02_control_calidad/Figura_sup_1B_Boxplot_PostNorm.png", 
    width = 2400, height = 1800, res = 300)
boxplot(counts_norm, main="Datos Normalizados (Log2)", las=2, outline=FALSE,
        names = F)
dev.off()


# Guardar datos -----------------------------------------------------------

save(ilmn_to_gene, counts_wg_filt, counts_norm, metadata_wg_filt, pval_humanWG_filt, file = paste0(results_folder, "datos_filtrados.RData"))

