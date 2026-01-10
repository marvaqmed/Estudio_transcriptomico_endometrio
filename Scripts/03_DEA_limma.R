# Libraries
library(readxl)
library(R.utils)
library(stringr)
library(dplyr)
library(DESeq2)
library(limma)
library(NOISeq)
library(tidyverse)
library(ggvenn)
library(UpSetR)

# Preparación escritorio 
rm(list=ls())
setwd("../Estudio_transcriptómico_endometrio/")
results_folder = "Results/03_DEA_limma/"
dir.create(results_folder)

# Carga de datos
load("Results/02_control_calidad/datos_filtrados.RData", verbose = T)

# Filtro por p-valor detección ------------------------------------------------

table(metadata_wg_filt$g_enfermedad)
metadata_wg_filt$g_enfermedad = gsub(" ", "_", metadata_wg_filt$g_enfermedad)
# Control_endometrium       Control_peritoneum       Deep_endoemtriosis    Ovarian_endometriosis 
# 34                        24                       69                    23 
# Patient_endometrium       Patient_peritoneum       Peritoneal_endometriosis 
# 60                        25                       70 

# ponemos como umbral que supere el p-valor en al menos la mitad de muestras del grupo más pequeño 
k = 12
umbral = 0.05

to_keep = pval_humanWG_filt < umbral
table(to_keep)
#   FALSE     TRUE 
# 3926441 3913889  
row_keep = rowSums(to_keep) >= k
table(row_keep)
# FALSE  TRUE 
#  8549 17157
counts_filt = counts_norm[row_keep, ]
dim(counts_filt)
# 17157   305

# Nos quedamos con 17.157 genes y 305 muestras


# Limma  ------------------------------------------------------------------

# Matriz de diseño
grupo = factor(metadata_wg_filt$g_enfermedad)
design = model.matrix(~ 0 + grupo)
colnames(design) = levels(grupo)

# Ajustar el modelo 
fit = lmFit(object = counts_filt, design = design)

# Definir contrastes 
cont_matrix = makeContrasts(
  pat_endom_vs_control_endom = Patient_endometrium - Control_endometrium,
  deep_vs_endometriosis = Deep_endoemtriosis - Patient_endometrium,
  ovario_vs_endom = Ovarian_endometriosis - Patient_endometrium, 
  peritoneo_vs_endom = Peritoneal_endometriosis - Patient_endometrium, 
  pat_peritoneo_vs_control = Patient_peritoneum - Control_peritoneum, 
  levels = design
)

# Cálculo contrastes 
fit2 = contrasts.fit(fit, cont_matrix)
fit2 = eBayes(fit2, trend = T)

# Sacamos resultados de cada comparación 
tt1 = limma::topTable(fit2, coef="pat_endom_vs_control_endom", adjust="BH", sort.by="none", number=nrow(fit2))
tt2 = limma::topTable(fit2, coef="deep_vs_endometriosis", adjust="BH", sort.by="none", number=nrow(fit2))
tt3 = limma::topTable(fit2, coef="ovario_vs_endom", adjust="BH", sort.by="none", number=nrow(fit2))
tt4 = limma::topTable(fit2, coef="peritoneo_vs_endom", adjust="BH", sort.by="none", number=nrow(fit2))
tt5 = limma::topTable(fit2, coef="pat_peritoneo_vs_control", adjust="BH", sort.by="none", number=nrow(fit2))

# Añadimos info de identificadores Entrez y Symbol y filtramos duplicados 
cambiar_ident = function(tt, ilmn_to_gene) {
  # Añadimos identificadores
  tt$Symbol = ilmn_to_gene[rownames(tt), "Symbol"]
  tt$Entrez = ilmn_to_gene[rownames(tt), "Entrez_Gene_ID"]
  # Nos quedamos con los genes significativos y sin NA
  # Eliminamos duplicados en base a Entrez y Symbol
  # Nos quedamos con la fila con menor p-valor
  tt_filt = tt %>%
    filter(adj.P.Val < 0.05, !(is.na(Entrez))) %>%
    arrange(adj.P.Val) %>%
    distinct(Entrez, .keep_all = T) %>% 
    distinct(Symbol, .keep_all = T)
  return(tt_filt)
}

tt1_filt = cambiar_ident(tt1, ilmn_to_gene)
# 5
tt2_filt = cambiar_ident(tt2, ilmn_to_gene)
# 8897
tt3_filt = cambiar_ident(tt3, ilmn_to_gene)
# 7489
tt4_filt = cambiar_ident(tt4, ilmn_to_gene)
# 8189
tt5_filt = cambiar_ident(tt5, ilmn_to_gene)
# 387

#### Diagrama de Venn #### 
genes1 = tt1_filt$Symbol
genes2 = tt2_filt$Symbol
genes3 = tt3_filt$Symbol
genes4 = tt4_filt$Symbol
genes5 = tt5_filt$Symbol

x = list("Endometrio (P vs C)" = genes1, "Lesión profunda" = genes2, "Lesión ovárica" = genes3, 
         "Lesión peritoneal" = genes4, "Peritoneo (P vs C)" = genes5)
ggvenn(x)


### UpSet ### 
png("Results/03_DEA_limma/Figura_2_UpSet.png", 
    width = 2400, height = 1800, res = 300)
upset(fromList(x), 
      nsets = 5,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.6, 0.4),
      main.bar.color = "black",
      sets.bar.color = "steelblue",
      text.scale = 2,
)
dev.off()

# Guardar datos  ----------------------------------------------------------

save(tt1, tt1_filt, tt2_filt, tt3_filt, tt4_filt, tt5_filt, genes1, genes2, genes3,
     genes4, genes5, ilmn_to_gene, file = paste0(results_folder, "result_dea.RData"))



