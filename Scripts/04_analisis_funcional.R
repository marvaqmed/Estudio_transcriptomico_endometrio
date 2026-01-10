# Libraries
library(R.utils)
library(dplyr)
library(stringr)
library(pheatmap)
library(limma)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(STRINGdb)

# Preparación escritorio 
rm(list=ls())
setwd("../Estudio_transcriptómico_endometrio/")
results_folder = "Results/04_analisis_funcional/"
dir.create(results_folder)

# Carga de datos
load("Results/03_DEA_limma/result_dea.RData", verbose = T)

# GSEA todos los genes ------------------------------------------------------

tt1$Symbol = ilmn_to_gene[rownames(tt1), "Symbol"]
tt1$Entrez = ilmn_to_gene[rownames(tt1), "Entrez_Gene_ID"]

tt1_gsea = tt1 %>%
  filter(!(is.na(Entrez)), !(is.na(Symbol))) %>%
  arrange(logFC) %>%
  distinct(Symbol, .keep_all = T) %>% 
  distinct(Entrez, .keep_all = T)

dim(tt1_gsea)
# 14499     8

geneList = tt1_gsea$logFC
names(geneList) = tt1_gsea$Symbol
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)
tail(geneList)

#### GO ####
gsea_endom = gseGO(geneList = geneList,
                   ont = "BP", 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "SYMBOL", 
                   verbose = T, 
                   minGSSize = 10, 
                   maxGSSize = 500, 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05, 
                   seed = 123)

# Dotplot
png(filename = paste0(results_folder, "Figura_3_dotplot_GO.png"),
    width = 2400, height = 1800, res = 300)
dotplot(gsea_endom, showCategory=8, orderBy = "p.adjust", split=".sign") + 
  facet_grid(.~.sign) + 
  theme(axis.text.y = element_text(size = 12))
dev.off()

# Gseaplot inmunidad mediada por leucocitos
png(filename = paste0(results_folder, "Figura_4A_gseaplot_inmune_leuco.png"),
    width = 2400, height = 1800, res = 300)
gseaplot(gsea_endom, geneSetID = "GO:0002443", title = "Inmunidad mediada por leucocitos")
dev.off()

# Gseaplot división celular
png(filename = paste0(results_folder, "Figura_4A_gseaplot_cell_div.png"),
    width = 2400, height = 1800, res = 300)
gseaplot(gsea_endom, geneSetID = "GO:0051301", title = "División celular")
dev.off()

# Cnetplot
png(filename = paste0(results_folder, "Figura_5_cnetplot.png"),
    width = 2400, height = 1800, res = 300)
cnetplot(gsea_endom, 
         showCategory = 10, 
         foldChange = geneList, 
         circular = TRUE,       
         color_category = "firebrick",
         node_label = "category")
dev.off()


#### KEGG #### 
geneList_kegg = tt1_gsea$logFC
names(geneList_kegg) = tt1_gsea$Entrez
geneList_kegg = sort(geneList_kegg, decreasing = TRUE)
gsea_kegg = gseKEGG(geneList = geneList_kegg,
                    organism = "hsa",
                    keyType = "ncbi-geneid", 
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH",
                    seed = 123)

png(filename = paste0(results_folder, "Figura_6_dotplot_kegg.png"),
    width = 2400, height = 1800, res = 300)
dotplot(gsea_kegg, showCategory=8, orderBy = "p.adjust", split = ".sign") + 
  facet_grid(.~.sign) + 
  theme(axis.text.y = element_text(size = 12))
dev.off()

pathview(gene.data  = geneList_kegg,     
         pathway.id = "hsa04110",  
         species    = "hsa",        
         limit      = list(gene=max(abs(geneList_kegg)), cpd=1))


# ORA ---------------------------------------------------------------------

# Seleccionamos los 4733 genes comunes y exclusivos de la intersección
sel_genes2 = Reduce(intersect, list(genes2, genes3, genes4))
sel_genes2 = sel_genes2[!sel_genes2 %in% c(genes1, genes5)]

#### Heatmap ####
# Comprobar que los DEGs tienen la misma direccionalidad en todas las lesiones
sel_tt2 = tt2_filt[tt2_filt$Symbol %in% sel_genes2, c("Symbol", "Entrez", "logFC")]
rownames(sel_tt2) = sel_tt2$Symbol
colnames(sel_tt2)[3] = "logFC_deep"
sel_tt3 = tt3_filt[tt3_filt$Symbol %in% sel_genes2, c("Symbol", "Entrez", "logFC")]
rownames(sel_tt3) = sel_tt3$Symbol
sel_tt3 = sel_tt3[rownames(sel_tt2),]
colnames(sel_tt3)[3] = "logFC_ovar"
sel_tt4 = tt4_filt[tt4_filt$Symbol %in% sel_genes2, c("Symbol", "Entrez", "logFC")]
rownames(sel_tt4) = sel_tt4$Symbol
sel_tt4 = sel_tt4[rownames(sel_tt2),]
colnames(sel_tt4)[3] = "logFC_perit"
all(rownames(sel_tt2) == rownames(sel_tt3) & 
      rownames(sel_tt3) == rownames(sel_tt4))
# TRUE 

toplot = merge(sel_tt2, sel_tt3, by="Symbol")
toplot = merge(toplot, sel_tt4, by= "Symbol")
rownames(toplot) = toplot$Symbol
toplot = toplot[, c("Symbol", "Entrez.x", "logFC_deep", "logFC_ovar", "logFC_perit")]
colnames(toplot)[2] = "Entrez"
pheatmap(toplot[-c(1,2)], 
         color = colorRampPalette(c("navy", "white", "red"))(50),
         main = "LogFC de las lesiones", 
         show_rownames = F, 
         angle_col = "0", 
         breaks = seq(-2, 2, length.out = 50))

png(filename = paste0(results_folder, "Figura_7_pheatmap_logFC.png"),
    width = 2400, height = 1000, res = 300)
pheatmap(toplot[-c(1,2)], 
         color = colorRampPalette(c("navy", "white", "red"))(50),
         main = "LogFC de las lesiones", 
         show_rownames = F, 
         angle_col = "0", 
         breaks = seq(-2, 2, length.out = 50))
dev.off()


# Hacemos un ORA con los genes significativos en los tres tipos de lesiones ectópicas
# Antes sacamos el logFC promedio de las tres lesiones 
df_logFC = cbind(toplot, 
                 rowMeans(toplot[-c(1,2)]))
colnames(df_logFC)[6] = "logFC_prom"

df_logFC_up = df_logFC[df_logFC$logFC_prom > 0, ]
df_logFC_down = df_logFC[df_logFC$logFC_prom < 0, ]


#### GO ####
# Up-regulated (activadas)
sel_genes2_up  = df_logFC_up$Entrez
ora_go_up = enrichGO(gene = sel_genes2_up, 
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID", 
                     ont = "BP")

png(filename = paste0(results_folder, "Figura_8A_ora_go_up.png"),
    width = 2400, height = 1800, res = 300)
barplot(ora_go_up, x = "GeneRatio", showCategory = 10) + 
  theme(axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("ORA GO Up-regulated")
dev.off()



# Down-regulated (suprimidas)
sel_genes2_down = df_logFC_down$Entrez
ora_go_down =  enrichGO(gene = sel_genes2_down, 
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID", 
                        ont = "BP")

png(filename = paste0(results_folder, "Figura_8B_ora_go_down.png"),
    width = 2400, height = 1800, res = 300)
barplot(ora_go_down, x = "GeneRatio", showCategory = 10)  + 
  theme(axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("ORA GO Down-regulated")
dev.off()


#### KEGG ####
# Up-regulated
ora_kegg_up = enrichKEGG(gene = sel_genes2_up, 
                         organism = "hsa", 
                         pvalueCutoff = 0.05)

png(filename = paste0(results_folder, "Figura_8C_ora_kegg_up.png"),
    width = 2400, height = 1800, res = 300)
barplot(ora_kegg_up, x = "GeneRatio", showCategory = 10) + 
  theme(axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("ORA KEGG Up-regulated")
dev.off()


# Down-regulated
ora_kegg_down = enrichKEGG(gene = sel_genes2_down, 
                           organism = "hsa", 
                           pvalueCutoff = 0.05)
png(filename = paste0(results_folder, "Figura_8D_ora_kegg_down.png"),
    width = 2400, height = 1800, res = 300)
barplot(ora_kegg_down, x = "GeneRatio", showCategory = 10) + 
  theme(axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("ORA KEGG Down-regulated")
dev.off()


# String ------------------------------------------------------------------

# Cogemos los 100 primeros genes up-regulated 
data_string = df_logFC_up %>% 
  arrange(desc(logFC_prom)) %>%
  dplyr::slice(1:100) %>% 
  dplyr::select(Symbol, logFC_prom)

string_db = STRINGdb$new(version = "12", species = 9606,
                         score_threshold = 400)
mapped = string_db$map(data_string, "Symbol", removeUnmappedRows = TRUE)
hits <- mapped$STRING_id
string_db$plot_network(hits)


# Guardar ----------------------------------------------------------------

save(gsea_endom, gsea_kegg, ora_go_up, ora_go_down, ora_kegg_up, ora_kegg_down,
     file = paste0(results_folder, "gsea.RData"))



