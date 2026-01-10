# Libraries
library(GEOquery)
library(readxl)
library(R.utils)
library(stringr)

# Preparación escritorio 
rm(list=ls())
setwd("../Estudio_transcriptómico_endometrio/")
results_folder = "Results/01_descarga_datos/"
dir.create(results_folder)

# Descarga de datos (GEO) -------------------------------------------------

gcel = getGEOSuppFiles("GSE141549")

# Vemos que ficheros se han descargado
list.files("./GSE141549/")
# [1] "GSE141549_batchCorrectednormalizedArrayscombined.xlsx"                          
# [2] "GSE141549_Illumina_HumanHT-12_V4.0_expression_beadchip_normalized.xlsx"         
# [3] "GSE141549_Illumina_HumanWG-6_v2.0_expression_beadchip_normalized.xlsx"          
# [4] "GSE141549_Non-normalized_data_GA_illumina_expression_HumanWG-6.xls.gz"          
# [5] "GSE141549_Non-normalized_data_GA_illumina_expression_platform_HumanHT-12.xls.gz"
# [6] "GSE141549_RAW.tar"                                                              
# [7] "GSE141549_Sample_link.xlsx"

# Nos interesan los datos no normalizados
# Descomprimimos los ficheros 
gunzip("GSE141549/GSE141549_Non-normalized_data_GA_illumina_expression_HumanWG-6.xls.gz", 
       destname = "GSE141549/GSE141549_Non-normalized_data_GA_illumina_expression_HumanWG-6.xls", 
       overwrite = TRUE, remove = FALSE)
gunzip("GSE141549/GSE141549_Non-normalized_data_GA_illumina_expression_platform_HumanHT-12.xls.gz", 
       destname = "GSE141549/GSE141549_Non-normalized_data_GA_illumina_expression_platform_HumanHT-12.xls", 
       overwrite = TRUE, remove = FALSE)


# Merge datos -------------------------------------------------------------

# Leemos ficheros 
datos_HumanHT = read_excel("GSE141549/GSE141549_Non-normalized_data_GA_illumina_expression_platform_HumanHT-12.xls", 
                           skip = 4)
datos_HumanHT = as.data.frame(datos_HumanHT)
dim(datos_HumanHT)
#  47323   173
ident = unique(datos_HumanHT$ID_REF)
length(ident)
# 40607
rep = which(duplicated(datos_HumanHT$ID_REF))

# Nos quedamos con la primera aparición 
datos_HumanHT_filt = datos_HumanHT[!(rownames(datos_HumanHT) %in% rep), ]
sum(duplicated(datos_HumanHT_filt$ID_REF))
# 0
rownames(datos_HumanHT_filt) = datos_HumanHT_filt$ID_REF

# Los datos de la plataforma HumanWG-6 vienen en 7 hojas excel
# Leemos cada hoja y los vamos juntando en un único objeto
# Pero seleccionando solo las filas con los mismos illumina ID 
join = data.frame("ID_REF" = ident)
rownames(join) = join$ID_REF
for (i in 1:7){
  print(i)
  datos = read_excel("GSE141549/GSE141549_Non-normalized_data_GA_illumina_expression_HumanWG-6.xls", 
                     skip = 4, sheet = i)
  print("Datos leídos")
  datos = as.data.frame(datos)
  rep = which(duplicated(datos$ID_REF))
  na = which(is.na(datos$ID_REF))
  datos_filt = datos[!(rownames(datos) %in% c(rep, na)), ]
  rownames(datos_filt) = datos_filt$ID_REF
  colnames(datos_filt) = paste(colnames(datos_filt), i, sep = "_")
  print("Datos filtrados")
  ids = intersect(rownames(join), datos_filt$ID_REF)
  sel_datos = datos_filt[ids, ]
  join = data.frame(join[ids, ])
  rownames(join) = ids
  if (all(rownames(join) == sel_datos$ID_REF)) {
    join = cbind(join, sel_datos[, -1])
  }
  print("Datos combinados")
}


# Guardar datos -----------------------------------------------------------

save(join, datos_HumanHT_filt, file = paste0(results_folder, "datos_crudos.RData"))
