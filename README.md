## Estudio transcriptómico en endometrio
Este repositorio contiene el código y flujo de trabajo asociado al trabajo **"Firmas transcriptómicas moleculares de la endometriosis reveladas mediante análisis in sílico"**.

### Descripción
El objetivo del estudio es  analizar, mediante herramientas bioinformáticas, los perfiles de expresión génica del endometrio eutópico y de las lesiones ectópicas
con el fin de identificar así firmas moleculares y funciones biológicas alteradas implicadas en la patogénesis de la endometriosis, a partir de los datos de la
serie GSE141549 de GEO.

### Estructura
- `scripts`: contiene los scripts de R para el flujo de análisis.
  - `01_descarga_datos`: descarga automática de los datos desde GEO y fusión de las plataformas (HT-12 y WG-6).
  - `02_control_calidad`: control de calidad, detección de efectos de tanda, normalización por cuantiles y curación metadata.
  - `03_DEA_limma`: ajuste de modelos lineales utilizando limma para realizar un análisis de expresión diferencial.
  - `04_analisis_funcional`: análisis de enriquecimiento funcional (GSEA y ORA) y construcción de redes de interacción proteína-proteína.

### Requisitos y librerías
Para el análisis se requiere R  y las siguientes librerías principales:
  ```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "NOISeq", "clusterProfiler", "orgs.Hs.eg.db", "STRINGdb"))

install.packages(c("tidyverse", "readxl", "R.utils", "pheatmap", "ggvenn", "UpsetR"))

library(readxl)
library(R.utils)
library(stringr)
library(ggplot2)
library(pheatmap)
library(limma)
library(dplyr)
library(NOISeq)
library(tidyverse)
library(pathview)
library(ggvenn)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
```
