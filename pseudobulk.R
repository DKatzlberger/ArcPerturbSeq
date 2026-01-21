## max. code width ============================================================

## ============================================================================
## Arc Institute in vivo Perturb-seq
## Author: dkatzlberger
## Date: 21.01.2026
## ============================================================================

## ============================================================================
## Libraries
## ============================================================================

library(anndataR)

## ============================================================================
## Datasets
## ============================================================================

## Single cell pseudobulks per organ
## Organs: Peritoneum, Liver, Spleen
dir_share <- "/vscratch/wes/arc/share"

dir_peritoneum <- file.path(dir_share, "JB22_PM")
dir_liver      <- file.path(dir_share, "JB22_Liver")
dir_spleen     <- file.path(dir_share, "JB22_Spleen")

## Available cell types per organ
## Peritoneum: Macrophages
## Liver:      Bcells, Monocytes, Neutrophiles
## Spleen:     Bcells

peritoneum_macrophages <- read_h5ad(file.path(dir_peritoneum, "pseudobulk_Macrophages.h5ad"))

