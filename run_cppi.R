#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Du Yang
# @Date: 2025-01-06 15:17:01
# @LastEditTime: 2025-01-09 18:20:37
# @LastEditors: Du Yang
# @Description: A single function to perform Context-specific PPI analysis
#

library(rio)
library(tidyverse)

# source functions
source("/mnt/dellfs/home/duyang/nsfc2021/00.scripts/cppi_helper.R")

# parameters
# 1. queried gene, e.g. "NDRG1"
# 2. cancer context, e.g. "LUAD"
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 7) {
    stop("Please provide 7 arguments: 
        queried_gene, cancer_context, data_source, mut_data, mrna_data, protein_data, phsopho_data")
}

queried_gene <- args[1]
cancer_context <- args[2]
data_source <- args[3]
mut_data <- args[4]
mrna_data <- args[5]
protein_data <- args[6]
phsopho_data <- args[7]

# demo
queried_gene <- "NDRG1"
cancer_context <- "LUAD"
data_source <- "TCGA-LUAD"
mut_data <- "/mnt/dellfs/pub/data/LinkedOmics/2023/TCGA-LUAD/Human__TCGA_LUAD__WUSM__Mutation__GAIIx__01_28_2016__BI__Gene__Firehose_MutSig2CV.cbt"
mrna_data <- "/mnt/dellfs/pub/data/LinkedOmics/2023/TCGA-LUAD/rna.exp_gene_revised.txt"
# protein_data <- "/mnt/dellfs/pub/data/LinkedOmics/2023/CPTAC-LUAD/protein.exp_gene_revised.txt"
protein_data <- NULL
phospho_data <- NULL


# config
output_dir <- paste0("/mnt/dellfs/home/duyang/nsfc2021/02.result/06.cppi_network/")

# export to json
net_final <- cppi(
    queried_gene = queried_gene,
    cancer_context = cancer_context,
    data_source = data_source,
    mut_data = mut_data,
    mrna_data = mrna_data,
    protein_data = protein_data,
    phospho_data = phospho_data
)


export(net_final, file = paste0(output_dir,"/",cancer_context, "/", queried_gene,".json"))
