#to draw venn diagram
library(VennDiagram)
library(readxl)
library(ggplot2)

#load your dataframes for significant loci in RA or unique to RA (compared to OA)
setwd("/Users/egeezen/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/three_of_them_separately///")
DMPs_Healthy_vs_RA_sig_all <- read.csv("after_x_y_removal/DMPs_Healthy_vs_RA_all_RA_with_deltaB.csv")
DMPs_Healthy_vs_RA_sig_unique_RA <- read.csv("after_x_y_removal/DMPs_Healthy_vs_RA_unique_to_RA_with_deltaB.csv")

# from the excel file 
DMPs_Healthy_vs_RA_sig_unique_RA <- read_excel("after_x_y_removal/delta beta > 0.1/attributes_genes_related_to_DMPs_sig_with_delta_B_unique_to_RA.xlsx")
# get unique hgnc_symbol entries
DMPs_Healthy_vs_RA_sig_unique_RA <- DMPs_Healthy_vs_RA_sig_unique_RA[!duplicated(DMPs_Healthy_vs_RA_sig_unique_RA$hgnc_symbol), ]

# load the risk loci attributes dataframe
risk_loci_att <- read_excel("~/Desktop/PhD/2024 First Term/risk loci genes/risk_loci_with_attributes_stimulation_and_expression.xlsx")

# active in SF risk RA genes
risk_loci_att <- risk_loci_att[(risk_loci_att$change_in_expression_upon_stimulation != "No") | (risk_loci_att$expression_in_fibroblasts != "No"),]
risk_loci_att <- risk_loci_att[(risk_loci_att$change_in_expression_upon_stimulation != "Not Known") | (risk_loci_att$expression_in_fibroblasts != "No"),]

set_1 <- unique(risk_loci_att$hgnc_symbol)
set_2 <- unique(DMPs_Healthy_vs_RA_sig_unique_RA$hgnc_symbol)

venn.plot <- venn.diagram(
  x = list(set1 = set_1, set2 = set_2),
  category.names = c("Risk Loci\nAssociated\nGenes", "(Unique to) RA vs Healthy\nDMPs related Genes"), filename=NULL,
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = c("#FFCBCB","#3C5B6F"),       
  cex = 1.0,
  fontface = "bold",
  fontfamily = "sans",       
  cat.cex = 0.9,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c( -170, 160 ),
  cat.dist = c(0.04,0.033),
  disable.logging = TRUE
)

# Plot Venn diagram
grid.draw(venn.plot)
ggsave(
  plot = venn.plot,
  filename = "vennplot_dmps_vs_activesf_risk.png",
  bg = "transparent",
  width = 5, height = 5
)

dev.off()
