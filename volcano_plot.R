#Code for Volcano Plot

#load the required dataframes
setwd("/Users/egeezen/Desktop/PhD/2024 First Term/DNA Methylation/Merged_samples/three_of_them_separately///")
DMPs_Healthy_vs_RA <- read.csv("after_x_y_removal/DMPs_Healthy_vs_RA.csv")
final_beta_values <- read.csv("beta_values_illumina_all_merged.csv")
rownames(final_beta_values) <- final_beta_values$Row.names
final_beta_values$Row.names <- NULL
final_beta_values <- final_beta_values[,-43] #getting rid of the outlier

#volcano plot for but instead of logFC in y axis, delta average beta ? 
#remove x and y chromosome loci from the final beta values df
final_beta_values_filtered <- final_beta_values[rownames(final_beta_values) %in% DMPs_Healthy_vs_RA$probeID,]
#getting healthy and ra from the beta values dataframe
healthy_volcano <- final_beta_values_filtered[, grepl("^Healthy", colnames(final_beta_values_filtered))]
ra_volcano <- final_beta_values_filtered[, grepl("^RA", colnames(final_beta_values_filtered))]
# Calculate the averages of healthy and ra, make sure that it is a dataframe
avg_all_ra <- rowMeans(as.data.frame(ra_volcano)) #rownames is a dataframe function, I don't know but sometimes it doesn't work(probably about tibbles)
avg_all_healthy <- rowMeans(as.data.frame(healthy_volcano))
#differences are calculated
difference_volcano <- avg_all_ra - avg_all_healthy
#now add delta b to the dmps_healthy_vs_ra
matching <- match(DMPs_Healthy_vs_RA$probeID, names(difference_volcano))
DMPs_Healthy_vs_RA$deltaB <- difference_volcano[matching]

#volcano plot drawing starts here:
#add a column to classify points as significant or not
DMPs_Healthy_vs_RA$Classification <- with(DMPs_Healthy_vs_RA, ifelse(adj.P.Val < 0.05 & deltaB > 0.2, "Hypermethylated", 
                                                                     ifelse(adj.P.Val < 0.05 & deltaB < -0.2, "Hypomethylated", 
                                                                            "Not Significant")))

#or for unique to RA, we should get sig_loci exclusive to the unique RA
#DMPs_Healthy_vs_RA$Classification <- with(DMPs_Healthy_vs_RA, ifelse(adj.P.Val < 0.05 & probeID %in% sig_loci & deltaB > 0.2, "Hypermethylated", 
#                                                                     ifelse(adj.P.Val < 0.05 & probeID %in% sig_loci & deltaB < 0, "Hypomethylated", 
#                                                                            "Not Significant")))

#calculating the -log10(padj)
DMPs_Healthy_vs_RA$log_adjP <- -log10(DMPs_Healthy_vs_RA$adj.P.Val)

#selecting points for geom text label
#label_data <- DMPs_Healthy_vs_RA[c(1,13,25,42,15,3,178,713,1210,12291,11233,37706,33675,23938,36695),,drop=FALSE]
label_data <- subset(DMPs_Healthy_vs_RA, (deltaB > 0.2 | deltaB < -0.2) & Classification != "Not Significant" )
rownames(label_data) <- NULL
label_data <- label_data[c(558,696,679,628,133,177,430,357,231,640,398,187),]

ggthemr('fresh')
ggplot(DMPs_Healthy_vs_RA, aes(x = deltaB, y = log_adjP)) +
  geom_point(aes(color = Classification,alpha=Classification), size = 1) + #alpha 0.5-0.8 looks "cool" if there is no alpha aes
  #if you want to use alpha as an aesthetic you have to put alpha=Classification so scale_alpha_manual could be used
  scale_alpha_manual(values=c(c("Hypermethylated" = 1, "Hypomethylated" = 1, "Not Significant" = 0.25))) + # if you want to add transparency
  scale_color_manual(values = c("Hypermethylated" = "#EE4E4E", "Hypomethylated" = "#2A629A", "Not Significant" = "lightgray")) +
  geom_vline(xintercept = c(-0.2,0.2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot of Differential Methylation Unique to RA",
    x = "Delta Beta (Δβ)",
    y = "-log10(Adjusted P-Value)",
    fill = "Classification"
  ) +   # Increase plot margins if necessary
  ylim(0,NA) +
  geom_label_repel(
    data = label_data,
    aes(label = genesUniq),
    box.padding = unit(0.1, "lines"),
    segment.color = "#2B2A4C",
    arrow=NULL,
    max.overlaps= 20
  ) +
  guides(color = guide_legend(title = "Methylation Level for \n RA Unique Loci", override.aes = list(size = 3)),alpha= FALSE)
