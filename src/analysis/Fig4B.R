df <- metadata %>% select(c(survivor, clusters, propCancer, aa_richness, cancerTCRsPerCluster))
df[is.na(df)] <- 0
df$clusters <- df$clusters/mean(df$clusters)
df$propCancer <- df$propCancer/mean(df$propCancer)
df$aa_richness <- df$aa_richness/mean(df$aa_richness)
df$cancerTCRsPerCluster <- df$cancerTCRsPerCluster/mean(df$cancerTCRsPerCluster)

df_long <- df %>%
  tidyr::gather(variable, value, -survivor)

df_long$variable = factor(df_long$variable, levels = c("aa_richness", "propCancer", 
                                                       "clusters", "cancerTCRsPerCluster"), ordered = TRUE)
ggplot(df_long, aes(variable, 
                    value, fill = survivor)) + 
  geom_boxplot(width = 0.4) + 
  scale_fill_brewer(palette = "Set1") +
  theme_classic2() +
  scale_y_log10() +
  labs(x = "", y = "Relative Magnitude") +
  scale_fill_manual(values = c("#FF9F80", "#C1E1C5"), name = "", 
                    labels = c("Short Term Surivors", "Long Term Survivors")) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face="bold"),
        axis.text.x = element_text()) +
  scale_x_discrete(labels = c("Richness", "Proportion TCRs\nCancer Associated", "Total\nClusters", "Cancer TCRs \nPer Cluster"))
