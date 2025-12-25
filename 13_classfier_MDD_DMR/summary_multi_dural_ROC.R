df <- read.csv("D:/MDD_result/Supp_Materials/AUC_summary_table.csv")
df <- df %>% mutate(FeatureModel = paste(Method, Features, sep = "_"))
anova_result <- aov(AUC ~ Group, data = df)
summary(anova_result)
anova_summary <- summary(anova_result)[[1]]
p_value <- anova_summary["Pr(>F)"][1, 1]
if (p_value < 0.05) {
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)
}
p <- ggplot(df, aes(x = FeatureModel, y = AUC, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 14) +
  labs(title = "AUC by Feature+Model and Group", x = "Feature + Model", y = "AUC") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 1))
print(p)