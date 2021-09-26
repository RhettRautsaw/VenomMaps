library(readxl)
library(dplyr)
library(ggpubr)

data<-read_xlsx("~/Desktop/PCA_vs_VIF.xlsx")

data<-data %>% filter(present_in_both=="yes")

ggboxplot(data, "variable_selection_method", "auc", fill="variable_selection_method", add="jitter") + stat_compare_means(method="t.test", label.x=1.5)

aggregate(data$auc, list(data$variable_selection_method), median)
