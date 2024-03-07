library(tidyverse)
library(magrittr)
library(randomForest)
library(caret)
library(caTools)
library(matrixStats)
library(pROC)
library(ggthemes)

setwd('/media/marcin/data2/Divya/')

df <- read_csv('/media/marcin/data2/Divya/divya_nmd.csv')

df<- art_rf


df_filtered <- df %>%
  dplyr::select(-c(1,2)) %>% #sequence is useless for classifier (every record has different sequence) and log2FC cannot be used because whole your classification is based on this value
  mutate(Solubility_AA1 = as.numeric(Solubility_AA1),
         Solubility_AA2 = as.numeric(Solubility_AA2)) ##--I'm not sure if I should do this since it changes freely to NA. 

nearZeroVar(df_filtered, saveMetrics = TRUE) #no features with near zero variance

na_vals <- data.frame(matrix(nrow = 0, ncol = 3))

for (i in seq_along(df_filtered)) {
  tmp = data.frame(col = colnames(df_filtered)[i],
                   no_na = table(is.na(df_filtered[,i]))[1],
                   na = table(is.na(df_filtered[,i]))[2])
  na_vals <- rbind(na_vals, tmp)
} ##-- This is determining number of NA vals to no NA vals 

na_vals <- na_vals %>%
  replace(is.na(.), 0) 
#pkR features have more NAs than values, may be confusing for classifier

df_filtered <- df_filtered %>%
  dplyr::select(-c(starts_with('pkR'))) %>% #deleting pkR values
  mutate_if(is.numeric , replace_na, replace = 0) %>% #replace NA with 0 in numeric columns
  mutate_if(is.character , replace_na, replace = 'None') #replace NA in categorical features as 'None' string
##--Can I do this (replace NA as 0) for the DNA FC too? or should I just delete it?

tmp_num <- df_filtered %>%
  select_if(., is.numeric) #subdataframe with numeric features

base_cor <- cor(tmp_num)
#we see that pkA values from carboxy and amino groups are highly correlated, 
#lets keep only carboxy
#also keep only the log RF rate
##--what is the cut-off for high correlation

df_filtered_final <- df_filtered %>%
  dplyr::select(-c(starts_with('pKaN'))) %>%
  dplyr::select(-c(starts_with('RF'))) %>%
  dplyr::select(c(ncol(.), 1:(ncol(.)-1))) %>%
  mutate_if(is.character, as.factor) #all categorical vars have to be transformed to factors

set.seed(109) # setting random seed
split <- sample.split(df_filtered_final$NMD, SplitRatio = 0.8) # splitting to train and test (validation) set (80% for training, 20% for testing)
train <- as.data.frame(subset(df_filtered_final, split == "TRUE"))
test <- as.data.frame(subset(df_filtered_final, split == "FALSE"))

##--Cross validation v bootstrap, grid method v. random method, how do you know which will perform best for the model?
trControl <- trainControl(method = "cv", #setting 5x train-test cross validation in training set 
                          number = 5,
                          search = "grid")
#hyperparameters tuning - adjusting 3 hyperparameters: mtry, maxnode number and tree number
tuneGrid <- expand.grid(.mtry = c(1:10)) 

rf_mtry <- train(NMD~.,
                 data = train,
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 nodesize = 14,
                 ntree = 300)
best_mtry <- rf_mtry$bestTune$mtry

store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(5:30)) {
  set.seed(109)
  rf_maxnode <- train(NMD~.,
                      data = train,
                      method = "rf",
                      metric = "Accuracy",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)

maxnode <- results_mtry$values %>%
  dplyr::select(grep('Accuracy', colnames(.))) %>%
  as.matrix() %>%
  colMedians(.)

best_maxnode <- which(maxnode == max(maxnode))[1]                

tree_n <- c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)
store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(109)
  rf_maxtrees <- train(NMD~.,
                       data = train,
                       method = "rf",
                       metric = "Accuracy",
                       tuneGrid = tuneGrid,
                       trControl = trControl,
                       importance = TRUE,
                       nodesize = 14,
                       maxnodes = best_maxnode,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
tree <- results_tree$values %>%
  dplyr::select(grep('Accuracy', colnames(.))) %>%
  as.matrix() %>%
  colMedians(.)

best_tree <- tree_n[which(tree == max(tree))[1]]

#training classifier with previously optimized hyperparameters
fit_rf <- train(NMD~.,
                train,
                method = "rf",
                metric = "Accuracy",
                tuneGrid = tuneGrid,
                trControl = trControl,
                importance = TRUE,
                nodesize = 14,
                ntree = best_tree,
                maxnodes = best_maxnode)
saveRDS(fit_rf, file = 'rf_model_art.Rds')
#prediction on test (validation) dataset
pred <- predict(fit_rf, test, type = 'prob', norm.votes = TRUE) #prediction with tree votes
pred_raw <- predict(fit_rf, test) #binary classification only
conf_mtx <- confusionMatrix(pred_raw, test$NMD, mode = 'everything') #confusion matrix

conf_mtx$table %>% 
  as_tibble() #%>%
  #write_csv('confusion_matrix_arti_noRF.csv')

conf_mtx$byClass %>%
  enframe()# %>%
  #write_csv('performance_arti_noRF.csv', col_names = FALSE)

importance <- varImp(fit_rf, useModel = TRUE) #feature importance
importance_df <- importance$importance %>%
  arrange(desc(`nmd escape`)) %>%
  dplyr::select(-c(2)) %>%
  rownames_to_column('feature') %>%
  dplyr::rename('importance' = 'nmd escape')
#write_csv(importance_df, 'feature_importance_arti_noRF.csv')

ROC_rf <- roc(test$NMD, pred[,2])
saveRDS(ROC_rf, file = "ROC_rf_arti.rds")

ROC_rf_arti<- readRDS("ROC_rf_arti.rds")
ROC_rf_art<- readRDS("ROC_rf_art.rds")

ggroc(list(art = ROC_rf_art, arti = ROC_rf_arti), legacy.axes = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "ROC Curves for ART and ARTi") +
  theme_classic()

feature_importance_art_withRF <- read_csv("feature_importance_art_withRF.csv")
feature_importance_arti_withRF <- read_csv("feature_importance_arti_withRF.csv")

feature_art<- ggplot(feature_importance_art_withRF[1:12,], aes(x = importance, y = reorder(feature, +importance))) +
  geom_bar(stat = "identity", fill = "#9f9f9f") +
  labs(x = "Importance", y = "Feature")

feature_arti<- ggplot(feature_importance_arti_withRF[1:12,], aes(x = importance, y = reorder(feature, +importance))) +
  geom_bar(stat = "identity", fill = "#9f9f9f") +
  labs(x = "Importance", y = "Feature")


#Making plots pretty
base_size <- 16
axis_text_rel_size = -1 
title_text_rel_size = +2

feature_arti+
  theme_foundation(
    base_size = base_size,
    base_family = "sans"
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )+
  theme(
    axis.line = element_line(colour="black", 
                             linewidth = 1),
    axis.ticks = element_line(colour="black", 
                              linewidth = 1)
  )+
  theme(
    text = element_text(colour = "black"),
    plot.title = element_text(
      face = "bold",
      size = rel((title_text_rel_size + base_size) / base_size),
      hjust = 0.5
    ),
    axis.title = element_text(face = "bold", lineheight = rel(1)),
    axis.title.y = element_text(angle = 90, vjust = 2),
    axis.title.x = element_text(vjust = -0.2),
    axis.text = element_text(face = "bold", 
                             size = rel((axis_text_rel_size +
                                           base_size) / base_size
                             )),
    axis.text.x = element_text(
      #angle = 45,
      hjust = 1,
      vjust = 1
    )
  )+
  scale_x_continuous(limits = c(0, 100), 
                     expand = expansion(mult = c(0, 0))
  )
  