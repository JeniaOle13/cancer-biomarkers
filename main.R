setwd("/home/acari/Documents/github/cancer-biomarkers/")

library(ggplot2)
library(ggpubr)
library(phyloseq)
library(MicrobiotaProcess)
library(vegan)
library(stringr)
library(Maaslin2)
library(RColorBrewer)
library(lmerTest)
library(multcomp)
library(clusterProfiler)

library(pROC)
library(caret)
library(ROSE)
library(tidyverse)
library(effectsize)

# import data
otu_table <- read.csv("data/mags_relab.tsv", sep = "\t", row.names = 1)

taxonomy <- read.csv("data/taxonomy.tsv", sep = "\t", row.names = 1)

kingdom <- sapply(str_split(taxonomy$Taxon, ";"), function(x) x[1])
phylum <- sapply(str_split(taxonomy$Taxon, ";"), function(x) x[2])
class_t <- sapply(str_split(taxonomy$Taxon, ";"), function(x) x[3])
order_t <- sapply(str_split(taxonomy$Taxon, ";"), function(x) x[4])
family_t <- sapply(str_split(taxonomy$Taxon, ";"), function(x) x[5])
genus <- sapply(str_split(taxonomy$Taxon, ";"), function(x) x[6])
genus[genus == "g__"] <- paste0(family_t[which(genus == "g__")], "_unknown")

phylum <- sub("p__", "", phylum)
genus <- sub("g__", "", genus)

tax_table <- data.frame(kingdom, phylum, class = class_t, order = order_t, family = family_t, genus)
rownames(tax_table) <- rownames(taxonomy)

meta_data <- read.csv("data/metadata.tsv", sep = "\t", row.names = 1)

distance_matrics <- read.csv("data/rpca/distance-matrix.tsv", sep = "\t", row.names = 1)
ordination <- read.csv("data/rpca/ordination.txt", sep = "\t")

mag.biomarkers <- read.csv("data/mag.biomarkers.tsv", sep = "\t")

songbird_Frankel_2017 <- read.csv("data/songbird/Frankel_2017/differentials.tsv", sep = "\t")
songbird_Gopalakrishnan_2019 <- read.csv("data/songbird/Gopalakrishnan_2019/differentials.tsv", sep = "\t")
songbird_Matson_2019 <- read.csv("data/songbird/Matson_2019/differentials.tsv", sep = "\t")

songbird_Spencer_2021 <- read.csv("data/songbird/Spencer_2021/differentials.tsv", sep = "\t")
songbird_Lee_2022 <- read.csv("data/songbird/Lee_2022/differentials.tsv", sep = "\t")

songbird_Liu_2022 <- read.csv("data/songbird/Liu_2022/differentials.tsv", sep = "\t")
songbird_McCulloch_2022 <- read.csv("data/songbird/McCulloch_2022/differentials.tsv", sep = "\t")
songbird_Peng_2020 <- read.csv("data/songbird/Peng_2020/differentials.tsv", sep = "\t")
songbird_Tsakmaklis_2023 <- read.csv("data/songbird/Tsakmaklis_2023/differentials.tsv", sep = "\t")
songbird_Gunjur_2024 <- read.csv("data/songbird/Gunjur_2024/differentials.tsv", sep = "\t")
songbird_Heshiki_2020 <- read.csv("data/songbird/Heshiki_2020/differentials.tsv", sep = "\t")

body_site <- read.csv("data/body_ogu.tsv", sep = "\t")[-2]
colnames(body_site)[2] <- "site"
body_site$genome <- sub(".fa", "", body_site$genome)

food <- read.csv("data/food_ogu.txt", sep = "\t", header = F)
colnames(food) <- "genome"
food$site <- "food"
food$genome <- sub(".fa", "", food$genome)

# Genus plot
meta_data_d <- meta_data
rownames(meta_data_d) <- meta_data_d$sampleid
meta_data_d <- meta_data_d[-1]

# Alpha diversity
shannon <- diversity(otu_table, index = "shannon")
shannon <- as.data.frame(shannon)
shannon <- cbind(rownames(shannon), shannon)

alpha.df <- merge(cbind(sampleid = rownames(meta_data), meta_data), shannon, by = 1)

shannon_boxplot <- ggplot(alpha.df)+
    geom_boxplot(mapping = aes(shannon, response), outlier.colour = "white")+
    geom_jitter(mapping = aes(shannon, response, col = response))+
    theme_classic()+
    theme(legend.position = "none")+
    scale_color_brewer(palette = "Set1")+
    xlab("Shannon index")+
    ylab("Dataset")

lmer.shannon <- lmer(shannon ~ response + (1|dataset), data=alpha.df)
summary(glht(lmer.shannon, linfct = mcp(response = "Tukey")), test = adjusted("bonferroni"))
shapiro.test(residuals(lmer.shannon))

# Beta diversity
ordination <- merge(cbind(sampleid = rownames(meta_data), meta_data), ordination, by = 1)

mds_plot <- ggplot(ordination)+
    geom_point(mapping = aes(PC1, PC2, shape = response), size = 3.5, stroke = 1)+
    geom_point(mapping = aes(PC1, PC2, col = response, shape = response), size = 3)+
    theme_classic()+
    stat_ellipse(mapping = aes(PC1, PC2, col = response))+
    scale_color_brewer(palette = "Set1")+
    theme(legend.position = "right")+
    xlab("PC1 (54%)")+
    ylab("PC2 (34%)")

adonis2(distance_matrics ~ dataset + response + Cancer.Type, data = meta_data, permutations = 999, by = "margin")

# Songbird biomarkers
songbird_Frankel_2017_sbs <- songbird_Frankel_2017[songbird_Frankel_2017$featureid %in% mag.biomarkers$featureid,]
songbird_Frankel_2017_sbs$dataset <- "Frankel_2017"

songbird_Gopalakrishnan_2019_sbs <- songbird_Gopalakrishnan_2019[songbird_Gopalakrishnan_2019$featureid %in% mag.biomarkers$featureid,]
songbird_Gopalakrishnan_2019_sbs$dataset <- "Gopalakrishnan_2019"

songbird_Matson_2019_sbs <- songbird_Matson_2019[songbird_Matson_2019$featureid %in% mag.biomarkers$featureid,]
songbird_Matson_2019_sbs$dataset <- "Matson_2019"

songbird_Spencer_2021_sbs <- songbird_Spencer_2021[songbird_Spencer_2021$featureid %in% mag.biomarkers$featureid,]
songbird_Spencer_2021_sbs$dataset <- "Spencer_2021"

songbird_Lee_2022_sbs <- songbird_Lee_2022[songbird_Lee_2022$featureid %in% mag.biomarkers$featureid,]
songbird_Lee_2022_sbs$dataset <- "Lee_2022"

songbird_Liu_2022_sbs <- songbird_Liu_2022[songbird_Liu_2022$featureid %in% mag.biomarkers$featureid,]
songbird_Liu_2022_sbs$dataset <- "Liu_2022"

songbird_McCulloch_2022_sbs <- songbird_McCulloch_2022[songbird_McCulloch_2022$featureid %in% mag.biomarkers$featureid,]
songbird_McCulloch_2022_sbs$dataset <- "Liu_2022"

songbird_Peng_2020_sbs <- songbird_Peng_2020[songbird_Peng_2020$featureid %in% mag.biomarkers$featureid,]
songbird_Peng_2020_sbs$dataset <- "Peng_2020"

songbird_Tsakmaklis_2023_sbs <- songbird_Tsakmaklis_2023[songbird_Tsakmaklis_2023$featureid %in% mag.biomarkers$featureid,]
songbird_Tsakmaklis_2023_sbs$dataset <- "Tsakmaklis_2023"

songbird_Gunjur_2024_sbs <- songbird_Gunjur_2024[songbird_Gunjur_2024$featureid %in% mag.biomarkers$featureid,]
songbird_Gunjur_2024_sbs$dataset <- "Gunjur_2024"

songbird_Heshiki_2020_sbs <- songbird_Heshiki_2020[songbird_Heshiki_2020$featureid %in% mag.biomarkers$featureid,]
songbird_Heshiki_2020_sbs$dataset <- "Heshiki_2020"

songbird_coef <- rbind(songbird_Frankel_2017_sbs, songbird_Gopalakrishnan_2019_sbs, songbird_Matson_2019_sbs,
      songbird_Spencer_2021_sbs, songbird_Lee_2022_sbs, songbird_Liu_2022_sbs, songbird_McCulloch_2022_sbs,
      songbird_Peng_2020_sbs, songbird_Tsakmaklis_2023_sbs, songbird_Gunjur_2024_sbs, songbird_Heshiki_2020_sbs)

songbird_coef <- songbird_coef %>% group_by(featureid) %>% summarise(mean = mean(response.T.R.), sd = sd(response.T.R.))
songbird_coef <- as.data.frame(songbird_coef)
songbird_coef <- merge(cbind(featureid = rownames(tax_table), tax_table), songbird_coef, by = 1)

songbird_coef$phylum <- sapply(str_split(songbird_coef$phylum, "_"), function(x) x[1])
songbird_coef$genus <- sapply(str_split(songbird_coef$genus, "_"), function(x) x[1])

songbird_coef <- songbird_coef[order(songbird_coef$mean, decreasing = T),]

coef_barplot <- ggplot(songbird_coef, aes(x = reorder(featureid, mean), y=mean)) + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
        position=position_dodge(0.05), col = "gray") +
    # geom_point(shape = 0, size = 0.25)+
    geom_hline(yintercept = 0, col = "red")+
    theme_bw()+
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA))+
    ylab("Mean R/NR + k")

BiocManager::install("biomformat")

df.log_ratio <- NULL
for (i in rownames(ogu_table)){
    df <- otu_table[rownames(ogu_table) %in% i,]
    df <- df[colSums(df) > 0]
    
    df_denominator <- df[colnames(df) %in% songbird_coef$featureid[songbird_coef$mean < 0]]
    df_numerator <- df[colnames(df) %in% songbird_coef$featureid[songbird_coef$mean > 0]]
    
    df_denominator <- df_denominator[df_denominator > 0]
    df_numerator <- df_numerator[df_numerator > 0]
    
    df_denominator <- df_denominator + 0.001
    df_numerator <- df_numerator + 0.001

    log_ratio <- log(sum(df_numerator)/sum(df_denominator))
    
    df.s <- data.frame(sampleid = i, log_ratio)
    df.log_ratio <- rbind(df.s, df.log_ratio)
}

df.log_ratio <- merge(df.log_ratio, cbind(rownames(meta_data), meta_data), by = 1)

lmer.shannon <- lmer(log_ratio ~ response + (1|dataset), data=df.log_ratio)
summary(glht(lmer.shannon, linfct = mcp(response = "Tukey")), test = adjusted("bonferroni"))
shapiro.test(residuals(lmer.shannon))
cohens_d(rank ~ response, data=df.log_ratio)

# Function to perform LOGO-CV
df.log_ratio <- df.log_ratio |>  mutate(response = as.factor(ifelse(response == "R", 1, 0)))
df.log_ratio$dataset <- as.factor(df.log_ratio$dataset)
df.log_ratio$cancer_type <- as.factor(df.log_ratio$cancer_type)

df.log_ratio_centered <- df.log_ratio %>%
    group_by(dataset) %>%
    mutate(log_ratio_centered = log_ratio - median(log_ratio, na.rm = TRUE))) %>%
    ungroup()
df.log_ratio_centered <- as.data.frame(df.log_ratio_centered)

unique_datasets <- unique(df.log_ratio$dataset)
roc_list <- list()
conf_matrix_list <- list()

class_counts <- table(df.log_ratio$response)

df.log_ratio$weights <- ifelse(df.log_ratio$response == "1", 
                  1/class_counts["1"], 
                  1/class_counts["0"])

df.log_ratio_centered$rank <- rank(df.log_ratio_centered$log_ratio_centered)

for (i in seq_along(unique_datasets)) {
    test_dataset <- unique_datasets[i]
    
    train_data <- df.log_ratio_centered %>% filter(dataset != test_dataset)
    test_data <- df.log_ratio_centered %>% filter(dataset == test_dataset)
    
    # Train logistic regression model
    model <- glm(response ~ log_ratio_centered, family = binomial, data = train_data, weights = weights)
    # model <- glm(response ~ log_ratio, family = binomial, data = train_data)
    
    # Predict on test data
    test_data$pred_prob <- predict(model, test_data, type = "response")
    roc_obj <- roc(test_data$response, test_data$pred_prob)
    roc_list[[i]] <- roc_obj$specificities
    
    # Confusion matrix (using 0.5 as threshold)
    test_data$pred_class <- ifelse(test_data$pred_prob >= 0.5, 1, 0)
    conf_matrix <- confusionMatrix(
        as.factor(test_data$pred_class), 
        test_data$response,
        positive = "1"
    )
    conf_matrix_list[[i]] <- conf_matrix
}

roc_combined <- roc_list[[1]]
for (i in 2:length(roc_list)) {
    roc_combined <- roc_combined + roc_list[[i]]
}
roc_mean <- roc_combined / length(roc_list)

conf_metrics <- map_dfr(conf_matrix_list, ~{
    data.frame(
        Accuracy = .x$overall["Accuracy"],
        Sensitivity = .x$byClass["Sensitivity"],
        Specificity = .x$byClass["Specificity"]
    )
}, .id = "Dataset") %>%
    mutate(Dataset = unique_datasets[as.numeric(Dataset)])

mean(conf_metrics$Accuracy)
mean(conf_metrics$Sensitivity)
mean(conf_metrics$Specificity)

logratio_boxplot <- ggplot(df.log_ratio, aes(response, log_ratio, col = response))+
    stat_halfeye()+
    facet_wrap(~dataset)+
    xlab("Response")+
    ylab("Log ratio")+
    theme_pubr()+
    theme(legend.position = "none")+
    scale_color_brewer(palette = "Set1")

# GSEA
vec <- ogu.markers$mean
names(vec) <- ogu.markers$featureid
vec <- sort(vec, decreasing = T)

## By taxonomy
fgsea_phylum <- GSEA(vec, TERM2GENE = ogu.markers[c(8,1)], eps = 0)
fgsea_phylum_df <- as.data.frame(fgsea_phylum)

fgsea_class <- GSEA(vec, TERM2GENE = ogu.markers[c(9,1)], eps = 0)
fgsea_class_df <- as.data.frame(fgsea_class)

fgsea_order <- GSEA(vec, TERM2GENE = ogu.markers[c(10,1)], eps = 0)
fgsea_order_df <- as.data.frame(fgsea_order)

fgsea_family <- GSEA(vec, TERM2GENE = ogu.markers[c(11,1)], eps = 0)
fgsea_family_df <- as.data.frame(fgsea_family)

fgsea_genus <- GSEA(vec, TERM2GENE = ogu.markers[c(12,1)], eps = 0)
fgsea_genus_df <- as.data.frame(fgsea_genus)

gs_species <- ogu.markers[c(13,1)]
gs_species$species <- sub(" ", "_", gs_species$species)

fgsea_species <- GSEA(vec, TERM2GENE = ogu.markers[c(13,1)], eps = 0)
fgsea_species_df <- as.data.frame(fgsea_species)

gseaNb(object = fgsea_phylum,
       geneSetID = fgsea_phylum_df$ID, 
       curveCol = brewer.pal(name = "Set1", n = 9))

gseaNb(object = fgsea_class,
       geneSetID = fgsea_class_df$ID, 
       curveCol = brewer.pal(name = "Set1", n = 9))

gseaNb(object = fgsea_order,
       geneSetID = fgsea_order_df$ID, 
       curveCol = brewer.pal(name = "Set1", n = 9)[-6])

gseaNb(object = fgsea_family,
       geneSetID = fgsea_family_df$ID, 
       curveCol = brewer.pal(name = "Set1", n = 9)[-6])

gseaNb(object = fgsea_genus,
       geneSetID = fgsea_genus_df$ID, 
       curveCol = c(brewer.pal(name = "Set1", n = 9)[-6], "cyan4"))

## By site
fgsea_food <- GSEA(vec, TERM2GENE = food[c(2,1)])
fgsea_bodysite <- GSEA(vec, TERM2GENE = body_site[c(2,1)])

gseaNb(object = fgsea_food, 
       geneSetID = 'food', 
       htCol = c("#A50F15", "#08519C"))

gseaNb(object = fgsea_bodysite,
       geneSetID = 'Oral cavity', 
       htCol = c("#A50F15", "#08519C"))

