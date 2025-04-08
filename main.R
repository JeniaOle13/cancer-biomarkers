setwd("/home/acari/Documents/github/cancer-biomarkers/")

library(ggplot2)
library(ggpubr)
library(phyloseq)
library(microViz)
library(vegan)
library(stringr)
library(Maaslin2)
library(RColorBrewer)
library(lmerTest)
library(multcomp)
library(clusterProfiler)
library(enrichplot)
library(fgsea)

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
phyobj <- phyloseq(otu_table(otu_table, taxa_are_rows=FALSE),
                        sample_data(meta_data),
                        tax_table(as.matrix(tax_table)))

phyobj_100 <- transform_sample_counts(phyobj, function(x) 100 * x/sum(x))

phylum_barplot <- phyobj_100 %>% comp_barplot("phylum", merge_other = FALSE, label = NULL)+
    facet_wrap(vars(response), scales = "free", nrow = 2)+
    theme_minimal()+
    theme(legend.position = "bottom")

genus_barplot <- phyobj_100 %>% comp_barplot("genus", merge_other = FALSE, label = NULL, n_taxa = 10)+
    facet_wrap(~response, scales = "free", nrow = 2)+
    theme_minimal()+
    theme(legend.position = "bottom")

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


df.log_ratio <- NULL
for (i in rownames(otu_table)){
    df <- otu_table[rownames(otu_table) %in% i,]
    df <- df[colSums(df) > 0]
    
    df_denominator <- df[colnames(df) %in% songbird_coef$featureid[songbird_coef$mean < 0]]
    df_numerator <- df[colnames(df) %in% songbird_coef$featureid[songbird_coef$mean > 0]]
    
    df_denominator <- df_denominator + 0.001
    df_numerator <- df_numerator + 0.001

    log_ratio <- log(sum(df_numerator)/sum(df_denominator))
    
    df.s <- data.frame(sampleid = i, log_ratio)
    df.log_ratio <- rbind(df.s, df.log_ratio)
}

df.log_ratio <- merge(df.log_ratio, cbind(rownames(meta_data), meta_data), by = 1)
df.log_ratio_sbs <- df.log_ratio[df.log_ratio$log_ratio != Inf,]
df.log_ratio_sbs <- df.log_ratio_sbs[df.log_ratio_sbs$log_ratio != -Inf,]

logratio_boxplot <- ggplot(df.log_ratio_sbs)+
    geom_boxplot(mapping = aes(log_ratio, response), outlier.color = "white")+
    geom_jitter(mapping = aes(log_ratio, response, col = response), 
                width = 0.25, alpha = 0.75, stroke = 1.25, shape = 1, size = 1)+
    coord_flip()+
    xlab("Response")+
    ylab("Log ratio")+
    theme_bw()+
    theme(legend.position = "none")+
    scale_color_brewer(palette = "Set1")

lmer.shannon <- lmer(log_ratio ~ response + (1|dataset), data=df.log_ratio_sbs)
summary(glht(lmer.shannon, linfct = mcp(response = "Tukey")), test = adjusted("bonferroni"))
shapiro.test(residuals(lmer.shannon))

# Function to perform LOGO-CV
df.log_ratio_sbs <- df.log_ratio_sbs |> mutate(response = factor(response, levels = c("NR", "R")))

logo_cv <- function(data, group_var, formula) {
    groups <- unique(data[[group_var]])
    predictions <- data.frame(
      obs = data$response,
      pred = NA,
      NR = NA,
      R = NA,
      group = data[[group_var]]
    )
    
    auc_values <- numeric(length(groups))
    
    for (i in seq_along(groups)) {
      # Split data
      test_index <- which(data[[group_var]] == groups[i])
      train_data <- data[-test_index, ]
      test_data <- data[test_index, ]
      
      # Train model
      model <- glm(formula, data = train_data, family = binomial)
      
      # Predict
      pred_probs <- predict(model, newdata = test_data, type = "response")
      predictions[test_index, "NR"] <- 1 - pred_probs
      predictions[test_index, "R"] <- pred_probs
      predictions[test_index, "pred"] <- ifelse(pred_probs > 0.5, "R", "NR")
      
      # Calculate AUC if both classes are present
      if (length(unique(test_data$response)) == 2) {
        roc_obj <- roc(test_data$response, pred_probs)
        auc_values[i] <- auc(roc_obj)
      } else {
        auc_values[i] <- NA
      }
    }
    
    return(list(predictions = predictions, auc_values = auc_values))
}

result_dataset <- logo_cv(df.log_ratio_sbs, "dataset", response ~ log_ratio)

conf_matrix_dataset <- confusionMatrix(
    factor(result_dataset$predictions$pred, levels = c("NR", "R")), 
    result_dataset$predictions$obs
)

roc_dataset <- roc(result_dataset$predictions$obs, result_dataset$predictions$R)

roc_plot_dataset <- ggroc(roc_dataset) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
    labs(title = "ROC Curve - LOGO-CV by Dataset",
    subtitle = paste0("AUC = ", round(auc(roc_dataset), 3), 
    " (95% CI: ", round(ci.auc(roc_dataset)[1], 3), 
    "-", round(ci.auc(roc_dataset)[3], 3), ")")) +
    theme_bw()+
    theme(panel.grid.major.y = element_blank())

# GSEA
## Taxonomic
vec <- songbird_coef$mean
names(vec) <- songbird_coef$featureid
vec <- sort(vec, decreasing = T)

gs_taxonomy <- cbind(featureid = rownames(tax_table), tax_table)
gs_taxonomy$phylum <- sapply(str_split(gs_taxonomy$phylum, "_"), function(x) x[1])
gs_taxonomy$genus <- sapply(str_split(gs_taxonomy$genus, "_"), function(x) x[1])

mag.biomarkers$species <- sub("_A|_D|_B|_E|_Q|_J|_H|_C|_F|_G|_K|_M", "", mag.biomarkers$species)

gsea_phylum <- GSEA(vec, TERM2GENE = gs_taxonomy[c(3,1)])
gsea_phylum_df  <- as.data.frame(gsea_phylum)

mag.biomarkers[mag.biomarkers$featureid %in% unlist(str_split(gsea_phylum_df$core_enrichment[[1]], "\\/")),]
mag.biomarkers[mag.biomarkers$featureid %in% unlist(str_split(gsea_phylum_df$core_enrichment[[2]], "\\/")),]

gsea_species <- GSEA(vec, TERM2GENE = mag.biomarkers[c(8,2)], eps = 0)
gsea_species_df  <- as.data.frame(gsea_species)

gsea_genus <- GSEA(vec, TERM2GENE = gs_taxonomy[c(7,1)], eps = 0)
gsea_genus_df  <- as.data.frame(gsea_genus)

## Body site-borne and food-borne
