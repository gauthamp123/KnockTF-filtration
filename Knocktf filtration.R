library(dplyr)
BiocManager::install("DT")
library(DT)
library(dorothea)
load("~/R/knocktf/KnockTF/R/generanking.Rdata")
data(dorothea_hs, package = "dorothea")
viper_regulon <- dorothea_hs %>%
  filter(confidence %in% c("A", "B","C"))



data = data.frame(KnockTF)

#filter by p-value
data = data %>% 
  filter(data$P_value < 0.05)

#zscore log2FC data and append to dataset
data$zscore <- (data$Log2FC - mean(data$Log2FC) / sd(data$Log2FC))

#inverse zscore value and add as new column 
data$inverseZ <- -data$zscore



#convert dataset into gmt file
gmt_data <- split(data$Gene, data$TF)
KnockTF.gmt <- as.data.frame(do.call(rbind, gmt_data))

#create adj file for regulon object (TF mode = zscore from -1 to 1 and likelyhood = 1)
adjfile <- data.frame(data$TF)
adjfile$Gene <- data$Gene
tfmode <- c(data$zscore[1])    #each dataset corresponds to a TF, so do we average the zscores?
list = list(tfmode)

for(i in 2:nrow(data)){
  tempVec <- c(data$zscore[i])
  list[[i]] <- tempVec
  
}
adjfile$TFmode <- list

frequency <- c(1)
for(j in 1:nrow(data)){
  frequency <- 1
}
adjfile$frequency <- frequency

#filter by inversez score and delete duplicate edges
knockTF_df <- data[,c(2,1,3,12)]
colnames(knockTF_df) <- c("tf", "confidence", "target", "mor")
knockTF_df <- as_tibble(knockTF_df)   #not sure if this should be there
knockTF_df %>% arrange(mor)   #zscore value is 18 should that also be a 1?


#keep the range of frequency between -1 and 1
for(k in 1:nrow(knockTF_df)){
  if(knockTF_df$mor[k] > 1){
    knockTF_df$mor[k] <- 1
  }
  else if(knockTF_df$mor[k] < -1){
    knockTF_df$mor[k] <- -1
  }
}

#extract tf and target columns
#run it through unique while arranging the table by highest mor
#find a way to merge the 2 tables together
temp_table <- knockTF_df[,c(1,2,3)]
temp_table <- unique(temp_table)  #how do I add the mor column?




knockTF_regulon <- df2regulon(knockTF_df)

mrs <- viper::msviper(ges = generanking, regulon = knockTF_regulon, minsize = 0, ges.filter = F)
summary(mrs)
TF_activities <- data.frame(Regulon = names(mrs$es$nes),
                            Size = mrs$es$size[ names(mrs$es$nes) ],
                            NES = mrs$es$nes,
                            p.value = mrs$es$p.value,
                            FDR = p.adjust(mrs$es$p.value, method = 'fdr'))

# Tidy up and order by FDR
TF_activitiesTidy <- TF_activities %>%
  as_tibble() %>%
  arrange()
# Show in a nice table
DT::datatable(TF_activitiesTidy)

#sample code
if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("mixtools")
BiocManager::install("bcellViper")
BiocManager::install("viper")

library(viper)

data(bcellViper, package="bcellViper")
adjfile <- system.file("aracne", "bcellaracne.adj", package = "bcellViper")
regul <- aracne2regulon(adjfile, dset, verbose = FALSE)
print(regul)


