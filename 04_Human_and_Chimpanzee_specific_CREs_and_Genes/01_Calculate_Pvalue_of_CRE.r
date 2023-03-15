library(tidyverse)

# Read data from literature
df <- read.table("human_chimp_macaque_DE_analysis/CREs_3Species_ReadCounts_Normalized_CRELength_ReadsInPeaks_NN-RS53875B.txt", header = T)

# Interation based on 8 brain region 
for (j in c("CaudateNucleus","Cerebellum","OccipitalPole","PrecentralGyrus","PrefrontalCortex","Putamen","ThalamicNuclei")){
	
	# Generate new csv file to save reslut
	sink(paste0(j,"_hu_ch_rh_DE.csv"))
	cat("CRE,Human_mean,Chimp_mean,Macaque_mean,DE_hu_ch,DE_hu_rh,DE_ch_rh\n")

	# Interation per line data(Total 60702 CREs)
	for(i in seq(1, 60702)){
		
		# Select data of per species
		hu_value <- select(df, contains(j))[i,] %>% select(contains("Human"))
		ch_value <- select(df, contains(j))[i,] %>% select(contains("Chimp"))
		rh_value <- select(df, contains(j))[i,] %>% select(contains("Rhesus"))

		# Calculate Pvalue of every pair species
		DE_hu_ch <- t.test(hu_value, ch_value)
		DE_hu_rh <- t.test(hu_value, rh_value)
		DE_ch_rh <- t.test(ch_value, rh_value)

		# Output result
		cat(toString(df[i,1]), rowMeans(hu_value), rowMeans(ch_value), rowMeans(rh_value), DE_hu_ch$p.value, DE_hu_rh$p.value, DE_ch_rh$p.value, sep=",", end="\n")
	}

	sink()
}