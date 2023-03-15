library(tidyverse)

# Read raw data
df <- read.table("01_simFiltered_rpkm_combat_nonParametric.txt", header =T )

# Read gene list
cdt_gene <- read.table("01_Overlaped_Gene_ID.bed", header = F)
cdt_gene <- as.vector(cdt_gene$V1)

# Brain region
cdt_region <- c("MFC","DFC","M1C","V1C","OFC","VFC","MD","STR","CBC")

# Loop
for(region in cdt_region){

        sink(paste0(region,"_DG_hu_ch_ma.csv"))
        cat("Gene_symble,Mean_Human_RPKM,Mean_Chimp_RPKM,Mean_Macaque_RPKM,DG_hu_ch,DG_hu_ma,DG_ch_ma\n")

        for(gene in cdt_gene){

                Human <- select(df, contains("HSB"))[grep(gene, df$ProbeID),] %>% select(contains(region))
                Chimp <- select(df, contains("PTB"))[grep(gene, df$ProbeID),] %>% select(contains(region))
                Macaque <- select(df, contains("RMB"))[grep(gene, df$ProbeID),] %>% select(contains(region))

				Human <- Human[ Human > 0 ]
				Chimp <- Chimp[ Chimp > 0 ]
				Macaque <- Macaque[ Macaque > 0 ]

				if (sum(!is.na(unique(Human))) > 1 && sum(!is.na(unique(Chimp))) > 1 && sum(!is.na(unique(Macaque))) > 1) {
					DG_hu_ch <- t.test(Human, Chimp)
					DG_hu_ma <- t.test(Human, Macaque)
					DG_ch_ma <- t.test(Chimp, Macaque)
				} else {
					DG_hu_ch <- data.frame(p.value=c("NA"))
					DG_hu_ma <- data.frame(p.value=c("NA"))
					DG_ch_ma <- data.frame(p.value=c("NA"))
				}

                cat(toString(df[grep(gene, df$ProbeID),1]), mean(Human), mean(Chimp), mean(Macaque), DG_hu_ch$p.value, DG_hu_ma$p.value, DG_ch_ma$p.value, sep=",", end="\n")

        }

        sink()
}