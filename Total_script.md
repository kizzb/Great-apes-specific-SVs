### Mapping and calling -F ( by Smartie-SV )
```
rule blasr:
	input:
		ref="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/NLE.fa",
		que="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/{que}.fa",
		bla="/gpfs/home/kizzhoub/zhoubin/tools/smartie-sv/smartie-sv/bin/blasr",		
	output:
		"mappings/NLE-{que}-aligned.sam",
		"unmappings/NLE-{que}-unaligned.fasta"
	shell:
		"""
		{input.bla} -clipping hard -alignContigs -sam -minMapQV 30 -nproc 6 -minPctIdentity 50 -unaligned {output[1]} {input.que} {input.ref} -out {output[0]}
		"""
		
rule callSV:
	input:
		sam="mappings/NLE-{que}-aligned.sam",
		pri="/gpfs/home/kizzhoub/zhoubin/tools/smartie-sv/smartie-sv/bin/printgaps",
		ref="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/NLE.fa"
	output:
		"variants/NLE-{que}.svs.bed"
	shell:
		"""
		cat {input.sam} | {input.pri} {input.ref} variants/NLE-{wildcards.que}
		"""

```

### Mapping and calling -R ( by Smartie-SV )
```
rule blasr:
	input:
		que="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/NLE.fa",
		ref="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/{ref}.fa",
		bla="/gpfs/home/kizzhoub/zhoubin/tools/smartie-sv/smartie-sv/bin/blasr",		
	output:
		"mappings/{ref}-NLE-aligned.sam",
		"unmappings/{ref}-NLE-unaligned.fasta"
	shell:
		"""
		{input.bla} -clipping hard -alignContigs -sam -minMapQV 30 -nproc 6 -minPctIdentity 50 -unaligned {output[1]} {input.que} {input.ref} -out {output[0]}
		"""
		
rule callSV:
	input:
		sam="mappings/{ref}-NLE-aligned.sam",
		pri="/gpfs/home/kizzhoub/zhoubin/tools/smartie-sv/smartie-sv/bin/printgaps",
		ref="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/{ref}.fa"
	output:
		"variants/{ref}-NLE.svs.bed"
	shell:
		"""
		cat {input.sam} | {input.pri} {input.ref} variants/{wildcards.ref}-NLE
		"""
```

### Mapping and calling - F ( by Minimap2 )
```
rule minimap2:
	input:
		ref="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/NLE.fa", 
		que="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/{sample}.fa"
	output:
		"F-mappings/NLE-{sample}.aligned.sam"
	shell:
		"""
		minimap2 -ax asm20 --cs {input.ref} {input.que} > {output}
		"""
		
rule sam2paf:
	input:
		"F-mappings/NLE-{sample}.aligned.sam"
	output:
		"F-paf_file/NLE-{sample}.aligned.paf"
	shell:
		"""
		paftools.js sam2paf {input} > {output}
		"""
		
rule callSV:
	input:
		"F-paf_file/NLE-{sample}.aligned.paf"
	output:
		"F-variants/NLE-{sample}.svs.txt"
	shell:
		"""
		sort -k6,6 -k8,8n {input} | paftools.js call - > {output}
		"""
```

### Mapping and calling - R ( by Minimap2 )
```
rule minimap2:
	input:
		que="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/NLE.fa", 
		ref="/gpfs/home/kizzhoub/zhoubin/data/diff_spe_assembly/{sample}.fa"
	output:
		"R-mappings/{sample}-NLE.aligned.sam"
	shell:
		"""
		minimap2 -ax asm20 --cs {input.ref} {input.que} > {output}
		"""
		
rule sam2paf:
	input:
		"R-mappings/{sample}-NLE.aligned.sam"
	output:
		"R-paf_file/{sample}-NLE.aligned.paf"
	shell:
		"""
		paftools.js sam2paf {input} > {output}
		"""
		
rule callSV:
	input:
		"R-paf_file/{sample}-NLE.aligned.paf"
	output:
		"R-variants/{sample}-NLE.svs.txt"
	shell:
		"""
		sort -k6,6 -k8,8n {input} | paftools.js call - > {output}
		"""
```
## Next analysis by Minimap2
### Extract DEL and INS
```
for i in V38 CCP GGO PAB RM10 CTJ
do 
awk 'BEGIN{OFS="\t"}{if($1=="V" && $4-$3==0 && $11-$10>49)print $2,$3,$4+$11-$10,$7,$8,$9,$10,$11,$12}' NLE-"$i"-SV.txt > NLE-"$i"-SV.txt-INS-plus
awk 'BEGIN{OFS="\t"}{if($1=="V" && $4-$3>49 && $11-$10==0)print $2,$3,$4,$7,$8,$9,$10,$11,$12}' NLE-"$i"-SV.txt > NLE-"$i"-SV.txt-DEL
done
```

### Intersect SV set of F and R
```
# for V38 CCP GGO PAB 
for i in V38 CCP GGO PAB 
do
bedtools intersect -a ../F-variants/NLE-"$i"-SV.txt-DEL -b ../R-variants/NLE-"$i"-SV.txt-DEL -wa -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > NLE-"$i"-SV.txt-DEL
bedtools intersect -a ../F-variants/NLE-"$i"-SV.txt-INS-plus -b ../R-variants/NLE-"$i"-SV.txt-INS-plus -wa -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > NLE-"$i"-SV.txt-INS-plus
done

# for RM10
cat ../F-variants/NLE-RM10-SV.txt-DEL ../R-variants/NLE-RM10-SV.txt-DEL |sort -k1,1 -s -V -k2n,2  > NLE-RM10-SV.txt-DEL
cat ../F-variants/NLE-RM10-SV.txt-INS-plus ../R-variants/NLE-RM10-SV.txt-INS-plus |sort -k1,1 -s -V -k2n,2 > NLE-RM10-SV.txt-INS-plus
```

### Cross species intersect
```
# for DEL
rule All:
	input: "N-V-C-G-P-v-R.bed-DEL"

rule itst_H_C:
	input: "../F-R-intsec/NLE-V38-SV.txt-DEL", "../F-R-intsec/NLE-CCP-SV.txt-DEL"
	output: "N-V-C.bed-DEL"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule itst_H_C_G:
	input: "N-V-C.bed-DEL", "../F-R-intsec/NLE-GGO-SV.txt-DEL"
	output: "N-V-C-G.bed-DEL"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule itst_H_C_G_P:
	input: "N-V-C-G.bed-DEL", "../F-R-intsec/NLE-PAB-SV.txt-DEL"
	output: "N-V-C-G-P.bed-DEL"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule filter_by_RM10:
	input: "N-V-C-G-P.bed-DEL", "../F-R-intsec/NLE-RM10-SV.txt-DEL"
	output: "N-V-C-G-P-v-R.bed-DEL"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -r -v |sort -k1,1 -s -V -k2n,2 -u > {output}
		"""

# for INS
rule All:
	input: "N-V-C-G-P-v-R.bed-INS-plus"

rule itst_H_C:
	input: "../F-R-intsec/NLE-V38-SV.txt-INS-plus", "../F-R-intsec/NLE-CCP-SV.txt-INS-plus"
	output: "N-V-C.bed-INS-plus"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule itst_H_C_G:
	input: "N-V-C.bed-INS-plus", "../F-R-intsec/NLE-GGO-SV.txt-INS-plus"
	output: "N-V-C-G.bed-INS-plus"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule itst_H_C_G_P:
	input: "N-V-C-G.bed-INS-plus", "../F-R-intsec/NLE-PAB-SV.txt-INS-plus"
	output: "N-V-C-G-P.bed-INS-plus"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule filter_by_RM10:
	input: "N-V-C-G-P.bed-INS-plus", "../F-R-intsec/NLE-RM10-SV.txt-INS-plus"
	output: "N-V-C-G-P-v-R.bed-INS-plus"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -v |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 
```

### Return the original INS coordinates
```
awk '{for(i=1;i<37;i++)
		{
			if(i==3)printf $2"\t";
			else if(i==12)printf $11"\t";
			else if(i==21)printf $20"\t";
			else if(i==30)printf $29"\t"; 
			else if(i==36)printf $i"\n";
			else printf $i"\t";
		}
	}' 01_N-V-C-G-P-v-R.bed-INS-plus > 01_N-V-C-G-P-v-R.bed-INS
```

### Add ID for SV file
```
# DEL
awk '{if(NR<10)print "NGSSV0000"NR"\t"$0; else if(NR<100)print "NGSSV000"NR"\t"$0; else if(NR<1000)print "NGSSV00"NR"\t"$0; else if(NR<10000)print "NGSSV0"NR"\t"$0; else print "NGSSV"NR"\t"$0}' 01_N-V-C-G-P-v-R.bed-DEL | awk '{printf $2"\t"$3"\t"$4"\t"$1"\t"; for(i=5;i<38;i++)printf $i"\t";printf "\n"}' > 02_N-V-C-G-P-v-R.bed-DEL_addNum
# INS
awk '{print "NGSSV"NR+14451"\t"$0}' 01_N-V-C-G-P-v-R.bed-INS-plus | awk '{printf $2"\t"$3"\t"$4"\t"$1"\t"; for(i=5;i<38;i++)printf $i"\t";printf "\n"}'> 02_N-V-C-G-P-v-R.bed-INS-plus_addNum

awk '{print "NGSSV"NR+14451"\t"$0}' 01_N-V-C-G-P-v-R.bed-INS | awk '{printf $2"\t"$3"\t"$4"\t"$1"\t"; for(i=5;i<38;i++)printf $i"\t";printf "\n"}'> 02_N-V-C-G-P-v-R.bed-INS_addNum

# Then move ID column to 3rd column
```

### Manual check for 23300 condidate GSSVs
```
# Get the coord for 4 great-ape species (V38 CCP GGO PAB) and extend them by 1k above and below
awk '{filename="./bed_file/"$4"_NLE_1k.bed"; print $1"\t"$2-1000"\t"$3+1000 >> filename; close(filename)}' 02_N-V-C-G-P-v-R.bed-All_addNum
awk '{filename="./bed_file/"$4"_V38_1k.bed"; print $7"\t"$8-1000"\t"$9+1000 >> filename; close(filename)}' 02_N-V-C-G-P-v-R.bed-All_addNum
awk '{filename="./bed_file/"$4"_CCP_1k.bed"; print $16"\t"$17-1000"\t"$18+1000 >> filename; close(filename)}' 02_N-V-C-G-P-v-R.bed-All_addNum
awk '{filename="./bed_file/"$4"_GGO_1k.bed"; print $25"\t"$26-1000"\t"$27+1000 >> filename; close(filename)}' 02_N-V-C-G-P-v-R.bed-All_addNum
awk '{filename="./bed_file/"$4"_PAB_1k.bed"; print $34"\t"$35-1000"\t"$36+1000 >> filename; close(filename)}' 02_N-V-C-G-P-v-R.bed-All_addNum

# Get fa by coord		
for i in `awk '{print $4}' 02_N-V-C-G-P-v-R.bed-All_addNum`
do 
for j in V38 CCP GGO PAB NLE
do 
bedtools getfasta -bed bed_file/"$i"_"$j"_1k.bed -fi ~/zhoubin/data/diff_spe_assembly/"$j".fa -fo fa_file/"$i"_"$j"_1k.fa
done  
done

# Find similar seq of RM10 by blastn
for i in `awk '{print $4}' 02_N-V-C-G-P-v-R.bed-All_addNum`
do 
blastn -max_target_seqs 5 -evalue 1e-100 -perc_identity 80 -query fa_file/"$i"_NLE_1k.fa -db ~/zhoubin/data/blast_db/Common_used_blstdb/RM10db -outfmt 6 -out asn_file/"$i"_NLE_RM10_1k.asn
done

# Calculate the math coord of RM10
sed -i -e 's/:/\t/g' -e 's/-/\t/g' asn_file/*asn 

for i in `awk '{print $4}' 02_N-V-C-G-P-v-R.bed-All_addNum`
do
head -1 asn_file/"$i"_NLE_RM10_1k.asn |awk '{if($12>$11)print $4"\t"$11-($9-1)"\t"$12+($3-$2-$10);else print $4"\t"$12-($3-$2-$10)"\t"$11+($9-1)}' > bed_file/"$i"_RM10_1k.bed
done

# Get fa of RM10
for i in `awk '{print $4}' 02_N-V-C-G-P-v-R.bed-All_addNum`
do
bedtools getfasta -bed bed_file/"$i"_RM10_1k.bed -fi ~/zhoubin/data/diff_spe_assembly/RM10.fa -fo fa_file/"$i"_RM10_1k.fa
done

# Alignment by Nucmer
for i in `awk '{print $4}' 02_N-V-C-G-P-v-R.bed-All_addNum`
do
for j in V38 CCP GGO PAB RM10 
do
nucmer -b 20 --prefix=delta_file/"$i"_NLE_"$j" fa_file/"$i"_NLE_1k.fa fa_file/"$i"_"$j"_1k.fa
done
done

# Plot by Mummerplot
for i in `awk '{print $4}' 02_N-V-C-G-P-v-R.bed-All_addNum`
do
for j in V38 CCP GGO PAB RM10 
do
mummerplot --png --small delta_file/"$i"_NLE_"$j".delta -p mummerplot_output/"$i"_NLE_"$j"
done
done

# Compress and package the file
mv mummerplot_output/*.png png_file

cd ./png_file
for i in `awk '{print $4}' ../02_N-V-C-G-P-v-R.bed-All_addNum`
do
mkdir "$i"
mv "$i"_*.png "$i"
done

tar --remove-files -czf 23300_Condidate_GSSV_Manual_Check.tar.gz  NGSSV*

# The last step: filter GSSVs by manual check
```

### Repeat analysis
```
# RepeatMasker annotation for V38 and NLE
RepeatMasker -parallel 10 -species primates -html -gff -dir Gibbon_NLE ~/zhoubin/data/diff_spe_assembly/NLE.fa
RepeatMasker -parallel 10 -species primates -html -gff -dir Human_V38 ~/zhoubin/data/diff_spe_assembly/V38.fa

# Overlap 15885 GSSVs and repeat region
bedtools intersect -a NLE_Coord.bed_DEL -b NLE.repeat.bed -wa -wb -f 0.5 |sort -k1,1 -s -V -k2n,2 -u > DEL_SVs_itsct_with_NLE_repeat.bed
bedtools intersect -a V38_Coord.bed_INS -b V38.repeat.bed -wa -wb -f 0.5 |sort -k1,1 -s -V -k2n,2 -u > INS_SVs_itsct_with_V38_repeat.bed
```

##### further classified the 15885 GSSVs into two sets which contain 2366 DELs and 4208 INSs

### VEP predict for 2366 DELs and 4208 INSs
Using the VEP online tool (http://www.ensembl.org/Tools/VEP) to predict the effect of SV, "Species" select Homo (GRCh38.p13) , "Restrict results" select "show most severe consequence per variant", and others were the default parameters.

### GO analysis
```
#####  get Genes which overlap with HC-GSSVs
# GRCh38_105_Coding_gene.bed download from ensemble (https://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/)
bedtools intersect -a DEL_hg38.bed -b GRCh38_105_Coding_gene.bed -wa -wb |awk '{print $7}' >> ENSG_Ovlp_with_6574_HCGSSV.bed
bedtools intersect -a INS_hg38.bed -b GRCh38_105_Coding_gene.bed -wa -wb |awk '{print $7}' >> ENSG_Ovlp_with_6574_HCGSSV.bed

#####  GO analysis by clusterProfiler
# Go analysis by clusterProfiler
library(clusterProfiler)

# data process
df <- read.table("ENSG_Ovlp_with_6574_HCGSSV.bed", header = F)
colnames(df) <- "ENSEMBL"
gene.df = bitr(df$ENSEMBL, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

# plot
barplot(ego, showCategory=20)
dotplot(ego, showCategory=20)
```

### Epigenomic and transcript data analysis
```
####################
### CRE analysis ###
####################

##### Calculate Pvalue between Human-Chimpanzee, Human-Macaque and Chimpanzee-Macaque
library(tidyverse)

# Calculate Pvalue of Human-Chimpanzee, Human-Macaque and Chimpanzee-Macaque
# Read data from literature
df <- read.table("CREs_3Species_ReadCounts_Normalized_CRELength_ReadsInPeaks_NN-RS53875B.txt", header = T)

# Interation based on 7 brain region 
for (j in c("CaudateNucleus","Cerebellum","OccipitalPole","PrecentralGyrus","PrefrontalCortex","Putamen","ThalamicNuclei")){
	
	# Generate new tab file to save reslut
	sink(paste0(j,"_hu_ch_rh_DE.tab"))
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
		cat(toString(df[i,1]), rowMeans(hu_value), rowMeans(ch_value), rowMeans(rh_value), DE_hu_ch$p.value, DE_hu_rh$p.value, DE_ch_rh$p.value, sep="\t", end="\n")
	}

	sink()
}

##### Filter CREs which are Human-Chimpanzee specific
for i in CaudateNucleus Cerebellum OccipitalPole PrecentralGyrus PrefrontalCortex Putamen ThalamicNuclei 
do 
awk '{if($5>0.05 && $6<0.05 && $7<0.05)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t'$i'"}' "$i"_hu_ch_rh_DE.tab >> Integrate_all_region_hu_ch_specific_CRE.bed 
done

##### Calculate Log2Fold
awk '{print $1"\t"log((($2+$3)/2)/$4)"\t"$8}' Integrate_all_region_hu_ch_specific_CRE.bed > Integrate_all_region_log_fold_CRE.bed
awk 'BEGIN{print "ID\tLog2_fold[mean(Hu+Ch)/Ma]\tRegion"}{print $1"\t"$2/log(2)"\t"$3}' Integrate_all_region_log_fold_CRE.bed > Integrate_all_region_Log2Fold_CRE.bed

##### Match CRE id with it's coord
awk 'NR==FNR{a[$4]=$0;next}{if($1 in a)print a[$1]"\t"$2"\t"$3}' CREs_3Species_hg38.sorted.bed Integrate_all_region_Log2Fold_CRE.bed |awk 'BEGIN{print "CRE_Chr\tCRE_Start\tCRE_End\tCRE_ID\tLog2Fold\tCRE_Region"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > Match_CRE_coord_Log2Fold.bed

##### Overlap SV with CRE
bedtools intersect -a DEL_coord_plus_hg38.bed -b Match_CRE_coord_Log2Fold.bed -wa -wb |sort -k1,1 -s -V -k2n,2 |uniq > DEL_ovlp_CRE_Log2Fold.bed
bedtools intersect -a INS_coord_hg38.bed -b Match_CRE_coord_Log2Fold.bed -wa -wb |sort -k1,1 -s -V -k2n,2 |uniq > INS_ovlp_CRE_Log2Fold.bed

##### Change the format
# Rscript
library(tidyverse)

for (i in c("DEL","INS")){
  df <- read.table(paste0(i,"_ovlp_CRE_Log2Fold.bed"), header=F)
  colnames(df) <- c("chr_SV","start_SV","end_SV","ID_SV","chr_CRE","start_CRE","end_CRE","ID_CRE","Log2Fold","Region")
  df <- spread(df, Region, Log2Fold, fill="NS")
  write.csv(df, paste0(i,"_ovlp_CRE_Log2Fold_ChangeFormat.csv"))
}

# csvtk trans
for i in DEL INS
do
csvtk csv2tab "$i"_ovlp_CRE_Log2Fold_ChangeFormat.csv |cut -f2- > "$i"_ovlp_CRE_Log2Fold_ChangeFormat.bed
done

##### Calculate CRE flanking 500k gene
# Generate CRE flank 500k coordinate and it's length
awk '{if($2-500000<=0)print $1"\t1\t"$3+500000"\t"$4"\t"$3-$2; else print $1"\t"$2-500000"\t"$3+500000"\t"$4"\t"$3-$2}' Match_CRE_coord_Log2Fold.bed|sort -k1,1 -s -V -k2n,2| uniq > Match_CRE_and_coord_fl500k.bed

# CRE flanking 500k region overlap with Gene
bedtools intersect -a Match_CRE_and_coord_fl500k.bed -b ../GRCh38_105_Coding_gene.bed -wa -wb |sort |uniq |awk 'BEGIN{OFS="\t"}{print $1,$3-500000-$5,$3-500000,$4,$6,$7,$8,$9,$10}'|awk '{print $0"\t"$6-$2}' |sort -k1,1 -s -V -k2n,2 |uniq > Result_CRE_fl500k_Gene.bed

##### Match SV, CRE and CRE flanking 500k gene
# Merge coord column for DEL and INS
sed '1d' DEL_ovlp_CRE_Log2Fold_ChangeFormat.bed |awk 'BEGIN{OFS="\t";print "SV_coord\tSV_ID\tCRE_coord\tCRE_ID\tCaudateNucleus\tCerebellum\tOccipitalPole\tPrecentralGyrus\tPrefrontalCortex\tPutamen\tThalamicNuclei"}{print $1"_"$2"_"$2,$4,$5"_"$6"_"$7,$8,$9,$10,$11,$12,$13,$14,$15}' > DEL_ovlp_CRE_Log2Fold_ChangeFormat_MergeCol.bed
sed '1d' INS_ovlp_CRE_Log2Fold_ChangeFormat.bed |awk 'BEGIN{OFS="\t";print "SV_coord\tSV_ID\tCRE_coord\tCRE_ID\tCaudateNucleus\tCerebellum\tOccipitalPole\tPrecentralGyrus\tPrefrontalCortex\tPutamen\tThalamicNuclei"}{print $1"_"$2"_"$3,$4,$5"_"$6"_"$7,$8,$9,$10,$11,$12,$13,$14,$15}' > INS_ovlp_CRE_Log2Fold_ChangeFormat_MergeCol.bed

# Merge coord column for CRE and CRE flanking 500k genes
awk 'BEGIN{OFS="\t";print "CRE_coord\tCRE_ID\tGene_coord\tENSG\tGene\tDistance"}{print $1"_"$2"_"$3,$4,$5"_"$6"_"$7,$8,$9,$10}' Result_CRE_fl500k_Gene.bed > Result_CRE_fl500k_Gene_MergeCol.bed

# Match
awk 'NR==FNR{a[$4]=$0;next}{if($2 in a)print a[$2]"\t"$3"\t"$4"\t"$5"\t"$6}' DEL_ovlp_CRE_Log2Fold_ChangeFormat_MergeCol.bed Result_CRE_fl500k_Gene_MergeCol.bed > Link_SV_CRE_Gene.bed_DEL
awk 'NR==FNR{a[$4]=$0;next}{if($2 in a)print a[$2]"\t"$3"\t"$4"\t"$5"\t"$6}' INS_ovlp_CRE_Log2Fold_ChangeFormat_MergeCol.bed Result_CRE_fl500k_Gene_MergeCol.bed > Link_SV_CRE_Gene.bed_INS

################################
### Gene experssion analysis ###
################################

# Extract Gene ID from simFiltered_rpkm_combat_nonParametric.txt
awk '{print $1}' simFiltered_rpkm_combat_nonParametric.txt |sed -e 's/"//g' -e 's/|/\t/g' |sed '1d' |awk '{print $1}'|sort |uniq > Gene_ID_of_simFiltered_RNA_data.bed

# Extract Gene ID from Link_SV_CRE_Gene.bed
awk '{print $13}' Link_SV_CRE_Gene.bed_DEL |sed '1d' > test.bed
awk '{print $13}' Link_SV_CRE_Gene.bed_INS |sed '1d' >> test.bed
sort test.bed |uniq > Linked_Gene_ID.bed 

# Match Gene ID from Gene_ID_of_simFiltered_RNA_data.bed and Linked_Gene_ID.bed 
awk 'NR==FNR{a[$1]=$1;next}{if($1 in a)print}' Gene_ID_of_simFiltered_RNA_data.bed Linked_Gene_ID.bed > Overlaped_Gene_ID.bed

# R package
library(tidyverse)

# Read raw data
df <- read.table("simFiltered_rpkm_combat_nonParametric.txt", header =T )

# Read gene list
cdt_gene <- read.table("Overlaped_Gene_ID.bed", header = F)
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

##### Filter genes which are Human-Chimpanzee specific
for i in MFC DFC M1C V1C OFC VFC MD STR CBC
do 
sed '1d' "$i"_DG_hu_ch_ma.csv |csvtk csv2tab|awk '{if($5>0.05 && $6<0.05 && $7<0.05)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t'$i'"}' >> Integrate_all_region_hu_ch_specific_gene.bed 
done

##### Calculate Log2Fold
awk '{print $1"\t"log((($2+$3)/2)/$4)"\t"$8}' Integrate_all_region_hu_ch_specific_gene.bed > Integrate_all_region_log_fold_gene.bed
awk 'BEGIN{print "ID\tLog2_fold[mean(Hu+Ch)/Ma]\tRegion"}{print $1"\t"$2/log(2)"\t"$3}' Integrate_all_region_log_fold_gene.bed > Integrate_all_region_Log2Fold_gene.bed

##### Change the format
# Rscript
library(tidyverse)

df <- read.table("Integrate_all_region_Log2Fold_gene.bed", header=T)
colnames(df) <- c("ENSG","Gene","Log2Fold","Region")
df <- spread(df, Region, Log2Fold, fill="NS")
write.csv(df, "Integrate_all_region_Log2Fold_gene_ChangeFormat.csv")

# csvtk trans
csvtk csv2tab Integrate_all_region_Log2Fold_gene_ChangeFormat.csv |cut -f2- > Integrate_all_region_Log2Fold_gene_ChangeFormat.bed

# Match SV-CRE-Gene
awk 'NR==FNR{a[$13]=$0;next}{if($1 in a)print a[$1]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' Link_SV_CRE_Gene.bed_DEL 07_DG_analysis/Integrate_all_region_Log2Fold_gene_ChangeFormat.bed > Link_SV_CRE_Gene_Log2Fold.bed_DEL
awk 'NR==FNR{a[$13]=$0;next}{if($1 in a)print a[$1]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' Link_SV_CRE_Gene.bed_INS 07_DG_analysis/Integrate_all_region_Log2Fold_gene_ChangeFormat.bed > Link_SV_CRE_Gene_Log2Fold.bed_INS
```

### TF analysis by FIMO for NOL3 (SV ID: NGSSV14990)
```
# Calculate TF of SV sequence of V38, CCP and RM10
for i in V38 CCP RM10
do
fimo --o NGSSV14990_TF_"$i"_result ~/zhoubin/data/JASPAR_TF_data/Vertebrates/all.meme NGSSV14990_CRE_"$i".fa
done

# Overlap between V38 and CCP
awk 'NR==FNR{a[$2]=$2;next}{if($2 in a)print a[$2]"\t"$2}' NGSSV14990_TF_V38_result/V38_fimo.tsv NGSSV14990_TF_CCP_result/CCP_fimo.tsv |sed '1d'|awk '{print $1}'|sort |uniq > TF_overlap_between_V38_CCP.bed
# Filter by RM10
awk 'NR==FNR{a[$2]=$2;next}{if($1 in a);else print $1}' RM10_fimo.tsv TF_overlap_between_V38_CCP.bed|sort |uniq > V38_CCP_specific_TF.bed
```