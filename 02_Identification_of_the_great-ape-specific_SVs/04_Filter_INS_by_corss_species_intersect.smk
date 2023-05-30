# The INS coordinate of gibbon is only 1 point, so, we add the SV length to the end of coordinate and then do filter

rule All:
	input: "N-V-C-G-P-v-R.bed-INS-plus"

rule step_1:
	input: "../F-R-intsec/NLE-V38-SV.txt-INS-plus", "../F-R-intsec/NLE-CCP-SV.txt-INS-plus"
	output: "N-V-C.bed-INS-plus"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule step_2:
	input: "N-V-C.bed-INS-plus", "../F-R-intsec/NLE-GGO-SV.txt-INS-plus"
	output: "N-V-C-G.bed-INS-plus"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule step_3:
	input: "N-V-C-G.bed-INS-plus", "../F-R-intsec/NLE-PAB-SV.txt-INS-plus"
	output: "N-V-C-G-P.bed-INS-plus"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule step_4:
	input: "N-V-C-G-P.bed-INS-plus", "../F-R-intsec/NLE-RM10-SV.txt-INS-plus"
	output: "N-V-C-G-P-v-R.bed-INS-plus"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -v |sort -k1,1 -s -V -k2n,2 -u > {output}
		"""