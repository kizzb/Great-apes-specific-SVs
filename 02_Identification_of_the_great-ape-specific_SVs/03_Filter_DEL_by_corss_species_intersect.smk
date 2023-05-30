rule All:
	input: "N-V-C-G-P-v-R.bed-DEL"

rule step_1:
	input: "../F-R-intsec/NLE-V38-SV.txt-DEL", "../F-R-intsec/NLE-CCP-SV.txt-DEL"
	output: "N-V-C.bed-DEL"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule step_2:
	input: "N-V-C.bed-DEL", "../F-R-intsec/NLE-GGO-SV.txt-DEL"
	output: "N-V-C-G.bed-DEL"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule step_3:
	input: "N-V-C-G.bed-DEL", "../F-R-intsec/NLE-PAB-SV.txt-DEL"
	output: "N-V-C-G-P.bed-DEL"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -f 0.5 -r |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 

rule step_4:
	input: "N-V-C-G-P.bed-DEL", "../F-R-intsec/NLE-RM10-SV.txt-DEL"
	output: "N-V-C-G-P-v-R.bed-DEL"
	shell:
		"""
		bedtools intersect -a {input[0]} -b {input[1]} -wa -wb -r -v |sort -k1,1 -s -V -k2n,2 -u > {output}
		""" 