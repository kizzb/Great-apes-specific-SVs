rule minimap2:
	input:
		ref="NLE.fa", 
		que="{sample}.fa"
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
