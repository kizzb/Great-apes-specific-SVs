rule blasr:
	input:
		ref="NLE.fa",
		que="{que}.fa",
		bla="tools/smartie-sv/smartie-sv/bin/blasr",		
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
		pri="tools/smartie-sv/smartie-sv/bin/printgaps",
		ref="/NLE.fa"
	output:
		"variants/NLE-{que}.svs.bed"
	shell:
		"""
		cat {input.sam} | {input.pri} {input.ref} variants/NLE-{wildcards.que}
		"""
