configfile:  "multiPrime-original.yaml"
#__version__ = "2.0.3"
#__date__ = "2022-7-28"
#__author__ = "Junbo Yang"
#__email__ = yang_junbo_hi@126.com; 1806389316@pku.edu.cn
#__description__ = "A Snakemake workflow to design multiPCR primers and get Max_primerset"

import os

virus = config["virus"]

def aggregate_input(wildcards):
	checkpoint_output = checkpoints.extract_cluster_fa.get(**wildcards).output[1]
	return expand(config["results_dir"] + "/Clusters_cprimer/{i}.candidate.primers.txt",
		i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fa")).i)
#Here a new directory will be created for each sample by the checkpoint. 
#After completion of the checkpoint. the aggregate_input function is re-evaluated as previously.
#The values of the wildcard i is this time used to expand the pattern "post/{sample}/{i}.txt",
#such that the rule intermediate is executed for each of the determined clusters.

rule all:
	## LOCAL ##
#	'''
#	Defines the target rule by specifying all final output files.
#	Additionally, the cluster_logs dir is deleted if
#	snakemake was run locally.
#	'''
	input: 
		config["results_dir"] + "/Primers_set/candidate_primers_sets.txt",
		config["results_dir"] + "/Primers_set/final_maxprimers_set.xls",
		config["results_dir"] + "/Primers_set/Coverage_stast.xls",
		config["results_dir"] + "/Primers_set/candidate_primers_sets.number",
		config["results_dir"] + "/Primers_set/final_maxprimers_set.fa.dimer",
		config["results_dir"] + "/Core_primers_set/core_Coverage_stast.xls",
		config["results_dir"] + "/Core_primers_set/core_final_maxprimers_set.fa.dimer",
		config["results_dir"] + "/Core_primers_set/core_candidate_primers_sets.number",
		config["results_dir"] + "/Core_primers_set/BWT_coverage/core_final_maxprimers_set.out"

#-------------------------------------------------------------------------------------------
# seq_format rule 1: Dependency packages - python
#-------------------------------------------------------------------------------------------
rule seq_format:
	input: 
		config["input_dir"] + "/{virus}.fa"
	output:
		config["results_dir"] + "/Total_fa/{virus}.format.fa"
	params:
		config["scripts_dir"]
	message:
		"Step1: format fasta file.."
	shell:
		'''
		python {params}/seq_format.py -i {input} -o {output}
		'''
#-------------------------------------------------------------------------------------------
# extract_ rule 2: Dependency packages - None
#-------------------------------------------------------------------------------------------
#rule extract_seq:
#	input:
#		config["results_dir"] + "/Total_fa/{virus}.format.fa"
#	output:
#		config["results_dir"] + "/Total_fa/{virus}.format.fa",
#	message:
#		"Step2: extract  from the raw fasta .."
#	shell:
#		'''
#		cat {input} | grep -A 1 "" > {output[0]}
#		'''
#-------------------------------------------------------------------------------------------
# rmdup rule 3: Dependency packages - cd-hit
#-------------------------------------------------------------------------------------------
rule rmdup:
	input:
		config["results_dir"] + "/Total_fa/{virus}.format.fa"
	output:
		config["results_dir"] + "/Total_fa/{virus}.format.rmdup.cluster.fa"
	message: 
		"Step3: Remove duplicated sequence .."
	shell:
		'''
		cd-hit -M 0 -T 0 -i {input} -o {output} -c 1
		'''
#-------------------------------------------------------------------------------------------
# bowtie2-index rule 2: Dependency packages - bowtie2
#-------------------------------------------------------------------------------------------
#rule makedb:
#	input:
#		config["results_dir"] + "/Total_fa/{virus}.format.rmdup.cluster.fa"
#	output:
#		touch(config["results_dir"] + "/Bowtie_db/done")
#	params:
#		output_prefix=config["results_dir"] + "/Bowtie_db/genome"
#	message:
#		"Step2: Build index for BWT (bowtie2) .."
#	shell:
#		'''
#		bowtie2-build {input} {params.output_prefix}
#		'''
#-------------------------------------------------------------------------------------------
# cluster_by_identity rule 4: Dependency packages - cd-hit || suggest: identity=0.8
#-------------------------------------------------------------------------------------------
rule cluster_by_identity:
	input:  
		config["results_dir"] + "/Total_fa/{virus}.format.rmdup.cluster.fa"
	output: 
		config["results_dir"] + "/Total_fa/{virus}.format.rmdup.cluster.uniq.fa",
		config["results_dir"] + "/Total_fa/{virus}.format.rmdup.cluster.uniq.fa.clstr"
	# It is allowed that output name not used in the shell.
	message: 
		"Step4: Cluster sequences by identity .."
	params:
		identity=config['identity']
	shell:
		'''
		cd-hit -M 0 -T 0 -i {input} -o {output[0]} -c {params.identity}
		'''
#-------------------------------------------------------------------------------------------
# extract_cluster_fa rule 5: Dependency packages - None
#-------------------------------------------------------------------------------------------
checkpoint extract_cluster_fa:
	input:
		expand(config["results_dir"] + "/Total_fa/{virus}.format.rmdup.cluster.fa",
			virus = virus),
		expand(config["results_dir"] + "/Total_fa/{virus}.format.rmdup.cluster.uniq.fa.clstr",
			virus = virus)
	output:
		config["results_dir"] + "/cluster.txt",
		directory(config["results_dir"] + "/Clusters_fa"),
		config["results_dir"] + "/cluster.identities.txt",
		config["results_dir"] + "/history.txt",
	params:
		script = config["scripts_dir"],
		max_seq = config["max_seq"],
		threshold = config["seq_number_ANI"],
		drop = config["drop"],
		ani = config["ani"]
	message:
		"Step5: Extract fasta from cd-hit results .."
	shell:
		'''
		python {params.script}/extract_cluster.py -i {input[0]} -c {input[1]} \
			 -m {params.max_seq} -o {output[0]} -y {output[2]} -d {output[1]};

		python {params.script}/merge_cluster_by_ANI.py -i {output[0]} -p 20 -t {params.threshold} \
			-o {output[3]} -d {params.drop} -a {params.ani}
		'''

#-------------------------------------------------------------------------------------------
# merge_cluster 6: Dependency packages - fastANI
#-------------------------------------------------------------------------------------------
#rule merge_cluster:
#	input:
#		config["results_dir"] + "/cluster.txt"
#	output:
#		config["results_dir"] + "/history.txt"
#	message:
#		"Step6: Merging cluster by ANI.."
#	shell:
#		'''
#		python {params.script}/merge_cluster_by_ANI.py -i {input} -p 20 -t 20 -o {output}
#		'''

#-------------------------------------------------------------------------------------------
# alignment_by_muscle rule 6: Dependency packages - None
#-------------------------------------------------------------------------------------------
rule alignment_by_muscle:
	input:
		config["results_dir"] + "/Clusters_fa/{i}.tfa"
	output:
		config["results_dir"] + "/Clusters_msa/{i}.tmsa"
	params:
		script = config["scripts_dir"]
	resources:
		mem_mb = 10000
	message:
		"Step6: Alignment by MAFFT/muscle .."
	shell:
		'''
		python {params.script}/run_mafft.py -i {input} -o {output}
		'''
#		muscle -in {input} -out {output}
#		for long input.
#		mafft --auto {input} > {output}
#-------------------------------------------------------------------------------------------
# multiPrime rule 7: Dependency packages - multiPrime-core
#-------------------------------------------------------------------------------------------
rule multiPrime:
	input:
		rules.alignment_by_muscle.output
	output:
		config["results_dir"] + "/Clusters_primer/{i}.top.primer.out"
	log:
		config["log_dir"] + "/multiPrime_{i}.log"
	resources:
		mem_mb = 10000
	params:
		script = config["scripts_dir"],
		dege_number = config["dege_number"],
		degeneracy = config["degeneracy"],
		primer_len = config["primer_len"],
		min_PCR_size = config["PRODUCT_size"].split(",")[0],
		variation = config["variation"],
		coordinate = config["coordinate"],
		GC = config["gc_content"],
		coverage = config["coverage"]
	message:
		"Step7: Design primers by multiPrime .."
	shell:
		'''
		python {params.script}/multiPrime-core_V15.py -i {input} -n {params.dege_number} \
			-d {params.degeneracy} -v {params.variation} -c {params.coordinate} \
			-g {params.GC} -s {params.min_PCR_size} -l {params.primer_len} \
			-o {output} -f {params.coverage} -p 1 2>&1 > {log}
		'''
#-------------------------------------------------------------------------------------------
# get_degePrimer rule 8: Dependency packages - pandas, biopython, math, operator,functools
#-------------------------------------------------------------------------------------------
rule get_multiPrime:
	input:
		primer = config["results_dir"] + "/Clusters_primer/{i}.top.primer.out",
		ref_fa = config["results_dir"] + "/Clusters_fa/{i}.tfa"
	output:
		config["results_dir"] + "/Clusters_cprimer/{i}.candidate.primers.txt"
	log:
		config["log_dir"] + "/get_multiPrime_{i}.log"
	resources:
		mem_mb = 10000
	params:
		script = config["scripts_dir"],
		fraction = config["coverage"],
		size = config["PRODUCT_size"],
		# maxseq=config["max_seq"],
		gc_content = config["gc_content"],
		distance = config["distance"],
		adaptor = config["adaptor"],
		end = config["end"]
	message:
		"Step8: Filter candidate primeri pairs for each cluster (hairpin, dimer (F-R) check) .."
	shell:
		'''
		python {params.script}/get_multiPrime.py -i {input.primer} -r {input.ref_fa} \
			-f {params.fraction} -s {params.size} -g {params.gc_content} -e {params.end} \
			-d {params.distance} -a {params.adaptor} -m 0\
			-o {output} -p 1 2>&1 > {log}
		'''

#-------------------------------------------------------------------------------------------
# aggregate_candidate_primers rule 9: Dependency packages - None
#-------------------------------------------------------------------------------------------
rule aggregate_candidate_primers:
	input:
		aggregate_input
	output:
		config["results_dir"] + "/Primers_set/candidate_primers_sets.txt"
	message:
		"Step9: Prepare all candidate primers for primer selection .."
	shell:
		'''
		cat {input} > {output}
		'''
#-------------------------------------------------------------------------------------------
# get_candidate_primer_fa rule 10: Dependency packages - None
#-------------------------------------------------------------------------------------------
rule get_candidate_primer_fa:
	input:
		config["results_dir"] + "/Primers_set/candidate_primers_sets.txt"
	output:
		config["results_dir"] + "/Primers_set/candidate_primers_sets.number",
		directory(config["results_dir"] + "/Primers_set/candidate_primers_sets")
	params:
		script = config["scripts_dir"],
		step = config["step"]
	message:
		"Step10: Get candidate primers .."
	shell:
		'''
		python {params.script}/candidate_primer_txt2fa.py -i {input} -s {params.step} \
			-n {output[0]} -o {output[1]}
		'''

#-------------------------------------------------------------------------------------------
# get_Maxprimerset rule 11: Dependency packages - python, pandas, biopython
#-------------------------------------------------------------------------------------------
rule get_Maxprimerset:
	input:
		config["results_dir"] + "/Primers_set/candidate_primers_sets.txt"
	output:
		config["results_dir"] + "/Primers_set/final_maxprimers_set.xls"
	log:
		config["log_dir"] + "/get_Maxprimerset.log"
	params:
		script = config["scripts_dir"],
		step = config["step"],
		method = config["method"]
	message:
		"Step11: Try to find Max_primer_set..."
	shell:
		'''
		python {params.script}/get_Maxprimerset.py -i {input} -s {params.step} \
			-m {params.method} -o {output} \
			 2>&1 > {log}
		'''
#-------------------------------------------------------------------------------------------
# get_core_primer_set rule 12: Dependency packages - python
#-------------------------------------------------------------------------------------------
rule get_core_primer_set:
	input:
		config["results_dir"] + "/Primers_set/candidate_primers_sets.txt"
	output:
		config["results_dir"] + "/Core_primers_set/core_candidate_primers_sets.txt"
	message:
		"step12: Extract core primer set..."
	params:
		script = config["scripts_dir"],
		number = config["core_number"]
	shell:
		"""
		python {params.script}/core_primerset_extraction.py -i {input} -o {output} \
			-n {params.number}
		"""
#-------------------------------------------------------------------------------------------
# get_core_Maxprimerset rule 13: Dependency packages - python
#-------------------------------------------------------------------------------------------
rule get_core_Maxprimerset:
	input:
		config["results_dir"] + "/Core_primers_set/core_candidate_primers_sets.txt"
	output:
		config["results_dir"] + "/Core_primers_set/core_final_maxprimers_set.xls"
	message:
		"Step13: Try to find core Max_primer_set .."
	params:
		script = config["scripts_dir"],
		step = config["step"],
		method = config["method"]
	log:
		config["log_dir"] + "/get_core_Maxprimerset.log"
	shell:
		"""
		python {params.script}/get_Maxprimerset.py -i {input} -s {params.step} \
			-m {params.method} -o {output} \
			 2>&1 > {log}
		"""
#-------------------------------------------------------------------------------------------
# format_transition rule 14: Dependency packages - python
#-------------------------------------------------------------------------------------------
rule format_transition:
	input:
		config["results_dir"] + "/Core_primers_set/core_candidate_primers_sets.txt"
	output:
		directory(config["results_dir"] + "/Core_primers_set/core_candidate_primers_sets"),
		config["results_dir"] + "/Core_primers_set/core_candidate_primers_sets.number"
	params:
		script = config["scripts_dir"],
		step = config["step"]
	message:
		"Step14: Format transition .."
	shell:
		"""
		python {params.script}/candidate_primer_txt2fa.py -i {input} -s {params.step} \
			-o {output[0]} -n {output[1]}
		"""
#-------------------------------------------------------------------------------------------
# get_all_PCR_product rule 15: Dependency packages - python
#-------------------------------------------------------------------------------------------
rule get_all_PCR_product:
	input:
		config["results_dir"] + "/Primers_set/final_maxprimers_set.xls",
		expand(config["results_dir"] + "/Total_fa/{virus}.format.fa",virus = virus)
	output:
		directory(config["results_dir"] + "/Primers_set/PCR_product"),
		config["results_dir"] + "/Primers_set/Coverage_stast.xls"
	params:
		config["scripts_dir"]
	message:
		"Step15: Extract PCR product from the input virus sequence (non-mismatch) .."
	shell:
		'''
		python {params}/extract_PCR_product.py -i {input[0]} -r {input[1]} -p 10 \
			-f xls -o {output[0]} -s {output[1]}
		'''
#-------------------------------------------------------------------------------------------
# get_core_PCR_product rule 16: Dependency packages - python
#-------------------------------------------------------------------------------------------
rule get_core_PCR_product:
	input:
		config["results_dir"] + "/Core_primers_set/core_final_maxprimers_set.xls",
		expand(config["results_dir"] + "/Total_fa/{virus}.format.fa",virus = virus)
	output:
		directory(config["results_dir"] + "/Core_primers_set/core_PCR_product"),
		config["results_dir"] + "/Core_primers_set/core_Coverage_stast.xls"
	params:
		config["scripts_dir"]
	message:
		"Step16: Extract core PCR product from the input virus sequence (non-mismatch) .."
	shell:
		'''
		python {params}/extract_PCR_product.py -i {input[0]} -r {input[1]} -p 10 \
			-f xls -o {output[0]} -s {output[1]}
		'''
#-------------------------------------------------------------------------------------------
# mfeprimer_check rule 17: Dependency packages - mfeprimer-3.2.6
#-------------------------------------------------------------------------------------------
rule all_mfeprimer_check:
	input:
		config["results_dir"] + "/Primers_set/final_maxprimers_set.xls"
	output:
		config["results_dir"] + "/Primers_set/final_maxprimers_set.fa",
		config["results_dir"] + "/Primers_set/final_maxprimers_set.fa.hairpin",
		config["results_dir"] + "/Primers_set/final_maxprimers_set.fa.dimer",
		config["results_dir"] + "/Primers_set/final_maxprimers_set.fa.findimer"
	params:
		config["scripts_dir"]
	message:
		"Step17: Hairpin and dimer check by mfeprimer .."
	
	shell:
		"""
		python {params}/primerset_format.py -i {input} -o {output[0]}
		{params}/mfeprimer-3.2.6 hairpin -i {output[0]} -o {output[1]}
		{params}/mfeprimer-3.2.6 dimer -i {output[0]} -o {output[2]}
		python {params}/finDimer.py -i {output[0]} -o {output[3]}
		"""
#-------------------------------------------------------------------------------------------
# core_mfeprimer_check rule 18: Dependency packages - mfeprimer-3.2.6
#-------------------------------------------------------------------------------------------
rule core_mfeprimer_check:
	input:
		config["results_dir"] + "/Core_primers_set/core_final_maxprimers_set.xls"
	output:
		config["results_dir"] + "/Core_primers_set/core_final_maxprimers_set.fa",
		config["results_dir"] + "/Core_primers_set/core_final_maxprimers_set.fa.hairpin",
		config["results_dir"] + "/Core_primers_set/core_final_maxprimers_set.fa.dimer",
		config["results_dir"] + "/Core_primers_set/core_final_maxprimers_set.fa.findimer"
	params:
		config["scripts_dir"]
	message:
		"Step18: Hairpin and dimer check by mfeprimer.. "
	shell:
		"""
		python {params}/primerset_format.py -i {input} -o {output[0]}
		{params}/mfeprimer-3.2.6 hairpin -i {output[0]} -o {output[1]}
		{params}/mfeprimer-3.2.6 dimer -i {output[0]} -o {output[2]}
		python {params}/finDimer.py -i {output[0]} -o {output[3]}
		"""
#-------------------------------------------------------------------------------------------
# core_primer_coverage rule 19: Dependency packages - bowtie2 
#-------------------------------------------------------------------------------------------
rule BWT_validation:
	input:
		config["results_dir"] + "/Core_primers_set/core_final_maxprimers_set.fa",
		expand(config["results_dir"] + "/Total_fa/{virus}.format.fa",virus = virus)
	output:
		config["results_dir"] + "/Core_primers_set/BWT_coverage/core_final_maxprimers_set.out"
	params:
		script = config["scripts_dir"],
		primer_len = config["primer_len"],
	message:
		"Step19: Primer coverage clculation .. "
	shell:
		"""
		python {params.script}/primer_coverage_validation_by_BWT.py -i {input[0]}  -r {input[1]} \
			-l {params.primer_len} -t 1 -s 50,2000 -o {output}
		"""

#-------------------------------------------------------------------------------------------
# Done!
#-------------------------------------------------------------------------------------------


