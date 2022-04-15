# import something (if necessary)

'''
workflow management for all post GutFunFind processes
'''

# definte wildcards
function = ["Azoreductase"]
mgx_dataset = ["hmp2", "prism"]

# NOTE: Snakefile is currently meant for archival purposes and is not fully functioning (yet :)

##### =========== index and align azoreductase genes using Salmon =========== #### 
rule index_azored_genes:
    input: "results/from-GutFunFind/{function}.gene_seqs.fna"
    output: directory("results/salmon_out/{function}.gene_seqs_index")
    shell: "code/cluster/index_ntseqs.sh {input} {output}"

rule quant_azored_genes:
    input: 
        paths = "processed/{mgx_dataset}_readfile_paths.txt",
        index = "processed/salmon_out/{function}.gene_seqs_index"
    output: "results/salmon_out/{dataset}/dummy/{function}.quantify.txt"
    shell: "sbatch code/cluster/quant_{mgx_dataset}_reads.sh {input.paths} {input.index} {output}"

rule process_salmon_out:
	input: 
	output:
		tpm = "results/salmon_processed/{dataset}/{function}_tpm.tsv",
		counts = "results/salmon_processed/{dataset}/{function}_counts.tsv"
	shell: "python code/cluster/process_salmon_out.py {dataset} {function}"

##### ============ index and align azoreductase genes using KMA ============ #### 

rule index_genes_kma:
    input: "results/from-GutFunFind/{function}.gene_seqs.fna"
    output: directory("processed/kma_out/index/{function}.gene_seqs_index")
    shell: "code/cluster/index_ntseqs_kma.sh {input} {output}"

rule align_reads_kma:
    input: 
        paths = "processed/{mgx_dataset}_readfile_paths.txt",
        index = "processed/kma_out/{function}.gene_seqs_index"
    output: "results/kma_out/{mgx_dataset}"
    shell: "code/cluster/align_{mgx_dataset}_kma.sh {input.paths} {input.index} {output}"
