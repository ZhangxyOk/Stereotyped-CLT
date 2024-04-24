import os
import sys
import glob
import time
## load config file
configfile: "/mnt/data/home/lzz/project/2021-4-9_lung_scRNA-seq/scripts/2021-4-9_pipline/PacBio_Analysis_byStrand_config_4.9.yaml"

# # input raw BAM file
raw_file_path = config["raw_subreads_path"]

# 10x cellranger output path
Seurat_path = config["Seurat_path"]

# filter BC base on 10x cellrangers
filterBC_dir = config["filterBC_dir"]
exp_mx_dir = config["exp_mx_dir"]

# reference 
bowtie_index = config["bowtie_index"]
ref_seq = config["ref_seq"]
tarPos = config["tarPos"]
iqtree_transM = config["iqtree_transMatrix"]

# some hand write scripts
# split_strand = config["split_strand"]
HigQulCCS = config["HigQulCCS"]
FilterCCS = config["FilterCCS"]
MergeR1R2CCS = config["MergeR1R2CCS"]
BCumi_con = config["BCumi_con"]
Align_toRef = config["Align_toRef"]
Call_Events = config["Call_Events"]
BC_conEvents = config["BC_conEvents"]
Events_toIqtreeInput = config["Events_toIqtreeInput"]
Run_Iqtree = config["Run_Iqtree"]
# optional
demultiplexing_sample = config["demultiplexing_sample"]
AdjustBC = config["AdjustBC"]

stainfo = config["stainfo"]
plotSta = config["plotSta"]
# R scripts
getPair_dist = config["getPair_dist"]
getExpMatrix = config["getExpMatrix"]
getTreeVSexpDist = config["getTreeVSexpDist"]
analysis_share = config["analysis_share"]

## some parameters
FP = config["P2"]
RP = config["P4"]
inFP = config["inFP"]
inRP = config["inRP"]
primerMis = config["primerMismatch"]
pAwinS = config["polyAwindowSize"]
pAmis = config["polyAmismatch"]
ref_len = config["ref_len"]

# demultiplexing sample bc
#Demul_BC = config["Demul_BC"]

# input raw data path list
def get_sample_name(sample_path, suffix=".subreads.bam"):
    return [file for file in os.listdir(sample_path) if file.endswith(suffix)]

raw_subreads_sample = get_sample_name(raw_file_path)
raw_subreads = os.path.join(raw_file_path, raw_subreads_sample[0])
# sample_name use to define outputs 
sample_name = ["A1-CBRAD5", "E5-HESC", "F11-HESC", "G11-CBRAD5", "G2-CBRAD5", "GS-HESC"]

# output directory
outDir = config["outDir"]
sample_out_paths = [os.path.join(outDir, j) for j in sample_name]
for each in sample_out_paths:
    if not os.path.exists(each):
        os.makedirs(each)

rule all:
    input: 
        expand(outDir + "/{sample}/{sample}.staPlots.pdf", sample = sample_name)
        #outDir + "/" + raw_subreads_sample[0].split(".subreads.bam")[0] + ".byStrand.ccs.raw.bam"

#rule runCCSbyStrand:
#    input:
#        raw_subreads
#    params:
#        raw_prefix = raw_subreads_sample[0].split(".subreads.bam")[0]
#    output:
#        outDir + "/" + raw_subreads_sample[0].split(".subreads.bam")[0] + ".byStrand.ccs.raw.bam",
#        outDir + "/" + raw_subreads_sample[0].split(".subreads.bam")[0] + ".byStrand.ccs.report.txt"
#    threads: 80
#    shell:
#        "ccs --report-file {output[1]} --min-passes 3 -j {threads} --by-strand {input[0]} {output[0]}"
#
#rule Demultiplexing_sample:
#    input:
#        rules.runCCSbyStrand.output[0]
#    output:
#        outDir + "/byStrand.ccs.raw.bam",
#        outDir + "/byStrand.ccs.raw.{sample}--{sample}.bam"
#    params:
#        out_prefix = "byStrand.ccs.raw.bam"
#    threads: 80
#    shell:
#        "lima --ccs --split-bam-named -j {threads} {input[0]} {Demul_BC} {output[0]}"

# rule get_filterBC_file:
#     input:
#         expand("{barcode_file}", barcode_file = Seurat_path.values())
#     output:
#         filterBC_dir + "/{sample}.filterBCs.txt"
#     run:
#         import os
#         import sys
#         for s, p in Seurat_path.items():
#             cmd = "zcat %s | awk -F '-' '{print $1}' > %s" % (p, "{}/{}.filterBCs.txt".format(filterBC_dir, s))
#             os.system(cmd)
    
rule ccsBAMtoFASTQ:
    input:
        outDir + "/byStrand.ccs.raw.{sample}--{sample}.bam"
        #filterBC_dir + "/{sample}.filterBCs.txt"
    output:
        outDir + "/{sample}/{sample}.ccs.raw.fastq"
    shell:
        "bedtools bamtofastq -i {input[0]} -fq {output[0]}"


rule getNoGenomeCCS_FASTQ:
    input:
        rules.ccsBAMtoFASTQ.output
    output:
        outDir + "/{sample}/{sample}.ccs.noG.fastq",
        outDir + "/{sample}/{sample}.ccs.unmap.sam",
        outDir + "/{sample}/{sample}.bowtie2_report.txt"
    threads: 20
    message: "Align to the human genome, exclude the nonspecific amplification"
    shell:
        """
        bowtie2 --threads {threads} -x {bowtie_index} --un {output[0]} -q {input[0]} -S {output[1]} > {output[2]} 2>&1
        """

# rule getHightestQualCCSforEachZwID:
#     input:
#         rules.getNoGenomeCCS_FASTQ.output[0]
#     output:
#         outDir + "/{sample}/{sample}.ccs.clean.fastq"
#     shell:
#         """
#         python {HigQulCCS} -i {input[0]} -o {output[0]}
#         """

rule filterCCSandExtractBCUMI:
    input:
        outDir + "/{sample}/{sample}.ccs.noG.fastq"
    output:
        outDir + "/{sample}/{sample}.filter.fa",
        outDir + "/{sample}/{sample}.staDF.txt",
        outDir + "/{sample}/{sample}.noPass.fa"
    threads: 20
    shell:
        """
        python {FilterCCS} -i {input[0]} -o {output[0]} -s {output[1]} --nopass {output[2]} \
            --format FASTQ --FP {FP} --RP {RP} --inFP {inFP} --inRP {inRP} --pAw {pAwinS} --pMs {pAmis} --mismatch {primerMis} --numcpu {threads}
        """

rule AdjustBCandUMI:
    input:
        rules.filterCCSandExtractBCUMI.output
    output:
        outDir + "/{sample}/{sample}.adjust.fa"
    params:
        filterBC_file=filterBC_dir + "/{sample}_filterBCs.txt"
    shell:
        "python {AdjustBC} -i {input[0]} -o {output[0]} -f {params.filterBC_file}"

rule getUMIconSeqEachBC:
    input:
        rules.AdjustBCandUMI.output
    output:
        outDir + "/{sample}/{sample}.UMIcon.fa"
    threads: 20
    shell:
        "python {BCumi_con} -i {input} -o {output} --numcpu {threads}"
    
rule AlignToRef:
    input:
        rules.getUMIconSeqEachBC.output
    output:
        outDir + "/{sample}/{sample}.UMIcon.align"
    threads: 40
    shell:
        "python {Align_toRef} -r {ref_seq} -s {input} -o {output} -c {threads}"

rule CallEditEventsFromAlignment:
    input:
        rules.AlignToRef.output
    output:
        outDir + "/{sample}/{sample}.Events.txt"
    threads: 40
    params: 
        output_dir = outDir + "/{sample}",
        name = "{sample}"
    shell:
        "python {Call_Events} -a {input} -t {tarPos} -n {params.name} --dir {params.output_dir} --cpu {threads}"

rule getBCconsensusEditEvents:
    input:
        rules.CallEditEventsFromAlignment.output
    output:
        outDir + "/{sample}/{sample}_comSCOREandUMIconEvents.txt",
        outDir + "/{sample}/{sample}_moreInfos.txt",
        outDir + "/{sample}/{sample}_BCstaInfos.txt"
    params:
        output_dir = outDir + "/{sample}",
        outprefix = "{sample}",
        filterBC_file=filterBC_dir + "/{sample}_filterBCs.txt"
    shell:
        "python {BC_conEvents} -i {input} -n {params.outprefix} --outDir {params.output_dir} -f {params.filterBC_file}"

rule eventsToIqtreeInput:
    input:
        rules.getBCconsensusEditEvents.output
    output:
        outDir + "/{sample}/{sample}.nwk",
        outDir + "/{sample}/{sample}.AllelesInfo.txt"
    params:
        output_dir = outDir + "/{sample}",
        outprefix = "{sample}",
        filterBC_file = filterBC_dir + "/{sample}_filterBCs.txt"
    threads:
        40
    shell:
        """
        python {Events_toIqtreeInput} -e {input[0]} -m {iqtree_transM} -d {params.output_dir} -n {params.outprefix} -l {ref_len} -f {params.filterBC_file}
        """

rule getTreeNodePairDist:
    input:
        rules.eventsToIqtreeInput.output
    output:
        outDir + "/{sample}/{sample}.oneCellNode.pair_Dist.txt"
    threads:
        40
    shell:
        """
        Rscript {getPair_dist} -i {input[0]} -o {output[0]}
        """

rule getExpMx:
    input:
        outDir + "/{sample}/{sample}.AllelesInfo.txt",
    output:
        outDir + "/{sample}/{sample}_oneCellexp.Rds",
        outDir + "/{sample}/{sample}_alltreeExp.Rds"
    params:
        filterBC_file = filterBC_dir + "/{sample}_filterBCs.txt",
        exp_mx = exp_mx_dir + "/{sample}_allCell_exp.Rds"
    shell:
        """
        Rscript {getExpMatrix} -i {params.exp_mx} -f {params.filterBC_file} -a {input[0]} --outOneCell {output[0]} --outAllTree {output[1]}
        """

rule getTreeVSexpDist:
    input:
        outDir + "/{sample}/{sample}_oneCellexp.Rds",
        outDir + "/{sample}/{sample}.oneCellNode.pair_Dist.txt"
    output:
        outDir + "/{sample}/{sample}_treeVSexpDist_result.Rds"
    shell:
        """
        Rscript {getTreeVSexpDist} --inTreeDist {input[1]} --inExpMx {input[0]} -o {output[0]}
        """

rule analysis_ANC_share:
    input:
        outDir + "/{sample}/{sample}.AllelesInfo.txt",
        outDir + "/{sample}/{sample}_moreInfos.txt",
        outDir + "/{sample}/{sample}.oneCellNode.pair_Dist.txt",
        outDir + "/{sample}/{sample}_treeVSexpDist_result.Rds"
    output:
        outDir + "/{sample}/{sample}_ANC_share_analysis_results.Rds"
    threads:
        60
    shell:
        "Rscript {analysis_share} --allele {input[0]} --moreInfo {input[1]} --pairDist {input[2]} -o {output}"

rule allele_staInfo:
    input:
        outDir + "/{sample}/{sample}.AllelesInfo.txt"
    output:
        outDir + "/{sample}/{sample}.Statistics.txt",
        outDir + "/{sample}/{sample}.EventInfo.txt"
    params: 
        out_name = "{sample}",
        output_dir = outDir + "/{sample}"
    shell:
        "python {stainfo} -i {input[0]} -n {params.out_name} -o {params.output_dir}"

rule plot_StaInfo:
    input:
        rules.allele_staInfo.output
    output:
        outDir + "/{sample}/{sample}.staPlots.pdf"
    shell:
        "Rscript {plotSta} -i {input[0]} -o {output[0]}"

## print some messages when success
onsuccess:
    print("Workflow finished, no error.")



