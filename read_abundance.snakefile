from os.path import join, abspath, expanduser
import sys

################################################################################
# specify project directories
PROJECT_DIR = config["output_directory"]
# convert PROJECT_DIR to absolute path
if PROJECT_DIR[0] == '~':
    PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)

# get input samples from sample table
with open(config["sample_table"]) as inf:
    insamps = [i for i in inf.readlines() if i != '\n']
    sample_dict = {sample: files.split(",") for sample, files in [l.strip().split("\t") for l in insamps]}
# ensure no comment lines
sample_dict = {k:sample_dict[k] for k in sample_dict.keys() if k[0] != '#'}

# get list of samples
sample_list = list(sample_dict.keys())

rule align_reads:
    input:
        #from config file, pulls out the value of align reads -> host genome
        index_amb = config['align_reads']['host_genome'] + '.amb',
        index_ann = config['align_reads']['host_genome'] + '.ann',
        index_bwt = config['align_reads']['host_genome'] + '.bwt',
        index_pac = config['align_reads']['host_genome'] + '.pac',
        index_sa = config['align_reads']['host_genome'] + '.sa',
        #figure out how to input this from samples table
        fwd       = sample_dict[wildcards.sample][0],
        rev       = sample_dict[wildcards.sample][1],
        orp       = sample_dict[wildcards.sample][2]
    #this will actually output mapped
    output:
        unmapped_1 = join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq.gz"),
        unmapped_2 = join(PROJECT_DIR, "01_processing/05_sync/{sample}_2.fq.gz"),
        unmapped_singletons = join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"),
    params:
        bwa_index_base = join(config['align_reads']['host_genome']),
        singelton_temp_1 = join(PROJECT_DIR, "01_processing/04_host_align/{sample}_rmHost_singletons1.fq.gz"),
        singelton_temp_2 = join(PROJECT_DIR, "01_processing/04_host_align/{sample}_rmHost_singletons2.fq.gz"),
    threads: 4
    resources:
        mem_mb=32000,
        mem=32,
        time=24
    singularity: "shub://bsiranosian/bens_1337_workflows:align"
    # conda: "envs/align.yaml"
    benchmark: join(PROJECT_DIR, "01_processing/04_host_align/{sample}_time.txt")
    shell: """
        mkdir -p {PROJECT_DIR}/01_processing/04_host_align/
        # if an index needs to be built, use bwa index ref.fa
        # run on paired reads
        bwa mem -t {threads} {params.bwa_index_base} {input.fwd} {input.rev} | \
            samtools view -F 4 -b > {output.unmapped_1} -2 {output.unmapped_2} -s {params.singelton_temp_1} -
        # run on unpaired reads
        bwa mem -t {threads} {params.bwa_index_base} {input.orp} | \
            samtools fastq -@ {threads} -t -T BX -f 4 - > {params.singelton_temp_2}
        # combine singletons
        zcat -f {params.singelton_temp_1} {params.singelton_temp_2} | pigz -p {threads} > {output.unmapped_singletons}
        rm {params.singelton_temp_1} {params.singelton_temp_2}
    """

    # #sorts read alignments
    # rule samtools_sort:
    #     input:
    #         "mapped_reads/{sample}.bam"
    #     output:
    #         "sorted_reads/{sample}.bam"
    #     shell:
    #         "samtools sort -T sorted_reads/{wildcards.sample} "
    #         "-O bam {input} > {output}"
    #
    # rule readcounts:
    # # htseq-count -r pos -t CDS -f bam $SAMPLE.map.markdup.bam $SAMPLE.map.gtf > $SAMPLE.count
    # input:
    #     raw = expand(join(DATA_DIR, "{sample}_" + READ_SUFFIX[0] + EXTENSION), sample=SAMPLE_PREFIX),
    #     dedup = expand(join(PROJECT_DIR, "01_processing/01_dedup/{sample}_1.fq.gz"), sample=SAMPLE_PREFIX),
    #     trimmed = expand(join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_1_val_1.fq" + gz_ext), sample=SAMPLE_PREFIX),
    #     rmhost = expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_1.fq.gz"), sample=SAMPLE_PREFIX),
    #     orphans = expand(join(PROJECT_DIR, "01_processing/05_sync/{sample}_orphans.fq.gz"), sample=SAMPLE_PREFIX)
    # output:
    #     join(PROJECT_DIR, "01_processing/readcounts.tsv")
    # resources:
    #     time = 24
    # run:
    #     outfile = str(output)
    #     if (os.path.exists(outfile)):
    #         os.remove(outfile)
    #     with open(outfile, 'w') as outf:
    #         outf.writelines('Sample\traw_reads\tdedup_reads\tdedup_frac\ttrimmed_reads\ttrimmed_frac\thost_removed_reads\thost_removed_frac\torphan_reads\torphan_frac\n')
    #         for sample in SAMPLE_PREFIX:
    #             raw_file = join(DATA_DIR, sample + "_" + READ_SUFFIX[0] + EXTENSION)
    #             dedup_file = join(PROJECT_DIR, "01_processing/01_dedup/" + sample + "_1.fq.gz")
    #             trimmed_file = join(PROJECT_DIR, "01_processing/02_trimmed/" + sample + "_1_val_1.fq" + gz_ext)
    #             rmhost_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_1.fq")
    #             orphans_file = join(PROJECT_DIR, "01_processing/05_sync/" + sample + "_orphans.fq")
    #
    #             raw_reads = int(file_len(raw_file) / 4)
    #             dedup_reads = int(file_len(dedup_file) / 4)
    #             trimmed_reads = int(file_len(trimmed_file) / 4)
    #             rmhost_reads = int(file_len(rmhost_file) / 4)
    #             orphans_reads = int(file_len(orphans_file) / 4)
    #
    #             dedup_frac = round(dedup_reads / float(raw_reads), 3)
    #             trimmed_frac = round(trimmed_reads / float(raw_reads), 3)
    #             rmhost_frac = round(rmhost_reads / float(raw_reads), 3)
    #             orphans_frac = round(orphans_reads / float(raw_reads), 3)
    #
    #             line = '\t'.join([sample, str(raw_reads),
    #                 str(dedup_reads), str(dedup_frac),
    #                 str(trimmed_reads), str(trimmed_frac),
    #                 str(rmhost_reads), str(rmhost_frac),
    #                 str(orphans_reads), str(orphans_frac)])
    #             outf.writelines(line+ '\n')
