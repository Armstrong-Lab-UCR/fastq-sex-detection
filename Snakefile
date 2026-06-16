import pandas as pd

OUTPUT_DIR = config["output_dir"]
SAMPLESHEET = config["samplesheet"]
REFERENCE = config["reference_genome"]
KMER_REF = config["kmer_reference"]
SDR_REGION = config["sdr_region"]

try:
    samples = pd.read_csv(SAMPLESHEET)
    samples.columns = samples.columns.str.lower()
    samples = samples.set_index("sample", drop=False)
except FileNotFoundError:
    raise FileNotFoundError(f"Samplesheet file not found at: {SAMPLESHEET}")
except KeyError:
    raise KeyError("The samplesheet must contain a 'sample' column.")

if not samples.index.is_unique:
    raise ValueError("Sample names in the samplesheet are not unique.")

if samples["read1"].duplicated().any() or samples["read2"].duplicated().any():
    raise ValueError("Read names (read1 or read2) in the samplesheet are not unique.")

SAMPLES = samples.index.tolist()


rule all:
    input:
        expand(f"{OUTPUT_DIR}/filtered/{{sample}}_R1_m.fastq", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/filtered/{{sample}}_R2_m.fastq", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/bam/{{sample}}.sorted.bam", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/bam/{{sample}}.sorted.bam.bai", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/coverage/{{sample}}.cov_out", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/coverage/{{sample}}.sdr.cov", sample=SAMPLES),


rule bwa_index:
    input:
        reference=REFERENCE
    output:
        amb=REFERENCE + ".amb",
        ann=REFERENCE + ".ann",
        bwt=REFERENCE + ".bwt",
        pac=REFERENCE + ".pac",
        sa=REFERENCE + ".sa"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=600
    shell:
        """
        bwa index -a bwtsw {input.reference}
        """


rule bbduk_filter:
    input:
        read1=lambda wildcards: samples.loc[wildcards.sample, "read1"],
        read2=lambda wildcards: samples.loc[wildcards.sample, "read2"],
        kmer_ref=KMER_REF
    output:
        read1=f"{OUTPUT_DIR}/filtered/{{sample}}_R1_m.fastq",
        read2=f"{OUTPUT_DIR}/filtered/{{sample}}_R2_m.fastq"
    threads: 4
    resources:
        mem_mb=32000,
        runtime=600
    shell:
        """
        mkdir -p {OUTPUT_DIR}/filtered

        bbduk.sh \
            in={input.read1} \
            in2={input.read2} \
            ref={input.kmer_ref} \
            outm={output.read1} \
            outm2={output.read2} \
            k=21
        """


rule bwa_mem_sort:
    input:
        read1=f"{OUTPUT_DIR}/filtered/{{sample}}_R1_m.fastq",
        read2=f"{OUTPUT_DIR}/filtered/{{sample}}_R2_m.fastq",
        reference=REFERENCE,
        bwa_index=rules.bwa_index.output
    output:
        bam=f"{OUTPUT_DIR}/bam/{{sample}}.sorted.bam"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=1200
    shell:
        """
        mkdir -p {OUTPUT_DIR}/bam

        bwa mem -M -t {threads} {input.reference} {input.read1} {input.read2} | \
            samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam} -
        """


rule samtools_index:
    input:
        bam=f"{OUTPUT_DIR}/bam/{{sample}}.sorted.bam"
    output:
        bai=f"{OUTPUT_DIR}/bam/{{sample}}.sorted.bam.bai"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=300
    shell:
        """
        samtools index {input.bam}
        """


rule whole_genome_coverage:
    input:
        bam=f"{OUTPUT_DIR}/bam/{{sample}}.sorted.bam",
        bai=f"{OUTPUT_DIR}/bam/{{sample}}.sorted.bam.bai"
    output:
        cov=f"{OUTPUT_DIR}/coverage/{{sample}}.cov_out"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=300
    shell:
        """
        mkdir -p {OUTPUT_DIR}/coverage

        samtools coverage {input.bam} \
            -o {output.cov} \
            -m \
            -w 100
        """


rule sdr_coverage:
    input:
        bam=f"{OUTPUT_DIR}/bam/{{sample}}.sorted.bam",
        bai=f"{OUTPUT_DIR}/bam/{{sample}}.sorted.bam.bai"
    output:
        cov=f"{OUTPUT_DIR}/coverage/{{sample}}.sdr.cov"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=300
    shell:
        """
        mkdir -p {OUTPUT_DIR}/coverage

        samtools coverage \
            -r {SDR_REGION} \
            {input.bam} \
            -o {output.cov}
        """
