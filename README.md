# Kmer Filter + Mapping Pipeline

This pipeline filters paired-end FASTQ files against a k-mer reference using `bbduk`, maps the retained reads to a reference genome using `bwa mem`, generates sorted/indexed BAM files, and calculates both whole-genome and region-specific coverage statistics using `samtools`.

## Input Files

You only need 2 input files: a `config.yaml` specifying paths to required files and a `samplesheet.csv` containing sample names and FASTQ file paths.

Your samplesheet can have any columns you want, but it **MUST** contain the following fields named *exactly* `sample`, `read1`, and `read2` (case insensitive).

### Example `config.yaml`

```yaml
reference_genome: "/path/to/reference.fa"

kmer_reference: "/path/to/male_kmers.txt.fasta"

sdr_region: "scaffold_35:1-24227604"

output_dir: "results"

samplesheet: "/path/to/samplesheet.csv"
```

### Example `samplesheet.csv`

```csv
sample,read1,read2,expected_sex
sample1,/full/path/sample1_R1.fastq.gz,/full/path/sample1_R2.fastq.gz,M
sample2,/full/path/sample2_R1.fastq.gz,/full/path/sample2_R2.fastq.gz,U
```

## Pipeline Steps

For each sample, the pipeline performs:

1. K-mer filtering using `bbduk`
2. Reference indexing using `bwa index`
3. Mapping filtered reads with `bwa mem`
4. Sorting BAM files with `samtools sort`
5. Indexing BAM files with `samtools index`
6. Whole-genome coverage calculation
7. Region-specific coverage calculation

## Output Files

The pipeline creates three subdirectories within the specified `output_dir`.

### Filtered Reads

```text
results/filtered/
```

Example outputs:

```text
results/filtered/sample1_R1_m.fastq
results/filtered/sample1_R2_m.fastq
```

### BAM Files

```text
results/bam/
```

Example outputs:

```text
results/bam/sample1.sorted.bam
results/bam/sample1.sorted.bam.bai

results/bam/sample2.sorted.bam
results/bam/sample2.sorted.bam.bai
```

### Coverage Files

```text
results/coverage/
```

Example outputs:

```text
results/coverage/sample1.cov_out
results/coverage/sample1.sdr.cov

results/coverage/sample2.cov_out
results/coverage/sample2.sdr.cov
```

Where:

* `.cov_out` contains whole-genome coverage statistics generated with:

```bash
samtools coverage -m -w 100
```

* `.sdr.cov` contains coverage statistics for the specified region in `sdr_region`.

## Multiple Libraries

If a sample was sequenced across multiple libraries, treat each library as a separate sample in the samplesheet.

For example:

```csv
sample,read1,read2,expected_sex
sample1,/path/sample1_R1.fastq.gz,/path/sample1_R2.fastq.gz,M
sample2_l1,/path/sample2_l1_R1.fastq.gz,/path/sample2_l1_R2.fastq.gz,M
sample2_l2,/path/sample2_l2_R1.fastq.gz,/path/sample2_l2_R2.fastq.gz,M
```

Each row must have a unique sample identifier.

## Software Requirements

The following software must be available in your environment:

```text
bbduk
bwa
samtools
snakemake
```

Recommended versions:

```text
bwa >= 0.7.17
samtools >= 1.15
bbmap >= 39
snakemake >= 7
```

## How to Run

### 1. Clone the repository

```bash
git clone https://github.com/Armstrong-Lab-UCR/kmer-filter-map.git
```

### 2. Create your `config.yaml` and `samplesheet.csv`

These can be stored either inside or outside of the repository directory.

### 3. Start a tmux session

This ensures the pipeline continues running if you disconnect.

```bash
module load tmux/3.3

tmux new -s [session_name]
```

### 4. Load Snakemake

```bash
module load snakemake/7.18
```

### 5. Navigate to the pipeline directory

```bash
cd kmer-filter-map
```

### 6. Dry run the workflow

```bash
snakemake --configfile /path/to/config.yaml --profile profiles/slurm -n
```

### 7. Run the workflow

```bash
snakemake --configfile /path/to/config.yaml --profile profiles/slurm
```

## tmux Notes

Detach from a session:

```bash
Ctrl+B, then D
```

Reconnect to a session:

```bash
tmux a -t [session_name]
```

List active sessions:

```bash
tmux list-sessions
```
