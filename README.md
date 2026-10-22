<!--
#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
-->

# Project

This is a package to help with the analysis of kallisto outputs.
One particularly useful application is to load and merge outputs of different samples.

From `path_to_kallisto_output`; the output path of `kallisto quant` (`-o` argument), you can merge all the samples outputs and read it as a pandas DataFrame as follows:

```python
from kallistomanager.loader import merge_kallisto_counts

merge_kallisto_counts(save_dir=path_to_kallisto_output, data_path=path_to_kallisto_output)
```

This writes `counts_matrix.csv`, `tpm_matrix.csv`, `gene_length_matrix.csv` and `gene_eff_length_matrix.csv` to
`save_dir`. `from_whole_transcriptome_quantification_to_coding_genes` then sums the counts per protein-coding gene.


## Installation of this package in another project

To install this library in another Python project, execute simply:

```bash
pip install git+https://github.com/AI-for-RNA-Biology/kallistomanager.git
```

## State of the project

Some functions have been developed by Christelle Schneuwly for analysis of kallisto outputs from SMART-seq data.
Note: the poetry environment has to be solved. In particular scanpy is not installed due to some conflicts.

# How to use the package

## Coverage

This is to visualize in IGV the coverage output from kallisto quant, when we used a transcriptome annotation at the transcript level and not at the chromosome level.

**Step 1:** Include the arrgument `--pseudobam` when running kallisto quant as follow:

```bash
kallisto quant --index path_to_kallisto_index -o path_to_output --pseudobam file1.fastq file2.fastq
```

**Step 2:** Sort the bam file using `samtools` as follow:

```bash
samtools sort pseudoalignments.bam -o pseudoalignments_sorted.bam
```

**Step 3:** Use `igvtools` to create a first coverage file, but based on transcripts instead of chromosomes. Then convert it into a bedfile:

```bash
igvtools count pseudoalignments_sorted.bam transcripts.tdf mygenome.chrom.sizes
igvtools tdftobedgraph transcripts.tdf transcripts.bed
```

**Step 4:** Sort the resulting bed file using `bedtools`:

```bash
bedtools sort transcripts.bed -g genome.transcripts.sizes > transcripts_sorted.bed
```

**Step 5:** Use the script `scripts/transcript_to_chromosome_coverage.py` of
[genomemanager](https://github.com/AI-for-RNA-Biology/genomemanager) to convert the coverage file from transcripts to chromosomes:

```bash
python scripts/transcript_to_chromosome_coverage.py --annotation_path path_to_gtf --bed_transcripts_path transcripts_sorted.bed --output chromosomes_coverage.bed
```

**Step 6:** Sort the resulting bedfile:

```bash
bedtools sort chromosomes_coverage.bed -g genome.chrom.sizes > chromosomes_coverage_sorted.bedgraph
```
