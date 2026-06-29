# Fun_ncRNA

Predicts core non-coding RNAs in fungal genomes using covariation models fine-tuned for fungi.

`Fun_ncRNA` is a command-line pipeline for fungal core ncRNA annotation. It searches fungal genome assemblies with curated/fine-tuned covariation models, integrates optional predictions from RNAmmer and tRNAscan-SE, removes overlapping/redundant predictions, and produces GFF files and a prediction summary.
Currently supported core ncRNAs are snRNA, snoRNA, RNase P, RNase MRP, SRP, rRNA and tRNA.

Simple use
```bash
perl Fungi_ncRNA.pl -fna genome.fna [options]
```

## Main features

* Predicts a set of fungal core non-coding RNAs from genome FASTA files.
* Uses fungal class-specific covariation models when requested.
* Supports optional comparison against an existing genome annotation in GFF format.
* Supports optional integration of RNAmmer rRNA predictions and tRNAscan-SE tRNA predictions.
* Produces GFF-formatted ncRNA annotations and a summary table suitable for publication.
* Can reanalyse a previous output directory to add optional RNAmmer/tRNAscan-SE results.
* Produces summary outputs for publication or submission

## Key outputs

* Machine readable gff3 formatted nc_RNA and combined files for downstream analysis (input_genome.ncrna.compatible.gff) 
* Simple summary (summary.txt). For example: Core ncRNAs were predicted with Fun_ncRNA (Chyou et al 2026). There were 11 snRNA, 35 snoRNA, 1 RNase P, 1 RNase MRP, 1 SRP, 268 tRNA and 23 rRNA predicted. Of the 23 rRNA, there were 2 5S rRNA, 12 LSU, 3 5.8S rRNA, and 6 SSU predicted.


## Repository structure

Important files and directories include:

```text
Fun_ncRNA/
├── Fungi_ncRNA.pl              # Main pipeline script
├── Fungi_ncRNA_search.sh       # Main search workflow
├── Fungi_ncRNA_recomp.sh       # Reanalysis workflow
├── Fungi_ncRNA_sep_merge.sh    # Optional RNAmmer/tRNAscan-aware merge workflow
├── GFFCompare_wrapper.pl       # Wrapper for optional GFFCompare analysis
├── PredictionSummary.pl        # Generates summary.txt
├── Rfam/                       # Default Rfam subset used by the pipeline
├── core_cm/                    # Core ncRNA covariation models
├── sno_u3_cm/                  # snoRNA/U3 covariation models
├── SRP_Basidio/                # SRP-related models
├── data/                       # Model lists and supporting metadata
└── example_output/             # Example output, if included
```

## Requirements

### Required dependencies

The following tools are required for the core pipeline:

* Perl
* Infernal
* bedtools
* gffread
* gffcompare

`Infernal` provides `cmsearch`, which is used for covariation-model searches. `bedtools`, `gffread`, and `gffcompare` are used for genomic interval and GFF/GTF processing.

### Optional dependencies

The following tools are optional:

* tRNAscan-SE
* RNAmmer
* HMMER/HMMER2, required by RNAmmer depending on the RNAmmer version/setup

`tRNAscan-SE` is used for additional tRNA prediction. RNAmmer is used for additional rRNA prediction. By default, these optional tools are not run unless `-run_rnammer_trnascan` is supplied.

## Installation and use

### 1. Clone the repository

```bash
git clone https://github.com/davidchyou/Fun_ncRNA.git
cd Fun_ncRNA
```

### 2. Install required dependencies with conda

Most dependencies can be installed from Bioconda. One convenient approach is to create a dedicated conda environment:

```bash
conda create -n fun_ncrna \
  -c conda-forge -c bioconda \
  perl infernal bedtools gffread gffcompare
```

Activate the environment:

```bash
conda activate fun_ncrna
```

Check that the required programs are working:

```bash
cmsearch -h
bedtools --version
gffread --version
gffcompare --version
perl -v
```

### 3. Install optional dependencies

To include tRNAscan-SE:

```bash
conda install -n fun_ncrna -c conda-forge -c bioconda trnascan-se
```

To install HMMER through Bioconda for RNAMMER:

```bash
conda install -n fun_ncrna -c conda-forge -c bioconda hmmer hmmer2 perl-xml-simple
```

RNAmmer is not distributed through Bioconda and must be downloaded and installed separately according to the developer’s instructions. RNAmmer may require manual configuration of its installation path and HMMER/HMMER2 path. After installation, make sure the RNAmmer executable is available in your `PATH`, or edit the relevant RNAmmer configuration files as instructed by the RNAmmer documentation.

Check optional tools with:

```bash
tRNAscan-SE -h
rnammer -h
hmmsearch -h
```

Depending on your RNAmmer installation, the required `hmmsearch` may be from HMMER2 rather than the current HMMER3 package. If RNAmmer fails while calling `hmmsearch`, check the RNAmmer configuration and install the HMMER version expected by your RNAmmer release.

### 4. Run Fun_ncRNA on a genome in fasta format.

The minimum required input is a fungal genome assembly in FASTA format:

```bash
perl Fungi_ncRNA.pl -fna genome.fna
```

This writes output to the default directory:

```text
Fungi_ncRNA/
```

A more typical run with a genome FASTA, existing annotation GFF, fungal class, and custom output directory is:

```bash
perl Fungi_ncRNA.pl \
  -fna genome.fna \
  -gff genome.gff \
  -fungal_class Sordariomycetes \
  -out Fun_ncRNA_output
```
(see section on the main output files below)

## Additional options
To also run RNAmmer and tRNAscan-SE if these options have been installed, summarise results by an assembly ID. 
Run GFFCompare and compare to an existing annotation:

```bash
perl Fungi_ncRNA.pl \
  -fna data/GCF_024613125.1_Pisori2_genomic.fna \
  -gff data/GCF_024613125.1_Pisori2_genomic.gff \
  -fungal_class Agaricomycetes \
  -gcf GCF_024613125.1 \
  -out Fun_ncRNA_Pisori2_output \
  -gffcmp  \
  -run_rnammer_trnascan
```

```bash
perl Fungi_ncRNA.pl \
  -fna genome.fna \
  -gff genome.gff \
  -out Fun_ncRNA_output \
  -gffcmp
```

## Command-line options

| Option                  |       Required?     | Default          | Description                                                                                                                                                         |
| ----------------------- | ------------------: | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-fna`                  |                 Yes | `NA`             | Input genome assembly in FASTA format.                                                                                                                              |
| `-gff`                  | No, but recommended | `NA`             | Input existing genome annotation in GFF format. Used for combined output and optional GFFCompare analysis.                                                                |
| `-fungal_class`         |                  No | `Fungi`          | Fungal taxonomic class (see below). If supplied, class-specific covariation models are used where available. Otherwise, fungi-wide models are used.                             |
| `-run_rnammer_trnascan` |                  No | off              | Run optional RNAmmer and tRNAscan-SE predictions.                                                                                                                   |
| `-gffcmp`               | No, but recommended | off              | Run GFFCompare. Skipped if no GFF annotation is provided.                                                                                                           |
| `-out`                  |                  No | `Fungi_ncRNA`    | Output directory. Existing directories with the same name will be removed and recreated (overwritten).                                                                            |
| `-gcf`                  |                  No | `NA`             | Assembly ID for multi-contig genome files.  Appears as the Assembly name in the summary output.                                                                                                                        |
| `-core_cm_cutoff`       |                  No | `60`             | `cmsearch` score cutoff for fungi fine-tuned or class-specific core ncRNA models.                                                                                   |
| `-sno_cm_cutoff`        |                  No | `30`             | `cmsearch` score cutoff for fungi fine-tuned snoRNA models.                                                                                                         |
| `-srp_cm_cutoff`        |                  No | `30`             | `cmsearch` score cutoff for fungal SRP covariation models.                                                                                                          |
| `-include_sn_variants`  |                  No | off              | Also predict snRNA variants. By default, these variants are excluded.                                                                                               |
| `-reanalyse`            |                  No | `NA`             | Reanalyse a previous output directory, for example to add RNAmmer and tRNAscan-SE predictions if they were not run previously.                                      |
| `-merge_separately`     | No, but recommended | off              | When overlaps occur, prefer RNAmmer rRNA predictions and tRNAscan-SE tRNA predictions where available; otherwise use the best-scoring covariation-model prediction. |
| `-rfam_path`            |                  No | `Rfam/Rfam_f.cm` | Path to an alternative Rfam covariance-model file. Useful when using an updated Rfam release or custom model set.                                                   |
| `-h`                    |                  No | off              | Display `help.txt`, the help manual.                                                                                                                                |

## Choosing `-fungal_class`

Use `-fungal_class` when the fungal class is known. The script supports class-specific model naming for classes such as:

* `Eurotiomycetes`
* `Sordariomycetes`
* `Saccharomycetes`
* `Dothideomycetes`
* `Agaricomycetes`
* `Leotiomycetes`
* `Tremellomycetes`
* `Microsporidia`

## Output files

For an output directory named `Fun_ncRNA_output`, important outputs include:

```text
Fun_ncRNA_output/
├── input_genome.ncrna.gff
├── input_genome.gff
├── input_genome.combined.gff
├── input_genome.ncrna.compatible.gff
├── input_genome.combined.compatible.gff
├── summary.txt
├── log.txt
├── log.err.txt
├── overlapped_ref_genes.id.txt
├── overlapped_ref_genes.txt
└── NCRNA/
    ├── input_genome.all.nr.gff
    ├── input_genome.all.nr.assm.gff
    └── input_genome.all.nr.assm.gffcmp.gff
```

Main output files:

| File                                        | Description                                                                                                                              |
| ------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `input_genome.ncrna.gff`                    | Final predicted ncRNA annotation in GFF-like gene/exon format.                                                                           |
| `input_genome.ncrna.bed`                    | Final predicted ncRNA annotation in BED format. Useful for generating an IGV report.                                                     |
| `input_genome.gff`                          | Copy of the input GFF annotation, if `-gff` was supplied.                                                                                |
| `input_genome.combined.gff`                 | Combined input annotation and predicted ncRNA annotation, sorted by contig and genomic start position. Produced when `-gff` is supplied. |
| `input_genome.ncrna.compatible.gff`         | Final predicted ncRNA annotation in GFF-like gene/mRNA/exon/CDS format. Useful when downstream software requires mRNA and CDS features.  |
| `input_genome.combined.compatible.gff`      | Combined input annotation and predicted ncRNA annotation with ncRNAs in gene/mRNA/exon/CDS format. Produced when `-gff` is supplied.     |
| `summary.txt`                               | Summary of predicted ncRNA features.                                                                                                     |
| `log.txt`                                   | Standard pipeline log.                                                                                                                   |
| `log.error.txt`                             | Error log.                                                                                                                               |
| `overlapped_ref_genes.id.txt`               | Gene IDs in the reference GFF annotation that are overlapped by predicted ncRNAs. Produced when `-gffcmp` is supplied.                   |
| `overlapped_ref_genes.txt`                  | As above, key GFFCompare outputs showing genes in the reference GFF annotation that are overlapped by predicted ncRNAs.                  |
| `NCRNA/input_genome.all.nr.gff`             | Non-redundant ncRNA predictions after merging/filtering.                                                                                 |
| `NCRNA/input_genome.all.nr.assm.gff`        | Non-redundant predictions with additional ID, assembly, and ncRNA class fields.                                                          |
| `NCRNA/input_genome.all.nr.assm.gffcmp.gff` | Predictions annotated with GFFCompare-style comparison codes, when available.                                                            |

## Reanalysing previous output

If you previously ran the pipeline without RNAmmer/tRNAscan-SE and later install these optional tools, you can reanalyse the existing output:

```bash
perl Fungi_ncRNA.pl \
  -reanalyse previous_output_directory \
  -out reanalysed_output \
  -run_rnammer_trnascan
```

The previous directory should contain an `NCRNA/` subdirectory generated by an earlier run.

## Examples

### Test example 

```bash
  perl Fungi_ncRNA.pl \
  -fna data/CopciAB_new_jgi_20220113.fasta \
  -gff data/CopciAB_new_jgi_20220113.gff \
  -fungal_class Agaricomycetes \
  --gffcmp \
  -out CociAB_ncRNA_output
```
The result should be the similar to 'example_output'

### Core ncRNA prediction only

```bash
perl Fungi_ncRNA.pl \
  -fna genome.fna \
  -out Fun_ncRNA_core
```

### ncRNA prediction with existing annotation

```bash
perl Fungi_ncRNA.pl \
  -fna genome.fna \
  -gff genome.gff \
  -out Fun_ncRNA_with_annotation
```

### Full workflow with optional rRNA/tRNA predictors and GFFCompare

```bash
perl Fungi_ncRNA.pl \
  -fna genome.fna \
  -gff genome.gff \
  -fungal_class Sordariomycetes \
  -gcf GCF_000000000.1 \
  -out Fun_ncRNA_full \
  -run_rnammer_trnascan \
  -gffcmp \
  -merge_separately
```
## Specific applications

### Preparing a GFF file for RNASeq analysis

If the focus of the RNASeq analysis is on protein coding genes, therefore ncRNA and spurious mRNA that overlap ncRNA. These overlapping and potentially spurious mRNA are provided a list of IDs of genes in the reference GFF annotation that overlap with predicted ncRNAs. Users can use the `grep` command to remove them from the reference GFF annotation.

```bash
grep -v -f overlapped_ref_genes.id.txt input_genome.gff | grep -v '^--' > input_genome.no_ncrna.gff
```

Then, to analyse or quantify ncRNAs, users can use one of the ncRNA-only GFF annotations. Because RNASeq pipelines often requires the gene/mRNA/exon/CDS formatting for all genes even if they are non-coding, the RNASeq software-compatible version, `input_genome.ncrna.compatible.gff`should be used.

### Mitochondrial genomes in the assembly

The initial assembly of a fungal genome will usually include contigs that are the mitochondrial genome (perhaps as a circular contig ~11-340kb). This should removed before eukaryotic genome analysis (e.g. using GetOrganelle) and annotated separetely (e.g. with Mitos or mfannot). If the mitogenome is not removed and the standard Rfam CM models are used then some ncRNAs from mitochondria will be detected e.g.`Intron_gpI` and `Intron_gpII` and prokaryotic rRNA and tRNA.

## Troubleshooting

### `cmsearch: command not found`

Install Infernal and make sure the conda environment is active:

```bash
conda activate fun_ncrna
which cmsearch
```

### `bedtools`, `gffread`, or `gffcompare` not found

Install the missing package through Bioconda:

```bash
conda install -n fun_ncrna -c conda-forge -c bioconda bedtools gffread gffcompare
```

### RNAmmer does not run

RNAmmer requires separate installation and may require manual configuration of installation paths and HMMER/HMMER2 paths. Check:

```bash
which rnammer
which hmmsearch
rnammer -h
```

If RNAmmer reports an HMMER-related error, verify that the HMMER version expected by your RNAmmer installation is installed and correctly configured.

### tRNAscan-SE does not run

Check that tRNAscan-SE is installed and available:

```bash
which tRNAscan-SE
tRNAscan-SE -h
```

If tRNAscan-SE reports Infernal-related errors, confirm that Infernal is installed in the same environment:

```bash
which cmsearch
cmsearch -h
```

### Existing output directory is deleted

The script removes and recreates the output directory if it already exists. To avoid losing previous results, always choose a new `-out` directory or back up the existing directory before rerunning.

### Reference annotation file cannot be processed, if provided:

A reference annotation file is required for merging ncRNA predictions with the referencing annotation, to excluding overlapping coding genes for RNASeq analyses. A reference annotation file need to be in GFF format. If it is in GFF format and it still cannot be processed, chances are that features are not presented following the NCBI-style gene/mRNA/exon/CDS formatting.

A common example is a JGI-formatted annotation file, which is in GFF format but features are recorded following a unique formatting. A Perl script `ReformatJGIGFF.pl` is provided, for converting JGI-formatted GFF annotation to the NCBI-style gene/mRNA/exon/CDS formatting.

```bash
perl ReformatJGIGFF.pl input_jgi.gff output_ncbi.gff
```

## Citation

If you use this pipeline, please cite the relevant tools used by your analysis, including Infernal, bedtools, gffread/gffcompare, tRNAscan-SE, RNAmmer, and any covariance-model resources used by the pipeline.

## License

Add license information here.

## Contact

For questions, issues, or feature requests, please use the GitHub Issues page for this repository.
