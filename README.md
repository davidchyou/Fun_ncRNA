# Fun_ncRNA

Predict non-coding RNAs in fungal genomes using covariation models fine-tuned for fungi.

`Fun_ncRNA` is a command-line pipeline for fungal ncRNA annotation. It searches fungal genome assemblies with curated/fine-tuned covariation models, integrates optional predictions from RNAmmer and tRNAscan-SE, removes overlapping/redundant predictions, and produces GFF files and a prediction summary.

The main entry point is:

```bash
perl Fungi_ncRNA.pl -fna genome.fna [options]
```

## Main features

* Predicts fungal non-coding RNAs from genome FASTA files.
* Uses fungal class-specific covariation models when requested.
* Supports optional comparison against an existing genome annotation in GFF format.
* Supports optional integration of RNAmmer rRNA predictions and tRNAscan-SE tRNA predictions.
* Produces GFF-formatted ncRNA annotations and a summary table.
* Can reanalyse a previous output directory to add optional RNAmmer/tRNAscan-SE results.

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

## Installation

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

Check that the required programs are available:

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

## Basic usage

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
  -gcf GCF_000000000.1 \
  -out Fun_ncRNA_output
```

To also run RNAmmer and tRNAscan-SE:

```bash
perl Fungi_ncRNA.pl \
  -fna genome.fna \
  -gff genome.gff \
  -fungal_class Sordariomycetes \
  -gcf GCF_000000000.1 \
  -out Fun_ncRNA_output \
  -run_rnammer_trnascan
```

To run GFFCompare against an existing annotation:

```bash
perl Fungi_ncRNA.pl \
  -fna genome.fna \
  -gff genome.gff \
  -out Fun_ncRNA_output \
  -gffcmp
```

## Command-line options

| Option                  |       Required? | Default          | Description                                                                                                                                                         |
| ----------------------- | --------------: | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-fna`                  |             Yes | `NA`             | Input genome assembly in FASTA format.                                                                                                                              |
| `-gff`                  | No, recommended | `NA`             | Existing genome annotation in GFF format. Used for combined output and optional GFFCompare analysis.                                                                |
| `-fungal_class`         |              No | `Fungi`          | Fungal taxonomic class. If supplied, class-specific covariation models are used where available. Otherwise, fungi-wide models are used.                             |
| `-run_rnammer_trnascan` |              No | off              | Run optional RNAmmer and tRNAscan-SE predictions.                                                                                                                   |
| `-gffcmp`               |              No | off              | Run GFFCompare. Skipped if no GFF annotation is provided.                                                                                                           |
| `-out`                  |              No | `Fungi_ncRNA`    | Output directory. Existing directories with the same name will be removed and recreated.                                                                            |
| `-gcf`                  |              No | `NA`             | Assembly ID for multi-contig genome files.                                                                                                                          |
| `-core_cm_cutoff`       |              No | `60`             | `cmsearch` score cutoff for fungi fine-tuned or class-specific core ncRNA models.                                                                                   |
| `-sno_cm_cutoff`        |              No | `30`             | `cmsearch` score cutoff for fungi fine-tuned snoRNA models.                                                                                                         |
| `-srp_cm_cutoff`        |              No | `30`             | `cmsearch` score cutoff for fungal SRP covariation models.                                                                                                          |
| `-include_sn_variants`  |              No | off              | Also predict snRNA variants. By default, these variants are excluded.                                                                                               |
| `-reanalyse`            |              No | `NA`             | Reanalyse a previous output directory, for example to add RNAmmer and tRNAscan-SE predictions if they were not run previously.                                      |
| `-merge_separately`     |              No | off              | When overlaps occur, prefer RNAmmer rRNA predictions and tRNAscan-SE tRNA predictions where available; otherwise use the best-scoring covariation-model prediction. |
| `-rfam_path`            |              No | `Rfam/Rfam_f.cm` | Path to an alternative Rfam covariance-model file. Useful when using an updated Rfam release or custom model set.                                                   |

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

If the class is unknown or not represented by a class-specific model set, omit this option or use the default:

```bash
-fungal_class Fungi
```

## Output files

For an output directory named `Fun_ncRNA_output`, important outputs include:

```text
Fun_ncRNA_output/
├── input_genome.ncrna.gff
├── input_genome.gff
├── input_genome.combined.gff
├── summary.txt
├── log.txt
├── log.err.txt
└── NCRNA/
    ├── input_genome.all.nr.gff
    ├── input_genome.all.nr.assm.gff
    └── input_genome.all.nr.assm.gffcmp.gff
```

Main output files:

| File                                        | Description                                                                                                                              |
| ------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `input_genome.ncrna.gff`                    | Final predicted ncRNA annotation in GFF-like gene/mRNA/exon format.                                                                      |
| `input_genome.gff`                          | Copy of the input GFF annotation, if `-gff` was supplied.                                                                                |
| `input_genome.combined.gff`                 | Combined input annotation and predicted ncRNA annotation, sorted by contig and genomic start position. Produced when `-gff` is supplied. |
| `summary.txt`                               | Summary of predicted ncRNA features.                                                                                                     |
| `log.txt`                                   | Standard pipeline log.                                                                                                                   |
| `log.err.txt`                               | Error log.                                                                                                                               |
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

## Suggested workflows

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

## Citation

If you use this pipeline, please cite the relevant tools used by your analysis, including Infernal, bedtools, gffread/gffcompare, tRNAscan-SE, RNAmmer, and any covariance-model resources used by the pipeline.

## License

Add license information here.

## Contact

For questions, issues, or feature requests, please use the GitHub Issues page for this repository.

