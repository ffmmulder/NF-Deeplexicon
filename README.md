# NF-Deeplexicon

Nextflow implementation of the bash wrapper (and as of now obsolete) pipeline script found here: https://github.com/UMCUGenetics/DeeplexiconDemuxAndSplit

Can be used to run demux and split on nanopore fast5 and fastq data using the Deeplexicon pipeline

See the Deeplexicon github for extended information:
https://github.com/Psy-Fer/deeplexicon

## Installing and setup

The following tools and steps are required for this pipeline:

1. [Nextflow v21+](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2. [Singularity](https://sylabs.io/guides/3.5/admin-guide/)
3. [Pull/Clone NF-PrepareGenome](#pull-or-clone)

## Pull or Clone (clone NextflowModules seperately)
```
git clone git@github.com:ffmmulder/NF-Deeplexicon.git
git -C ./NF-Deeplexicon clone git@github.com:UMCUGenetics/NextflowModules.git
```

## Summary

The workflow basically has 3 main steps.

1. Predict the barcodes based on the fast5 files in the nanopre run folder, this is done for each fast5 file in parallel for fastest processing
2. Concatenate the resulting files from 1) as well as all the seperate fastq.gz files in nanopore run folder for final processing
3. Split the fastq file based on the barcodes found and recorded in the tsv files

## Execution

The script is ran in the following way, the supplied parameters are obligatory:

```
# Minimal example
nextflow main.nf \
    --out_dir /path/to/output/ \
    --fast5_path /path/to/fast5_pass/ \
    -profile slurm \
    -c /path/to/NF_Deeplexicon/nextflow.config \
    -resume
```

At least the path to the passed fast5 files as well as the output folder must be specified. Optionally a folder to the corresponding passed fastq files can be supplied, if not it will look in the parent folder of the fast5_pass folder and look for a fastq_pass folder in there.

## Optional parameters

--fastq_path  Path containign the passed fastq files (if not found in the default location mentioned above)
--tsv_path  If barcode detection was already done and tsv files were already generated, this can be used to specify path where tsv files can be found for further concatenation and processing



