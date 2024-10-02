# namlab-mapper
Little workflow which can download and map multiple RNA sequencing files from the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) as well as any local FASTQ files to a common reference using kallisto. 
Because it is written in [Nextflow](https://www.nextflow.io/), it can automatically parallelize steps across CPUs or nodes, if you are running it on a cluster (see [this page](https://www.nextflow.io/docs/latest/executor.html) for more details).
It is also built to be economical with disk space by removing large intermediary files when they are no longer needed.
The output is a combined table containing abundance quantifications as well as [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports for each of sequence files.

## Prerequisites
rnaseq-mapper will try to load the following [modules](http://modules.sourceforge.net/): `sratoolkit`, `kallisto`, `R`, `fastqc`.
If your system doesn't use modules, make sure the execs are available in your PATH.

## Usage
1. Set up nextflow (if not installed already):
```
curl -s https://get.nextflow.io | bash
```
2. Create a file called `nextflow.config` (exactly this name) by using the `example_nextflow.config` from this directory as a template and adapting it to your use case.
3. Create an input file with the sequences you want to map in the format of `example_input.csv` and make sure it is referred to in your config file.
4. If desired, place any FASTQ files in the directories referenced in your `nextflow.config` (if you don't have any, make sure the folders still exist and just leave them empty).
5. Run the pipeline:
```
./nextflow run NAMlab/rnaseq-mapper
```

### Singularity Container
If you prefer, you can also make use of the [Singularity](https://sylabs.io/) container that packages all the required software.
This requires Singularity or [Apptainer](https://apptainer.org/) to be installed in your system.
You can then simply execute the pipeline (step 5 above, the other steps stay the same) via `./nextflow run NAMlab/rnaseq-mapper -with-singularity library://merlin/default/rnaseq-mapper:latest` or `./nextflow run NAMlab/rnaseq-mapper -with-apptainer library://merlin/default/rnaseq-mapper:latest` respectively.

## Output
You will get out a TSV file with the combined kallisto outputs for all your sequence files like this one (by default in the `work/out` folder):

| target_id | length | SRR1805735_eff_length | SRR1805735_est_counts | SRR1805735_tpm | SRR1805737_eff_length | SRR1805737_est_counts | SRR1805737_tpm | SRR6512869_eff_length | SRR6512869_est_counts | SRR6512869_tpm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Solyc00g005280.1.1 | 411  | 252.224 | 0 | 0 | 241.253 | 0 | 0 | 212     | 0 | 0 |
| Solyc00g005285.2.1 | 216  | 68.6464 | 0 | 0 | 63.7937 | 0 | 0 | 31.5146 | 0 | 0 | 
| Solyc00g006483.2.1 | 390  | 231.296 | 0 | 0 | 220.691 | 0 | 0 | 191     | 0 | 0 | 
| Solyc00g006487.2.1 | 276  | 120.525 | 0 | 0 | 114.108 | 0 | 0 | 77.4659 | 2 | 22.2662 |
| Solyc00g006560.2.1 | 1317 | 1158    | 0 | 0 | 1145.76 | 0 | 0 | 1118    | 0 | 0 | 
| Solyc00g006890.2.1 | 300  | 143.123 | 0 | 0 | 135.795 | 0 | 0 | 101.044 | 0 | 0 | 
| Solyc00g006900.2.1 | 576  | 416.999 | 0 | 0 | 404.931 | 0 | 0 | 377     | 0 | 0 | 
| Solyc00g007225.2.1 | 1275 | 1116    | 0 | 0 | 1103.76 | 0 | 0 | 1076    | 0 | 0 | 
| Solyc00g007330.1.1 | 516  | 356.999 | 0 | 0 | 345.082 | 0 | 0 | 317     | 0 | 0 | 

You will also get [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports for each of sequence files in the same folder.
