# namlab-mapper
Convenience script in [Nextflow](https://www.nextflow.io/) to download and map multiple RNA sequences from SRA to the same reference at once using kallisto and combine all the abundance quantifications into one table.

## Prerequisites
rnaseq-mapper will try to load the following [modules](http://modules.sourceforge.net/): `sratoolkit`, `kallisto`, `R`, `fastqc`.
If your system doesn't use modules, make sure the execs are available in your PATH.
I'm currently working on a [Singularity](https://sylabs.io/) container that packages all the required software (see `singularity_container.def`).

## Usage
1. Set up nextflow (if not installed already):
```
curl -s https://get.nextflow.io | bash
```
2. Create a file `nextflow.config` by using the `example_nextflow.config` from this directory as a template and adapting it to your use case.
3. Create an input file with the sequences you want to map in the format of `example_input.csv` and make sure it is referred to in your config file.
4. Run the pipeline:
```
./nextflow run NAMlab/rnaseq-mapper
```

## Output
You will get out a TSV file with the combined kallisto outputs for all your sequences like this one:

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

You will also get [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports for each of the downloaded sequences.
