sequences_ch = Channel.from(file(params.input_file).text).splitCsv(header: true)

process indexCDS {
  input:
    path "cds.fa" from file(params.ref_cds)
  output:
    path "cds_index" into cds_ch

  script:
  """
   module load kallisto
   kallisto index -i cds_index cds.fa
  """
}

process GetnMapSequence {
  scratch params.scratch_dir

  input:
    val sequence from sequences_ch
    path "cds_index" from cds_ch
  output:
    path "${sequence.sra_run_id}.tsv" into abundances_ch

  script:
    if (sequence.layout == "paired")
      """
       module load kallisto sratoolkit
       fasterq-dump --split-files ${sequence.sra_run_id}
       kallisto quant -i cds_index -o ./ ${sequence.sra_run_id}_1.fastq ${sequence.sra_run_id}_2.fastq
       mv abundance.tsv ${sequence.sra_run_id}.tsv
      """
    else if (sequence.layout == "single")
      """
       module load kallisto sratoolkit
       fasterq-dump ${sequence.sra_run_id}
       kallisto quant -i cds_index -o ./ --single -l ${sequence.fragment_length_average} -s ${sequence.fragment_length_sd} ${sequence.sra_run_id}.fastq
       mv abundance.tsv ${sequence.sra_run_id}.tsv
      """
    else
      error "Invalid sequence layout: ${sequence.layout}"
}

process combineAll {
  publishDir "work/out/", mode: 'copy'

  input:
    path "*" from abundances_ch.collect()
  output:
    path "combined_abundance.tsv" into combined_ch

  script:
  """
  #!/usr/bin/env Rscript

  abundance.files = list.files(pattern=("*.tsv"))
  accession.ids = sub(".tsv", "", abundance.files)
  combined_df = read.csv(abundance.files[1], sep="\t")
  names(combined_df) = c("target_id", "length", 
                         paste0(accession.ids[1], "_eff_length"), 
                         paste0(accession.ids[1], "_est_counts"), 
                         paste0(accession.ids[1], "_tpm"))
  for(i in 2:length(abundance.files)) {
    additional_df = read.csv(abundance.files[i], sep="\t")
    names(additional_df) = c("target_id", "length", 
                             paste0(accession.ids[i], "_eff_length"), 
                             paste0(accession.ids[i], "_est_counts"), 
                             paste0(accession.ids[i], "_tpm"))
    combined_df = merge(combined_df, additional_df, all=T)
  }
  write.table(combined_df, "combined_abundance.tsv", row.names=F, sep="\t", quote=F)
  """
}
