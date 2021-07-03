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
    """
     module load kallisto sratoolkit
     fasterq-dump --split-files ${sequence.sra_run_id}
     kallisto quant -i cds_index -o ./ ${sequence.sra_run_id}_1.fastq ${sequence.sra_run_id}_2.fastq
     mv abundance.tsv ${sequence.sra_run_id}.tsv
    """

}

process combineAll {
  publishDir "work/out/", mode: 'copy'

  input:
    path "*" from abundances_ch.collect()
  output:
    path "combined_abundance.tsv" into combined_ch

  // @TODO: this needs to be properly combined (e.g. in R or something)
  script:
  """
  cat * > combined_abundance.tsv
  """
}
