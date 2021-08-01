sequences_ch = Channel.from(file(params.input_file).text).splitCsv(header: true)

process indexReference {
  input:
    path "ref.fa" from file(params.ref_fasta)
  output:
    path "ref_index" into reference_ch

  script:
  """
   module load kallisto
   kallisto index -i ref_index ref.fa
  """
}

process GetnMapSequence {
  publishDir "work/out/fastqc-reports", mode: 'move', pattern: '*_fastqc.zip'
  scratch params.scratch_dir
  errorStrategy 'retry'
  maxRetries 4

  input:
    val sequence from sequences_ch
    path "ref_index" from reference_ch
  output:
    path "${sequence.sra_run_id}.tsv" into abundances_ch
    path "${sequence.sra_run_id}*_fastqc.zip" into fastqc_ch

  script:
    if (sequence.layout == "paired")
      """
       module load kallisto sratoolkit fastqc
       fasterq-dump --split-files ${sequence.sra_run_id}
       fastqc -t 2 ${sequence.sra_run_id}_1.fastq ${sequence.sra_run_id}_2.fastq
       kallisto quant -i ref_index -o ./ ${sequence.sra_run_id}_1.fastq ${sequence.sra_run_id}_2.fastq
       mv abundance.tsv ${sequence.sra_run_id}.tsv
      """
    else if (sequence.layout == "single")
      """
       module load kallisto sratoolkit fastqc
       fasterq-dump ${sequence.sra_run_id}
       fastqc ${sequence.sra_run_id}.fastq
       kallisto quant -i ref_index -o ./ --single -l ${sequence.fragment_length_average} -s ${sequence.fragment_length_sd} ${sequence.sra_run_id}.fastq
       mv abundance.tsv ${sequence.sra_run_id}.tsv
      """
    else
      error "Invalid sequence layout: ${sequence.layout}"
}

process combineAll {
  publishDir "work/out/", mode: 'move'

  input:
    path "*" from abundances_ch.collect()
  output:
    path "combined_abundance.tsv" into combined_ch
    
  module 'R'

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
