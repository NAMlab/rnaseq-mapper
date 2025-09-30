process indexReference {
  module 'kallisto'
  input:
    path "ref.fa"
  output:
    path "ref_index" 

  script:
  """
   kallisto index -i ref_index ref.fa
  """
}

process SRA {
  publishDir "work/out/fastqc-reports", mode: 'move', pattern: '*_fastqc.zip'
  scratch params.scratch_dir
  errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }
  maxRetries 4
  maxForks params.max_sra_forks
  module 'kallisto'
  module 'sratoolkit'
  module 'fastqc'

  input:
    val sequence 
    path "ref_index" 
  output:
    path "${sequence.sra_run_id}.tsv" 
    path "${sequence.sra_run_id}*_fastqc.zip"

  script:
    if (sequence.layout == "paired")
      """
       fasterq-dump --split-files ${sequence.sra_run_id}
       fastqc -t 2 ${sequence.sra_run_id}_1.fastq ${sequence.sra_run_id}_2.fastq
       kallisto quant -i ref_index -o ./ ${sequence.sra_run_id}_1.fastq ${sequence.sra_run_id}_2.fastq
       mv abundance.tsv ${sequence.sra_run_id}.tsv
      """
    else if (sequence.layout == "single")
      """
       fasterq-dump ${sequence.sra_run_id}
       fastqc ${sequence.sra_run_id}.fastq
       kallisto quant -i ref_index -o ./ --single -l ${sequence.fragment_length_average} -s ${sequence.fragment_length_sd} ${sequence.sra_run_id}.fastq
       mv abundance.tsv ${sequence.sra_run_id}.tsv
      """
    else
      error "Invalid sequence layout: ${sequence.layout}"
}

process LocalPaired {
  publishDir "work/out/fastqc-reports", mode: 'move', pattern: '*_fastqc.zip'
  scratch params.scratch_dir
  module 'kallisto'
  module 'sratoolkit'
  module 'fastqc'

  input:
    tuple val(pair_id), path(reads)
    path "ref_index"
  output:
    path "${pair_id}.tsv" 
    path "${pair_id}*_fastqc.zip" 

  script:
    """
     module load kallisto sratoolkit fastqc
     gunzip -c ${reads[0]} > ${pair_id}_1.fastq
     gunzip -c ${reads[1]} > ${pair_id}_2.fastq
     fastqc -t 2 ${pair_id}_1.fastq ${pair_id}_2.fastq
     kallisto quant -i ref_index -o ./ ${pair_id}_1.fastq ${pair_id}_2.fastq
     mv abundance.tsv ${pair_id}.tsv
    """
}

process LocalSingle {
  publishDir "work/out/fastqc-reports", mode: 'move', pattern: '*_fastqc.zip'
  scratch params.scratch_dir
  module 'kallisto'
  module 'sratoolkit'
  module 'fastqc'

  input:
    path fastq_file 
    path "ref_index"

  output:
    path "${fastq_file.simpleName}.tsv"
    path "${fastq_file.simpleName}*_fastqc.zip" 

  script:
    """
     gunzip -c ${fastq_file} > ${fastq_file.simpleName}.fastq
     fastqc ${fastq_file.simpleName}.fastq
     kallisto quant -i ref_index -o ./ --single -l ${params.local_reads_single_fraglength_mean} -s ${params.local_reads_single_fraglength_sd} ${fastq_file.simpleName}.fastq
     mv abundance.tsv ${fastq_file.simpleName}.tsv
    """
}

process combineAll {
  publishDir "work/out/", mode: 'move'
  module 'R'

  input:
    path "*" 
  output:
    path "combined_abundance.tsv"
    

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

workflow {
  sequences_local_single_ch = Channel.fromPath(params.local_reads_single)
  sequences_local_paired_ch = Channel.fromFilePairs(params.local_reads_paired)
  sequences_sra_ch = Channel.from(file(params.input_file).text).splitCsv(header: true)
  ref_index = indexReference(params.ref_fasta)
  (abundances_sra, fastq_sra) = SRA(sequences_sra_ch, ref_index)
  (abundances_single, fastq_single) = LocalSingle(sequences_local_single_ch, ref_index)
  (abundances_paired, fastq_paired) = LocalPaired(sequences_local_paired_ch, ref_index)
  combineAll(abundances_single.mix(abundances_paired, abundances_sra).collect())
}
