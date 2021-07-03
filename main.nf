key_vals_ch = Channel.from(file("input.csv").text).splitCsv(header: true)
scratch_dir = "/Users/psaroudakis/nextflow_test/tmp"

process printOut {
  scratch "$scratch_dir"

  input:
    val keyVal from key_vals_ch
  output:
    path "outfile" into outfiles_ch

  script:
    """
     echo ${keyVal.KEY.toUpperCase()} > outfile
     echo ${keyVal.VAL.toUpperCase()} >> outfile
     touch garbage_file
     echo '-----'
    """

}

process combineAll {
  publishDir "work/out/"

  input:
    path "outfile*" from outfiles_ch.collect()
  output:
    path "combined_out" into combined_ch

  script:
  """
  cat outfile* > combined_out
  """
}
