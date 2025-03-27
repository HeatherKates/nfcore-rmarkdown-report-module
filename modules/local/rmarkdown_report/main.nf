process RENDER_RNASEQ_REPORT {
  tag "${params.sample_id}"                  // <-- must use params.
  publishDir "${params.rmarkdown_outdir}", mode: 'copy' // <-- must use params.

  input:
  path "RNAseq_report.Rmd"
  path "HRK_funcs.R"
  path "analysis.R"
  path "params_file"
  val rmarkdown_container
  val sample_id
  val rmarkdown_outdir

  output:
  path "*.html"

  script:
  """
  apptainer exec ${rmarkdown_container} Rscript -e "rmarkdown::render(
    'RNAseq_report.Rmd',
    output_file = '${sample_id}.Report.html',
    params = list(params_file = '${params_file}')
  )"
  """
}

