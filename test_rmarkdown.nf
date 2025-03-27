nextflow.enable.dsl=2

include { RENDER_RNASEQ_REPORT } from './modules/local/rmarkdown_report'


workflow {
  RENDER_RNASEQ_REPORT(
    RNAseq_report       = file("modules/local/rmarkdown_report/RNAseq_report.Rmd"),
    HRK_funcs           = file("modules/local/rmarkdown_report/HRK_funcs.R"),
    analysis            = file("modules/local/rmarkdown_report/analysis.R"),
    params_file         = file(params.params_file),
    rmarkdown_container = params.rmarkdown_container,
    sample_id           = params.sample_id,
    rmarkdown_outdir    = params.rmarkdown_outdir
  )
}

