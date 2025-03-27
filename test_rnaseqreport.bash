module load nextflow
module load apptainer

nextflow run test_rmarkdown.nf \
  -profile standard \
  --sample_id CX_lung \
  --params_file /blue/cancercenter-dept/hkates/RNAseq_reporting/CX_lung_report_params.txt \
  --rmarkdown_container /blue/cancercenter-dept/hkates/RNAseq_reporting/rmarkdown_report/rmarkdown_report.sif \
  --rmarkdown_outdir results/reports
