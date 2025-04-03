module load nextflow
module load apptainer

nextflow run test_rmarkdown.nf \
  -profile singularity \
  --sample_id CX_lung_test_contrast \
  --params_file CX_lung_report_params.txt  \
  --rmarkdown_container /blue/cancercenter-dept/PIPELINES/rmarkdown_report_module/rmarkdown_report.sif \
  --rmarkdown_outdir results/reports \
  --log CX_test.log
