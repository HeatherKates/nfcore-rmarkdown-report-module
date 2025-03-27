# Local Nextflow Module: `rmarkdown_report`

This module renders a custom RMarkdown report (`RNAseq_report.Rmd`) using provided analysis code and parameter files. It's designed to run as a local module within a Nextflow pipeline, such as `nf-core/rnaseq`.

---

## 📁 Included Files

```
modules/local/rmarkdown_report/
├── main.nf              # Nextflow process for rendering report
├── RNAseq_report.Rmd    # Main RMarkdown report template
├── analysis.R           # Custom analysis code
├── HRK_funcs.R          # Custom utility functions
```

---

## 🚀 How to Run the Module Standalone

### 1. Create a test script (e.g. `test_rnaseqreport.bash`)

```bash
module load nextflow
module load apptainer

nextflow run test_rmarkdown.nf \
  -profile standard \
  --sample_id CX_lung \
  --params_file /full/path/to/CX_lung_report_params.txt \
  --rmarkdown_container /full/path/to/rmarkdown_report.sif \
  --rmarkdown_outdir results/reports
```

> 🔁 Replace paths with your actual file locations.

---

### 2. Create a minimal `test_rmarkdown.nf`

```nextflow
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
```

---

## 🐳 Container Notes

This module uses a custom Singularity container (`.sif`) that includes `downloadthis` and many CRAN/Bioconductor packages.

**Best practice:** Store the container in a stable location, e.g.:

```
/blue/cancercenter-dept/containers/rmarkdown_report.sif
```

Then reference it in your `test_rnaseqreport.bash` or main pipeline.

---

## 📝 Parameters

| Parameter              | Description                          |
|------------------------|--------------------------------------|
| `sample_id`            | Sample name for output HTML file     |
| `params_file`          | Plaintext file of parameters          |
| `rmarkdown_container`  | Path to `.sif` container file         |
| `rmarkdown_outdir`     | Output directory for rendered report |

---

## ✅ Usage in Pipeline

To integrate into a full Nextflow project (like `nf-core/rnaseq`), place this module in:

```
modules/local/rmarkdown_report/
```

Then call it from any workflow with the appropriate inputs.

---

