# ASCAT Shiny App

An interactive Shiny application for ASCAT (Allele-Specific Copy number Analysis of Tumours) analysis with real-time visualization of copy number profiles.

## Features

- **Interactive Analysis**: Upload ASCAT .Rdata files and run analysis with custom parameters
- **Manual Control**: Option to manually set rho (tumour purity) and psi (tumour ploidy)
- **Real-time Visualization**: Four synchronized plots:
  - Segmented data
  - Raw/non-rounded profile
  - ASCAT copy-number profile
  - Sunrise/goodness-of-fit plot
- **Modern UI**: Clean, responsive interface built with bslib

## Requirements

### R Version
- R >= 4.0.0

### R Packages
```r
# CRAN packages
install.packages(c("shiny", "bslib", "shinycssloaders", "ggplot2", "dplyr", "tidyr", "png"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ASCAT", "GenomicRanges", "IRanges"))
```

## Installation

1. Clone this repository:
```bash
git clone https://github.com/evitjohn/ASCATShinyApp.git
cd ASCATShinyApp
```

2. Ensure you have the required R packages installed (see Requirements above)

3. Make sure you have the `ascat_runAscat_ForShiny.R` file in the same directory as the main app file

## Usage

Run the app from R or RStudio:

```r
shiny::runApp("ascat_shiny_app_V9.R")
```

Or from the command line:

```bash
R -e "shiny::runApp('ascat_shiny_app_V9.R')"
```

### Input

Upload an ASCAT .Rdata file containing an ASCAT object (typically named `ascat.bc`, `ASCATobj`, or `ascat_bc`).

### Parameters

- **Gamma**: Segmentation penalty parameter (default: 1)
- **Manual rho/psi**: Optional manual override of tumour purity and ploidy
  - **Rho**: Tumour purity percentage (1-100%)
  - **Psi**: Tumour ploidy (1.0-10.0)

### Output

The app displays four plots:
1. **Segmented Data**: Shows the segmented genomic data
2. **Raw Profile**: Non-rounded copy number profile
3. **ASCAT Profile**: Final copy number profile with allele-specific information
4. **Sunrise Plot**: Goodness-of-fit visualization

Results including purity, ploidy, and goodness-of-fit are displayed in the status box.

## File Structure

```
.
├── ascat_shiny_app_V9.R          # Main Shiny app
├── ascat_runAscat_ForShiny.R     # Custom ASCAT function (required)
└── README.md                      # This file
```

## Notes

- Maximum upload file size: 100 MB
- Temporary files are created in the system temp directory during analysis
- If automatic rho/psi optimization fails, try enabling manual selection

## Acknowledgments

Built using the [ASCAT](https://github.com/VanLoo-lab/ascat) package for copy number analysis.
