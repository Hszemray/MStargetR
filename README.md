MStargetR
<img src='man/figures/logo.png' style='float: right; height:50px;' />
================

[![license](https://img.shields.io/github/license/vimalnathnambiar/exfilms.svg?color=%23A31F34)](LICENSE)
[![semantic-release](https://img.shields.io/badge/semantic--release-angular-E10079.svg?logo=semantic-release&logoColor=%23E10079)](https://github.com/semantic-release/semantic-release)
[![doi-zenodo](https://img.shields.io/badge/zenodo-10.5281/zenodo.17889838-blue.svg?logo=doi&logoColor=blue)](https://doi.org/10.5281/zenodo.17889838)

## Overview

MStargetR provides a pipeline for processing and quality control of
targeted mass spectrometry vendor files to concentration outputs ready
for further analysis. The package allows for cross platform usage
through Docker.

Currently Apple silicon based products(M1, M2, M3, M4) are not supported
due to system architecture limitations. User’s with these devices will
be unable to use msConvertR or PeakForgeR with method set to skyline.

The toolkit consists of three main functions:

**`msConvertR()`** Processes mass spectrometry vendor files into open
source mzml files.

**`PeakForgeR()`** Performs retention time and peak boundary
optimisation, peak picking and integration.

**`qcCheckR()`** Performs quality control, batch/signal drift
correction, and generates analytical reports.

## Installation

**`Docker Desktop`** Users are required to install Docker in order to
run MStargetR::msConvertR and MStargetR::PeakForgeR.  
Please visit <https://www.docker.com/get-started/> for installation
process.

Note: Please ensure Docker Application is running in the background when
using MStargetR if you are using a Mac or Windows OS.

**`MStargetR`** Installation is effortless with the below helper
function. Once the package is attached you will be asked for a username
and a default default workflow will open for your convenience.

``` r
source("https://raw.githubusercontent.com/Hszemray/MStargetR/master/R/install.R")
install_MStargetR()
rm(install_MStargetR)
library(MStargetR)
```

## Project Setup

Create this project folder structure before execution:

    User/specified/directory
    └── raw_data #Store vendor files here

For example: “C:/User/Desktop/Project1/raw_data”

- Project1 is the name of the project

- raw_data is a subfolder to keep vendor files for the project

## Main User Functions

### 1. msConvertR() Workflow:

1.  Converts vendor files → .mzML via ProteoWizard msConvert.
2.  Organises mzml outputs into file structure for each vendor file

### 2. PeakForgeR() Workflow:

Process iterates by plate to:

1.  Optimise retention times based on QC samples
2.  Optimise peak boundaries
3.  Generate PeakForgeR input files
4.  Execute PeakForgeR peak picking and integration
5.  Organise output into plate-specific folders

### 3. qcCheckR() Workflow:

1.  Imputation of missing values
2.  Batch / signal drift correction
3.  Sample and feature filtering  
4.  Quality metric calculation (RSD, missing values)  
5.  Generation of:
    - Interactive PCA and control chart plots  
    - Quality summary tables  
    - HTML/Excel reports

## Output Structure

Post-execution directory of msConvertR, PeakForgeR, and qcCheckR:

    User Project Directory
    ├── Archive #Archive of wiff files and MStargetR logs by plate
    ├── ALL
    │   ├── data
    │   │   ├── batch_correction # batch  and singnal drift correction outputs
    │   │   ├── rda # saved rda from qcCheckR 
    │   └── html_report #Summary of project, PCA, and control charts
    │   └── xlsx_report # Data for further user analysis
    ├── PLATE_NAME # This structure will be created per plate
        ├── data
            ├── mzml #converted mzml files
            ├── rda  #saved rda from PeakForgeR processing
            ├── raw_data #original vendor files
            └── PeakForgeR #input and output files of PeakForgeR processing

## User Guide

For a more detailed guide please visit the following link
[here](https://Hszemray.github.io/MStargetR/MStargetR-vignette.html)

## Technical Notes

1.  **Dependencies**:
    - Docker Desktop (<https://www.docker.com/get-started/>)
2.  **Memory Management**:
    - Process large datasets plate-by-plate
    - Offload processed plate data before processing next plate
3.  **Troubleshooting**:
    - Run R as systems administrator
    - Docker Application is installed and up to date
    - Ensure your R, Rstudio and package versions are all up to date
    - Check vendor files/mzml conversions for corruption
    - Ensure sufficient disk space
    - Minimum of 8GB RAM required

## Support

For technical assistance, please open an issue at:  
<https://github.com/Hszemray/MStargetR>

## License

Please refer to license information [here](LICENSE)

## Citation

If you use this software in your work, please cite it using the
following metadata:

Szemray, H., Nambiar, V., Lawler, N., Lodge, S., Wist, J.,& Whiley, L.
(2025). MStargetR (Version 1.0.0) \[Computer software\].
Zenodo.(zenodo link)

<!-- Links -->
