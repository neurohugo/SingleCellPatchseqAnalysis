# SingleCellPatchseqAnalysis 

<!-- badges: start -->
[![](https://img.shields.io/badge/devel%20version-1.1.1-blue.svg)](https://github.com/ludvigla/semla/releases) ![](https://img.shields.io/github/last-commit/neurohugo/RegenOrNoRegen.svg) [![License: UCSD](https://img.shields.io/badge/License-UCSD-yellow.svg)](https://opensource.org/license/ucsd/)


<!-- badges: end -->

<br>


`RegenOrNoRegen` is an R package that classifies Regenerating or Non-regenerating neurons from scRNAseq data.

Please Download the package from this link [website]((https://github.com/neurohugo/RegenOrNoRegen)). 
<br>

## Package Link: RegenOrNoRegen [website]((https://github.com/neurohugo/RegenOrNoRegen))

### Prerequisite: Seurat, Garnett and its requirements

You should install Garnett [website]((https://cole-trapnell-lab.github.io/garnett/docs/)). 

````
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")

BiocManager::install(c('monocle','DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))
install.packages("devtools")
devtools::install_github("cole-trapnell-lab/garnett")

````


## Installation

The dev version of the package can be installed through GitHub using;

````
devtools::install_github("neurohugo/RegenorNoRegen")
````

## Usage

After load the library you can use it this way.
````

ResultSeurat=RegenorNoRegen(OriginalSeurat)

````

<br>


