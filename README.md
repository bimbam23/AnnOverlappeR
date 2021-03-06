AnnOverlappeR
================

### What is AnnOverlappeR?

This is the AnnOverlappeR a tool to find the best overlaps between NCBI EntrezGeneIDs and EnsemblGeneIDs from the Ensembl database using GFF and GTF files from each databases. The tool it self is performing an overlap based on the genomic positions.

### Install

1. Clone the [GitHub repo](https://github.com/bimbam23/AnnOverlappeR), e.g. with `git clone https://github.com/bimbam23/AnnOverlappeR.git`

### Usage

    Rscript AnnOverlappeR.R
    
    Usage: Rscript AnnOverlappeR.R [options]
    
      --species character
        Species in scientific name (e.g.: homo_sapiens), string, mandatory option
    
      --output_path character
        output_path, string, optional option
    
      --numcores integer
        Number of cpus , integer, optional (default 4)
    
      --download 
        --no-download for no downloads, boolean, mandatory option (default TRUE)
    
      --latest_ncbi 
        --no-latest_NCBI for no downloads, boolean, mandatory option (default
        FALSE)
    
      --latest_name character
        latest name found for ncbi latest assembly, string, optional option
    
      --verbose 
        print messages
    
      --help 
        Print help message and exit.
    
      --version 
        Print version information and exit.
    
    Please contact jochen.bick@usys.ethz.ch for comments

### Example

    Rscript AnnOverlappeR.R -s homo_sapiens --numcores 12 -o ./2018-06-13/ --download --latest_ncbi >hsa_new.out &
    
### Contact
jochen.bick@usys.ethz.ch

### Tool dependencies:

``` {r}
library(GetoptLong)
library(GenomicRanges)
library(rtracklayer)
library(XML)
library(pbmcapply)
```

##### more dependencies
perl, gzip, wget

#### Galaxy app
The Galaxy platform is available at https://galaxyproject.org/admin/get-galaxy/.
All files have to be copied to the listed folders.
The file content of all xml files in the config folder have to be added to the corresponding xml files. 

