# chromunity

![chromunity](inst/extdata/chromunity_logo.jpg)

### <font color=black> Discovery of communities in Pore-C concatemers.</font>

## <font color=black> Installations </font>

Install devtools from CRAN


```R
install.packages('devtools')
```

Install dependent packages and latest Bioconductor (if you haven't already)


```R
source('https://bioconductor.org/biocLite.R')
biocLite('GenomicRanges')
```

Install mskilab R dependencies (gUtils)


```R
devtools::install_github('mskilab/gUtils')
```

Install chromunity


```R
devtools::install_github('mskilab/chromunity')
```

## <font color=black> Tutorial </font>

This is a graph based community detection method for Pore-C concatemers. It takes a GRanges with each row as the Pore-C monomer. For each monomer, the corresponding concatemer is under a column called "read_idx". For a genomic window of interest, divided into bins, chromunity creates a binary matrix of concatemer by bin overlap. From this matrix, a weighted concatemer by concatemer similarity matrix is computed by quantifying bin-wise intersections between concatemer pairs (i.e. Jaccard similarity). A second similarity matrix is built from the first matrix using the fraction of k-shared nearest neighbor. This matrix, masked on a minimum of k shared nearest neighbors, is used to populate a graph whose concatemers are nodes. Finally, the Louvain method for community detection is used to detect communities of frequently interacting concatemers and these communities are reported as combinatorial chromatin states. 

###  <font color=black> 1. The GRanges of Pore-C output </font>

Here is an exapmle of GRanges output from Pore-C. To generate GRanges, a good starting point are the parquet files created by ONT pore-c pipeline. Details are found here: https://github.com/nanoporetech/pore-c . Using the Apache arrow API found here: https://arrow.apache.org/docs/r/ , these parquet files can be loaded in R in tibble format, that can be converted to GRanges. THe following is the end product that becomes the input for chromunity. The file may have several metadata columns, but shown below are the minimum data needed to run chromunity.

```R
example_gr = readRDS("~/git/chromunity/inst/extdata/example_gr.rds")
example_gr
```
```R
GRanges object with 104730 ranges and 1 metadata column:                                                                                                                                                                                                    
            seqnames          ranges strand |  read_idx                                                                                
               <Rle>       <IRanges>  <Rle> | <integer>                                                                                 
        [1]    chr12 4313520-4313978      * |      5024                                                                                 
        [2]    chr12 4488762-4489153      * |      5024                                                                                 
        [3]    chr12 3979507-3980071      * |      6148                                                                                 
        [4]    chr12 4319746-4320318      * |      7064                                                                                 
        [5]    chr12 4065933-4066241      * |      7886                                                                                 
        ...      ...             ...    ... .       ...                                                                                 
   [104726]    chr12 4193349-4193989      * |  85343061                                                                                 
   [104727]    chr12 4048222-4049016      * |  85343838                                                                                 
   [104728]    chr12 3974695-3975095      * |  85343979                                                                                 
   [104729]    chr12 4464156-4465082      * |  85345063                                                                                 
   [104730]    chr12 4448631-4448896      * |  85345063                                                                                 
   -------                                                                                                                               
   seqinfo: 1 sequence from an unspecified genome
```

###  <font color=black> 2. Running chromunity </font>

chromunity takes in pore-c GRanges object, a GRanges object defining the window of interest, tiled GRanges object spanning the window of interest. k.nn and k.min parameter are set to 25 and 3 respectively. User can change these based on necessity. We found them to be stable so far. 

```R
example_gr = readRDS("~/git/chromunity/inst/extdata/example_gr.rds")
window_gr = readRDS("~/git/chromunity/inst/extdata/window_gr.rds")
tiles_gr = readRDS("~/git/chromunity/inst/extdata/tiles_gr.rds")

chromunity_out_gr = chromunity(example_gr, which.gr = window_gr, tiles = tiles_gr, k.knn = 25, k.min = 3)    
chromunity_out_gr
```
```R
GRanges object with 32262 ranges and 5 metadata columns:                                                                                 
           seqnames          ranges strand |  read_idx       tix     count                                                               
              <Rle>       <IRanges>  <Rle> | <integer> <integer> <integer>                                                               
       [1]    chr12 4313520-4313978      * |      5024        16         2                                                               
       [2]    chr12 4488762-4489153      * |      5024        23         2                                                               
       [3]    chr12 4065933-4066241      * |      7886         6         3                                                               
       [4]    chr12 4103452-4103794      * |      7886         8         3                                                               
       [5]    chr12 4190469-4190667      * |      7886        11         3                                                               
       ...      ...             ...    ... .       ...       ...       ...                                                               
   [32258]    chr12 4123666-4124130      * |  85340541         8         5                                                               
   [32259]    chr12 4233600-4233808      * |  85340541        13         5                                                               
   [32260]    chr12 4234321-4234488      * |  85340541        13         5                                                               
   [32261]    chr12 4464156-4465082      * |  85345063        22         2                                                               
   [32262]    chr12 4448631-4448896      * |  85345063        21         2                                                               
           community  num.memb                                                                                                           
           <numeric> <integer>                                                                                                           
       [1]         9       223                                                                                                           
       [2]         9       223                                                                                                           
       [3]         6       732                                                                                                           
       [4]         6       732                                                                                                           
       [5]         6       732                                                                                                           
       ...       ...       ...                                                                                                           
   [32258]        37       555                                                                                                           
   [32259]        37       555                                                                                                           
   [32260]        37       555                                                                                                           
   [32261]         1       132                                                                                                           
   [32262]         1       132                                                                                                           
   -------                                                                                                                               
   seqinfo: 1 sequence from an unspecified genome 
```
The algorithm adds community annotation and number of members in each community. Other metadata relate to concatemers used internally and can be ignored.

###  <font color=black> 3. Annotating communities that as simple (unimodal) or complex (multimodal) </font>

This is the last step, where the community interactions can be classified as simple or complex based on the frquency distribution of interactions on linear genome. Loess smoothing is used to find peaks. Th smoothing window "w" is defined as 5% of the window of interest. Peaks are called in each "w" and the community is called multimodal if there are 2 or more peaks at least 10% * (width of window of interest) apart. 

```R
chromunity_out_gr = readRDS("~/git/chromunity/inst/extdata/chromunity_out_gr.rds")
annotate_c = annotate_multimodal_communities(granges = chromunity_out_gr, which.gr = window_gr)
```

```R
GRanges object with 32262 ranges and 6 metadata columns:                                                                                 
           seqnames          ranges strand | community  read_idx       tix                                                               
              <Rle>       <IRanges>  <Rle> | <numeric> <integer> <integer>                                                               
       [1]    chr12 4470337-4470845      * |         1    117294        22                                                               
       [2]    chr12 4432256-4432527      * |         1    117294        21                                                               
       [3]    chr12 4482745-4486317      * |         1   1204030        23                                                               
       [4]    chr12 4456416-4459046      * |         1   1204030        22                                                               
       [5]    chr12 4434671-4435764      * |         1   1204030        21                                                               
       ...      ...             ...    ... .       ...       ...       ...                                                               
   [32258]    chr12 4452177-4452452      * |        40  84815407        22                                                               
   [32259]    chr12 4522903-4523488      * |        40  84849078        24                                                               
   [32260]    chr12 4294488-4294958      * |        40  84849078        15                                                               
   [32261]    chr12 4475485-4475860      * |        40  84849078        23                                                               
   [32262]    chr12 4472102-4472414      * |        40  84849078        22                                                               
               count  num.memb multimodal                                                                                               
           <integer> <integer>  <logical>                                                                                               
       [1]         2       132      FALSE                                                                                               
       [2]         2       132      FALSE                                                                                               
       [3]         4       132      FALSE                                                                                               
       [4]         4       132      FALSE                                                                                               
       [5]         4       132      FALSE                                                                                               
       ...       ...       ...        ...                                                                                               
   [32258]         2       135      FALSE                                                                                               
   [32259]         4       135      FALSE                                                                                               
   [32260]         4       135      FALSE                                                                                               
   [32261]         4       135      FALSE                                                                                               
   [32262]         4       135      FALSE                                                                                               

   seqinfo: 1 sequence from an unspecified genome 
```

This adds a column called "multimode" that is boolean if a community is unimodal (simple) or multimodal (complex). Important to note, the output is still pore-c GRanges format we started with where each row is a monomer. Metadata is duplicated for each concatemer and each community. 
