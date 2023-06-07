#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table set
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table setkeyv
#' @importFrom data.table setnames
#' @importFrom data.table transpose
#' @importFrom Matrix sparseMatrix
#' @importFrom plyr round_any
#' @import zoo
#' @import arrow
#' @importFrom pbmcapply pbmclapply 
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom MASS ginv
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#' @import GenomicRanges
#' @import gUtils
#' @importFrom igraph graph.adjacency
#' @importFrom igraph cluster_louvain
#' @importFrom igraph cluster_fast_greedy
#' @importFrom BiocGenerics t

globalVariables(c("::", ":::", "num.memb", "community", "max.local.dist", "read_idx", "as.data.table", "count", "dcast.data.table", "combn", "tot", "copy", "nfrac", "nmat", ".", "bx2", "bx1", "as", "seqlevels", "loess", ".N", ".SD", ":="))


#' @name parquet2gr
#' @description
#' Converts parquet format from Pore-C snakemake pipeline to a single gRanges object
#'
#' @param path string Path to a directory containing all pore_c.parquet files to be read in.
#' @param col_names vector Column names from parquet files to read in.
#' @param save_path string Path to save the gRanges to.
#' @param prefix string File prefix
#' @param mc.cores integer Number of cores in mclapply (default = 5)
#' @param verbose boolean "verbose" flag (default = FALSE)
#' @return GRanges format of the parquet files combined
#' @author Aditya Deshpande, Marcin Imielinski
#' 
#' @export
parquet2gr = function(path = NULL, col_names = NULL, save_path = NULL, prefix = "NlaIII_this_sample", mc.cores = 5, verbose = TRUE){

    if(is.null(path)){
        stop("Need a valid path to all Pore-C parquets.")
    }

    all.paths = data.table(file_path = dir(path, all.files = TRUE, recursive = TRUE, full = TRUE))[grepl("*pore_c.parquet*", file_path)]

    if(nrow(all.paths) == 0){
        stop("No valid files files with suffix pore_c.parquet found.")
    } 

    if(verbose){"Beginning to read parquet files"}

    if(is.null(col_names)){col_names = c("read_name", "chrom", "start", "end", "pass_filter")}
    
    parq.list = pbmclapply(1:nrow(all.paths), function(k){
        parq.al = read_parquet(all.paths[k]$file_path, col_select = col_names)
        parq.al = as.data.table(parq.al)
        return(parq.al)
    }, mc.cores = mc.cores)

    parq.dt = rbindlist(parq.list, fill = TRUE)

    parq.dt[, read_idx := .GRP, by = read_name]
    uni.dt = unique(parq.dt[, .(read_name, read_idx)])
    uni.dt[, fr.read_name := .N, by = read_name]
    uni.dt[, fr.read_idx := .N, by = read_idx]
    if (!any(uni.dt$fr.read_name > 1)){
        if(!any(uni.dt$fr.read_idx > 1)){
            rm(uni.dt)
        } else {
            message("Duplicate read_idx found, check this field in the output")
        }
    } else {
        message("Duplicate read_name found, make all files belong to a single, unique sample")
    }

    parq.gr = dt2gr(parq.dt)

    if (!is.null(save_path)){
        saveRDS(parq.gr, paste0(save_path, "/", prefix, ".rds"))
    }
    return(parq.gr)
}



#' @name csv2gr
#' @description
#' Converts csv format loaded to GEO to a single gRanges object
#'
#' @param path string Path to a directory containing all csv files to be read in.
#' @param col_names vector Column names from csv files to read in.
#' @param save_path string Path to save the gRanges to.
#' @param prefix string File prefix
#' @param mc.cores integer Number of cores in mclapply (default = 5)
#' @param verbose boolean "verbose" flag (default = FALSE)
#' @return GRanges format of the csv files combined
#' @author Aditya Deshpande, Marcin Imielinski
#' 
#' @export
csv2gr = function(path = NULL, col_names = NULL, save_path = NULL, prefix = "NlaIII_this_sample", mc.cores = 5, verbose = TRUE){

    if(is.null(path)){
        stop("Need a valid path to all Pore-C csvs.")
    }

    all.paths = data.table(file_path = dir(path, all.files = TRUE, recursive = TRUE, full = TRUE))[grepl("*fragment_alignments.csv.gz*", file_path)]

    if(nrow(all.paths) == 0){
        stop("No valid files files with suffix pore_c.csv found.")
    } 

    if(verbose){"Beginning to read csv files"}

    if(is.null(col_names)){col_names = c("read_name", "chrom", "start", "end", "pass_filter")}
    
    csv.list = pbmclapply(1:nrow(all.paths), function(k){
        csv.al = fread(all.paths[k]$file_path)[, .(read_name, chrom, start, end, pass_filter)]
        csv.al = as.data.table(csv.al)
        csv.al = csv.al[pass_filter ==  TRUE]  
        return(csv.al)
    }, mc.cores = mc.cores)

    csv.dt = rbindlist(csv.list, fill = TRUE)

    gc()
    
    csv.dt[, read_idx := .GRP, by = read_name]
    uni.dt = unique(csv.dt[, .(read_name, read_idx)])
    uni.dt[, fr.read_name := .N, by = read_name]
    uni.dt[, fr.read_idx := .N, by = read_idx]
    if (!any(uni.dt$fr.read_name > 1)){
        if(!any(uni.dt$fr.read_idx > 1)){
            rm(uni.dt)
        } else {
            message("Duplicate read_idx found, check this field in the output")
        }
    } else {
        message("Duplicate read_name found, make all files belong to a single, unique sample")
    }

    csv.gr = dt2gr(csv.dt)

    if (!is.null(save_path)){
        saveRDS(csv.gr, paste0(save_path, "/", prefix, ".rds"))
    }
    return(csv.gr)
}

    


#' @name re_background
#' @description
#' Given n binsets generates random "background" binsets that mirrors the input binset characteristics with respect to chromosome, width, and distance.
#'
#' Generates (chromosome specific) kernels for width, cardinality, distance and models interchromosomal binsets by allowing
#' a chromosome transition matrix which then lands you on a coordinate that is also chosen from a kernel. 
#' x
#' @param binsets GRanges of bins with fields seqnames, start, end, and $bid specifying binset id
#' @param resolution the resolution to use for distance simulation
#' @param n integer scalar specifying number of sets to sample (default = nrow(binsets))
#' @param pseudocount integer used to populate chr-chr transition matrix
#' @param resolution integer resolution to use for simulating bin-sets
#' @param gg gGnome object to be used for simulating bin-sets in the presence of structural variation
#' @param interchromosomal.dist integer inferred average inter-chr distance to be used in case of whole genome simulation
#' @param interchromosomal.table data.table optional table that contains specific chr-chr inferred distances.
#' @param verbose boolean "verbose" flag (default = FALSE)
#' @param mc.cores integer Number of cores in mclapply (default = 10)
#' @author Aditya Deshpande, Marcin Imielinski
#' @export

re_background = function(binsets, n = length(binsets), pseudocount = 1, resolution = 5e4, gg=NULL, interchromosomal.dist = 1e8, interchromosomal.table = NULL, verbose = TRUE, mc.cores = 5)
{
  if (!length(binsets))
    stop('empty binsets')

  ## sort from left to right
  ## binsets = gr2dt(binsets)
  ## setkeyv(binsets, c('seqnames','start', 'end'))

  if (verbose) smessage('Making kernels')

  ## populate data for kernels
####
  ## Approach to put all dostances in one bag
  binsets$binid = 1:length(binsets)
  if (is.null(gg))
  {
    distance.kernel = gr2dt(binsets)[, as.data.table(expand.grid(i = binid, j = binid))[i<j, ], by = bid] %>% setkeyv(c('i', 'j'))
    distance.kernel[, val := GenomicRanges::distance(binsets[i], binsets[j])]
    distance.kernel = distance.kernel[j-i==1]
  } else {
    if (verbose) smessage('Using graph distance')
######
    distance.kernel = gr2dt(binsets)[, as.data.table(expand.grid(i = binid, j = binid))[i<j, ], by = bid] %>% setkeyv(c('i', 'j'))
    distance.kernel[, val := GenomicRanges::distance(binsets[i], binsets[j])]
    distance.kernel.intra = distance.kernel[!is.na(val)]    
####
    gg = gg$copy$disjoin(disjoin(binsets))
    
    binsetd = data.table(binid = binsets$binid, bid = binsets$bid) %>% setkey('binid')
    binsetd[, gid := gr.match(binsets, gr.chr(gg$nodes$gr))] ## will streamline distance computation to get index
    distance.kernel.g = pbmclapply(unique(binsetd$bid), function(this.bid){
        this.dist = tryCatch(gg$dist(binsetd[bid == this.bid, gid]) %>% melt %>% as.data.table %>% setnames(., c('gi', 'gj', 'val')),
                             error = function(e) NULL)
        if(!is.null(this.dist)){
            this.dist[, bid := as.factor(this.bid)]
            return(this.dist)
         }
    }, mc.cores = mc.cores) %>% rbindlist
    distance.kernel = merge(distance.kernel.g, binsetd, by.x = c("gi", "bid"), by.y = c("gid", "bid"), allow.cartesian=TRUE)
    setnames(distance.kernel, "binid", "i")
    distance.kernel = merge(distance.kernel, binsetd, by.x = c("gj", "bid"), by.y = c("gid", "bid"), allow.cartesian=TRUE)
    setnames(distance.kernel, "binid", "j")
    distance.kernel = unique(distance.kernel[, .(bid, i, j, val)][i<j])
    distance.kernel = merge(distance.kernel, distance.kernel.intra, by = c("bid", "i", "j"), all.x = T)
    distance.kernel[, val := ifelse(val.x <= 1, val.y, val.x)]  
    distance.kernel[, val := round_any(val, resolution, f = ceiling)]
    distance.kernel = distance.kernel[, .(bid, i , j, val)]  
    setkeyv(distance.kernel, c('i', 'j'))
  }

  if (!is.null(interchromosomal.table)){
      interchromosomal.table[, dist := round_any(dist, resolution)]
      distance.kernel = merge(distance.kernel, gr2dt(binsets)[, .(seqnames, binid)], by.x = "i", by.y = "binid")
      distance.kernel = merge(distance.kernel, gr2dt(binsets)[, .(seqnames, binid)], by.x = "j", by.y = "binid")
      setnames(distance.kernel, c("seqnames.x", "seqnames.y"), c("V1", "V2"))
      distance.kernel = merge(distance.kernel, interchromosomal.table, by = c("V1", "V2"), all.x = T)
      distance.kernel[, val := ifelse(is.na(val), dist, val)]
      distance.kernel = distance.kernel[, .(j,   i,  bid, val)]
      distance.kernel[, val := ifelse(is.na(val), max(interchromosomal.table$dist), val)]
      setkeyv(distance.kernel, c('i', 'j'))
  } else {
      distance.kernel[is.infinite(val), val := interchromosomal.dist]
      distance.kernel[is.na(val), val := interchromosomal.dist]
      #distance.kernel = distance.kernel[!is.infinite(val)]
      #distance.kernel = distance.kernel[!is.na(val)] 
  }

  ## distance.kernel[is.na(val), val := interchromosomal.dist]
  
  ####
  binsets = gr2dt(binsets)
  setkeyv(binsets, c('seqnames','start', 'end'))

  ####
  ends.kernel = binsets[, .(val = max(end)), by = .(seqnames, bid)] %>% setkey('seqnames')
  ends.kernel[, freq := .N, by = seqnames]
  if (any(ends.kernel[, freq] == 1)){
      ends.kernel = rbind(ends.kernel[freq == 1], ends.kernel) %>% setkey('seqnames')
  }

####
  ## distance.kernel = binsets[, .(val = start - (data.table::shift(end))), by = .(seqnames, bid)][!is.na(val), ] %>% setkey('seqnames')
  ## distance.kernel[, val := ifelse(is.na(val), 0, val)] %>% setkey('seqnames')
  ## distance.kernel[, freq := .N, by = .(seqnames, bid)]
  ## distance.kernel =  rbind(distance.kernel[freq == 1], distance.kernel) %>% setkey('seqnames')
  
  ####
  width.kernel = binsets[, .(seqnames, val = width)] ## %>% setkey('seqnames')
  ## width.kernel[, freq := .N, by = seqnames]
  ## width.kernel = rbind(width.kernel[freq == 1], width.kernel)%>% setkey('seqnames')
  
  ends.bw =  ends.kernel[, .(bw = density(val)$bw), keyby = seqnames]

  ## distance.bw =  distance.kernel[, .(bw = density(val)$bw), keyby = seqnames]
  ## width.bw =  width.kernel[, .(bw = density(val)$bw), keyby = seqnames]

  distance.bw =  distance.kernel[, .(bw = density(val)$bw)]
  width.bw =  width.kernel[, .(bw = density(val)$bw)]
  
  if (verbose) smessage('Making transition matrices')

  ## populate cardinality 
  cardinality.multinomial = binsets[, .(cardinality = .N), by = .(bid)][, .N, by = cardinality][, .(cardinality, prob = N/sum(N))] %>% setkey(cardinality)

  ## populate chrom multi-nomial and chr-chr transition matrix
  chrom.multinomial = binsets[, .(count = .N + pseudocount), by = seqnames][, .(seqnames, prob = count / sum(count))] %>% setkey('seqnames')

  ## populate chrom multi-nomial and chr-chr transition matrix

  chrom.transition = merge(binsets, binsets, by = 'bid', allow.cartesian = TRUE)[, .(count = .N + pseudocount), by = .(seqnames.x, seqnames.y)][, .(seqnames.y, prob = count / sum(count)), keyby = .(seqnames.x)] 
  
  ## chrom.transition = merge(binsets, binsets, by = 'bid', allow.cartesian = TRUE)[, .(count = .N + pseudocount), by = .(seqnames.x, seqnames.y)]
  ## chrom.transition[, exclude := ifelse(seqnames.x %in% unique(distance.kernel$seqnames), FALSE, TRUE)]
  ## chrom.transition[, count := ifelse(exclude & seqnames.x == seqnames.y, 0, count)]
  ## chrom.transition = chrom.transition[, .(seqnames.y, prob = count / sum(count)), keyby = .(seqnames.x)]
  ## chrom.transition[, prob := ifelse(is.nan(prob), 0, prob)]
  ## chrom.transition[, sum.p := sum(prob), by = seqnames.x]
  ## chrom.transition = chrom.transition[sum.p > 0]
  
  
  if (verbose) smessage('Generating random sets')
  out = pbmclapply(1:n, mc.cores = mc.cores, function(i)
  {
    cardinality = cardinality.multinomial[, sample(cardinality, 1, prob = prob)]
    binset = data.table(bid = i, seqnames = chrom.multinomial[, sample(seqnames, prob = prob, 1)])
    binset$end = rdens(1, ends.kernel[.(binset$seqnames), val], width = ends.bw[.(binset$seqnames), bw])
    binset$start = binset$end - round_any(rdens(1, width.kernel[, val], width = width.bw[, bw]),resolution)
    ## binset$start = binset$end - rdens(1, width.kernel[, val], width = width.bw[, bw])
    for (k in setdiff(1:cardinality, 1))
    {
      lastchrom = tail(binset$seqnames, 1)
      laststart = tail(binset$start, 1)
      newchrom = chrom.transition[.(lastchrom), sample(seqnames.y, 1, prob = prob)]
      ##print(as.character(newchrom))
      if (newchrom == lastchrom)
      {
        ## newbin = data.table(bid = i,
        ##                     seqnames = newchrom,
        ##                     end = laststart - rdens(1, distance.kernel[, val], resolution, width = distance.bw[, bw]))
        ## newbin[, start := end - rdens(1, width.kernel[, val], width = width.bw[, bw])]
        ##
        newbin = data.table(bid = i,
                           seqnames = newchrom,
                           end = laststart - round_any(rdens(1, distance.kernel[, val], resolution, width = distance.bw[, bw]), resolution))
        newbin[, start := end - round_any(rdens(1, width.kernel[, val], width = width.bw[, bw]), resolution)]
        ##print(newbin)
      }
      else
      {
        newbin = data.table(bid = i,
                            seqnames = newchrom,
                            end = rdens(1, ends.kernel[.(newchrom), val], width = ends.bw[.(newchrom), bw]))
        newbin[, start := end - round_any(rdens(1, width.kernel[, val], width = width.bw[, bw]), resolution)]
        ##newbin[, start := end - rdens(1, width.kernel[, val], width = width.bw[, bw])]
        ##print(newbin)
      }
      binset = rbind(binset, newbin, fill = TRUE)
    }
    binset
  }) %>% rbindlist(fill = TRUE)

  ## fix coordinates
  out[, start := pmax(1, ceiling(start))]
  out[, end := pmax(1, ceiling(end))]
  out = gr2dt(dt2gr(out))
  return(out)
}



#' @name rdens
#' @description
#'
#' Function to sample kernel density from a distribution adapted from
#' https://stats.stackexchange.com/questions/321542/how-can-i-draw-a-value-randomly-from-a-kernel-density-estimate
#' 
#' @param n integer number of samples
#' @param values vector data to be used for building the kernel
#' @param width integer kernel width
#' @param kernel string kernel type, use Gaussian
#' @return data.able with values drawn from the kernel

rdens = function(n, values, width = NULL, kernel="gaussian") {
  if (is.null(width)) width = density(values)$bw
  rkernel <- function(n) rnorm(n, sd=width)  # Kernel sampler
  ret = data.table(value = sample(values, n, replace=TRUE) + rkernel(n))
  ret[, value := ifelse(value < 0, value*-1, value)]
  return(ret$value)
}



#' @name rdens_sliding
#' @description
#'
#' Function to sample kernel density from a distribution adapted from
#' https://stats.stackexchange.com/questions/321542/how-can-i-draw-a-value-randomly-from-a-kernel-density-estimate
#' for sliding_window implementation
#' 
#' @param n integer number of samples
#' @param den vector data to be used for building the kernel
#' @param dat integer kernel width
#' @param kernel string kernel type, use Gaussian
#' @return data.able with values drawn from the kernel


rdens_sliding = function(n, den, dat, kernel="gaussian") {
    width <- den$bw                              # Kernel width
    rkernel <- function(n) rnorm(n, sd=width)  # Kernel sampler
    ret = data.table(value = sample(dat, n, replace=TRUE) + rkernel(n))
    ret[, value := ifelse(value < 0, value*-1, value)]
    return(ret$value)
}



#' @name .chr2str
#' @description
#'
#' Internal function to convert chromosome string to int
#' 
#' @param chr_str chromosome string
#' @return int chr


.chr2str = function(chr_str){
    chr.int.ch = gsub("chr", "", chr_str)
    if (any(chr.int.ch == "X")){
        chr.int.ch[which(chr.int.ch == "X")] = 23
    } else if (any(chr.int.ch == "Y")){
        chr.int.ch[which(chr.int.ch == "Y")] = 24
    }
        chr.int = as.integer(chr.int.ch)
    return(chr.int)
}


#' @name extract_dw
#' @description
#'
#' Internal function to extract covariates
#' 
#' @param chromunity.outbinsets from chromunity algorithms
#' @return data.table with covariates



extract_dw = function(chromunity.out, num.cores = 10){
####
    unique.chrom = as.character(unique(chromunity.out$bid))
    message("Generating distributions")
    list.dw.vec = pbmclapply(1:length(unique.chrom), function(i){
        this.uc = unique.chrom[i]
        this.chrom = chromunity.out %Q% (bid == this.uc)
#####
        ## modified for immidiate dists
        dist.vec = as.data.table(gr.dist(this.chrom)[lower.tri(gr.dist(this.chrom))])[, type := "dist"]
        width.vec = as.data.table(width(this.chrom))[, type := "width"]
        this.card = as.data.table(length(this.chrom))[, type := "cardinality"]
        all.vec = rbind(dist.vec, width.vec, this.card)
        all.vec[, bid := this.uc]
        if (this.card$V1 > 2){
            return(all.vec)
        }
        }, mc.cores = num.cores)
    all.dt = rbindlist(list.dw.vec, fill = T)
    return(all.dt)
}




#' @name annotate
#' @description
#'
#' Given n binsets annotates their k-power set (i.e. all sets in the power set with cardinality <= k) with
#' basic (min median max distance, cardinality, width) and optional numerical and interval covariates provided as input.
#' 
#' @param binsets GRanges of bins with fields $bid specifying binset id
#' @param concatemers GRanges of monomers with fields seqnames, start, end, and $cid specifying concatemer id, which will be counted across each binset
#' @param covariates Covariate object specifying numeric or interval covariates to merge with this binset
#' @param k the max cardinality of sets to annotate from the power set of each binmer
#' @param gg optional gGraph input specifying alternate distance function
#' @param interchromosomal.dist numeric scalar of "effective" distance for inter chromosomal bins [1e8]
#' @param verbose logical flag
#' @param mc.cores integer how many cores to parallelize
#' @param threads used to set number of data table threads to use with setDTthreads function, segfaults may occur if >1
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
#' @return data.table of sub-binsets i.e. k-power set of binsets annotated with $count field representing covariates, ready for fitting, **one row per binset

annotate = function(binsets, concatemers, covariates = NULL, k = 5, interchromosomal.dist = 1e8, interchromosomal.table = NULL, gg = NULL, mc.cores = 5, numchunks = 200*mc.cores-1, seed = 42, verbose = TRUE, unique.per.setid = TRUE, resolution = 5e4, threads = 1)
{
  setDTthreads(threads)
  set.seed(seed)
  if (!inherits(binsets, 'GRanges'))
    binsets = dt2gr(binsets)

  binsets = gr.stripstrand(binsets)
  binsets$bid = as.factor(binsets$bid)
  ## each row of binsets is a bin
  binsets$binid = 1:length(binsets)
  #####
  ## frontload / batch some useful calculations
  #####
  ## bin vs concatemer overlaps
  if (verbose) smessage('Overlapping ', length(binsets), ' bins with ', length(concatemers), ' monomers')

  binsets
  concatemers
  ov = binsets %*% concatemers[, 'cid'] %>% as.data.table
  ##
  if (verbose) smessage('Computing bin by bin pairwise distance')
  ## bin x bin pairwise distance within each set
  if (is.null(gg))
  {
    bindist = gr2dt(binsets)[, as.data.table(expand.grid(i = binid, j = binid))[i<j, ], by = bid] %>% setkeyv(c('i', 'j'))
    bindist[, distance := GenomicRanges::distance(binsets[i], binsets[j])]
  } else {
    if (verbose) smessage('Using graph distance')
######
    bindist = gr2dt(binsets)[, as.data.table(expand.grid(i = binid, j = binid))[i<j, ], by = bid] %>% setkeyv(c('i', 'j'))
    bindist[, distance := GenomicRanges::distance(binsets[i], binsets[j])]
    bindist.intra = bindist[!is.na(distance)]    
######
    gg = gg$copy$disjoin(disjoin(binsets))
    binsetd = data.table(binid = binsets$binid, bid = binsets$bid) %>% setkey('binid')
    binsetd[, sun.bin.id := 1:.N, by = bid]

    ## binsetd[, gid := gr.match(binsets, gr.chr(gg$nodes$gr))]
    ## will streamline distance computation to get index
    ## bindist.g = pbmclapply(unique(binsetd$bid), function(this.bid){
    ##     this.dist = tryCatch(gg$dist(binsetd[bid == this.bid, gid]) %>% melt %>% as.data.table %>% setnames(., c('gi', 'gj', 'distance')),
    ##                          error = function(e) NULL)
    ##     if(!is.null(this.dist)){
    ##         this.dist[, bid := as.factor(this.bid)]
    ##         return(this.dist)
    ##      }
    ## }, mc.cores = mc.cores) %>% rbindlist
    ## bindist = merge(bindist.g, binsetd, by.x = c("gi", "bid"), by.y = c("gid", "bid"), allow.cartesian=TRUE)
    ## setnames(bindist, "binid", "i")
    ## bindist = merge(bindist, binsetd, by.x = c("gj", "bid"), by.y = c("gid", "bid"), allow.cartesian=TRUE)
    ## setnames(bindist, "binid", "j")
    ## bindist = unique(bindist[, .(bid, i, j, distance)][i<j])
    ## bindist = merge(bindist, bindist.intra, by = c("bid", "i", "j"), all.x = T)
    ## bindist[, distance := ifelse(distance.x <= 1, distance.y, distance.x)]  
    ## bindist[, distance := round_any(distance, resolution)]
    ## bindist = bindist[, .(bid, i , j, distance)]  
    ## ## bindist[, distance := ifelse(distance < resolution, resolution, distance)]
    ## setkeyv(bindist, c('i', 'j')) 

    ## Better approach to get graph distances, need further testing
    bindist.g = pbmclapply(unique(binsetd$bid), function(this.bid){
        this.dist = tryCatch(gg$dist(gr.nochr(binsets[binsetd[bid == this.bid, binid]])) %>% melt %>% as.data.table %>% setnames(., c('gi', 'gj', 'distance')), error = function(e) NULL)
        if(!is.null(this.dist)){
            this.dist[, bid := as.factor(this.bid)]
            return(this.dist)
        }
    }, mc.cores = mc.cores) %>% rbindlist
    bindist = merge(bindist.g, binsetd, by.x = c("gi", "bid"), by.y = c("sun.bin.id", "bid"), allow.cartesian=TRUE)
    setnames(bindist, "binid", "i")
    bindist = merge(bindist, binsetd, by.x = c("gj", "bid"), by.y = c("sun.bin.id", "bid"), allow.cartesian=TRUE)
    setnames(bindist, "binid", "j")
    bindist = unique(bindist[, .(bid, i, j, distance)][i<j])
    bindist = merge(bindist, bindist.intra, by = c("bid", "i", "j"), all.x = T)
    bindist[, distance := ifelse(distance.x <= 1, distance.y, distance.x)]  
    bindist[, distance := round_any(distance, resolution)]
    bindist = bindist[, .(bid, i , j, distance)]  
    ## bindist[, distance := ifelse(distance < resolution, resolution, distance)]
    setkeyv(bindist, c('i', 'j'))
  }
  if (!is.null(interchromosomal.table)){
      interchromosomal.table[, dist := round_any(dist, resolution)]
      bindist = merge(bindist, gr2dt(binsets)[, .(seqnames, binid)], by.x = "i", by.y = "binid")
      bindist = merge(bindist, gr2dt(binsets)[, .(seqnames, binid)], by.x = "j", by.y = "binid")
      setnames(bindist, c("seqnames.x", "seqnames.y"), c("V1", "V2"))
      bindist = merge(bindist, interchromosomal.table, by = c("V1", "V2"), all.x = T)
      bindist[, distance := ifelse(is.na(distance), dist, distance)]
      bindist = bindist[, .(j,   i,  bid,    distance)]
      bindist[, distance := ifelse(is.na(distance), max(interchromosomal.table$dist), distance)]
      setkeyv(bindist, c('i', 'j'))
  } else {
      bindist[is.infinite(distance), distance := interchromosomal.dist]
      bindist[is.na(distance), distance := interchromosomal.dist]
  }
#####
  ## now start annotating sub.binsets aka sets
  #####
  if (verbose) smessage('Making sub binsets')
  ## make all sub power sets up to k for all bids
  ## each row is now a binid with a setid and bid
  ## adding ones for sum.cov to be removed later
  sub.binsets = gr2dt(binsets)[, powerset(binid, 1, k), by = bid] %>% setnames(c('bid', 'setid', 'binid')) %>% setkey(bid)
  sub.binsets[, ":="(iid = 1:.N, tot = .N), by = .(setid, bid)] ## label each item in each sub-binset, and total count will be useful below

  if (verbose) smessage('Made ', nrow(sub.binsets), ' sub-binsets')

  ## first use ov to count how many concatemers fully overlap all the bins in the subbinset
  if (verbose) smessage('Counting concatemers across sub-binsets across ', mc.cores, ' cores')
  ref.counts = unique(sub.binsets[, .(bid, setid)])[, count := 0]
  if (nrow(ov))
    {
      ubid = unique(sub.binsets$bid) ## split up to lists to leverage pbmclapply
      ubidl = split(ubid, ceiling(runif(length(ubid))*numchunks)) ## randomly chop up ubid into twice the number of mc.coreso
      counts = pbmclapply(ubidl, mc.cores = mc.cores, function(bids)
      {
          out = tryCatch(merge(sub.binsets[.(as.factor(bids)), ], ov[bid %in% bids], by = c('binid', 'bid'), allow.cartesian = TRUE), error = function(e) NULL)
        if (!is.null(out) & nrow(out) > 0){
            if (unique.per.setid){
                out = out[, .(binid, bid, setid, iid, tot, cid)]
                this.step1 = out[, .(hit = all(1:tot[1] %in% iid)), by = .(cid, setid, bid, tot)][hit == TRUE]
                setkeyv(this.step1, c("bid", "cid", "tot"))
                this.step1 = this.step1[, tail(.SD, 1), by = .(cid, bid)]
                this.counts = this.step1[, .(count = sum(hit, na.rm = TRUE)), by = .(setid, bid)]
                this.counts = merge(ref.counts[bid %in% bids], this.counts, by = c("bid", "setid"), all.x = T)
                this.counts[, count := sum(count.x, count.y, na.rm = T), by = .(bid, setid)][, .(bid, setid, count)]
            } else {
                out[, .(hit = all(1:tot[1] %in% iid)), by = .(cid, setid, bid)][, .(count = sum(hit, na.rm = TRUE)), by = .(setid, bid)]
            }
        } else {
            NULL
        }
      })  %>% rbindlist
    }
  
  ## other goodies
  ## changing to mean instead of median
  if (verbose) smessage('Computing min median max distances per setid')
  ## dists = sub.binsets[, bindist[as.data.table(expand.grid(i = binid, j = binid))[i<j, ], .(dist = c('min.dist', 'mean.dist', 'max.dist'), value = as.numeric(summary(distance+1, na.rm = T)[c(1, 4, 6)]))], by = .(setid, bid)] %>% dcast(bid + setid ~ dist, value.var = 'value')

  sub.binsets = sub.binsets[tot > 1]
  ubid = unique(sub.binsets$bid) ## split up to lists to leverage pbmclapply
  ubidl = split(ubid, ceiling(runif(length(ubid))*numchunks)) ## randomly chop up ubid into twice the number of mc.coreso
  dists = pbmclapply(ubidl, mc.cores = mc.cores, function(bids)
  {
      this.sub.binsets = sub.binsets[bid %in% bids]
      this.dists = this.sub.binsets[, bindist[as.data.table(expand.grid(i = binid, j = binid))[i<j, ], .(dist = c('min.dist',  'max.dist'), value = quantile(distance+1,  c(0, 1)))], by = .(setid, bid)] %>% dcast(bid + setid ~ dist, value.var = 'value')
      })  %>% rbindlist


  if (verbose) smessage('Computing marginal sum per bid')
  margs = counts[, .(sum.counts = sum(count)), by = .(bid)]

  if (verbose) smessage('Computing total width and cardinality per setid')
  ubid = unique(sub.binsets$bid) ## split up to lists to leverage pbmclapply
  ubidl = split(ubid, ceiling(runif(length(ubid))*numchunks)) ## randomly chop up ubid into twice the number of mc.coreso
  widths = pbmclapply(ubidl, mc.cores = mc.cores, function(bids)
  {
      this.sub.binsets = sub.binsets[bid %in% bids]
      this.widths = this.sub.binsets[, .(width = sum(width(binsets)[binid])+1, cardinality = .N), by = .(bid, setid)]
  })  %>% rbindlist


  if (verbose) smessage('Merging counts, distance, and width')
  annotated.binsets = unique(sub.binsets[, .(bid, setid)]) %>%
    merge(counts, by = c('bid', 'setid')) %>%
    merge(dists, all.x = TRUE, by = c('bid', 'setid')) %>%
    merge(widths, all.x = TRUE, by = c('bid', 'setid')) %>%
    merge(margs, all.x = TRUE, by = c('bid'))  

  annotated.binsets = annotated.binsets[cardinality > 1]
  
  ## now cycle through interval / numeric covariates
  if (!is.null(covariates))
  {
    setkeyv(sub.binsets, c('bid', 'binid'))
    if (verbose) smessage('Adding covariates')
    for (i in 1:length(covariates))
    {
      if (covariates$type[i] == 'numeric')
        {
          dat = covariates$data[[i]][, covariates$field[[i]]]
          names(values(dat)) = 'val'
        }
      else
        dat = covariates$data[[i]][, c(0)]

      if (verbose) smessage('Overlapping ', length(dat), ' ranges from covariate ', covariates$names[i])
      ov = binsets %*% dat %>% as.data.table

      if (verbose) smessage('Tallying values from ', covariates$type[i], ' covariate ', covariates$names[i])
      covdata = data.table(val = numeric(), bid = factor(), setid = numeric()) %>% setnames('val', covariates$names[i])

      if (nrow(ov))
        {
          sov = merge(sub.binsets, ov, by = c('bid', 'binid'), allow.cartesian = TRUE)
          
          if (covariates[i]$type == 'interval') ## count vs ?sum width TODO: add flag on covariates specifying whether to count
          {
            covdata = sov[ , .(val = .N), by = .(bid, setid)]
            ## covdata = ov[sub.binsets, , allow.cartesian = TRUE][, .(val = sum(width, na.rm = TRUE)), by = .(bid, setid)]          
          }
          else if (covariates[i]$type == 'numeric') ## average value
            covdata = sov[, .(val = sum(width*val, na.rm = TRUE)/sum(width + 0*val, na.rm = TRUE)), by = .(bid, setid)]

          ## eps is small amount to add to covariate to make it loggable
          eps = range(covdata$val, na.rm = TRUE) %>% diff * 1e-4
          covdata$val = covdata$val + eps

          if (any(is.na(covdata$val)) || any(covdata$val<=0))
            warning(paste('Covariate ', covariates$names[[i]], ' has NA or <= values, check or fix before fitting '))
          setnames(covdata, 'val', covariates$names[i])
        }
          
      if (verbose) smessage('Merging data from ', covariates$type[i], ' covariate ', covariates$names[i])
      annotated.binsets = merge(annotated.binsets, covdata, all.x = TRUE, by = c('bid', 'setid'))
    }
  }
  return(annotated.binsets)
}


#' @name powerset
#' @description
#'
#' 
#' @param set vector
#' @param min.k minimum cardinality subset to include in power set
#' @param max.k maximum cardinality subset to include in power set
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
powerset = function(set, min.k = 1, max.k = 5)
  (lapply(pmin(length(set), min.k):pmin(length(set), max.k), function(k) combn(set, k, simplify = FALSE)) %>% do.call('c', .) %>% dunlist)[, .(setid = listid, item = V1)]



#' @name fit
#' @description
#'
#' Fit model to annotated binsets (outputs of annotate)
#'
#' @param annotated.binsets data.table of annotated binsets with at least the fields $bid (binset id) $setid (subset id) $count $min.dist $med.dist $max.dist $width $cardinality $width +/- additional covariates as well, note: every column in this data.table must be a covariate 
#' @param nb logical flag specifying whether to run negative binomial vs. Poisson model
#' @param return.model logical flag whether to return the model or the scored input
#' @param verbose logical flag
#' @return fitted model that can be applied to new cases
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
fit = function(annotated.binsets, nb = TRUE, return.model = TRUE, verbose = TRUE, maxit = 50)
{
  if (!nb) stop('not yet supported')
  ## added sumcounts as cov and width as only offset
  covariates = setdiff(names(annotated.binsets), c('bid', 'setid', 'mean.dist', 'count'))
  fmstring = paste('count ~', paste(paste0('log(', covariates, ')', collapse = ' + ')))
  ## fmstring = paste0(fmstring, " + ", "offset(log(width))")
  fm = formula(fmstring)
##
  model = tryCatch(glm.nb(formula = fm, data = annotated.binsets, control = glm.control(maxit = maxit)), error = function(e) NULL)

  return(list(model = model, covariates = covariates))
}


#' @name sscore
#' @description
#'
#' Takes as input annotated binsets and a model and returns each binset scored according to the model
#' Meaning it annotates the columns with $p, $p.left, p.right, $count.predicted
#'
#' @param annotated.binsets data.table of annotated binsets with field $count $min.width $med.width $max.width $cardinality $width +/- additional covariates as well
#' @param model model output of fit(), the fields of model must match the columns of annotated.binsets
#' @param verbose logical flag 
#' @return annotated.binsets with additional fields $p and $count.predicted specifying predicted count at each binset
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
sscore = function(annotated.binsets, model, verbose = TRUE)
{
  if (is.null(annotated.binsets$count))
    stop('annotsynated.binsets need to have $count column, did you forget to annotate?')

  if (length(missing <- setdiff(model$covariates, names(annotated.binsets))))
    stop('annotated.binsets missing covariates: ', paste(missing, collapse = ', '))

  annotated.binsets$count.predicted = (predict(model$model, type = "response", newdata = annotated.binsets))

  ## randomized p value procedure for negative binomial
  pval = annotated.binsets[, pnbinom(count -1, mu = count.predicted, size = model$model$theta, lower.tail = F)]
  pval.right = annotated.binsets[, pnbinom(count, mu = count.predicted, size = model$model$theta, lower.tail = F)]
  pval.right = ifelse(is.na(pval.right), 1, pval.right)
  pval = ifelse(is.na(pval), 1, pval)
  annotated.binsets[, enrichment := count / count.predicted]
  annotated.binsets$pval = runif(nrow(annotated.binsets), min = pval.right, max = pval)
  return(annotated.binsets)
}

#' @name synergy
#' @description
#'
#' Applies synergy model by either taking as input (or training) a background model using a set of providved (or generated) background.binsets
#' and testing a set of binmers against a set of concatemers.  The model can be generated using a set of covariates. 
#' @export
#' @param concatemers data.table of monomers with fields seqnames, start, end, and $cid specifying concatemer id, which will be counted across each binset
#' @param binsets data.table of bins with fields seqnames, start, end, and $bid specifying binset id
#' @param background.binsets background binsets drawn from output of background or similar, if NULL and model is NULL will be computed [NULL]
#' @param model model to fit to binsets / concatemers, if NULL then model is fit after annotating background.binsets, model features must be compatible with provided covariates 
#' @param covariates covariates to use when annotating binsets and background.binsets
#' @param k cardinality k to which to limit k-power set for each binset or background.binset
#' @param gg gGraph around which to compute bin to bin distances in binsets and also extract copy numbers using node / edge feature $cn
#' @param maxit glmnb max iterations
#' @param mc.cores threads to parallelize
#' @param verbose logical flag whether to fit
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
synergy = function(binsets, concatemers, background.binsets = NULL, model = NULL, covariates = NULL, annotated.binsets = NULL, k = 5, gg = NULL, mc.cores = 5, verbose = TRUE, resolution = 5e4, maxit = 50)
{
  if (!inherits(binsets, 'GRanges'))
    {
      binsets = muffle(dt2gr(binsets))
      if (is.null(binsets))
        stop('binsets failed conversion to GRanges, please provide valid GRnges')
    }

  if (is.null(background.binsets) & is.null(model))
  {
    if (verbose) smessage('Computing random background binsets using features of provided binsets')
    background.binsets = background(binsets, n = 1500, mc.cores = mc.cores, resolution = resolution)
    background.binsets = background.binsets[!bid %in% background.binsets[width < resolution]$bid]
  }
  
  if (is.null(model))
  {
    if (verbose) smessage('Annotating k-power sets of background binsets using features and covariates')
    annotated.background.binsets = annotate(binsets = background.binsets, concatemers = concatemers, gg = gg, covariates = covariates, k = k, mc.cores = mc.cores, verbose = verbose)

    if (verbose) smessage('Fitting model to k-power sets of annotated background binsets')
    model = fit(annotated.background.binsets, nb = TRUE, return.model = TRUE, verbose = verbose, maxit = maxit)
  }

  if (is.null(annotated.binsets))
    annotated.binsets = annotate(binsets = binsets, concatemers = concatemers, covariates = covariates, k = k, gg = gg, mc.cores = mc.cores, verbose = verbose)

  ## browser()
  if (verbose) smessage('Scoring binsets')
  scored.binsets = sscore(annotated.binsets, model, verbose = verbose)

  scored.binsets[, multiway := cardinality > 2]
  setkey(scored.binsets, bid)
  ubid = unique(scored.binsets$bid)

  res = pbmclapply(ubid, function(this.bid) muffle(dflm(glm.nb.fh(data = scored.binsets[.(this.bid),], count ~ multiway + offset(log(count.predicted)), theta = model$model$theta, control = glm.control(maxit = maxit)))[2, ][, name := this.bid])) %>% rbindlist
  setnames(res, 'name', 'bid')

  return(res)
}


#' @name re_chromunity
#' @description
#'
#' Runs genome-wide chromunity detection across a sliding or provided genomic window
#'
#' @param concatemers GRanges with $cid
#' @param resolution bin size for community detection [5e4]
#' @param region region to run on [si2gr(concatemers)]
#' @param windows GRanges or GRangesList of windows to test, more appropriate for testing specific windows, e.g. targets such as promoters/enhancers windows. If sliding window approach is desired, keep this parameter as NULL and use piecewise argument.
#' @param shave  logical flag specifying whether to iteratively "shave" concatemers and bins to a subset C and B where every bin in B has at leat bthresh concatemer support across C and every concatemer in C has order / cardinality of at least cthresh across B, this is done by iteratively removing bins and concatemers until you reach a fixed point
#' @param window.size window size to do community detection within
#' @param max.size max size of problem (concatemers x bins) to consider (default is 2^31-1).  If we hit this problem size, will subsample concatemers so that concatemers * bins is < max.dim, this step is done downstream of shaving 
#' @param tiles.k.knn KNN parameter specifying how many nearest neighbors to sample when building KNN graph
#' @param peak.thresh peak threshold with which to call a peak
#' @param k.min minimal number of nearest neighbors an edge in KNN graph needs to have before community detection
#' @param pad integer pad to use when computing the footprint of each chromunity and finding peak regions which become binsets
#' @return list with items $binset,  $support, $params: $binsets is GRanges of bins with field $bid corresponding to binset id and $support which is the concatemer community supporting the binset which are GRanges with $bid
#' @author Aditya Deshpande, Marcin Imielinski
#' @export

re_chromunity = function(concatemers, resolution = 5e4, region = si2gr(concatemers), windows = NULL, piecewise = TRUE, shave = FALSE, bthresh = 3, cthresh = 3, max.size = 2^31-1, subsample.frac = NULL, window.size = 2e6, max.slice = 1e6, min.support = 5, stride = window.size/2, mc.cores = 5, k.knn = 25, k.min = 5, pad = 1e3, peak.thresh = 0.85, seed = 42, verbose = TRUE)
{
  if (is.null(windows))
      windows = gr.start(gr.tile(region, stride))+window.size/2
      windows = dt2gr(gr2dt(windows)[, start := ifelse(start < 0, 1, start)])
  

  if (inherits(windows, 'GRanges'))
    windows = split(windows, 1:length(windows))

  if (is.null(concatemers$cid))
  {
      if ('read_idx' %in% names(values(concatemers)))
        names(values(concatemers))[match('read_idx', names(values(concatemers)))] = 'cid'
      else
        stop("concatemer GRanges must have metadata column $cid or $read_idx specifying the concatemer id")
  }
  
  params = data.table(k.knn = k.knn, k.min = k.min, seed = seed)

  if (!is.null(resolution)){
      bins = gr.tile(reduce(gr.stripstrand(unlist(windows))), resolution)[, c()]
  } else {
      bins = windows %>% unlist %>% gr.stripstrand %>% disjoin
  }
    
  if (verbose) cmessage('Generated ', length(bins), ' bins across ', length(windows), ' windows')

  if (verbose) cmessage('Matching concatemers with bins, and bins with windows using gr.match with max.slice ', max.slice, ' and ', mc.cores, ' cores')

  ## (batch) match up concatemers with binids
  concatemers$binid = gr.match(concatemers, bins, max.slice = max.slice, mc.cores =  mc.cores, verbose = verbose)

  ## maybe NA need to be removed
  concatemers = concatemers %Q% (!is.na(binid))

  ## match window ids and bins 
  binmap = bins %*% grl.unlist(windows)[, c('grl.ix')] %>% as.data.table %>% setnames('query.id', 'binid') %>% setnames('grl.ix', 'winid') %>% setkeyv('winid')

  ## cycle through (possibly complex) windows call cluster_concatemers and convert to gr.sums
  ## winids = unique(binmap$winid)
  winids = unique(binmap[binid %in% unique(concatemers$binid)]$winid)

  if (shave)
    {
      if (verbose)
        cmessage('Shaving concatemers with bthresh = ', bthresh, ' and cthresh = ', cthresh)
      concatemers = shave_concatemers(concatemers, bthresh = bthresh, cthresh = cthresh, verbose = verbose) 
    }

    ncat = concatemers$cid %>% unique %>% length
    nbin = concatemers$binid %>% unique %>% length


  if (verbose)
      cmessage(sprintf('Running concatemer communities with %s concatemers and %s bins', ncat, nbin))
  
    cc = concatemer_communities(concatemers, k.knn = k.knn, k.min = k.min, seed = seed, max.size = max.size, verbose = verbose, subsample.frac = subsample.frac)

    if (length(cc))
      cc = cc %Q% (support>=min.support)
      
  if (!length(cc))
    return(Chromunity(concatemers = GRanges(), binsets = GRanges(), meta = params))

  uchid = unique((cc %Q% (support >= min.support))$chid)

  if (verbose) cmessage('Analyzing gr.sums associated with ', length(uchid), ' concatemer communities to generate binsets')


  binsets = pbmclapply(uchid, mc.cores = mc.cores, function(this.chid)
  {
    suppressWarnings({
      this.cc = cc %Q% (chid == this.chid)
      peaks = gr.sum(this.cc + pad) %>% gr.peaks('score')
      binset = bins[, c()] %&% (peaks[peaks$score > quantile(peaks$score, peak.thresh)])
      if (length(binset))
      {
        binset$chid = this.chid
        binset$bid = this.chid
        binset$winid = this.cc$winid[1]
        binset = gr.reduce(binset)
      }
    })
    binset
  })  %>% do.call(grbind, .)

  return(Chromunity(concatemers = cc[cc$chid %in% binsets$chid], binsets = binsets, meta = params))
}


#' @name shave_concatemers
#' @description
#'
#' "Shaves" concatemers and bins to a set C and B, respectively such that every bin in B
#' has at least bthresh concatemer support and every concatemer in C has at least cthresh
#' bin-wise order.  This is done by iteratively removing concatemers and bins until a fixed point is reached.
#'
#' @param concatemers GRanges of concatemers with $cid and $binid
#' @param cthresh integer minimum concatemer order across binset B
#' @param bthresh integer minimum bin support across set C of concatemers
#' @param verbose logical flag whether to print diff statements cmessages for each step of iteration
#' @param threads used to set number of data table threads to use with setDTthreads function, segfaults may occur if >1
shave_concatemers = function(concatemers, cthresh = 3, bthresh = 2, verbose = TRUE, threads = 1)
{
  setDTthreads(threads)
  .shave = function(concatemers, bthresh = 2, cthresh = 2)
  {
    dt = unique(gr2dt(concatemers), by = c('cid', 'binid'))
    ccount = dt[, .N, by = cid]
    bcount = dt[, .N, by = binid]
    coolc = ccount[N>=bthresh, cid]
    coolb = bcount[N>=cthresh, binid]
    concatemers %Q% (cid %in% coolc) %Q% (binid %in% coolb) 
  }
  old = concatemers
  new = .shave(concatemers, bthresh = bthresh, cthresh = cthresh)
  while (length(new)<length(old))
  {
    if (verbose)
      cmessage('shaving --> Diff: ', length(old)-length(new), ', Old: ', length(old), ', New:', length(new))
    old = new;
    new = .shave(old, bthresh = bthresh, cthresh = cthresh)
  }
  new
}


#' @name concatemer_communities
#' @description
#'
#' Low level function that labels concatemers with chromunity ids $chid using community detection on a graph. 
#'
#' Given a GRanges of monomers labeled by concatemer id $cid
#'
#' @param concatemers GRanges of monomers with field $cid indicating concatemer id and $binid represent bin id
#' @param tiles.k.knn KNN parameter specifying how many nearest neighbors to sample when building KNN graph
#' @param k.min minimal number of nearest neighbors an edge in KNN graph needs to have before community detection
#' @param drop.small logical flag specifying whether to remove "small" concatemers ie those with a footprint <= small argument [FALSE]
#' @param small integer threshold for bases that define small concatemers, only relevant if drop.small = TRUE
#' @param subsample.frac optional arg specifying fraction of concatemers to subsample [NULL]
#' @param seed seed for subsampling
#' @param threads used to set number of data table threads to use with setDTthreads function, segfaults may occur if >1
#' @return GRanges of concatemers labeled by $c mmunity which specifies community id
concatemer_communities = function (concatemers, k.knn = 25, k.min = 5,
    drop.small = FALSE, small = 1e4, max.size = 2^31-1,
    subsample.frac = NULL, seed = 42, verbose = TRUE, debug = FALSE, threads = 1)  
{
  setDTthreads(threads) #horrible segfaults occur if you don't include this
  reads = concatemers

  if (is.null(reads$cid))
  {
    if ('read_idx' %in% names(values(reads)))
      names(values(reads))[match('read_idx', names(values(reads)))] = 'cid'
    else
      stop("concatemer GRanges must have metadata column $cid or $read_idx specifying the concatemer id")
  }
  
  if (drop.small) {
    if (verbose) cmessage(paste0("Filtering out reads < ", small))
    reads = gr2dt(reads)
    setkeyv(reads, c("seqnames", "start"))
    reads[, `:=`(max.local.dist, end[.N] - start[1]), by = cid]
    reads = reads[max.local.dist > small]
    reads = dt2gr(reads)
  }

  if (debug)
    browser()

  if (verbose) cmessage("Matching reads to tiles")
  reads = as.data.table(reads)[, `:=`(count, .N), by = cid]
  reads2 = reads[count > 2, ] 

  if (!nrow(reads2))
    {
      warning('no high order concatemers, returning empty result')
      return(reads[, chid := NA][c(), ])
    }

  reads2$cid = factor(reads2$cid)
  ucid = levels(reads2$cid)

  if (!is.null(subsample.frac))
  {
      if (verbose) cmessage("Using fraction subsampling")
    setkey(reads2, cid)
    reads2 = reads2[.(sample(ucid, length(ucid)*subsample.frac)), ]
    reads2$cid = factor(reads2$cid)
    ucid = levels(reads2$cid)
  }

    
  if (verbose) cmessage("Matrices made")
  ## gc()
 
  ## remove bins hit only by one concatemer
  reads2$binid = factor(reads2$binid)
  setkey(reads2, binid)
  ubid = reads2[, .N, by = binid][N>1 , binid]
 
  if (!length(ubid))
  {
    warning('No bins hit by two concatemers, returning empty result')
    return(reads[, chid := NA][c(), ])
  }

  reads2 = reads2[.(ubid), ]
  ## added for subsamp
  reads2$binid = factor(reads2$binid) 

  ## refactor may be necessary
  reads2$cid = factor(reads2$cid)
  ucid = levels(reads2$cid)
  
  ## size is the number of pairs which can't exceed the max.size (or integer max)
  size = data.table(cid = reads2$cid, binid = reads2$binid)[, choose(.N,2), by = binid][, sum(as.numeric(V1))]

  ncat = reads2$cid %>% unique %>% length
  nbin = reads2$binid %>% unique %>% length
  
  if (size > max.size)
    stop(sprintf('size %s of the problem %s with %s concatemers and %s bins  exceeds max.size %s, please considering subsampling concatemers, using fewer windows or bins, or shaving concatemers with more aggressive parameters (bthresh, cthresh)', size, ifelse(shave, 'after shaving', ''), ncat, nbin, max.size))


  ## all pairs of concatemers that share a bin
  reads2[, cidi := as.integer(cid)]

  if (verbose) cmessage("Making Pairs object") 
  pairs = reads2[, t(combn(cidi, 2)) %>% as.data.table, by = binid]
  if (verbose) cmessage("Pairs object made") 

  setkey(reads2, cid)
  concatm = reads2[.(ucid),  sparseMatrix(factor(cid, ucid) %>% as.integer, binid %>% as.integer, x = 1, dimnames = list(ucid, levels(binid)))]

  p1 = concatm[pairs[, V1], -1, drop = FALSE]
  p2 = concatm[pairs[, V2], -1, drop = FALSE]
  matching = Matrix::rowSums(p1 & p2)
  total = Matrix::rowSums(p1 | p2)
  dt = data.table(bx1 = pairs[, V1], bx2 = pairs[, V2], mat = matching, 
                  tot = total)[, `:=`(frac, mat/tot)]
  dt2 = copy(dt)
  dt2$bx2 = dt$bx1
  dt2$bx1 = dt$bx2
  dt3 = rbind(dt, dt2)
  dt3$nmat = dt3$mat
  dt3$nfrac = dt3$frac
  setkeyv(dt3, c("nfrac", "nmat"))
  dt3 = unique(dt3)
  dt3.2 = dt3[order(nfrac, nmat, decreasing = T)]
  if (verbose) cmessage("Pairs made")
  gc()
  k = k.knn
  knn.dt = dt3.2[mat > 2 & tot > 2, .(knn = bx2[1:k]), by = bx1][!is.na(knn), ]
  if (!nrow(knn.dt))
  {
    warning('no concatemers with neighbors, perhaps data too sparse? returning empty result')
    return(reads[, chid := NA][c(), ])
  }
  setkey(knn.dt)
  knn = sparseMatrix(knn.dt$bx1, knn.dt$knn, x = 1)
  knn.shared = knn %*% knn
  if (verbose) cmessage("KNN done")
  KMIN = k.min
  A = knn.shared * sign(knn.shared > KMIN)
  A[cbind(1:nrow(A), 1:nrow(A))] = 0
#  A <- as(A, "sparseMatrix")
  A = A + t(A)
  G = graph.adjacency(A, weighted = TRUE, mode = "undirected")
  cl.l = cluster_fast_greedy(G)
  cl = cl.l$membership
  ## rename so chid has most support
  cls = 1:max(cl)
  names(cls) = cl %>% table %>% sort %>% rev %>% names
  cl = cls[as.character(cl)]
  if (verbose) cmessage("Communities made")
  memb.dt = data.table(cid = ucid[1:nrow(A)], chid = cl)
  reads[, cid := as.character(cid)]
  reads = merge(reads, memb.dt, by = "cid")
  reads[, `:=`(support, length(unique(cid))), by = chid]
  reads = dt2gr(reads)
  return(reads)
}



##############

#' @name sliding_window_chromunity
#' @description
#'
#' Runs genome-wide chromunity detection across a sliding or provided genomic window
#'
#' @param concatemers GRanges with $cid
#' @param resolution bin size for community detection [5e4]
#' @param region region to run on [si2gr(concatemers)]
#' @param windows GRanges or GRangesList of windows to test, more appropriate for testing specific windows, e.g. targets such as promoters/enhancers windows. If sliding window approach is desired, keep this parameter as NULL and use piecewise argument.
#' @param window.size window size to do community detection within
#' @param max.size max size of problem (concatemers x bins) to consider (default is 2^31-1).  If we hit this problem size, will subsample concatemers so that concatemers * bins is < max.dim, this step is done downstream of shaving 
#' @param k.knn KNN parameter specifying how many nearest neighbors to sample when building KNN graph
#' @param peak.thresh peak threshold with which to call a peak
#' @param k.min minimal number of nearest neighbors an edge in KNN graph needs to have before community detection
#' @param pad integer pad to use when computing the footprint of each chromunity and finding peak regions which become binsets
#' @param take_sub_sample take subsample to be used for training
#' @param subsample.frac proprtion of concatemers to subsample
#' @return list with items $binset,  $support, $params: $binsets is GRanges of bins with field $bid corresponding to binset id and $support which is the concatemer community supporting the binset which are GRanges with $bid
#' @author Aditya Deshpande, Marcin Imielinski
#' @export

sliding_window_chromunity = function(concatemers, resolution = 5e4, region = si2gr(concatemers), windows = NULL, chr = NULL, take_sub_sample = TRUE, subsample.frac = 0.5, window.size = 2e6, max.slice = 1e6, min.support = 3, stride = window.size/2, mc.cores = 5, k.knn = 25, k.min = 5, pad = 1e3, peak.thresh = 0.85, seed = 145, verbose = TRUE, genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens")
{
    
  set.seed(seed)
    
  if (is.null(chr)){
      stop("Provide chr as this is sliding window implementation")
  }
    
  if (is.null(windows))
    windows = gr.tile(hg_seqlengths(genome = genome), window.size/2)+window.size/4
    windows = dt2gr(gr2dt(windows)[, start := ifelse(start < 0, 1, start)])
    windows = windows %Q% (seqnames == chr)
    windows = sortSeqlevels(windows)
    windows = sort(windows)
   
  if (is.null(concatemers$cid))
  {
      if ('read_idx' %in% names(values(concatemers)))
        names(values(concatemers))[match('read_idx', names(values(concatemers)))] = 'cid'
      else
        stop("concatemer GRanges must have metadata column $cid or $read_idx specifying the concatemer id")
  }
  
  params = data.table(k.knn = k.knn, k.min = k.min, seed = seed)

  if (!is.null(resolution)){
      bins = gr.tile(hg_seqlengths(genome = genome), resolution)
      #bins = bins %Q% (seqnames == chr) 
      #bins = gr.tile(region, resolution)
  } else {
      bins = windows %>% unlist %>% gr.stripstrand %>% disjoin
  }

  if (verbose) cmessage('Generated ', length(bins), ' bins across ', length(windows), ' windows')

  ## (batch) match up concatemers with binids
  concatemers$binid = gr.match(concatemers, bins, max.slice = max.slice, mc.cores =  mc.cores, verbose = verbose)

  ## maybe NA need to be removed
  concatemers = concatemers %Q% (!is.na(binid))

  if (!is.null(subsample.frac)){
      take_sub_sample = TRUE
      if (verbose) cmessage('Taking sub-sample with fraction ', subsample.frac)
  }

    chrom.comm = mclapply(1:length(windows), mc.cores = mc.cores, function(j){
        suppressWarnings({
            which.gr = windows[j]
            these.cids = (concatemers %&% which.gr)$cid
            this.pc.gr =  concatemers %Q% (cid %in% these.cids)
            this.bins = bins %&&% which.gr
            this.chromunity.out = tryCatch(gr2dt(concatemer_chromunity_sliding(this.pc.gr,
                                                                               k.knn = k.knn,
                                                                               k.min = k.min,
                                                                               tiles = this.bins,
                                                                               take_sub_sample = take_sub_sample,
                                                                               frac = subsample.frac,
                                                                               verbose = FALSE)), error = function(e) NULL)  
            if (!is.null(this.chromunity.out)){
                this.chromunity.out[, winid := j]
            } else {
                this.chromunity.out = data.table(NA)
            }
        })
        return(this.chromunity.out)
    })
    chrom.comm = rbindlist(chrom.comm, fill = T)

    chrom.comm[, tix := ifelse(is.na(tix), 0, tix)]
    chrom.comm[, V1 := NULL]
    chrom.comm = na.omit(chrom.comm)
    chrom.comm[, chid := .GRP, by = .(chid, winid)]
    chrom.comm.filt = chrom.comm[support >= min.support]
    chrom.comm.filt = dt2gr(chrom.comm.filt)
      
  if (!length(chrom.comm.filt))
    return(Chromunity(concatemers = GRanges(), binsets = GRanges(), meta = params))

  uchid = unique((chrom.comm.filt %Q% (support >= min.support))$chid)

  if (verbose) cmessage('Analyzing gr.sums associated with ', length(uchid), ' concatemer communities to generate binsets')

  binsets = mclapply(uchid, mc.cores = mc.cores, function(this.chid)
  {
    suppressWarnings({
      this.chrom.comm = chrom.comm.filt %Q% (chid == this.chid)
      winid = unique(this.chrom.comm$winid)
      peaks = gr.sum(this.chrom.comm + pad) 
      binset = bins[, c()] %&% (peaks[peaks$score > quantile(peaks$score, peak.thresh)])
      binset = binset %&% windows[winid]
      if (length(binset))
      {
        binset$chid = this.chid
        binset$bid = this.chid
        binset$winid = winid
        binset = gr.reduce(binset)
      }
    })
    binset
  })  %>% do.call(grbind, .)

  chrom.comm = dt2gr(chrom.comm)
  
  return(Chromunity(concatemers = chrom.comm, binsets = binsets, meta = params))
}



#' @name concatemer_chromunity_sliding
#' @description
#'
#' Low level function that labels concatemers with chromunity ids $chid using community detection on a graph. A faster implementation  for sliding window 
#'
#' Given a GRanges of monomers labeled by concatemer id $cid
#'
#' @param concatemers GRanges of monomers with field $cid indicating concatemer id and $binid represent bin id
#' @param tiles.k.knn KNN parameter specifying how many nearest neighbors to sample when building KNN graph
#' @param k.min minimal number of nearest neighbors an edge in KNN graph needs to have before community detection
#' @param drop.small logical flag specifying whether to remove "small" concatemers ie those with a footprint <= small argument [FALSE]
#' @param small integer threshold for bases that define small concatemers, only relevant if drop.small = TRUE
#' @param take_sub_sample optional arg specifying if fraction of concatemers to subsample [FALSE]  
#' @param subsample.frac optional arg specifying fraction of concatemers to subsample [0.5]
#' @param seed seed for subsampling
#' @param threads used to set number of data table threads to use with setDTthreads function, segfaults may occur if >1
#' @return GRanges of concatemers labeled by $c mmunity which specifies community id

concatemer_chromunity_sliding <- function (concatemers, k.knn = 10, k.min = 1, tiles,  
    drop.small = FALSE, small = NULL, 
    take_sub_sample = FALSE, frac = 0.5, seed.n = 154, verbose = FALSE, threads = 1) 
{
    setDTthreads(threads) #horrible segfaults occur if you don't include this
    reads = concatemers
    
    if (drop.small) {
        if (verbose) cmessage(paste0("Filtering out reads < ", small))
        reads = gr2dt(reads)
        setkeyv(reads, c("seqnames", "start"))
        reads[, `:=`(max.local.dist, end[.N] - start[1]), by = cid]
        reads = reads[max.local.dist > small]
        reads = dt2gr(reads)
    }
    reads$tix = gr.match(reads, tiles)
    reads = as.data.table(reads)[, `:=`(count, .N), by = cid] 
    mat = suppressMessages(dcast.data.table(reads[count > 2, ] %>% gr2dt, cid ~ 
        tix, value.var = "strand", fill = 0))
    mat2 = mat[, c(list(cid = cid), lapply(.SD, function(x) x >= 
        1)), .SDcols = names(mat)[-1]]
    mat2 = suppressWarnings(mat2[, `:=`("NA", NULL)])
    reads.ids = mat2$cid
    mat2 = as.data.table(lapply(mat2, as.numeric))
    if (take_sub_sample) {
        tot.num = nrow(mat2[rowSums(mat2[, -1]) > 1, ])
        if (verbose) cmessage(paste0("Total number of rows are: ", tot.num))
        if (verbose) cmessage("taking a subsample")
        number.to.subsample = pmax(round(tot.num * frac), 1000)
        if (verbose) cmessage(paste0("Number sampled: ", number.to.subsample))
        set.seed(seed.n)
        gt = mat2[rowSums(mat2[, -1]) > 1, ][sample(.N, number.to.subsample), 
            ]
    }
    else {
        gt = mat2[rowSums(mat2[, -1]) > 1, ]
    }
    ubx = gt$cid
    if (verbose) cmessage("Matrices made")
    gc()
    pairs = t(do.call(cbind, apply(gt[, setdiff(which(colSums(gt) > 
        1), 1), with = FALSE] %>% as.matrix, 2, function(x) combn(which(x != 0), 2)))) 
    gt = as(as.matrix(as.data.frame(gt)), "sparseMatrix")    
    p1 = gt[pairs[, 1], -1]
    p2 = gt[pairs[, 2], -1]
    matching = rowSums(as.array(p1 & p2))
    total = rowSums(as.array(p1 | p2))    
    dt = data.table(bx1 = pairs[, 1], bx2 = pairs[, 2], mat = matching, 
        tot = total)[, `:=`(frac, mat/tot)]
    dt2 = copy(dt)
    dt2$bx2 = dt$bx1
    dt2$bx1 = dt$bx2
    dt3 = rbind(dt, dt2)
    dt3$nmat = dt3$mat
    dt3$nfrac = dt3$frac
    setkeyv(dt3, c("nfrac", "nmat"))
    dt3 = unique(dt3)
    dt3.2 = dt3[order(nfrac, nmat, decreasing = T)]
    if (verbose) cmessage("Pairs made")
    gc()
    k = k.knn
    knn.dt = dt3.2[mat > 2 & tot > 2, .(knn = bx2[1:k]), by = bx1][!is.na(knn), 
        ]
    setkey(knn.dt)
    knn = sparseMatrix(knn.dt$bx1, knn.dt$knn, x = 1)
    knn.shared = knn %*% knn
    if (verbose) cmessage("KNN done")
    KMIN = k.min
    A = knn.shared * sign(knn.shared > KMIN)
    A[cbind(1:nrow(A), 1:nrow(A))] = 0
    A <- as(A, "matrix")
    A <- as(A, "sparseMatrix")
    A = A + t(A)
    G = graph.adjacency(A, weighted = TRUE, mode = "undirected")
    cl.l = cluster_fast_greedy(G)
    cl = cl.l$membership
    if (verbose) cmessage("Communities made")
    memb.dt = data.table(cid = ubx[1:nrow(A)], chid = cl)
    reads = merge(reads, memb.dt, by = "cid")
    reads[, `:=`(support, length(unique(cid))), by = chid]
    reads = dt2gr(reads)
    return(reads)
}



#' @name sliding_window_background
#' @description
#'
#' Given n binsets generates random "background" binsets that mirrors the input binset characteristics with respect to chromosome, width, and distance.
#'
#' @param chrtomosome the chromosome to work on, string.
#' @param binsets GRanges of bins with fields seqnames, start, end, and $bid specifying bin-set id
#' @param n number of binsets to generate [1000]
#' @param resolution to use for the simulation [5e4]
#' @param seed seed for subsampling
#' @param genome.to.use genome build to use for the simulation [hg38]
#' @return GRanges of simulated binsets


sliding_window_background = function(chromosome, binsets, seed = 145, n = 1000, resolution = 5e4, num.cores = 10, genome.to.use = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"){
    set.seed(seed)
    this.final.dt = data.table()
    chr.int = .chr2str(chromosome)
    upper.bound = hg_seqlengths(genome = genome.to.use)[chr.int]
#### get the relevant distributions
    dist.pdf.dt = extract_dw(binsets, num.cores = num.cores)
    this.num = 0
    this.tot = 0
    i = 0
    message("Generating GRanges")
    this.list = pbmclapply(1:n, function(i){
        this.iter = i
        this.card = round(rdens_sliding(1, den = density(dist.pdf.dt[type == "cardinality"]$V1),
                                dat = dist.pdf.dt[type == "cardinality"]$V1))
        if (this.card > 0){
            this.dists = tryCatch(round_any(c(0, (rdens_sliding(this.card-1, den = density(dist.pdf.dt[type == "dist"]$V1), dat = dist.pdf.dt[type== "dist"]$V1))), resolution), error = function(e) NULL)
            this.width = tryCatch(round_any(rdens_sliding(this.card, den = density(dist.pdf.dt[type == "width"]$V1), dat = dist.pdf.dt[type == "width"]$V1), resolution),  error = function(e) NULL)
            if (!is.null(this.dists) & !is.null(this.width)){
                this.width = this.width-1
                this.dists = this.dists+1
####
                this.loc.dt = data.table()
                for (j in 1:this.card){
                    if (j == 1){
                        anchor.pt = start(gr.sample(hg_seqlengths(genome = genome.to.use)[chr.int], 1, wid = 1))
                        sts = anchor.pt + this.dists[j]
                        this.gr = GRanges(seqnames = Rle(c(chromosome), c(1)),
                                          ranges = IRanges(c(sts)))
                        this.gr = this.gr + (this.width[j]/2)
                        this.dt = gr2dt(this.gr)
                        this.loc.dt = rbind(this.loc.dt, this.dt)
                    } else {
                        sts = tail(this.loc.dt, n = 1)$end + this.dists[j]
                        ends = sts + this.width[j]
                        this.gr = GRanges(seqnames = Rle(c(chromosome), c(1)),
                                          ranges = IRanges(sts, end = ends))
                        this.dt = gr2dt(this.gr)
                        this.loc.dt = rbind(this.loc.dt, this.dt)
                    }
                }
                this.loc.dt[, bid := paste0("rand_", i)]
                this.loc.gr = dt2gr(this.loc.dt)
                if (!any(start(this.loc.gr) > upper.bound)){
                    this.num = this.num + 1
                    this.final.dt = rbind(this.final.dt, this.loc.dt)
                } else {this.final.dt = data.table(NA)}
            } else {this.final.dt = data.table(NA)}
        } else {this.final.dt = data.table(NA)}
        return(this.final.dt)
    }, mc.cores = num.cores)
    this.dt = rbindlist(this.list, fill = TRUE)
    return(this.dt)
}



#' @name smessage
#' @description
#'
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
#' @private 
smessage = function(..., pre = 'Synergy')
  message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)

#' @name cmessage
#' @description
#'
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
#' @private 
cmessage = function(..., pre = 'Chromunity')
  message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)



##########################################################################
## CHROMUNITY.CLASS
#########################################################################

#' @name Chromunity
#' @title Chromunity
#'
#' 
Chromunity = function(binsets = GRanges(), concatemers = GRanges(), meta = data.table(), verbose = TRUE)
{
  ChromunityObj$new(binsets = binsets, concatemers = concatemers, meta = meta, verbose = verbose)
}


#' @name ChromunityObj
#' @title Chromunity object
#' @description
#'
#' Vectorized object that stores the output of chromunity analysis and can be inputted into annotate / synergy functions.
#' The main accessors are $binsets and $concatemers which each return GRanges linked through a chromunity id $chid field
#'
#' Chromunities can be subsetted, concatenated.  Concatenation will result in deduping any $chid that are present in two inputs. 
ChromunityObj = R6::R6Class("Chromunity", 
  public = list(
    initialize = function(binsets = GRanges(), concatemers = GRanges(), meta = data.table(), verbose = TRUE)
    {
      if (!inherits(binsets, 'GRanges'))
        stop('binsets must non-NULL and a GRanges')
      
      if (!inherits(concatemers, 'GRanges'))
        stop('concatemers must be a GRanges')

      private$pseqlengths = data.table(len = c(seqlengths(binsets), seqlengths(concatemers)),
                                       nm = c(names(seqlengths(binsets)), names(seqlengths(concatemers))))[, max(len, na.rm = TRUE), keyby = nm][, structure(V1, names = nm)]

      if (!is.null(meta))
        {
          if (!is.data.table(meta))
            stop('meta should be a data.table')

          private$pmeta = copy(meta)
        }
      
      ## empty binset return base object 
      if (!length(binsets))
        return(self)

      if (is.null(binsets$chid))
        stop('binset must have chromunity $chid field defined')
      
      private$pbinsets = as.data.table(binsets) %>% setkey(chid)

      private$pchid = unique(binsets$chid)

      ## empty concatemers, nothing to do 
      if (!length(concatemers))
        return(self)

      if (is.null(concatemers$cid))
      {
        if ('read_idx' %in% names(values(concatemers)))
          names(values(concatemers))[match('read_idx', names(values(concatemers)))] = 'cid'
        else
          stop("concatemer GRanges must have metadata column $cid or $read_idx specifying the concatemer id")
      }


      
      if (is.null(concatemers$chid))
      {
        warning('concatemers do not have $chid field pre-defined, doing overlap query to match up ')
        concatemers = concatemers %*% binsets[, 'chid']

      }
      concatemers = as.data.table(concatemers) %>% setkey(chid)

      ## if (!retain){
      ##     if (length(setdiff(concatemers$chid, binsets$chid)))
      ##     {
      ##         warning("concatemers defined that don't map to provided binsets, removing")
      ##         concatemers = concatemers[.(unique(binsets$chid)), ]
      ##     }
      ## }

      ## tally support as keyed vector
      support = merge(private$pbinsets, concatemers, by = 'chid', allow.cartesian = TRUE)[, .(support = length(unique(cid))), keyby = chid]
      private$psupport = support[.(private$pchid), ][is.na(support), support := 0][, chid := private$pchid][, structure(support, names = private$pchid)]

      private$pconcatemers = concatemers %>% copy

      ## all done .. 
      return(invisible(self))
    },
    subset = function(ix, ...){
      if (is.integer(ix) | is.numeric(ix))
        ix = private$pchid[ix]
      
      if (any(is.na(ix)) || length(setdiff(ix, private$pchid)))
        stop('indices outside of ix of object, check query against chromunity id $chid')
      


      out = data.table(chid = ix, chid.new = dedup(ix))
      private$pchid = out$chid.new
      private$psupport = structure(private$psupport[out$chid], names = out$chid.new)

      binsets = merge(private$pbinsets, out, by = 'chid', allow.cartesian = TRUE);
      binsets$chid = binsets$chid.new
      binsets$chid.new = NULL
      private$pbinsets = binsets;

      concatemers = merge(private$pconcatemers, out, by = 'chid', allow.cartesian = TRUE);
      concatemers$chid = concatemers$chid.new
      concatemers$chid.new = NULL
      private$pconcatemers = concatemers;

      return(invisible(self))
    },
    print = function(...)
    {
      quants = quantile(private$psupport, c(0, 0.5, 1))
      message('Chromunity object with ', self$length, ' chromunities (', paste(head(private$pchid, 3), collapse = ', '), ifelse(self$length>3, ', ...', ''),  ') spanning ', nrow(private$pconcatemers), ' concatemers');
      message('\t with median concatemer support:', quants[2], ', range: [', quants[1], '-', quants[3], ']')
    },
    gtrack = function(name = '', binsets = TRUE, concatemers = TRUE, heatmap = FALSE, binsize = 1e4, clim = NA, ...)
    {
      out = gTrack()
      if (binsets)
      {
        binsets = self$binsets
        out = c(out, gTrack(split(binsets, binsets$chid) %>% unname, height = 5, name = paste(name, 'binsets')));
      }
      if (concatemers)
      {
        concatemers = self$concatemers
        out = c(out, gTrack(split(concatemers, concatemers$cid) %>% unname, height = 10, name = paste(name, 'concatemers')));
      }
      if (heatmap)
      {
        concatemers = self$concatemers
        gm = self$gm(binsize = binsize)
        out = c(out, gm$gtrack(clim = clim, ...))
      }
      return(out)
    },
    gm = function(binsize = 5e3) GxG::cocount(self$concatemers, bins = gr.tile(self$footprint, binsize), by = 'chid')
  ),
  private = list(
    pchid = c(),
    pseqlengths = c(),
    pmeta = data.table(),
    psupport = c(),
    pbinsets = data.table(seqnames = factor(), start = integer(), end =  integer(), chid = factor(), support = integer()) %>% setkey(chid),
    pconcatemers = data.table(seqnames = factor(), start = integer(), end =  integer(), chid = factor()) %>% setkey(chid)
    ),
  active = list(
    chid = function() private$pchid,
    meta = function() copy(private$pmeta),
    length = function() length(private$pchid),
    names = function() private$pchid,
    gt = function(x) self$gtrack(),
    footprint = function(x) self$binsets %>% gr.stripstrand %>% sort %>% reduce, 
    binsets = function() dt2gr(private$pbinsets, seqlengths = private$pseqlengths),
    concatemers = function() dt2gr(private$pconcatemers, seqlengths = private$pseqlengths),
    support = function() private$psupport,
    seqlengths = function() private$pseqlengths
  )
)


#' @name [.Chromunity
#' @title Subset chromunity
#' @description
#'
#' Overrides the subset operator x[] for use with Chromunity
#'
#' @param obj Covariate This is the Covariate to be subset
#' @param range vector This is the range of Covariates to return, like subsetting a vector. e.g. c(1,2,3,4,5)[3:4] == c(3,4)
#' @return A new Covariate object that contains only the Covs within the given range
#' @author Marcin Imielinski
#' @export
'c.Chromunity' = function(...){

  ##Ensure that all params are of type Covariate
    cc = list(...)
    isc = sapply(cc, function(x)  class(x)[1] == 'Chromunity')

    if(any(!isc)){
        stop('Error: All inputs must be of class Chromunity.')
    }


  chids = lapply(cc, function(x) x$chid)
  ## this will map old chids to new ones, preventing collisions from chid in the input lists
  chmap = dunlist(chids)[, .(chid.old = V1, listid = listid)][, chid.new := dedup(chid.old)] %>% setkeyv(c('listid', 'chid.old'))

  
  binsets = lapply(1:length(cc), function(x) {tmp = cc[[x]]$binsets; tmp$chid = chmap[.(x, tmp$chid), chid.new];tmp %>% as.data.table})
  concatemers = lapply(1:length(cc), function(x) {tmp = cc[[x]]$concatemers; tmp$chid = chmap[.(x, tmp$chid), chid.new]; tmp %>% as.data.table})
  metas = lapply(cc, function(x) x$meta)

  sl = (lapply(cc, function(x) data.table(nm = names(x$seqlengths), len = x$seqlengths)) %>% rbindlist)[, .(len = max(len)), keyby = nm][, structure(len, names = nm)]
  binsets = rbindlist(binsets) %>% dt2gr(seqlengths = sl)
  concatemers = rbindlist(concatemers) %>% dt2gr(seqlengths = sl)
  return(ChromunityObj$new(binsets = binsets, concatemers = concatemers, meta = rbindlist(metas, fill = TRUE)))
}


#' @name [.Chromunity
#' @title Subset chromunity
#' @description
#'
#' Overrides the subset operator x[] for use with Chromunity
#'
#' @param obj Covariate This is the Covariate to be subset
#' @param range vector This is the range of Covariates to return, like subsetting a vector. e.g. c(1,2,3,4,5)[3:4] == c(3,4)
#' @return A new Covariate object that contains only the Covs within the given range
#' @author Marcin Imielinski
#' @export
'[.Chromunity' = function(obj, range){
  if (any(is.na(range)))
    stop('NA in subscripts not allowed')

  if (any(range>length(obj)))
    stop('Subscript out of bounds')

  ##Clone the object so we don't mess with the original
  ret = obj$clone()
  ##Call the subset function of the Covariate class that will modify the cloned Covariate
  ret$subset(range)

  return (ret)
}

#' @name names.Chromunity
#' @title title
#' @description
#'
#' Overrides the names function for Covariate object
#'
#' @param obj Covariate This is the Covariate whose names we are extracting'
#' @return names of covariates
#' @author Zoran Z. Gajic
#' @export
'names.Chromunity' = function(x) return (x$chid)

#' @name length.Chromunity
#' @title length.Chromunity
#' @description
#' Number of binsets in chromunity object 
#'
#' @param obj Covariate object that is passed to the length function
#' @return number of covariates contained in the Covariate object as defined by length(Covariate$data)
#' @author Zoran Z. Gajic
#' @export
'length.Chromunity' = function(obj,...) return(obj$length)


##########################################################################
## COVARIATE.CLASS
#########################################################################


#' @name Covariate
#' @title title
#' @description
#'
#' Stores Covariates for passing to synergy or annotate function
#'
#' These are GRanges of interval or numeric type.  The numeric can be optionally log transformed while the interval covariates can be
#' count (ie count the total number of intervals in the binset) or width based (count the total width in the binset)
#'
#' 
#'
#' @param name character vector Contains names of the covariates to be created, this should not include the names of any Cov objects passed
#' @param pad numeric vector Indicates the width to extend each item in the covarite. e.g. if you have a GRanges covariate with two ranges (5:10) and (20:30) with a pad of 5,
#' These ranges wil become (0:15) and (15:35)
#' @param type character vector Contains the types of each covariate (numeric, interval, sequencing)
#' @param signature, see ffTrack, a vector of signatures for use with ffTrack sequence covariates
#' fftab signature: signatures is a named list that specify what is to be tallied.  Each signature (ie list element)
#' consist of an arbitrary length character vector specifying strings to %in% (grep = FALSE)
#' or length 1 character vector to grepl (if grep = TRUE)
#' or a length 1 or 2 numeric vector specifying exact value or interval to match (for numeric data)
#' Every list element of signature will become a metadata column in the output GRanges
#' specifying how many positions in the given interval match the given query
#' @param field, a chracter vector for use with numeric covariates (NA otherwise) the indicates the column containing the values of that covarites.
#' For example, if you have a covariate for replication timing and the timings are in the column 'value', the parameter field should be set to the character 'Value'
#' @param na.rm, logical vector that indicates whether or not to remove nas in the covariates
#' @param grep, a chracter vector of  grep for use with sequence covariates of class ffTrack
#' The function fftab is called during the processing of ffTrack sequence covariates grep is used to specify inexact matches (see fftab)
#' @param data, a list of covariate data that can include any of the covariate classes (GRanges, ffTrack, RleList, character)
#' @param log logical flag specifying whether to log the numeric covariate, only applicable to numeric covariate
#' @param count logical flag specifying whether to count the number of intervals 
#' @return Covariate object that can be passed directly to the Chromunity object constructor
#' @author Marcin Imielinski
#' @import R6
#' @export
Covariate = R6::R6Class('Covariate',
    public = list(

      ## See the class documentation
        initialize = function(data = NULL,
                              field = as.character(NA),
                              name = as.character(NA),
                              log = FALSE,
                              count = TRUE, 
                              pad = 0,
                              type = 'numeric',
                              signature = as.character(NA),
                              na.rm = as.logical(NA),
                              grep = as.logical(NA)){

        ##If data are valid and are a list of tracks concatenate with any premade covs
        if(is.null(data))
        {
          self$data = NULL
          self$names = name
          self$type = type
          self$signature = signature
          self$field = field
          self$pad = pad
          self$log = log
          self$count = count
          self$na.rm = na.rm
          self$grep = grep
          return()
        }

        if(class(data) != 'list'){
          data = list(data)
        }

        ## replicate params and data if necessary
        params = data.table(id = 1:length(data), field = field, name = name, pad = pad, type = type, signature = signature, na.rm = na.rm, grep = grep, log = log, count = count)

        if (length(data)==1 & nrow(params)>1)
          data = rep(data, nrow(params))

        self$data = data

        params$dclasses = sapply(self$data, class)

        if (any(ix <<- params$dclasses == 'character'))
        {
          if (any(iix <<- !file.exists(unlist(self$data[ix]))))
          {
            stop(sprintf('Some covariate files not found:\n%s', paste(unlist(self$data[ix][iix]), collapse = ',')))
          }
        }

        ## for any GRanges data that are provided where there is more than one metadata
        ## column, there should be a field given, otherwise we complain
        dmeta = NULL
        if (any(ix <<- params$dclasses != 'character'))
        {
          dmeta = lapply(self$data[ix], function(x) names(values(x)))
        }

        ## we require field to be specified if GRanges have more than one metadata column, otherwise
        ## ambiguous
        if (length(dmeta)>0)
        {
          ## check to make sure that fields actually exist in the provided GRanges arguments
          found = mapply(params$field[ix], dmeta, FUN = function(x,y) ifelse(is.na(x), NA, x %in% y))

          if (any(!found, na.rm = TRUE))
          {
            stop('Some provided Covariate fields not found in their corresponding GRanges metadata, please check arguments')
          }
        }

        if (na.ix <<- any(is.na(params$name)))
        {
          params[na.ix, name := ifelse(!is.na(field), field, paste0('Covariate', id))]
          params[, name := dedup(name)]
        }

        ## label any type = NA covariates for which a field has not been specified
        ## as NA by default
        if (any(na.ix <<- !is.na(params$field) & is.na(params$type)))
        {
          params[na.ix, type := 'numeric']
        }

        if (any(na.ix <<- is.na(params$field) & is.na(params$type)))
        {
          params[na.ix, type := 'interval']
        }

        ## check names to make sure not malformed, i.e. shouldn't start with number or contain spaces or special
        ## characters

        if (any(iix <<- grepl('\\W', params$name)))
        {
          warning('Replacing spaces and special characters with "_" in Covariate names')
          params$names[iix] = gsub('\\W+', '_', params$name[iix])
        }
        
        if (!is.null(params$name))
          {
            if (any(iix <<- duplicated(params$name)))
            {
              warning('Deduping covariate names')
              params$name = dedup(params$name)
            }
          }

        self$names = params$name
        self$type = params$type
        self$signature = params$signature
        self$field = params$field
        self$pad = params$pad
        self$log = params$log
        self$count = params$count
        self$na.rm = params$na.rm
    },

    ## Params:
    ## ... Other Covariates to be merged into this array, note that it can be any number of Covariates
    ## Return:
    ## A single Covariate object that contains the contents of self and all passed Covariates
    ## UI:
    ## None
    ## Notes:
    ## This is linked to the c.Covariate override and the c.Covariate override should be used preferentially over this
    merge = function(...){
        return (c(self,...))
    },


    ## Params:
    ## No params required, included arguments will be ignored.
    ## Return:
    ## returns a list of character vectors. If the respective covariate is of class GRanges, the vector will contain all of the chromosome names,
    ## if it is not of class GRanges, will return NA
    ## UI:
    ## None
    seqlevels = function(...){
        if(length(private$pCovs) == 0){
            return(NULL)
        }
        seqs = lapply(c(1:length(private$pCovs)), function(x){
            cov = private$pCovs[[x]]
            if(class(cov) == 'GRanges'){
                return(GenomeInfoDb::seqlevels(cov))
            } else{
                return(NA)
            }
        })
        return(seqs)
    },

    ## Params:
    ## range, a numeric vector of the covariates to include. e.g. if the Covariate contains the covariates (A,B,C) and the range is c(2:3),
    ## this indicates you wish to get a Covariate containing (B,C). NOTE THAT THIS DOES NOT RETURN A NEW COV_ARR, IT MODIFIES THE CURRENT.
    ## Return:
    ## None, this modifies the Covariate on which it was called
    ## UI:
    ## None
    ## Notes:
    ## If you want to create a new Covariate containing certain covariates, use the '[' operator, e.g. Covariate[2:3]
    subset = function(range, ...){

      private$pCovs = private$pCovs[range]
      private$pnames = private$pnames[range]
      private$ptype = private$ptype[range]
      private$psignature = private$psignature[range]
      private$pfield = private$pfield[range]
      private$ppad = private$ppad[range]
      private$plog = private$plog[range]
      private$pcount = private$pcount[range]
      private$pna.rm = private$pna.rm[range]
    },

    ## Params:
    ## No params required, included arguments will be ignored.
    ## Return:
    ## A list of lists where each internal list corresponds to the covariate and is for use internally in the annotate.hypotheses function
    ## The list representation of the covariate will contain the following variables: type, signature, pad, na.rm, field, grep
    ## UI:
    ## None
    toList = function(...){
        if(length(private$pCovs) == 0){
            return(list())
        }
        out = lapply(c(1:length(private$pCovs)), function(x){
            return (list(track = private$pCovs[[x]],
                         type = private$ptype[x],
                         signature = private$psignature[x],
                         pad = private$ppad[x],
                         na.rm = private$pna.rm[x],
                         field = private$pfield[x],
                         log = private$plog[x],
                         count = private$pcount[x]
                         ))
        })
        names(out) = private$pnames
        return(out)

        },

    ## Params:
    ## No params required, included arguments will be ignored.
    ## Return:
    ## Nothing
    ## UI:
    ## Prints information about the Covariate to the console with all of covariates printed in order with variables printed alongside each covariate
    print = function(...){
        if(length(private$pCovs) == 0){
            fmessage('Empty Covariate Object')
            return(NULL)
        }

        message('', length(private$pCovs), ' Covariates with features:')
        print(data.table(
          name = private$pnames,
          type = private$ptype,
          class = sapply(private$pCovs, class),
          field = private$pfield,
          signature = private$psignature,
          na.rm = private$pna.rm,
          pad = private$ppad,
          log = private$plog,
          count = private$pcount
        ))

        ## out= sapply(c(1:length(private$pCovs)),
        ##     function(x){
        ##         cat(c('Covariate Number: ' , x, '\nName: ', grepprivate$pnames[x],
        ##         '\ntype: ',private$ptype[x], '\tsignature: ', private$psignature[x],
        ##         '\nfield: ',private$pfield[x], '\tpad: ', private$ppad[x],
        ##         '\nna.rm: ', private$pna.rm[x], '\tgrep: ', private$pgrep[x],
        ##         '\nCovariate Class: ', class(private$pCovs[[x]]), '\n\n'), collapse = '', sep = '')
    }
    ),

    ## Private variables are internal variables that cannot be accessed by the user
    ## These variables will have active representations that the user can interact with the update
    ## and view these variables, all internal manipulations will be done with these private variables
    private = list(
      ## The list of covariates, each element can be of class: 'GRanges', 'character', 'RleList', 'ffTrack'
      pCovs = list(),
      ## A string vector containing the names of the covariates, the covariate will be refered to by its name in the final table
      pnames = c(),
      ## Type is a string vector of types for each covariate, can be: 'numeric','sequence', or 'interval'
      ptype = c(),
      ## A vector of signatures for use with ffTrack, se fftab
      psignature = c(),
      ## A character vector of field names for use with numeric covariates, see the Covariate class definition for more info
      pfield = c(),
      ## A numeric vector of paddings for each covariate, see the 'pad' param in Covariate class definition for more info
      ppad = c(),
      ## A logical vector for each covariate, see the 'na.rm' param in Covariate class definition for more info
      pna.rm = c(),
      ## logical vector specifying whether to log
      plog = c(),
      ## logical vector specifying whether to count
      pcount = c(),
      ##  Valid Covariate Types
      COV.TYPES = c('numeric', 'interval'),
      ##  Valid Covariate Classes
      COV.CLASSES = c('GRanges', 'RleList', 'ffTrack', 'character')
    ),

    ## The active list contains a variable for each private variable.
    ## Active variables are for user interaction,
    ## Interactions can be as such
    ## class$active will call the active variable function with the value missing
    ## class$active = value will call the active variable function with the value = value
    active = list(

            ## Covariate Names
            ## Here we check to make sure that all names are of class chracter and that they are the same length as pCovs -> the internal covariate list
            ## If the names vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the names vector, will will replicate the names vector such that it matches in length
            ## to the pCovs list.
            names = function(value) {

                if(!missing(value)){

                  if(!is.character(value) && !all(is.na(value)) ){
                    stop('Error: names must be of class character')
                  }

                  if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                    stop('Error: Length of names must be of length equal to the number of Covariates or a divisor of number of covariates.')
                  }

                  if(length(private$pCovs) / length(value) != 1){
                    private$pnames = rep(value, length(private$pCovs)/length(value))
                    return(private$pnames)
                  }

                  if (any(iix <<- grepl('\\W', value)))
                  {
                    warning('Replacing spaces and special characters with "_" in provided Covariate names')
                    value[iix] = gsub('\\W+', '_', value[iix])
                  }

                  private$pnames = value
                  return(private$pnames)

                } else{
                  return(private$pnames)
                }
            },

            ## Covariate type
            ## Here we check to make sure that all types are of class chracter and that they are the same length as pCovs -> the internal covariate list
            ## We then check to make sure that each type is a valid type as defined by the COV.TYPES private parameter
            ## If the types vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the types vector, will will replicate the types vector such that it matches in length
            ## to the pCovs list.
            type = function(value) {

                if(!missing(value)){
                    if(!is.character(value) && !all(is.na(value))){
                        stop('Error: type must be of class character')
                    }

                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop('Error: Length of type must be of length equal to the number of Covariates or a divisor of number of covariates.')
                    }

                    if(!all(value %in% private$COV.TYPES)){
                        stop('"type" must be "numeric", "sequence", or "interval"')
                    }

                    if(length(private$pCovs) / length(value) != 1){
                        private$ptype = rep(value, length(private$pCovs)/length(value))
                        return(private$ptype)
                    }

                    private$ptype = value
                    return(private$ptype)

                } else{
                    return(private$ptype)
                }
            },

            ##Covariate Signature
            ## Here we check to make sure that all signatures are list within lists  and that they are the same length as pCovs -> the internal covariate list
            ## If the signature vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the signature vector, will will replicate the signature vector such that it matches in length
            ## to the pCovs list.
            signature = function(value) {
                if(!missing(value)){

                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop('Error: Length of signature must be of length equal to the number of Covariates or a divisor of number of covariates.')
                    }
                    if(length(private$pCovs) / length(value) != 1){
                        private$psignature = rep(value, length(private$pCovs)/length(value))
                        return(private$psignature)
                    }

                    private$psignature = value
                    return(private$psignature)

                } else{
                    return(private$psignature)
                }
            },

            ##Covariate Field
            ## Here we check to make sure that all fields are of class chracter and that they are the same length as pCovs -> the internal covariate list
            ## If the fields vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the fields vector, will will replicate the fields vector such that it matches in length
            ## to the pCovs list.
            field = function(value) {
                if(!missing(value)){
                    if(!is.character(value) && !all(is.na(value))){
                       stop('Error: field must be of class character')
                    }
                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop('Error: Length of field must be of length equal to the number of Covariates or a divisor of number of covariates.')
                    }

                    if(length(private$pCovs) / length(value) != 1){
                        private$pfield = rep(value, length(private$pCovs)/length(value))
                        return(private$pfield)
                    }

                    private$pfield = value
                    return(private$pfield)

                } else{
                    return(private$pfield)
                }
            },

          ## Covariate log
          ## Here we check to make sure that all pad are of class numeric and that they are the same length as pCovs -> the internal covariate list
          ## If the pad vector is equal in length to the pCovs list length we will allow the assignment
          ## If the pCovs list length is a multiple of the pad vector, will will replicate the pad vector such that it matches in length
          ## to the pCovs list.
          log = function(value) {
            if(!missing(value)){
              if(!is.logical(value) && !all(is.na(value))){
                stop("Error: log must be of class logical")
              }
              if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                stop("Error: Length of pad must be of length equal to the number of Covariates or a divisor of number of covariates.")
              }
              if(length(private$pCovs) / length(value) != 1){
                private$plog = rep(value, length(private$pCovs)/length(value))
                return(private$plog)
              }
              
              private$plog = value
                    return(private$plog)
              
            } else{
              return(private$plog)
            }
          },
          
          ## Covariate Paddinig
          ## Here we check to make sure that all pad are of class numeric and that they are the same length as pCovs -> the internal covariate list
          ## If the pad vector is equal in length to the pCovs list length we will allow the assignment
          ## If the pCovs list length is a multiple of the pad vector, will will replicate the pad vector such that it matches in length
          ## to the pCovs list.
          count = function(value) {
            if(!missing(value)){
              if(!is.logical(value) && !all(is.na(value))){
                stop("Error: count must be of class logical")
              }
              if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                stop("Error: Length of count must be of length equal to the number of Covariates or a divisor of number of covariates.")
              }
              if(length(private$pCovs) / length(value) != 1){
                private$pcount = rep(value, length(private$pCovs)/length(value))
                return(private$pcount)
              }

              private$pcount = value
              return(private$pcount)

            } else{
              return(private$pcount)
            }
          },



            ## Covariate Padding
            ## Here we check to make sure that all pad are of class numeric and that they are the same length as pCovs -> the internal covariate list
            ## If the pad vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the pad vector, will will replicate the pad vector such that it matches in length
            ## to the pCovs list.
            pad = function(value) {
                if(!missing(value)){
                    if(!is.numeric(value) && !all(is.na(value))){
                        stop("Error: pad must be of class numeric")
                    }
                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop("Error: Length of pad must be of length equal to the number of Covariates or a divisor of number of covariates.")
                    }
                    if(length(private$pCovs) / length(value) != 1){
                        private$ppad = rep(value, length(private$pCovs)/length(value))
                        return(private$ppad)
                    }

                    private$ppad = value
                    return(private$ppad)

                } else{
                    return(private$ppad)
                }
            },




            ##Covariate na.rm
            ## Here we check to make sure that all na.rms are of class logical and that they are the same length as pCovs -> the internal covariate list
            ## If the na.rms vector is equal in length to the pCovs list length we will allow the assignment
            ## If the pCovs list length is a multiple of the na.rms vector, will will replicate the na.rms vector such that it matches in length
            ## to the pCovs list.
            na.rm = function(value) {
                if(!missing(value)){
                    if(!is.logical(value) && !all(is.na(value))){
                        stop("Error: na.rm must be of class logical")
                    }
                    if(length(value) != length(private$pCovs) & length(private$pCovs) %% length(value) != 0){
                        stop("Error: Length of na.rm must be of length equal to the number of Covariates or a divisor of number of covariates.")
                    }

                    if(length(private$pCovs) / length(value) != 1){
                        private$pna.rm = rep(value, length(private$pCovs)/length(value))
                        return(private$pna.rm)
                    }

                    private$pna.rm = value
                    return(private$pna.rm)

                } else{
                    return(private$pna.rm)
                }
            },

           

            ##Covariate Covs
            data = function(value) {
                if(!missing(value)){
                    private$pCovs = value
                    return(private$pCovs)
                } else{
                    return(private$pCovs)
                }
            }
    ),
)

#' @name c.Covariate
#' @title title
#' @description
#'
#' Override the c operator for covariates so that you can merge them like a vector
#'
#' @param ... A series of Covariates, note all objects must be of type Covariate
#' @return Covariate object that can be passed directly into the Chromunity object constructor that contains all of the Covariate covariates
#' Passed in the ... param
#' @author Zoran Z. Gajic
#' @export
'c.Covariate' = function(...){

  ##Ensure that all params are of type Covariate
    Covariates = list(...)
    isc = sapply(Covariates, function(x)  class(x)[1] == 'Covariate')

    if(any(!isc)){
        stop('Error: All inputs must be of class Covariate.')
    }

    ## Merging vars of the covariates
  names  = unlist(lapply(Covariates, function(x) x$names))
  type  = unlist(lapply(Covariates, function(x) x$type))
  signature  = unlist(lapply(Covariates, function(x) x$signature))
  field  = unlist(lapply(Covariates, function(x) x$field))
  pad  = unlist(lapply(Covariates, function(x) x$pad))
  na.rm  = unlist(lapply(Covariates, function(x) x$na.rm))
  log  = unlist(lapply(Covariates, function(x) x$log))
  count  = unlist(lapply(Covariates, function(x) x$count))

  ## Merging Covariates
  covs = lapply(Covariates, function(x) x$data)
  Covs = unlist(covs, recursive = F)

  ##Creating a new Covariate and assigning all of the merged variables to it
  ret = Covariate$new(data = Covs, name = names, type = type, signature = signature, field = field, pad = pad, na.rm = na.rm, grep = grep)
  
  return(ret)
}



#' @name [.Covariate
#' @title title
#' @description
#'
#' Overrides the subset operator x[] for use with Covariate to allow for vector like subsetting
#'
#' @param obj Covariate This is the Covariate to be subset
#' @param range vector This is the range of Covariates to return, like subsetting a vector. e.g. c(1,2,3,4,5)[3:4] == c(3,4)
#' @return A new Covariate object that contains only the Covs within the given range
#' @author Zoran Z. Gajic
#' @export
'[.Covariate' = function(obj, range){
  if (any(is.na(range)))
    stop('NA in subscripts not allowed')

  if (any(range>length(obj)))
    stop('Subscript out of bounds')

  ##Clone the object so we don't mess with the original
  ret = obj$clone()
  ##Call the subset function of the Covariate class that will modify the cloned Covariate
  ret$subset(range)
  return (ret)
}


#' @name names.Covariate
#' @title title
#' @description
#'
#' Overrides the names function for Covariate object
#'
#' @param obj Covariate This is the Covariate whose names we are extracting'
#' @return names of covariates
#' @author Zoran Z. Gajic
#' @export
'names.Covariate' = function(x){
    ##Call the subset function of the Covariate class that will modify the cloned Covariate
    return (x$names)
}


#' @name length.Covariate
#' @title title
#' @description
#'
#' Overrides the length function 'length(Covariate)' for use with Covariate
#'
#' @param obj Covariate object that is passed to the length function
#' @return number of covariates contained in the Covariate object as defined by length(Covariate$data)
#' @author Zoran Z. Gajic
#' @export
'length.Covariate' = function(obj,...){
    return(length(obj$data))
}


#' @name covariate
#' @title covariate
#' @description
#'
#' function to initialize Covariates for passing to Chromunity object constructor.
#'
#' Can also be initiated by passing a vector of multiple vectors of equal length, each representing one of the internal variable names
#' You must also include a list containg all of the covariates (Granges, chracters, RLELists, ffTracks)
#'
#' Covariate serves to mask the underlieing list implemenations of Covariates in the Chromunity Object.
#' This class attempts to mimic a vector in terms of subsetting and in the future will add more vector like operations.
#'
#'
#' @param name character vector Contains names of the covariates to be created, this should not include the names of any Cov objects passed
#' @param log  logical flag whether to log output for numeric covariate after averaging [TRUE]
#' @param count logical flag whether to count the number of intervals (count = TRUE) or aggregate the width (count = FALSE) [TRUE]
#' @param pad numeric vector Indicates the width to extend each item in the covarite. e.g. if you have a GRanges covariate with two ranges (5:10) and (20:30) with a pad of 5,
#' These ranges wil become (0:15) and (15:35)
#' @param type character vector Contains the types of each covariate (numeric, interval, sequencing)
#' @param signature, see ffTrack, a vector of signatures for use with ffTrack sequence covariates
#' fftab signature: signatures is a named list that specify what is to be tallied.  Each signature (ie list element)
#' consist of an arbitrary length character vector specifying strings to %in% (grep = FALSE)
#' or length 1 character vector to grepl (if grep = TRUE)
#' or a length 1 or 2 numeric vector specifying exact value or interval to match (for numeric data)
#' Every list element of signature will become a metadata column in the output GRanges
#' specifying how many positions in the given interval match the given query
#' @param field, a chracter vector for use with numeric covariates (NA otherwise) the indicates the column containing the values of that covarites.
#' For example, if you have a covariate for replication timing and the timings are in the column 'value', the parameter field should be set to the character 'Value'
#' @param na.rm, logical vector that indicates whether or not to remove nas in the covariates
#' @param grep, a chracter vector of  grep for use with sequence covariates of class ffTrack
#' The function fftab is called during the processing of ffTrack sequence covariates grep is used to specify inexact matches (see fftab)
#' @param data, a list of covariate data that can include any of the covariate classes (GRanges, ffTrack, RleList, character)
#' @return Covariate object that can be passed directly to the Chromunity object constructor
#' @author Zoran Z. Gajic
#' @import R6
#' @export
covariate = function(data = NULL, field = as.character(NA), name = as.character(NA), log = TRUE, count = TRUE, pad = 0, type = as.character(NA),
               signature = as.character(NA), na.rm = NA, grep = NA){
  Covariate$new(name = name, data = data, pad = pad, type = type, count = count, log = log, signature = signature,
                field = field, na.rm = na.rm, grep = grep)
}



#' @name dunlist
#' @title dunlist
#' @description
#'
#' Utility function to unlist a list into a data.table
#'
#' @param x list
dunlist = function (x) 
{
    listid = rep(1:length(x), elementNROWS(x))
    if (!is.null(names(x))) 
        listid = names(x)[listid]
    xu = unlist(x, use.names = FALSE)
    if (is.null(xu)) {
        return(as.data.table(list(listid = c(), V1 = c())))
    }
    if (!(inherits(xu, "data.frame")) | inherits(xu, "data.table")) 
        xu = data.table(V1 = xu)
    out = cbind(data.table(listid = listid), xu)
    setkey(out, listid)
    return(out)
}


################################
#' @name dedup
#' @title dedup
#'
#' @description
#' relabels duplicates in a character vector with .1, .2, .3
#' (where "." can be replaced by any user specified suffix)
#'
#' @param x input vector to dedup
#' @param suffix suffix separator to use before adding integer for dups in x
#' @return length(x) vector of input + suffix separator + integer for dups and no suffix for "originals"
#' @author Marcin Imielinski
################################
dedup = function(x, suffix = '.')
{
  dup = duplicated(x);
  udup = setdiff(unique(x[dup]), NA)
  udup.ix = lapply(udup, function(y) which(x==y))
  udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
  out = x;
  out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
  return(out)
}


################################
#' @name glm.nb.fh
#' @title glm.nb.fh
#'
#' @description
#' modified glm.nb that takes a pre-defined theta
#' 
#' 
#' @author Marcin Imielinski
################################
glm.nb.fh = function (formula, data, weights, subset, na.action, start = NULL,
                      etastart, mustart, control = glm.control(...), method = "glm.fit",
                      model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,  theta = NULL,
                      init.theta, link = log)
{
    loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th +
                                                        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
                                                 log(mu + (y == 0)) - (th + y) * log(th + mu)))
    link <- substitute(link)
    fam0 <- if (missing(init.theta))
                do.call("poisson", list(link = link))
            else do.call("negative.binomial", list(theta = init.theta,
                                                   link = link))
    mf <- Call <- match.call()
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(mf)
    Terms <- attr(mf, "terms")
    if (method == "model.frame")
        return(mf)
    Y <- model.response(mf, "numeric")
    X <- if (!is.empty.model(Terms))
             model.matrix(Terms, mf, contrasts)
         else matrix(, NROW(Y), 0)
    w <- model.weights(mf)
    if (!length(w))
        w <- rep(1, nrow(mf))
    else if (any(w < 0))
        stop("negative weights not allowed")
    offset <- model.offset(mf)
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    n <- length(Y)
    if (!missing(method)) {
        if (!exists(method, mode = "function"))
            stop(gettextf("unimplemented method: %s", sQuote(method)),
                 domain = NA)
        glm.fitter <- get(method)
    }
    else {
        method <- "glm.fit"
        glm.fitter <- stats::glm.fit
    }
    if (control$trace > 1)
        message("Initial fit:")

    fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart,
                      mustart = mustart, offset = offset, family = fam0, control = list(maxit = control$maxit,
                                                                                        epsilon = control$epsilon, trace = control$trace >
                                                                                                                       1), intercept = attr(Terms, "intercept") > 0)
    class(fit) <- c("glm", "lm")
    mu <- fit$fitted.values
    if (is.null(theta))
    {
        th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                                 trace = control$trace > 2))
    }
    else
    {
        th = theta

        if (control$trace > 1)
            message(gettextf("Fixing theta value to 'theta': %f", signif(th)),
                    domain = NA)
    }

    if (control$trace > 1)
        message(gettextf("Initial value for 'theta': %f", signif(th)),
                domain = NA)
    fam <- do.call("negative.binomial", list(theta = th[1], link = link))
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    g <- fam$linkfun
    Lm <- loglik(n, th, mu, Y, w)
    Lm0 <- Lm + 2 * d1
    while (
    ((iter <- iter + 1) <= control$maxit) &
    ((abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon)
    ){
        eta <- g(mu)
        fit <- glm.fitter(x = X, y = Y, w = w, etastart = eta,
                          offset = offset, family = fam, control = list(maxit = control$maxit,
                                                                        epsilon = control$epsilon, trace = control$trace >
                                                                                                       1), intercept = attr(Terms, "intercept") >
                                                                                                               0)
        t0 <- th
        if (is.null(theta))
        {
            th <- theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                           trace = control$trace > 2)
        } else
        {
            th = theta
        }

        fam <- do.call("negative.binomial", list(theta = th[1],  ## we don't need all the thetas here if theta is vectorized
                                                 link = link))

        mu <- fit$fitted.values
        del <- t0 - th ## this is where the vectorized theta makes a difference
        Lm0 <- Lm
        Lm <- loglik(n, th, mu, Y, w) ## and here - log likelihood computation
        if (control$trace) {
            Ls <- loglik(n, th, Y, Y, w)
            Dev <- 2 * (Ls - Lm)
            message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f",
                            iter, signif(th), signif(Dev)), domain = NA)
        }
    }
    if (!is.null(attr(th, "warn")))
        fit$th.warn <- attr(th, "warn")
    if (iter > control$maxit) {
        warning("alternation limit reached")
        fit$th.warn <- gettext("alternation limit reached")
    }
    if (length(offset) && attr(Terms, "intercept")) {
        null.deviance <- if (length(Terms))
                             glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w,
                                        offset = offset, family = fam, control = list(maxit = control$maxit,
                                                                                      epsilon = control$epsilon, trace = control$trace >
                                                                                                                     1), intercept = TRUE)$deviance
        else fit$deviance
        fit$null.deviance <- null.deviance
    }
    class(fit) <- c("negbin", "glm", "lm")
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    Call$init.theta <- signif(as.vector(th), 10)
    Call$link <- link
    fit$call <- Call
    if (model)
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x)
        fit$x <- X
    if (!y)
        fit$y <- NULL
    fit$theta <- as.vector(th)
    fit$SE.theta <- attr(th, "SE")
    fit$twologlik <- as.vector(2 * Lm)
    fit$aic <- -fit$twologlik + 2 * fit$rank + 2
    fit$contrasts <- attr(X, "contrasts")
    fit$xlevels <- .getXlevels(Terms, mf)
    fit$method <- method
    fit$control <- control
    fit$offset <- offset
    fit
}



################
## Legacy code
################


## ##############################
## ## chromunity
## ##############################
## #' @name chromunity
## #'
## #' @title Discovery of communities in Pore-C concatemers
## #' 
## #' @description This function takes in a GRanges with each row as Pore-C monomer with a metadata of the corresponding concatemer id 
## #' 
## #' @export
## #' @param this.pc.gr is the GRanges with each row as Pore-C monomer and must have a cmetadata column with corresponding concatemer annotation. The column should be called "read_idx"
## #' 
## #' @param k.knn numeric threshoold on k nearest concatemer neighbors to be used to create the graph. 
## #' 
## #' @param k.min numeric the threshold to number of concatemers pairs to be considered "similar"
## #' 
## #' @param tiles GRanges object dividing the genome in a fixed size bin to be defined by user
## #' 
## #' @param which.gr the GRanges for the window of interest
## #'
## #' @param filter_local_number boolean (default == FALSE) filters out concatemers in the window that fall below filter_local_thresh length
## #' 
## #' @param filter_local_thresh (default == NULL) numeric minimum length of concatemer to be considered for downstream analyses. To be set if filter_local_number == TRUE.
## #' @param take_sub_sample take a sub sample of concatemers. A random sample of fraction frac is taken
## #' 
## #' @param frac fraction of concatemer to be subsampled. To be set if take_sub_sample == TRUE.
## #'
## #' @param seed.n numeric set a seed when doing random subsampling
## #' 
## #' @return \code{chromunity} returns input GRanges with additional columns as follows:
## #' 
## #'    \item{community}{  numeric; \cr
## #'              the community annotation for each concatemer
## #'    }
## #'    \item{num.memb}{  numeric; \cr
## #'              number of members in each community
## #' 
## #' @author Aditya Deshpande


## chromunity <- function(this.pc.gr, k.knn = 10, k.min = 1, tiles, which.gr = which.gr, filter_local_number = FALSE, filter_local_thresh = NULL, take_sub_sample = FALSE, frac = 0.25, seed.n = 154){
    
##     reads = this.pc.gr 
##     if (filter_local_number){
##         message(paste0("Filtering out reads < ", filter_local_thresh))
##         reads = gr2dt(reads)
##         setkeyv(reads, c("seqnames", "start"))
##         reads[, max.local.dist := end[.N]-start[1], by = read_idx]
##         reads = reads[max.local.dist > filter_local_thresh]
##         reads = dt2gr(reads)
##     }
    
##     reads$tix = gr.match(reads, tiles)
##     reads = as.data.table(reads)[, count := .N, by = read_idx]
##     mat = dcast.data.table(reads[count > 2 ,]  %>% gr2dt, read_idx ~ tix, value.var = "strand", fill = 0)
##     mat2 = mat[, c(list(read_idx = read_idx), lapply(.SD, function(x) x >= 1)),.SDcols = names(mat)[-1]]
##     mat2 = suppressWarnings(mat2[, "NA" := NULL])
##     reads.ids = mat2$read_idx
##     mat2 = as.data.table(lapply(mat2, as.numeric))

##     if (take_sub_sample){
##         tot.num = nrow(mat2[rowSums(mat2[, -1]) > 1, ])
##         message(paste0("Total number of rows are: ", tot.num))
##         message("taking a subsample")
##         number.to.subsample = pmax(round(tot.num*frac), 1000)
##         message(paste0("Number sampled: ", number.to.subsample))
##         set.seed(seed.n)
##         gt = mat2[rowSums(mat2[, -1]) > 1, ][sample(.N, number.to.subsample), ] 
##     }

##     else {
##         gt = mat2[rowSums(mat2[, -1]) > 1, ]
##     }
    
##     ubx = gt$read_idx
##     message("Matrices made")
##     gc()

##     ## Prepare pairs for KNN
##     pairs = t(do.call(cbind, apply(gt[,setdiff(which(colSums(gt) > 1),1), with = FALSE] %>% as.matrix, 2, function(x) combn(which(x!=0),2))))
##     p1 = gt[pairs[,1], -1]
##     p2 = gt[pairs[,2], -1]
##     matching = rowSums(p1 & p2)
##     total = rowSums(p1 | p2)
##     dt = data.table(bx1 = pairs[,1], bx2 = pairs[,2], mat = matching, tot = total)[, frac := mat/tot]
##     dt2 = copy(dt)
##     dt2$bx2 = dt$bx1
##     dt2$bx1 = dt$bx2
##     dt3 = rbind(dt, dt2)
##     dt3$nmat = dt3$mat
##     dt3$nfrac = dt3$frac
##     setkeyv(dt3, c('nfrac', 'nmat'))
##     dt3 = unique(dt3)
##     dt3.2 = dt3[order(nfrac, nmat, decreasing = T)]
##     message("Pairs made")
##     gc()

##     ## Clustering    
##     k = k.knn
##     knn.dt = dt3.2[mat > 2 & tot > 2, .(knn = bx2[1:k]), by = bx1][!is.na(knn), ]
##     setkey(knn.dt)
##     knn = sparseMatrix(knn.dt$bx1, knn.dt$knn, x = 1)
##     knn.shared = knn %*% knn
##     message("KNN done")

##     ## community find in graph where each edge weight is # nearest neighbors (up to k) shared between the two nodes
##     KMIN = k.min
##     A = knn.shared*sign(knn.shared > KMIN)
##     A[cbind(1:nrow(A), 1:nrow(A))] = 0
##     A <- as(A, "matrix")
##     A <- as(A, "sparseMatrix")
##     A = A + t(A)
##     G = graph.adjacency(A, weighted = TRUE, mode = 'undirected')
##     cl.l = cluster_fast_greedy(G)
##     cl = cl.l$membership
##     message("Communities made")
##     memb.dt = data.table(read_idx = ubx[1:nrow(A)], community = cl)
##     reads = merge(reads, memb.dt, by = "read_idx")
##     reads[, num.memb := length(unique(read_idx)), by = community]
##     reads = dt2gr(reads)
##     return(reads)
## }


## ##############################
## ## annotate_multimodal_communities
## ##############################
## #' @name  annotate_multimodal_communities 
## #' 
## #'
## #' @title Annotates communities that are very dense with respect to genomic coordinates. 
## #' 
## #' @description Using the nature of distribution of contacts on genomic coordinates, annotates community that are local  
## #' 
## #' @export
## #' @param granges GRanges output from chromunity function
## #' 
## #' @param which.gr the GRanges for the window of interest
## #'
## #' @param min.memb numeric minimum number of members needed to be in a community to be considered for further analyses
## #' 
## #' @return \code{annotate_local_communities} returns input GRanges with additional columns as follows:
## #' 
## #'    \item{multimodal}{  boolean; \cr
## #'              whether a community had multimodal contacts based on parameters set by user
## #'    }
## #'  
## #' @author Aditya Deshpande

## annotate_multimodal_communities <- function(granges, which.gr, min.memb = 50){
##     this.dt = gr2dt(granges)[num.memb > min.memb]
##     ind.int = unique(this.dt$community)
##     dt.stat = data.table()
##     for (i in 1:length(ind.int)){
##         which.int = (ind.int)[i]  
##         message(which.int)
##         dt.int = this.dt[community == which.int]
##         dt.int.tmp = gr2dt(dt2gr(dt.int) %&&% which.gr)
##         setkeyv(dt.int.tmp, c("seqnames", "start"))
##         dt.sum = suppressWarnings(gr.sum(dt2gr(dt.int.tmp)+1e3))
##         y.max.val = max(dt.sum$score)
##         peak.gr = find_multi_modes( dt.sum[-1], w =  round(0.05*length(dt.sum[-1])), distance = 0.1*width(which.gr))
##         status = unique(peak.gr$bimodal)
##         dt.stat = rbind(dt.stat, data.table(community = which.int, multimodal = status))
##     }

##     this.dt = merge(this.dt, dt.stat, by = "community", allow.cartesian = T)
##     return(dt2gr(this.dt))
## }
        
        

## ##############################
## ## find_multi_modes
## ##############################
## #' @name  find_multi_modes 
## #' 
## #'
## #' @title Determines if a community has multimodal ditribution of contacts along genome
## #' 
## #' @description Using loess smoothing, determines if a community has multimodal ditribution of contacts along genome  
## #' 
## #' @export
## #' @param granges GRanges output from chromunity function for one community
## #' 
## #' @param x.field This is the X axis along which smoothing is done
## #' 
## #' @param y.field values to be smoothed
## #'
## #' @param w window size over which to smooth the distribution
## #' 
## #' @param distance numeric genomic distance beyond which a peak is not considered local 
## #' 
## #' @return \code{find_multi_modes} returns boolean value if more than one peak is present
## #'  
## #' @author Aditya Deshpande


## find_multi_modes <- function(granges, x.field = "start", y.field = "score", w = 1,  distance = 1e5) {
##     which.chr = seqlevels(granges)
##     x = x = gr2dt(granges)$start
##     y = values(granges)[, y.field]
##     n = length(y)
##     y.smooth <- loess(y ~ x, span = 0.1)$fitted
##     y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
##     delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
##     i.max <- which(delta <= 0) + w
##     max.gr <- granges[i.max]
##     if (any(gr.dist(max.gr) > distance)){
##         max.gr$bimodal = TRUE
##     } else {
##         max.gr$bimodal = FALSE
##     }
##     return(max.gr)
## }

## #' @name chromunity
## #' @description
## #'
## #' Runs genome-wide chromunity detection across a sliding or provided genomic window
## #'
## #' @param concatemers GRanges with $cid
## #' @param resolution bin size for community detection [5e4]
## #' @param region region to run on [si2gr(concatemers)]
## #' @param windows GRanges or GRangesList of windows to test, if empty will be generated by splitting region GRanges into window.size tiles with stride stride
## #' @param window.size window size to do community detection within
## #' @param tiles.k.knn KNN parameter specifying how many nearest neighbors to sample when building KNN graph
## #' @param peak.thresh peak threshold with which to call a peak
## #' @param k.min minimal number of nearest neighbors an edge in KNN graph needs to have before community detection
## #' @param pad integer pad to use when computing the footprint of each chromunity and finding peak regions which become binsets
## #' @return list with items $binset,  $support, $params: $binsets is GRanges of bins with field $bid corresponding to binset id and $support which is the concatemer community supporting the binset which are GRanges with $bid
## #' @author Aditya Deshpande, Marcin Imielinski
## #' @export
## chromunity = function(concatemers, resolution = 5e4, region = si2gr(concatemers), windows = NULL, window.size = 2e6, max.slice = 1e6, min.support = 5, stride = window.size/2, mc.cores = 20, k.knn = 25, k.min = 5, pad = 1e3, peak.thresh = 0.85, seed = 42, verbose = TRUE)
## {
##   if (is.null(windows))
##     windows = gr.start(gr.tile(region, stride))+window.size/2

##   if (inherits(windows, 'GRanges'))
##     windows = split(windows, 1:length(windows))

##   if (is.null(concatemers$cid))
##   {
##       if ('read_idx' %in% names(values(concatemers)))
##         names(values(concatemers))[match('read_idx', names(values(concatemers)))] = 'cid'
##       else
##         stop("concatemer GRanges must have metadata column $cid or $read_idx specifying the concatemer id")
##   }
  
##   params = data.table(k.knn = k.knn, k.min = k.min, seed = seed)

##   bins = gr.tile(reduce(gr.stripstrand(unlist(windows))), 5e4)[, c()]

##   if (verbose) cmessage('Generated ', length(bins), ' bins across ', length(windows), ' windows')

##   if (verbose) cmessage('Matching concatemers with bins, and bins with windows using gr.match with max.slice ', max.slice, ' and ', mc.cores, ' cores')

##   ## (batch) match up concatemers with binids
##   concatemers$binid = gr.match(concatemers, bins, max.slice = max.slice, mc.cores =  mc.cores, verbose = verbose)

##   ## match window ids and bins 
##   binmap = bins %*% grl.unlist(windows)[, c('grl.ix')] %>% as.data.table %>% setnames('query.id', 'binid') %>% setnames('grl.ix', 'winid') %>% setkeyv('winid')

##   ## cycle through (possibly complex) windows call cluster_concatemers and convert to gr.sums
##   winids = unique(binmap$winid)

##   if (verbose) cmessage('Starting concatemer community detection across ', length(winids), ' windows')

##   cc = pbmclapply(winids, mc.cores = mc.cores, function(win)
##   {
##     suppressWarnings({
##       these.bins = binmap[.(win), ]
##       cc = concatemer_communities(concatemers %Q% (binid %in% these.bins$binid), k.knn = k.knn, k.min = k.min, seed = seed, verbose = verbose>1)
##       if (length(cc))
##       {
##         cc = cc[cc$support >= min.support]
##         cc$winid = win
##       }
##     })
##     cc
##   }) %>% do.call(grbind, .)

##   if (!length(cc))
##     return(Chromunity(concatemers = cc, chromunities = GRanges(), params = params))

##   uchid = unique(cc$chid)

##   if (verbose) cmessage('Analyzing gr.sums associated with ', length(uchid), ' concatemer communities to generate binsets')

##   binsets = pbmclapply(uchid, mc.cores = mc.cores, function(this.chid)
##   {
##     suppressWarnings({
##       this.cc = cc %Q% (chid == this.chid)
##       peaks = gr.sum(this.cc + pad) %>% gr.peaks('score')
##       binset = bins[, c()] %&% (peaks[peaks$score > quantile(peaks$score, peak.thresh)])
##       if (length(binset))
##       {
##         binset$chid = this.chid
##         binset$winid = this.cc$winid[1]
##       }
##     })
##     binset
##   })  %>% do.call(grbind, .)

##   return(Chromunity(concatemers = cc[cc$chid %in% binsets$chid], binsets = binsets, meta = params))
## }


## #' @name concatemer_communities
## #' @description
## #'
## #' Low level function that labels concatemers with chromunity ids $chid using community detection on a graph. 
## #'
## #' Given a GRanges of monomers labeled by concatemer id $cid
## #'
## #' @param concatemers GRanges of monomers with field $cid indicating concatemer id and $binid represent bin id
## #' @param tiles.k.knn KNN parameter specifying how many nearest neighbors to sample when building KNN graph
## #' @param k.min minimal number of nearest neighbors an edge in KNN graph needs to have before community detection
## #' @param drop.small logical flag specifying whether to remove "small" concatemers ie those with a footprint <= small argument [FALSE]
## #' @param small integer threshold for bases that define small concatemers, only relevant if drop.small = TRUE
## #' @param subsample.frac optional arg specifying fraction of concatemers to subsample [NULL]
## #' @param seed seed for subsampling
## #' @return GRanges of concatemers labeled by $c mmunity which specifies community id
## concatemer_communities = function (concatemers, k.knn = 25, k.min = 5, 
##     drop.small = FALSE,small = 1e4, 
##     subsample.frac = NULL, seed = 42, verbose = TRUE) 
## {
##   reads = concatemers

##   if (is.null(reads$cid))
##   {
##     if ('read_idx' %in% names(values(reads)))
##       names(values(reads))[match('read_idx', names(values(reads)))] = 'cid'
##     else
##       stop("concatemer GRanges must have metadata column $cid or $read_idx specifying the concatemer id")
##   }
  
##   if (drop.small) {
##     if (verbose) cmessage(paste0("Filtering out reads < ", small))
##     reads = gr2dt(reads)
##     setkeyv(reads, c("seqnames", "start"))
##     reads[, `:=`(max.local.dist, end[.N] - start[1]), by = cid]
##     reads = reads[max.local.dist > small]
##     reads = dt2gr(reads)
##   }

##   if (verbose) cmessage("Matching reads to tiles")
##   reads = as.data.table(reads)[, `:=`(count, .N), by = cid]
##   mat = dcast.data.table(reads[count > 2, ] %>% gr2dt, cid ~  binid, value.var = "strand", fun.aggregate = length, fill = 0)                                                        
##   mat2 = mat[, c(list(cid = cid), lapply(.SD, function(x) x >= 1)), .SDcols = names(mat)[-1]]                                                          
##   mat2 = suppressWarnings(mat2[, `:=`("NA", NULL)])
##   reads.ids = mat2$cid
##   mat2 = as.data.table(lapply(mat2, as.numeric))
##   if (!is.null(subsample.frac)) {
##     if (verbose) cmessage("Subsampling concatemers")
##     tot.num = nrow(mat2[rowSums(mat2[, -1]) > 1, ])
##     if (verbose) cmessage(paste0("Total number of rows are: ", tot.num))
##     if (verbose) cmessage("taking a subsample")
##     number.to.subsample = pmax(round(tot.num * subsample.frac), 1000)
##     if (verbose) cmessage(paste0("Number sampled: ", number.to.subsample))
##     set.seed(seed)
##     concatm = mat2[rowSums(mat2[, -1]) > 1, ][sample(.N, min(c(.N, number.to.subsample))), ]
##   }
##   else {
##     concatm = mat2[rowSums(mat2[, -1]) > 1, ]
##   }
##   ubx = concatm$cid

##   if (verbose) cmessage("Matrices made")
##   gc()

##   concatm = concatm[, setdiff(which(colSums(concatm) > 1), 1), with = FALSE]

##   if (!ncol(concatm))
##   {
##     warning('No concatemers found hitting two bins, returning empty result')
##     return(reads[, chid := NA][c(), ])
##   }

##   pairs = t(do.call(cbind, apply(concatm[, setdiff(which(colSums(concatm) > 1), 1), with = FALSE] %>% as.matrix, 2, function(x) combn(which(x != 0), 2))))
                                                    
##   concatm = as(as.matrix(as.data.frame(concatm)), "sparseMatrix")    
##   p1 = concatm[pairs[, 1], -1]
##   p2 = concatm[pairs[, 2], -1]
##   matching = Matrix::rowSums(p1 & p2)
##   total = Matrix::rowSums(p1 | p2)
##   dt = data.table(bx1 = pairs[, 1], bx2 = pairs[, 2], mat = matching, 
##                   tot = total)[, `:=`(frac, mat/tot)]
##   dt2 = copy(dt)
##   dt2$bx2 = dt$bx1
##   dt2$bx1 = dt$bx2
##   dt3 = rbind(dt, dt2)
##   dt3$nmat = dt3$mat
##   dt3$nfrac = dt3$frac
##   setkeyv(dt3, c("nfrac", "nmat"))
##   dt3 = unique(dt3)
##   dt3.2 = dt3[order(nfrac, nmat, decreasing = T)]
##   if (verbose) cmessage("Pairs made")
##   gc()
##   k = k.knn
##   knn.dt = dt3.2[mat > 2 & tot > 2, .(knn = bx2[1:k]), by = bx1][!is.na(knn), ]
##   setkey(knn.dt)
##   knn = sparseMatrix(knn.dt$bx1, knn.dt$knn, x = 1)
##   knn.shared = knn %*% knn
##   if (verbose) cmessage("KNN done")
##   KMIN = k.min
##   A = knn.shared * sign(knn.shared > KMIN)
##   A[cbind(1:nrow(A), 1:nrow(A))] = 0
##   A <- as(A, "matrix")
##   A <- as(A, "sparseMatrix")
##   A = A + t(A)
##   G = graph.adjacency(A, weighted = TRUE, mode = "undirected")
##   cl.l = cluster_fast_greedy(G)
##   cl = cl.l$membership
##   if (verbose) cmessage("Communities made")
##   memb.dt = data.table(cid = ubx[1:nrow(A)], chid = cl)
##   reads = merge(reads, memb.dt, by = "cid")
##   reads[, `:=`(support, length(unique(cid))), by = chid]
##   reads = dt2gr(reads)
##   return(reads)
## }

## #' @name smessage
## #' @description
## #'
## #' @author Aditya Deshpande, Marcin Imielinski
## #' @export
## #' @private 
## smessage = function(..., pre = 'Synergy')
##   message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)

## #' @name cmessage
## #' @description
## #'
## #' @author Aditya Deshpande, Marcin Imielinski
## #' @export
## #' @private 
## cmessage = function(..., pre = 'Chromunity')
##   message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)


