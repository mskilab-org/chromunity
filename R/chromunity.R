#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table set
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table setkeyv
#' @importFrom data.table setnames
#' @importFrom data.table transpose
#' @importFrom Matrix sparseMatrix 
#' @import zoo
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom MASS ginv
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#' @import GenomicRanges
#' @import gUtils
#' @importFrom igraph graph.adjacency
#' @importFrom igraph cluster_louvain
#' @importFrom BiocGenerics t

globalVariables(c("::", ":::", "num.memb", "community", "max.local.dist", "read_idx", "as.data.table", "count", "dcast.data.table", "combn", "tot", "copy", "nfrac", "nmat", ".", "bx2", "bx1", "as", "seqlevels", "loess", ".N", ".SD", ":="))



##############################
## chromunity
##############################
#' @name chromunity
#'
#' @title Discovery of communities in Pore-C concatemers
#' 
#' @description This function takes in a GRanges with each row as Pore-C monomer with a metadata of the corresponding concatemer id 
#' 
#' @export
#' @param this.pc.gr is the GRanges with each row as Pore-C monomer and must have a cmetadata column with corresponding concatemer annotation. The column should be called "read_idx"
#' 
#' @param k.knn numeric threshoold on k nearest concatemer neighbors to be used to create the graph. 
#' 
#' @param k.min numeric the threshold to number of concatemers pairs to be considered "similar"
#' 
#' @param tiles GRanges object dividing the genome in a fixed size bin to be defined by user
#' 
#' @param which.gr the GRanges for the window of interest
#'
#' @param filter_local_number boolean (default == FALSE) filters out concatemers in the window that fall below filter_local_thresh length
#' 
#' @param filter_local_thresh (default == NULL) numeric minimum length of concatemer to be considered for downstream analyses. To be set if filter_local_number == TRUE.
#' @param take_sub_sample take a sub sample of concatemers. A random sample of fraction frac is taken
#' 
#' @param frac fraction of concatemer to be subsampled. To be set if take_sub_sample == TRUE.
#'
#' @param seed.n numeric set a seed when doing random subsampling
#' 
#' @return \code{chromunity} returns input GRanges with additional columns as follows:
#' 
#'    \item{community}{  numeric; \cr
#'              the community annotation for each concatemer
#'    }
#'    \item{num.memb}{  numeric; \cr
#'              number of members in each community
#' 
#' @author Aditya Deshpande


chromunity <- function(this.pc.gr, k.knn = 10, k.min = 1, tiles, which.gr = which.gr, filter_local_number = FALSE, filter_local_thresh = NULL, take_sub_sample = FALSE, frac = 0.25, seed.n = 154){
    
    reads = this.pc.gr 
    if (filter_local_number){
        message(paste0("Filtering out reads < ", filter_local_thresh))
        reads = gr2dt(reads)
        setkeyv(reads, c("seqnames", "start"))
        reads[, max.local.dist := end[.N]-start[1], by = read_idx]
        reads = reads[max.local.dist > filter_local_thresh]
        reads = dt2gr(reads)
    }
    
    reads$tix = gr.match(reads, tiles)
    reads = as.data.table(reads)[, count := .N, by = read_idx]
    mat = dcast.data.table(reads[count > 2 ,]  %>% gr2dt, read_idx ~ tix, value.var = "strand", fill = 0)
    mat2 = mat[, c(list(read_idx = read_idx), lapply(.SD, function(x) x >= 1)),.SDcols = names(mat)[-1]]
    mat2 = suppressWarnings(mat2[, "NA" := NULL])
    reads.ids = mat2$read_idx
    mat2 = as.data.table(lapply(mat2, as.numeric))

    if (take_sub_sample){
        tot.num = nrow(mat2[rowSums(mat2[, -1]) > 1, ])
        message(paste0("Total number of rows are: ", tot.num))
        message("taking a subsample")
        number.to.subsample = pmax(round(tot.num*frac), 1000)
        message(paste0("Number sampled: ", number.to.subsample))
        set.seed(seed.n)
        gt = mat2[rowSums(mat2[, -1]) > 1, ][sample(.N, number.to.subsample), ] 
    }

    else {
        gt = mat2[rowSums(mat2[, -1]) > 1, ]
    }
    
    ubx = gt$read_idx
    message("Matrices made")
    gc()

    ## Prepare pairs for KNN
    pairs = t(do.call(cbind, apply(gt[,setdiff(which(colSums(gt) > 1),1), with = FALSE] %>% as.matrix, 2, function(x) combn(which(x!=0),2))))
    p1 = gt[pairs[,1], -1]
    p2 = gt[pairs[,2], -1]
    matching = rowSums(p1 & p2)
    total = rowSums(p1 | p2)
    dt = data.table(bx1 = pairs[,1], bx2 = pairs[,2], mat = matching, tot = total)[, frac := mat/tot]
    dt2 = copy(dt)
    dt2$bx2 = dt$bx1
    dt2$bx1 = dt$bx2
    dt3 = rbind(dt, dt2)
    dt3$nmat = dt3$mat
    dt3$nfrac = dt3$frac
    setkeyv(dt3, c('nfrac', 'nmat'))
    dt3 = unique(dt3)
    dt3.2 = dt3[order(nfrac, nmat, decreasing = T)]
    message("Pairs made")
    gc()

    ## Clustering    
    k = k.knn
    knn.dt = dt3.2[mat > 2 & tot > 2, .(knn = bx2[1:k]), by = bx1][!is.na(knn), ]
    setkey(knn.dt)
    knn = sparseMatrix(knn.dt$bx1, knn.dt$knn, x = 1)
    knn.shared = knn %*% knn
    message("KNN done")

    ## community find in graph where each edge weight is # nearest neighbors (up to k) shared between the two nodes
    KMIN = k.min
    A = knn.shared*sign(knn.shared > KMIN)
    A[cbind(1:nrow(A), 1:nrow(A))] = 0
    A <- as(A, "matrix")
    A <- as(A, "sparseMatrix")
    A = A + t(A)
    G = graph.adjacency(A, weighted = TRUE, mode = 'undirected')
    cl.l = cluster_fast_greedy(G)
    cl = cl.l$membership
    message("Communities made")
    memb.dt = data.table(read_idx = ubx[1:nrow(A)], community = cl)
    reads = merge(reads, memb.dt, by = "read_idx")
    reads[, num.memb := length(unique(read_idx)), by = community]
    reads = dt2gr(reads)
    return(reads)
}


##############################
## annotate_multimodal_communities
##############################
#' @name  annotate_multimodal_communities 
#' 
#'
#' @title Annotates communities that are very dense with respect to genomic coordinates. 
#' 
#' @description Using the nature of distribution of contacts on genomic coordinates, annotates community that are local  
#' 
#' @export
#' @param granges GRanges output from chromunity function
#' 
#' @param which.gr the GRanges for the window of interest
#'
#' @param min.memb numeric minimum number of members needed to be in a community to be considered for further analyses
#' 
#' @return \code{annotate_local_communities} returns input GRanges with additional columns as follows:
#' 
#'    \item{multimodal}{  boolean; \cr
#'              whether a community had multimodal contacts based on parameters set by user
#'    }
#'  
#' @author Aditya Deshpande

annotate_multimodal_communities <- function(granges, which.gr, min.memb = 50){
    this.dt = gr2dt(granges)[num.memb > min.memb]
    ind.int = unique(this.dt$community)
    dt.stat = data.table()
    for (i in 1:length(ind.int)){
        which.int = (ind.int)[i]  
        message(which.int)
        dt.int = this.dt[community == which.int]
        dt.int.tmp = gr2dt(dt2gr(dt.int) %&&% which.gr)
        setkeyv(dt.int.tmp, c("seqnames", "start"))
        dt.sum = suppressWarnings(gr.sum(dt2gr(dt.int.tmp)+1e3))
        y.max.val = max(dt.sum$score)
        peak.gr = find_multi_modes( dt.sum[-1], w =  round(0.05*length(dt.sum[-1])), distance = 0.1*width(which.gr))
        status = unique(peak.gr$bimodal)
        dt.stat = rbind(dt.stat, data.table(community = which.int, multimodal = status))
    }

    this.dt = merge(this.dt, dt.stat, by = "community", allow.cartesian = T)
    return(dt2gr(this.dt))
}
        
        


##############################
## find_multi_modes
##############################
#' @name  find_multi_modes 
#' 
#'
#' @title Determines if a community has multimodal ditribution of contacts along genome
#' 
#' @description Using loess smoothing, determines if a community has multimodal ditribution of contacts along genome  
#' 
#' @export
#' @param granges GRanges output from chromunity function for one community
#' 
#' @param x.field This is the X axis along which smoothing is done
#' 
#' @param y.field values to be smoothed
#'
#' @param w window size over which to smooth the distribution
#' 
#' @param distance numeric genomic distance beyond which a peak is not considered local 
#' 
#' @return \code{find_multi_modes} returns boolean value if more than one peak is present
#'  
#' @author Aditya Deshpande


find_multi_modes <- function(granges, x.field = "start", y.field = "score", w = 1,  distance = 1e5) {
    which.chr = seqlevels(granges)
    x = x = gr2dt(granges)$start
    y = values(granges)[, y.field]
    n <- length(y)
    y.smooth <- loess(y ~ x, span = 0.1)$fitted
    y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
    delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
    i.max <- which(delta <= 0) + w
    max.gr <- granges[i.max]
    if (any(gr.dist(max.gr) > distance)){
        max.gr$bimodal = TRUE
    } else {
        max.gr$bimodal = FALSE
    }
    return(max.gr)
}




#' @name background
#' @description
#' Given n binsets generates random "background" binsets that mirrors the input binset characteristics with respect to chromosome, width, and distance.
#'
#' Generates (chromosome specific) kernels for width, cardinality, distance and models interchromosomal binsets by allowing
#' a chromosome transition matrix which then lands you on a coordinate that is also chosen from a kernel. 
#' x
#' @param binsets GRanges of bins with fields seqnames, start, end, and $bid specifying binset id
#' @param n integer scalar specifying number of sets to sample (default = nrow(binsets))
#' @export
background = function(binsets, n = nrow(binsets), pseudocount = 1, verbose = TRUE, mc.cores = 20)
{
  if (!length(binsets))
    stop('empty binsets')

  ## sort from left to right
  binsets = gr2dt(binsets)
  setkeyv(binsets, c('seqnames','start', 'end'))

  if (verbose) smessage('Making kernels')

  ## populate data for kernels
  ends.kernel = binsets[, .(val = max(end)), by = .(seqnames, bid)] %>% setkey('seqnames')
  distance.kernel = binsets[, .(val = start-shift(end)), by = .(seqnames, bid)][!is.na(val), ] %>% setkey('seqnames')
  width.kernel = binsets[, .(seqnames, val = width)] %>% setkey('seqnames')

  ends.bw =  ends.kernel[, .(bw = density(val)$bw), keyby = seqnames]
  distance.bw =  distance.kernel[, .(bw = density(val)$bw), keyby = seqnames]
  width.bw =  width.kernel[, .(bw = density(val)$bw), keyby = seqnames]


  if (verbose) smessage('Making transition matrices')

  ## populate cardinality 
  cardinality.multinomial = binsets[, .(cardinality = .N), by = .(bid)][, .N, by = cardinality][, .(cardinality, prob = N/sum(N))] %>% setkey(cardinality)

  ## populate chrom multi-nomial and chr-chr transition matrix
  chrom.multinomial = binsets[, .(count = .N + pseudocount), by = seqnames][, .(seqnames, prob = count / sum(count))] %>% setkey('seqnames')

  ## populate chrom multi-nomial and chr-chr transition matrix
  chrom.transition = merge(binsets, binsets, by = 'bid', allow.cartesian = TRUE)[, .(count = .N + pseudocount), by = .(seqnames.x, seqnames.y)][, .(seqnames.y, prob = count / sum(count)), keyby = .(seqnames.x)] 

  if (verbose) smessage('Generating random sets')
  out = pbmclapply(1:n, mc.cores = mc.cores, function(i)
  {
    cardinality = cardinality.multinomial[, sample(cardinality, 1, prob = prob)]
    binset = data.table(bid = i, seqnames = chrom.multinomial[, sample(seqnames, prob = prob, 1)])
    binset$end = rdens(1, ends.kernel[.(binset$seqnames), val], width = ends.bw[.(binset$seqnames), bw])
    binset$start = binset$end - rdens(1, width.kernel[.(binset$seqnames), val], width = width.bw[.(binset$seqnames), bw])
    for (k in setdiff(1:cardinality, 1))
    {
      lastchrom = tail(binset$seqnames, 1)
      laststart = tail(binset$start, 1)
      newchrom = chrom.transition[.(lastchrom), sample(seqnames.y, 1, prob = prob)]
      if (newchrom == lastchrom)
      {
        newbin = data.table(bid = i,
                            seqnames = newchrom,
                            end = laststart - rdens(1, distance.kernel[.(newchrom), val], width = distance.bw[.(newchrom), bw]))
        newbin[, start := end - rdens(1, width.kernel[.(newchrom), val], width = width.bw[.(newchrom), bw])]
      }
      else
      {
        newbin = data.table(bid = i,
                            seqnames = newchrom,
                            end = rdens(1, ends.kernel[.(newchrom), val], width = ends.bw[.(newchrom), bw]))
        newbin[, start := end - rdens(1, width.kernel[.(newchrom), val], width = width.bw[.(newchrom), bw])]
      }
      binset = rbind(binset, newbin, fill = TRUE)
    }
    binset
  }) %>% rbindlist(fill = TRUE)

  ## fix coordinates
  out[, start := pmax(1, ceiling(start))]
  out[, end := pmax(1, ceiling(end))]
  return(out)
}



#' @name rdens
#' @description
#'
#' Function to sample kernel density from a distribution adapted from
#' https://stats.stackexchange.com/questions/321542/how-can-i-draw-a-value-randomly-from-a-kernel-density-estimate
#' 
#' @param n integer number of samples
#' @param density
rdens = function(n, values, width = NULL, kernel="gaussian") {
  if (is.null(width)) width = density(values)$bw
  rkernel <- function(n) rnorm(n, sd=width)  # Kernel sampler
  ret = data.table(value = sample(values, n, replace=TRUE) + rkernel(n))
  ret[, value := ifelse(value < 0, value*-1, value)]
  return(ret$value)
}



#' @name annotate
#' @description
#'
#' Given n binsets annotates their k-power set (i.e. all sets in the power set with cardinality <= k) with
#' basic (min median max distance, cardinality, width) and optional numerical and interval covariates provided as input.
#' 
#' @param binsets GRanges of bins with fields $bsid specifying binset id
#' @param concatemers GRanges of monomers with fields seqnames, start, end, and $cid specifying concatemer id, which will be counted across each binset
#' @param covariates Covariate object specifying numeric or interval covariates to merge with this binset
#' @param k the max cardinality of sets to annotate from the power set of each binmer
#' @param gg optional gGraph input specifying alternate distance function
#' @param interchromosomal.dist numeric scalar of "effective" distance for inter chromosomal bins [1e8]
#' @param verbose logical flag
#' @param mc.cores integer how many cores to parallelize
#' @export
#' @return data.table of sub-binsets i.e. k-power set of binsets annotated with $count field representing covariates, ready for fitting, **one row per binset
annotate = function(binsets, concatemers, covariates = NULL, k = 5, interchromosomal.dist = 1e8, gg = NULL, mc.cores = 20, numchunks = 2*mc.cores-1, seed = 42, verbose = TRUE)
{
  set.seed(seed)
  if (!inherits(binsets, 'GRanges'))
    binsets = dt2gr(binsets)

  binsets = gr.stripstrand(binsets)

  ## each row of binsets is a bin
  binsets$binid = 1:length(binsets)

  #####
  ## frontload / batch some useful calculations
  #####

  ## bin vs concatemer overlaps
  if (verbose) smessage('Overlapping ', length(binsets), ' bins with ', length(concatemers), ' monomers')
  ov = binsets %*% concatemers[, 'cid'] %>% as.data.table

  if (verbose) smessage('Computing bin by bin pairwise distance')
  ## bin x bin pairwise distance within each set
  if (is.null(gg))
  {
    bindist = gr2dt(binsets)[, as.data.table(expand.grid(i = binid, j = binid))[i<j, ], by = bid] %>% setkeyv(c('i', 'j'))
    bindist[, distance := GenomicRanges::distance(binsets[i], binsets[j])]
  }
  else
  {
    if (verbose) smessage('Using graph distance')
    gg = gg$copy$disjoin(disjoin(binsets))
    
    binsetd = data.table(binid = binsets$binid, bid = binsets$bid) %>% setkey('bid')
    binsetd[, gid := gr.match(binsets, gg$nodes)] ## will streamline distance computation to get index
    bindist = pbmclapply(unique(bindsetd$bid), function(bid)
      gg$dist(binsetl[.(bid), gid]) %>% melt %>% as.data.table %>% setnames('i', 'j', 'distance')[i<j, ][, bid := bid], mc.cores = mc.cores) %>% rbindlist      
  }
  bindist[is.infinite(distance), distance := interchromosomal.dist]
  
  #####
  ## now start annotating sub.binsets aka sets
  #####

  if (verbose) smessage('Making sub binsets')
  ## make all sub power sets up to k for all bids
  ## each row is now a binid with a setid and bid
  sub.binsets = gr2dt(binsets)[, powerset(binid, 2, k), by = bid] %>% setnames(c('bid', 'setid', 'binid')) %>% setkey(bid)
  sub.binsets[, ":="(iid = 1:.N, tot = .N), by = .(setid, bid)] ## label each item in each sub-binset, and total count will be useful below

  if (verbose) smessage('Made ', nrow(sub.binsets), ' sub-binsets')

  ## first use ov to count how many concatemers fully overlap all the bins in the subbinset
  if (verbose) smessage('Counting concatemers across sub-binsets across ', mc.cores, ' cores')
  counts = unique(sub.binsets[, .(bid, setid)])[, count := 0]
  if (nrow(ov))
    {
      ubid = unique(sub.binsets$bid) ## split up to lists to leverage pbmclapply
      ubidl = split(ubid, ceiling(runif(length(ubid))*numchunks)) ## randomly chop up ubid into twice the number of mc.coreso
      counts = pbmclapply(ubidl, mc.cores = mc.cores, function(bids)
      {
        out = merge(sub.binsets[.(bids), ], ov, by = c('binid', 'bid'), allow.cartesian = TRUE)
        if (nrow(out))
          out[, .(hit = all(1:tot[1] %in% iid)), by = .(cid, setid, bid)][, .(count = sum(hit, na.rm = TRUE)), by = .(setid, bid)]
        else
          NULL
      })  %>% rbindlist
    }
   
  ## other goodies
  if (verbose) smessage('Computing min median max distances per setid')
  dists = sub.binsets[, bindist[as.data.table(expand.grid(i = binid, j = binid))[i<j, ], .(dist = c('min.dist', 'median.dist', 'max.dist'), value = quantile(distance+1,  c(0, 0.5, 1)))], by = .(setid, bid)] %>% dcast(bid + setid ~ dist, value.var = 'value')

  if (verbose) smessage('Computing total width and cardinality per setid')
  widths = sub.binsets[, .(width = sum(width(binsets)[binid])+1, cardinality = .N), by = .(bid, setid)]

  if (verbose) smessage('Merging counts, distance, and width')
  annotated.binsets = unique(sub.binsets[, .(bid, setid)]) %>%
    merge(counts, by = c('bid', 'setid')) %>%
    merge(dists, all.x = TRUE, by = c('bid', 'setid')) %>%
    merge(widths, all.x = TRUE, by = c('bid', 'setid'))

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
 
  covariates = setdiff(names(annotated.binsets), c('bid', 'setid', 'width', 'count'))
  fmstring = paste('count ~', paste(paste0('log(', covariates, ')', collapse = ' + ')))
  fmstring = paste0(fmstring, " + ", "offset(log(width))")
  fm = formula(fmstring)

  model = tryCatch(glm.nb(formula = fm, data = annotated.binsets, control = glm.control(maxit = 50)), error = function(e) NULL)

  return(list(model = model, covariates = covariates))
}

#' @name score
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
score = function(annotated.binsets, model, verbose = TRUE)
{
  if (is.null(annotated.binsets$count))
    stop('annotated.binsets need to have $count column, did you forget to annotate?')

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
#' @param mc.cores threads to parallelize
#' @param verbose logical flag whether to fit
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
synergy = function(binsets, concatemers, background.binsets = NULL, model = NULL, covariates = NULL, annotated.binsets = NULL, k = 5, gg = NULL, mc.cores = 20, verbose = TRUE)
{
  if (is.null(background.binsets) & is.null(model))
  {
    if (verbose) fmessage('Computing random background binsets using features of provided binsets')
    background.binsets = background(binsets, mc.cores = mc.cores)
  }
  
  if (is.null(model))
  {
    if (verbose) fmessage('Annotating k-power sets of background binsets using features and covariates')
    annotated.backround.binsets = annotate(binsets = background.binsets, concatemers = concatemers, covariates = covariates, k = k, mc.cores = mc.cores, verbose = verbose)

    if (verbose) fmessage('Fitting model to k-power sets of annotated background binsets')
    model = fit(annotated.background.binsets, nb = TRUE, return.model = TRUE, verbose = verbose)
  }

  if (is.null(annotated.binsets))
    annotated.binsets = annotate(binsets = binsets, concatemers = concatemers, covariates = covariates, k = k, mc.cores = mc.cores, verbose = verbose)

  scored.binsets = score(annotated.binsets, model, verbose = verbose)

  scored.binsets[, multiway := cardinality > 2]
  setkey(scored.binsets, bid)
  ubid = unique(scored.binsets$bid)

  res = pbmclapply(ubid, function(this.bid) muffle(dflm(glm.nb(data = scored.binsets[.(this.bid),], count ~ offset(log(count.predicted))))[1, name := this.bid])) %>% rbindlist
  setnames(res, 'name', 'bid')

  return(res)
}

#' @name chromunity
#' @description
#'
#' Runs genome-wide chromunity detection across a sliding or provided genomic window
#'
#' @param concatemers GRanges with $cid
#' @param resolution bin size for community detection [5e4]
#' @param region region to run on [si2gr(concatemers)]
#' @param windows GRanges or GRangesList of windows to test, if empty will be generated by splitting region GRanges into window.size tiles with stride stride
#' @param window.size window size to do community detection within
#' @param tiles.k.knn KNN parameter specifying how many nearest neighbors to sample when building KNN graph
#' @param peak.thresh peak threshold with which to call a peak
#' @param k.min minimal number of nearest neighbors an edge in KNN graph needs to have before community detection
#' @param pad integer pad to use when computing the footprint of each chromunity and finding peak regions which become binsets
#' @return list with items $binset,  $support, $params: $binsets is GRanges of bins with field $bid corresponding to binset id and $support which is the concatemer community supporting the binset which are GRanges with $bid
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
chromunity = function(concatemers, resolution = 5e4, region = seqinfo2gr(concatemers), windows = NULL, window.size = 2e6, max.slice = 1e6, min.support = 5, stride = window.size/2, mc.cores = 20, k.knn = 25, k.min = 5, pad = 1e3, peak.thresh = 0.85, seed = 42, verbose = TRUE)
{
  if (is.null(windows))
    windows = gr.start(gr.tile(region, stride))+window.size/2

  if (inherits(windows, 'GRanges'))
    windows = split(windows, 1:length(windows))

  params = list(k.knn = k.knn, k.min = k.min, seed = seed)

  bins = gr.tile(reduce(gr.stripstrand(windows)), 5e4)[, c()]

  if (verbose) cmessage('Generated ', length(bins), ' bins across ', length(windows), ' windows')

  if (verbose) cmessage('Matching concatemers with bins, and bins with windows using gr.match with max.slice ', max.slice, ' and ', mc.cores, ' cores')

  ## (batch) match up concatemers with binids
  concatemers$binid = gr.match(concatemers, bins, max.slice = max.slice, mc.cores =  mc.cores, verbose = verbose)

  ## match window ids and bins 
  binmap = bins %*% grl.unlist(windows)[, c('grl.ix')] %>% as.data.table %>% setnames('query.id', 'binid') %>% setnames('grl.ix', 'winid') %>% setkeyv('winid')

  ## cycle through (possibly complex) windows call cluster_concatemers and convert to gr.sums
  winids = unique(binmap$winid)

  if (verbose) cmessage('Starting concatemer community detection across ', length(winids), ' windows')

  cc = pbmclapply(winids, mc.cores = mc.cores, function(win)
  {
    these.bins = binmap[.(win), ]
    cc = concatemer_communities(concatemers %Q% (binid %in% these.bins$binid), k.knn = k.knn, k.min = k.min, seed = seed, verbose = verbose>1)
    if (length(cc))
      {
        cc = cc %Q% (support >= min.support)
        cc$winid = win
      }
    cc
  }) %>% do.call(grbind, .)

  if (!length(cc))
    return(list(concatemers = cc, chromunities = GRanges(), params = params))

  ubid = unique(cc$bid)

  if (verbose) cmessage('Analyzing gr.sums associated with ', length(ubid), ' concatemer communities to generate binsets')

  binsets = pbmclapply(ubid, mc.cores = mc.cores, function(this.bid)
  {
    this.cc = cc %Q% (bid == this.bid)
    peaks = gr.sum(this.cc + pad) %>% gr.peaks('score')
    binset = bins[, c()] %&% (peaks %Q% (score > quantile(footprint$score, peak.thresh)))
    if (length(binset))
      {
        binset$bid = this.bid
        binset$winid = this.cc$winid[1]
      }
    binset
  })  %>% do.call(grbind, .)

  return(list(concatemers = cc %Q% (bid %in% binset$bid), binsets = binsets, params = params))
}


#' @name concatemer_communities
#' @description
#'
#' Low level function that labels concatemers with community ids $bid using community detection on a graph. 
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
#' @return GRanges of concatemers labeled by $c mmunity which specifies community id
concatemer_communities = function (concatemers, k.knn = 25, k.min = 5, 
    drop.small = FALSE,small = 1e4, 
    subsample.frac = NULL, seed = 42, verbose = TRUE) 
{
  reads = concatemers

  if (is.null(reads$cid))
  {
    if ('read_idx' %in% names(values(reads)))
      names(values)[match('read_idx', names(values(reads)))] = 'cid'
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

  if (verbose) cmessage("Matching reads to tiles")
  reads = as.data.table(reads)[, `:=`(count, .N), by = cid]
  mat = dcast.data.table(reads[count > 2, ] %>% gr2dt, cid ~  binid, value.var = "strand", fun.aggregate = length, fill = 0)                                                        
  mat2 = mat[, c(list(cid = cid), lapply(.SD, function(x) x >= 1)), .SDcols = names(mat)[-1]]                                                          
  mat2 = suppressWarnings(mat2[, `:=`("NA", NULL)])
  reads.ids = mat2$cid
  mat2 = as.data.table(lapply(mat2, as.numeric))
  if (!is.null(subsample.frac)) {
    if (verbose) cmessage("Subsampling concatemers")
    tot.num = nrow(mat2[rowSums(mat2[, -1]) > 1, ])
    if (verbose) cmessage(paste0("Total number of rows are: ", tot.num))
    if (verbose) cmessage("taking a subsample")
    number.to.subsample = pmax(round(tot.num * subsample.frac), 1000)
    if (verbose) cmessage(paste0("Number sampled: ", number.to.subsample))
    set.seed(seed)
    concatm = mat2[rowSums(mat2[, -1]) > 1, ][sample(.N, min(c(.N, number.to.subsample))), ]
  }
  else {
    concatm = mat2[rowSums(mat2[, -1]) > 1, ]
  }
  ubx = concatm$cid

  if (verbose) cmessage("Matrices made")
  gc()

  concatm = concatm[, setdiff(which(colSums(concatm) > 1), 1), with = FALSE]

  if (!ncol(concatm))
  {
    warning('No concatemers found hitting two bins, returning empty result')
    return(reads[, bid := NA][c(), ])
  }

  pairs = t(do.call(cbind, apply(concatm[, setdiff(which(colSums(concatm) > 1), 1), with = FALSE] %>% as.matrix, 2, function(x) combn(which(x != 0), 2))))
                                                    
  concatm = as(as.matrix(as.data.frame(concatm)), "sparseMatrix")    
  p1 = concatm[pairs[, 1], -1]
  p2 = concatm[pairs[, 2], -1]
  matching = rowSums(p1 & p2)
  total = rowSums(p1 | p2)
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
  memb.dt = data.table(cid = ubx[1:nrow(A)], bid = cl)
  reads = merge(reads, memb.dt, by = "cid")
  reads[, `:=`(support, length(unique(cid))), by = bid]
  reads = dt2gr(reads)
  return(reads)
}

#' @name smessage
#' @description
#'
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
#' @private 
;smessage = function(..., pre = 'Synergy')
  message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)

#' @name cmessage
#' @description
#'
#' @author Aditya Deshpande, Marcin Imielinski
#' @export
#' @private 
cmessage = function(..., pre = 'Chromunity')
  message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)
