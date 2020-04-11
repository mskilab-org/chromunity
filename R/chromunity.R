#' @import GenomicRanges
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom data.table set
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table setkeyv
#' @importFrom data.table setnames
#' @importFrom data.table transpose
#' @import Matrix
#' @import zoo
#' @importFrom gUtils gr2dt
#' @importFrom gUtils dt2gr
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom MASS ginv
#' @importFrom utils globalVariables
#' @import dplyr




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
#' @param verbose boolean (default == TRUE). Outputs progress.
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
    
    require(Matrix)
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
    mat = dcast.data.table(reads[count > 1 ,]  %>% gr2dt, read_idx ~ tix, value.var = "strand", fill = 0)
    mat2 = mat[, c(list(read_idx = read_idx), lapply(.SD, function(x) x >= 1)),.SDcols = names(mat)[-1]]
    mat2 = mat2[, "NA" := NULL]
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
    knn.dt = dt3.2[mat >= 2 & tot >= 2, .(knn = bx2[1:k]), by = bx1][!is.na(knn), ]
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
    A = A+t(A)
    G = graph.adjacency(A, weighted = TRUE, mode = 'undirected')
    cl.l = cluster_louvain(G)
    cl = cl.l$membership
    message("Communities made")
    memb.dt = data.table(read_idx = ubx[1:nrow(A)], community = cl)
    reads = merge(reads, memb.dt, by = "read_idx")
    reads[, max.clust.dist := max(max.dist), by = community]
    reads[, min.clust.dist := min(max.dist), by = community]
    reads[, num.memb := length(unique(read_idx)), by = community]
    reads = dt2gr(reads)
    return(reads)
}


##############################
## annotate_multimodal_clusters
##############################
#' @name  annotate_multimodal_clusters 
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
#' @param k.min numeric the threshold to number of concatemers pairs to be considered "similar"
#' 
#' @param tiles GRanges object dividing the genome in a fixed size bin to be defined by user
#' 
#' @param which.gr the GRanges for the window of interest
#'
#' @param min.memb numeric minimum number of members needed to be in a community to be considered for further analyses
#' 
#' @return \code{annotate_local_clusters} returns input GRanges with additional columns as follows:
#' 
#'    \item{multimodal}{  boolean; \cr
#'              whether a community had multimodal contacts based on parameters set by user
#'    }
#'  
#' @author Aditya Deshpande

annotate_multimodal_clusters <- function(granges, which.gr, min.memb = 50){
    this.dt = gr2dt(granges)[num.memb > min.memb]
    ind.int = unique(this.dt$community)
    dt.stat = data.table()
    for (i in 1:length(ind.int)){
        which.int = (ind.int)[i]  
        message(which.int)
        dt.int = this.dt[community == which.int]
        dt.int.tmp = gr2dt(dt2gr(dt.int) %&&% which.gr)
        setkeyv(dt.int.tmp, c("seqnames", "start"))
        dt.sum = gr.sum(dt2gr(dt.int.tmp)+1e3)
        y.max.val = max(dt.sum$score)
        peak.gr = find_multi_modes( dt.sum[-1], w =  round(0.05*length(dt.sum[-1])), distance = 0.1*width(which.gr))
        status = unique(peak.gr$bimodal)
        dt.stat = rbind(dt.stat, data.table(community = which.int, multimodal = status))
    }

    this.dt = merge(this.dt, dt.stat, by = "community", allow.cartesian = T)
    return(this.dt)
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
#' @param which.gr the GRanges for the window of interest
#' 
#' @param x.field This is the X axis along which smoothing is done
#' 
#' @param y.field values to be smoothed
#' 
#' @param distance numeric genomic distance beyond which a peak is not considered local 
#' 
#' @return \code{find_multi_modes} returns boolean value if more than one peak is present
#'  
#' @author Aditya Deshpande


find_multi_modes <- function(granges, x.field = "start", y.field = "score", w = 1,  distance = 1e5) {
    which.chr = seqlevels(granges)
    x = start(granges)
    y = values(granges)[, y.field]
    require(zoo)
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

