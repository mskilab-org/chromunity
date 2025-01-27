#for now, dumping all chromunity functions into one file.

rebin_community = function(concatemers, this.chrom.w, resolution = 5e4, rebin_thresh=0.85) {
    tiles = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), resolution)
    this.chrom = gr2dt(concatemers %Q% (chid %in% this.chrom.w))
    this.list.chrom = pbmclapply(1:length(this.chrom.w), function(j){
        this.pr = dt2gr(this.chrom[chid %in% this.chrom.w[j]])
        sum.this.com = gr.sum((this.pr)+1e4)
        sum.this.com = gr2dt(sum.this.com)
        sum.this.com[, q := quantile(score, rebin_thresh), by = seqnames]
        sum.this.com[, q := ifelse(q < 5, 5, q)]
        active.cont = tryCatch((tiles %&% dt2gr(sum.this.com[score > q])), error = function(e) NULL)
        this.clust = gr2dt(gr.reduce(active.cont))
        this.clust[, chid := this.chrom.w[j]]
        return(this.clust)
    }, mc.cores  = 10)
    if(length(this.chrom.w) == 1){
        this.list.chrom = this.list.chrom$value
    }
    this.chrom.dt = rbindlist(this.list.chrom, fill = TRUE)
    return(this.chrom.dt)
}


load_bad_regions = function(chromosome, genome.to.use = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens") {
    this.chr = chromosome
    if(genome.to.use == "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"){
        bands.td = gTrack::karyogram(file = "/gpfs/commons/groups/imielinski_lab/DB/UCSC/hg38.cytoband.txt")
    } else {
        bands.td = gTrack::karyogram(file = "/gpfs/commons/groups/imielinski_lab/DB/UCSC/hg19.cytoband.txt")
    }
    bands = bands.td@data
    bands = grl.unlist(do.call(`GRangesList`, bands))
    cen = bands %Q% (stain=="acen")
    if (!(this.chr %in% c('chrX', 'chrY'))) {
        chr.ind = as.numeric(sub("chr*","",this.chr))
    } else if (this.chr == 'chrX'){
        chr.ind = 23
    } else {
        chr.ind = 24
    }
    
    this.max = GRanges(paste0(this.chr, ":", hg_seqlengths(genome = genome.to.use)[chr.ind]-1e6, "-",  hg_seqlengths(genome = genome.to.use)[chr.ind]))         
    this.min = GRanges(paste0(this.chr, ":", "1-1e6"))                                                                                                                       
    this.cen = (cen %Q% (seqnames == this.chr))+1e6
    this.bad = c(this.min, this.cen, this.max) 
    return(this.bad)
}
    





##let's try to run giga chromunity interchromosomally
evaluate_synergy_interchr = function(concatemers, leave_out_concatemers, chid.to.test, chromosome = NULL, filter_binsets = TRUE, folder = NULL, rebin_thresh = 0.85, mc.cores = 20, numchunks = mc.cores*200 + 1) {
    this.chr = chromosome
    bands.td = gTrack::karyogram(file = "/gpfs/commons/groups/imielinski_lab/DB/UCSC/hg38.cytoband.txt")
    bands = bands.td@data
    bands = grl.unlist(do.call(`GRangesList`, bands))
    cen = bands %Q% (stain=="acen")
    this.max = GRanges(paste0(this.chr, ":", hg_seqlengths()[sub("chr*","",this.chr)]-1e6, "-",  hg_seqlengths()[sub("chr*","",this.chr)]))         
    this.min = GRanges(paste0(this.chr, ":", "1-1e6"))                                                                                                                       
    this.cen = (cen %Q% (seqnames %in% this.chr))+1e6
    this.bad = c(this.min, this.cen, this.max) 
    if (!dir.exists(folder)) {
        stop("output folder does not exist")
    }
    if(is.null(chromosome)){
        chromosome = c(paste0("chr", c(as.character(1:22), "X")))
    }
    resolution = 1e4
    tiles = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), resolution)
    gc5b = readRDS("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/gc.38.rds")
    ## create a list of covariates
    cov_list = list(gc5b)
    ## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
    names(cov_list) <- c("score:gc.cov")
    ## Make the covariate object
    gc_cov = covariate(name = c("gc"), type = c("numeric"), field = c("score"), data = cov_list)
    gc.cov = gc_cov
    this.dat.chrom = data.table()
    this.chrom = gr2dt(concatemers)
    if (filter_binsets) {
        this.chrom = this.chrom[support > summary(unique(this.chrom[, .(support, chid)])$support)[3]]####filtering out all bins below median support
    }
    this.chrom.w = unique(chid.to.test)
    this.chrom.dt = rebin_community(concatemers, this.chrom.w, resolution=resolution)
    this.chrom.dt[, cardinality := .N, by = chid]
    this.chrom.dt = na.omit(this.chrom.dt)
    ##
    this.all.dat = copy(this.chrom.dt)
    this.all.dat = this.all.dat[cardinality > 2]
    this.all.dat[, bid := chid]
    this.all.dat = this.all.dat[seqnames %in% chromosome]
    this.all.dat[, overall.cardinality := cardinality, by = bid]
    this.chrom.card = unique(this.all.dat[, .(overall.cardinality, bid)])
    ######HELLL NAH
    #this.all.dat = this.all.dat[cardinality < 100]
    ####
    this.sub.parq = leave_out_concatemers
    this.sub.parq$cid = this.sub.parq$read_idx
    this.all.dat[, binid := .I]
    ##
    ###filtering out the bad regions
    ###this makes sense why it wouldn't be here for RE chromunity cause you're only looking at annotated regions anyway
    this.all.dat$bid = this.all.dat$chid

    ##only drop the binids in bad regions, keep the binsets as a whole perhaps
    this.bad.binids = unique(as.character((dt2gr(this.all.dat) %&% (this.bad))$binid))                                                 
    this.all.dat = this.all.dat[!binid %in% this.bad.binids]  
    #browser()
    debug(annotate)
    this.chrom.dat = annotate(binsets = dt2gr(this.all.dat[, bid := chid]),
                              k = 3,
                              concatemers = this.sub.parq,
                              covariates = gc.cov, resolution = resolution,
                              mc.cores = mc.cores, numchunks = numchunks)
    this.chrom.dat.2 = this.chrom.dat
    this.chrom.dat = merge(this.chrom.dat, this.chrom.card[, bid := as.factor(bid)], by = "bid")
    this.chrom.dat[, annotation := "chromunity"]
    set.seed(198)
    message("generating background binsets for the model")
##
    #back.dt = sliding_window_background(chromosome = chromosome, binsets = dt2gr(this.all.dat), n = 1000, resolution = resolution)#, mc.cores=mc.cores)
    ##browser()
    back.dt = re_background(binsets = dt2gr(this.all.dat), resolution = resolution, n=dim(this.all.dat)[1])#, mc.cores=mc.cores)
    back.dt = back.dt[start != 1]
    upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
    setkeyv(back.dt, c("seqnames", "start"))
    back.dt[, V1 := NULL]
    back.dt = na.omit(back.dt)
    back.dt = back.dt[!bid %in% back.dt[width < (resolution-1)]$bid]
    back.dt = gr2dt(gr.reduce(dt2gr(back.dt), by = "bid"))
    back.dt$bid <- as.factor(back.dt$bid)
    back.dt = merge(back.dt, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.dt = back.dt[end < V2][start < V2]
    back.dt[, overall.cardinality := .N, by = bid]
    back.dt = back.dt[overall.cardinality > 1]
    ##
    this.card = unique(back.dt[, .(bid, overall.cardinality)])
    back_gr = dt2gr(back.dt)
    message("extracting background binsets distances")
    this.back.train.dat = annotate(binsets = dt2gr(back.dt),
                                   interchromosomal.table = NULL, #all.hr.dt.mean,
                                   gg = NULL,
                                   concatemers = this.sub.parq, k = 3,
                                   covariates = gc.cov, resolution = resolution,
                                   mc.cores = mc.cores, numchunks = numchunks)
    this.back.train.dat = merge(this.back.train.dat, this.card[, bid := as.factor(bid)], by = "bid")
    this.back.train.dat[, annotation := "random"]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][width <= resolution]$bid)]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][min.dist < resolution]$bid)]
    back.dt[, binid := .I]
    this.bad.train = unique(as.character((back_gr %&% (this.bad))$bid))
    this.back.train.dat = this.back.train.dat[!bid %in% this.bad.train]
    this.back.train.dat = this.back.train.dat[!bid %in% this.back.train.dat[, .(sum(count)), by = bid][V1 == 0]$bid]
####    
    back.model = fit(na.omit(this.back.train.dat)[sum.counts > 0][, setdiff(names(this.back.train.dat), c('overall.cardinality', 'chr', 'annotation')), with = F])

    this.chrom.dat = sscore(this.chrom.dat, model = back.model.ep)
    message("generating random binsets for testing")
    n.chrom = length(unique(this.chrom.dat$bid))
    this.all.dat = this.all.dat[seqnames %in% chromosome]
    back.test = re_background(binsets = dt2gr(this.all.dat), resolution = resolution, n=dim(this.all.dat)[1]*3)#, mc.cores=mc.cores)
    #back.test = sliding_window_background(binsets = dt2gr(this.all.dat), chromosome = chromosome,
    #                                      n = n.chrom*3,
    #                                      resolution = resolution)
    back.test[, V1 := NULL]
    back.test = na.omit(back.test)
    setkeyv(back.test, c("seqnames", "start"))
    back.test[, start := ifelse(start < 0, 0, start)]
    back.test = back.test[!bid %in% back.test[width < (resolution-1)]$bid]
    back.test = gr2dt(gr.reduce(dt2gr(back.test), by = "bid"))
    back.test = merge(back.test, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.test = back.test[end < V2][start < V2]
    back.test[, overall.cardinality := .N, by = bid]
    back.test = back.test[overall.cardinality > 1]
    back_test = dt2gr(back.test)
###
    back.test$bid <- as.factor(back.test$bid)
    this.back.card = unique(back.test[, .(bid, overall.cardinality)])
    bid.back = unique(back.test$bid)
    message("extracting random binsets distances") 
    this.back.test.dat = annotate(binsets = dt2gr(back.test), k=sub.binset.order,
                                  concatemers = this.sub.parq,
                                  covariates = gc.cov, resolution = resolution,
                                  mc.cores = mc.cores, numchunks=numchunks)
    this.back.test.dat = merge(this.back.test.dat, this.back.card, by = "bid")
    this.back.test.dat[, annotation := "random"]
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][width <= resolution]$bid)]
###
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][min.dist < resolution]$bid)]
####
####
    this.bad.dat = unique(as.character((back_test %&% (this.bad))$bid))
    this.back.test.dat = this.back.test.dat[!bid %in% this.bad.dat]
    back_test = gr2dt(back_test)[bid %in% unique(this.back.test.dat$bid)]
    #this.back.test.dat = sscore(this.back.test.dat, model = back.model)
####
    back.test[, binid := .I]
#####
    set.seed(178)
    all.bid = unique(this.all.dat$bid)
###
    #browser()
    message("starting random walks")
    message("extracting shuffled binsets distances")
    this.all.dat.shuff5 = pbmclapply(1:length(all.bid), function(nr){
         this.clust = dt2gr(this.all.dat[bid == all.bid[nr]])
         this.chrs = .chr2str(as.character(unique(seqnames(this.clust))))
         this.clust.wind = gr.reduce(this.clust+2e6)
         upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T)
         this.clust.wind = (gr2dt(this.clust.wind)[, start := ifelse(start < 0, 1, start)])
         this.clust.wind = merge(this.clust.wind, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
         this.clust.wind[, end := ifelse(end > V2, V2, end)]
         this.clust.wind = dt2gr(this.clust.wind)
         this.sub.win = gr2dt(this.sub.parq %&% this.clust.wind)
         this.sub.win[, new.count := .N, by = read_idx]
         ##
         card = unique(this.sub.win[, .(read_idx, new.count)])
         this.steps = sum(card$new.count)
         this.tiles.orig = gr.tile(this.clust.wind, resolution)
         if (nrow(this.sub.win) > 0){
             this.tgm = cocount(dt2gr(this.sub.win), bins = (this.tiles.orig), by = 'read_idx', full = T)
             A = this.tgm$mat %>% as.matrix
             rownames(A) <- NULL
             colnames(A) <- NULL
             A[cbind(1:nrow(A), 1:nrow(A))] = 0
             A = A + t(A)
             An = A 
             An = round(10*An/min(An[An>0]))  ##can you remove this background 1 value? does this change anything. peculiar
             An = round(1+10*An/min(An[An>0]))
             edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
             ##
             G = graph.edgelist(edges[, cbind(row, col)])
             RW = random_walk(G, start = 200, steps = sum(card$new.count)) %>% as.numeric
             ##
             rm(G)
             rm(edges)
             gc()
             out = this.tgm$gr[RW]%>% gr2dt()
             out$read_idx = card[, rep(read_idx, new.count)] 
             out[, bid := all.bid[nr]]
             out[, cid := read_idx]
             
             this.chrom.sh.dat = annotate(binsets = this.clust, verbose = F,
                                          k = 3,
                                          concatemers = dt2gr(out.old),
                                          covariates = gc.cov, resolution = resolution, 
                                          mc.cores = 2, numchunks = numchunks)
         } else {
             this.chrom.sh.dat = data.table(NA)
         }
         return(this.chrom.sh.dat)
    }, mc.cores = 5, mc.preschedule = T)
####
##     bid.chrom = unique(this.all.dat$bid)
##     chrom.this.chr.shuff5 = pbmclapply(1:length(bid.chrom), function(nr){
##         this.clust = dt2gr(this.all.dat[bid == bid.chrom[nr]])
##         #this.clust.wind = .collapse_gr(this.clust) + ((window.size))
##         #this.sub.win = gr2dt(this_gr_testing %&% this.clust.wind)
##         this.sub.win[, new.count := .N, by = read_idx]
##         ##
##         card = unique(this.sub.win[, .(read_idx, new.count)])
##         this.steps = sum(card$new.count)
##         this.tiles.orig = gr.tile(this.clust.wind, resolution)
##         if (nrow(this.sub.win) > 0){
##             this.tgm = cocount(dt2gr(this.sub.win), bins = (this.tiles.orig), by = 'read_idx', full = T)
##             A = this.tgm$mat %>% as.matrix
##             rownames(A) <- NULL
##             colnames(A) <- NULL
##             A[cbind(1:nrow(A), 1:nrow(A))] = 0
##             A = A + t(A)
##             An = A
##             An = round(1+10*An/min(An[An>0]))
##             edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
##             ##
##             G = graph.edgelist(edges[, cbind(row, col)])
##             RW = random_walk(G, start = 1, steps = sum(card$new.count)) %>% as.numeric
##             ##
##             out = this.tgm$gr[RW]%>% gr2dt()
##             out$read_idx = card[, rep(read_idx, new.count)]
##             out[, bid := bid.chrom[nr]]
##         } else {
##             out = data.table(NA)
##         }
##         return(out)
##     }, mc.cores = 20)
##     kill.zombies()
##     message("running synergy")
## ###
##     chrom.this.chr.shuff5.ne = chrom.this.chr.shuff5[sapply(chrom.this.chr.shuff5, function(x) !inherits(x, "try-error"))]
##     chrom.this.chr.shuff.dt = rbindlist(chrom.this.chr.shuff5.ne, fill = T)
##     chrom.this.chr.shuff.dt[, `:=`(V1 = NULL)]
##     chrom.this.chr.shuff.dt = na.omit(chrom.this.chr.shuff.dt)
##     chrom.this.chr.shuff.dt[, unique.chr.comm := paste0(seqnames, "_", bid)]
##     chrom.this.chr.shuff.dt[, cid := read_idx]
##     set.seed(198)
## ############
##     message("extracting shuffled binsets distances")
##     the.shuff.list = pbmclapply(1:length(bid.chrom), function(x){
##         this.clust = dt2gr(this.all.dat[bid == bid.chrom[x]])
##         chromosome = unique(as.character(seqnames(this.clust)))
##         if (length(this.clust) %in% c(3:10)){
##             this.clust.wind = (.collapse_gr(this.clust)+((window.size)))
##             if (length(this.clust.wind %&% this.bad)  == 0){
##                 this.chrom = tryCatch(annotate(concatemers = dt2gr(chrom.this.chr.shuff.dt[bid == bid.chrom[x]]), binsets = this.clust, covariates = gc_cov, verbose = F, resolution = resolution, k = 5, mc.cores = 1), error = function(e) NULL)
##                 this.chrom = sscore(this.chrom, model = back_model)
##                 }
##             }
##         },  mc.cores = 20)
##         ##kill.zombies()
## ###
    #this.shuff.chrom = tryCatch(rbindlist(the.shuff.list, fill = T), error = function(e) NULL)
    #browser()
    this.all.dat.shuff5.ne = this.all.dat.shuff5[sapply(this.all.dat.shuff5, function(x) !inherits(x, "try-error"))]
    this.all.sh = rbindlist(this.all.dat.shuff5.ne, fill =T)
###
    this.all.sh = sscore(this.all.sh, model = back.model)
    sh.chrom = tryCatch(synergy(binsets = dt2gr(this.all.dat), annotated.binsets = this.all.sh, model = back.model), error = function(e) NULL)
    sh.chrom$fdr = signif(p.adjust(sh.chrom$p, "BH"), 2)
    ##
#####
    this.back.test.dat = sscore(this.back.test.dat, model = back.model)
    theta = back.model$model$theta
    s.chrom = synergy(binsets = dt2gr(this.all.dat), #theta = back.model$model$theta,
                      annotated.binsets = this.chrom.dat, model = back.model)
    s.chrom$fdr = signif(p.adjust(s.chrom$p, "BH"), 2)
    b.chrom = synergy(binsets = dt2gr(back.test), #theta = back.model$model$theta,
                      annotated.binsets = na.omit(this.back.test.dat), model = back.model)
    b.chrom$fdr = signif(p.adjust(b.chrom$p, "BH"), 2)
    ##
    synergy.inter.EP = rbind(s.chrom[, annotation := "chromunity"],
                             b.chrom[, annotation := "random"], 
                             sh.chrom[, annotation := "shuffled"], fill = T)
    if (!is.null(folder)) {
        saveRDS(this.chrom.dat, paste0(folder,'chrom_annotate.rds'))
        saveRDS(this.back.test.dat, paste0(folder,'back_annotate.rds'))
        saveRDS(this.all.sh, paste0(folder,'shuffled_annotate.rds'))
        saveRDS(this.all.dat, paste0(folder,'binsets.rds'))
        saveRDS(back.model, paste0(folder,'back_model.rds'))
        saveRDS(synergy.inter.EP, paste0(folder,'synergy_results.rds'))
    }
    return(synergy.inter.EP)
}




shuffle_concatemers = function(concatemers, contact_matrix) {
    A = contact_matrix$mat %>% as.matrix
    rownames(A) <- NULL
    colnames(A) <- NULL
    A[cbind(1:nrow(A), 1:nrow(A))] = 0
    A = A + t(A)
    An = A 
    An = round(1+10*An/min(An[An>0]))  
    edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]

    G = graph.edgelist(edges[, cbind(row, col)])
    
    concats.binned = bin_concatemers(concatemers, contact_matrix$gr)
    concats.dt = concats.binned
    concats.dt = unique(concats.dt, by=c('binid','cidi')) ###dedupe!!

    concat.counts = concats.dt[, new.count := .N, by='read_idx']


    card = unique(concat.counts[, .(read_idx, new.count)])

    this.steps = sum(card$new.count)

    row.labels = 1:dim(An)[1]
    start.index = row.labels[rowSums(An)>1]
    
    RW = random_walk(G, start = start.index, steps = sum(card$new.count)) %>% as.numeric
    out = contact_matrix$gr[RW] %>% gr2dt()
    out$read_idx = card[, rep(read_idx, new.count)] 

}




    






derive_binsets_from_network = function(G.kant, pairwise, binned.concats, bins, rr.thresh = 0, dist.decay.test=NULL, all.pairs.tested=NULL, num.members=10, pairwise.trimmed=pairwise, expansion.cutoff = 0.5, fdr.expansion.thresh = 0.25, fdr.thresh=0.1, rebin.resolution=10000) {
##num.members=2
    ##G.kant = bin.pair.network
    
    pairwise.trimmed$agg = do.call(Map, c(f = c, pairwise.trimmed[, c('i','j')]))
    G.kant.10 = subgraph.edges(G.kant, E(G.kant)[E(G.kant)$weight > 0])
    G.kant.undir = as.undirected(G.kant.10)
    cl.l = cluster_leiden(G.kant.undir, objective_function='modularity')
    cl = cl.l$membership
    dt.membership = data.table(cluster = cl)
    dt.membership$pair = V(G.kant.10)$name
    ##dt.membership$classification = V(G.kant.10)$classification
    dt.membership[, num.memb := .N, by='cluster']

    print(dt.membership)
    unique.clusters = dt.membership[num.memb > num.members]$cluster %>% unique
    print(unique.clusters)
    print(dist.decay.test)

    dt.eig.sub = pbmclapply(unique.clusters, mc.cores = 5, function(cluster.id) {
        G.sub = subgraph(G.kant.10, dt.membership[cluster==cluster.id]$pair)
        sub.eig = eigen_centrality(G.sub, directed=TRUE)##, cutoff=-1)##, options)
        dt.subgraph = data.table(pairs = V(G.sub)$name)
        dt.subgraph$eigen = sub.eig[[1]]
        if(class(pairwise.trimmed$pair.hashes)=='integer'){
            dt.subgraph$pairs = strtoi(dt.subgraph$pairs)
        }
        setorder(dt.subgraph, -eigen)
        dt.eig.sub = merge(dt.subgraph, pairwise.trimmed[, c('pair.hashes','agg')], by.x='pairs', by.y='pair.hashes')
        dt.eig.sub$cluster = cluster.id
        return(dt.eig.sub)
    })

    
    if(length(unique.clusters) == 1) {
        dt.eig.sub.agg = rbindlist(dt.eig.sub[[1]])
        dt.eig.sub.agg$eigen = 1
    } else {
        dt.eig.sub.agg = rbindlist(dt.eig.sub)
    }
        
    dt.eig.sub.agg[, max.eig := max(eigen), by='cluster']
    seeds = dt.eig.sub.agg[eigen == max.eig]

    colnames(seeds)[1] = 'pair.hashes'
    seeds = unique(seeds[max.eig>0], by='cluster')

    print(seeds)
    print("FDR EXPANSION THRESH")
    print(fdr.expansion.thresh)
    binsets.gr = expand_seeds(seeds, dist.decay.test, pairwise, all.pairs.tested=NULL, bins, fdr.thresh=fdr.thresh, fdr.expansion.thresh, expansion.cutoff=expansion.cutoff)
    binsets.dt = gr2dt(binsets.gr)
    binsets.dt$chid = binsets.dt$bid
    ##browser()

    ppdf(plot(plot_binsets(dt2gr(binsets.dt))), 'debug_plot.pdf')
    ## bins = dist.decay.test[pair.hashes=='c(2547, 2552)']

    binned.concats[, bincount := .N, by='binid']
    ##binned.concats = binned.concats[mapping_quality > 5]
    thresh = quantile(unique(binned.concats, by='binid')$bincount, 0.9999) ###Remove very highly mapped bins
    standard.dev = sd(unique(binned.concats, by='binid')$bincount)
    binned.concats = binned.concats[bincount < (thresh + standard.dev)]
    
    
    binned.concats[, bid := NULL]
    binned.concats[, chid := NULL]
    concat.merge = merge(binned.concats, binsets.dt[, c('binid','bid')], by='binid', allow.cartesian=TRUE)
    concat.count = concat.merge[, .N, by=c('cid','bid')]

##    binsets.dt
    concat.contacts = concat.count[N>2]
##    tpe1.dt = gr2dt(tpe1_chr8)
    ##concat.contacts
    
    concat.contacts$cid = concat.contacts$cid %>% as.integer
##    tpe1.chr8.binned$cid = tpe1.chr8.binned$cid %>% as.integer
    binned.concats$cid = binned.concats$cid %>% as.integer
    ##binned.concats = unique(binned.concats, by=c('cid','binid')) ##we are so BACK


    ##binned.concats.2 = binned.concats[, c('min.start','max.end','seqnames','strand','cidi','cid','binid')]
    ##binned.concats.2[, start := min.start]
    ##binned.concats.2[, end := max.end]

    binned.concats = unique(binned.concats, by=c('cidi','binid'))

    ###browser()

    contacted.concatemers = merge(binned.concats, concat.contacts, by='cid', all.y=TRUE, allow.cartesian=TRUE)
    contacted.concatemers = dt2gr(contacted.concatemers)
    contacted.concatemers$chid = contacted.concatemers$bid

    
    new.binsets = rebin_community(contacted.concatemers, unique(concat.contacts$bid), resolution=rebin.resolution, rebin_thresh=0.85)

    
    new.binsets.dt = gr2dt(new.binsets)
    new.binsets.dt[, numbins := .N, by='chid']
    new.binsets.dt = new.binsets.dt[width<200000 & numbins<500] ##I will hardcode this as limit: something has gone horribly wrong if you reach this

    new.binsets = dt2gr(new.binsets.dt[numbins > 2])
    chrom = Chromunity(binsets=new.binsets, concatemers=contacted.concatemers, meta=data.table())
    return(chrom)

}



######PIPELINE





expand_seeds = function(seeds, dist.decay.test, pairwise, all.pairs.tested, bins, expansion.cutoff = 0.5, fdr.expansion.thresh=0.25, fdr.thresh=0.1, plot_expansion=FALSE, debug_mode=FALSE) {
    seeds.archive = seeds %>% copy

    seeds = unique(seeds, by=c('cluster'))
    seeds.loop = seeds[, c('pair.hashes','cluster')]
    
    
    memb.dt = seeds[, .(unlist(agg)), by='cluster']
    memb.dt$memb = TRUE
   
    i=1
    bins=gr.stripstrand(bins)
    if(class(dist.decay.test$pair.hashes)=='integer'){
        seeds.loop$pair.hashes = strtoi(seeds.loop$pair.hashes)
        all.pairs.tested$pair.hashes = strtoi(all.pairs.tested$pair.hashes)
    }

    pairwise.2 = pairwise %>% copy
    pairwise.2$i = pairwise$j
    pairwise.2$j = pairwise$i
    pairwise.sym = rbind(pairwise, pairwise.2)
    ##colnames(pairwise.sym)[colnames(pairwise.sym)=='id'] = 'pair.hashes'

    ##brwoser()

    memb.dt.aggregate = data.table()
    
    for(i in 1:50){
        if(class(dist.decay.test$pair.hashes)=='integer'){
            seeds.loop$pair.hashes = strtoi(seeds.loop$pair.hashes)
        }

        seed.decays = merge(dist.decay.test, seeds.loop, by.x='pair.hashes', by.y='pair.hashes', allow.cartesian=TRUE)
        ##seed.decays = dist.decay.test[pair.hashes %in% seeds.loop$pair.hashes]

        if(i==1) {
            seed.decays = seed.decays[fdr<fdr.thresh]
        }

        ##seed.decays = merge(seed.decays, seeds.loop[, c('pair.hashes','cluster')], by='pair.hashes', allow.cartesian=TRUE)
        print(fdr.expansion.thresh)
        seed.decays[fdr < fdr.expansion.thresh, sum.relative.risk := sum(relative.risk), by=c('binterrogate','cluster')]
        seed.decays[, num.concats := sum(num.concats), by=c('binterrogate','cluster')]
        seed.decays[, num.concats.pred := sum(num.concats.pred), by=c('binterrogate','cluster')]
        
        
        if (plot_expansion == TRUE) {
            plot_expansion_increment(seed.decays, 1, bins, folder='MYC_binset_increment_old/', filename=i)
        }

        if (debug_mode==TRUE) {
            seed.decays$iteration = i
            memb.dt.aggregate = rbind(memb.dt.aggregate, seed.decays, fill=TRUE)
        }
            ##seed.decays = merge(seed.decays, all.pairs.tested, by=c('pair.hashes','binterrogate'))

        ####if we want to be verbose with it

        ##to.drop = merge(memb.dt, seed.decays, by.x=c('cluster','V1'), by.y=c('cluster','binterrogate'), all.y=TRUE, allow.cartesian=TRUE)
        
###also: DROP ANY BIN THAT IS DIRECTLY ADJACENT TO ANOTHER MEMBER BIN
        
        
        adj.bins = memb.dt[, .(adj.bins = list((V1-1):(V1+1))), by=c('cluster','V1')]
        adj.bins[, V1 := NULL]
        adj.bins = adj.bins[, .(unlist(adj.bins)), by='cluster']
        adj.bins = unique(adj.bins, by=c('cluster','V1'))
        adj.bins$memb = TRUE
        
        to.drop = merge(adj.bins, seed.decays, by.x=c('cluster','V1'), by.y=c('cluster','binterrogate'), all.y=TRUE, allow.cartesian=TRUE)
        seed.decays = to.drop[is.na(memb)][, memb := NULL]
        seed.decays[, binterrogate := V1]
        
        
        seed.decays[!is.na(sum.relative.risk), max.rr := max(sum.relative.risk), by='cluster']
        
        bin.to.add = seed.decays[sum.relative.risk == max.rr, bin.to.add := V1][!is.na(bin.to.add), c('cluster','bin.to.add')] %>% unique(by='cluster')
        
        colnames(bin.to.add)[2] = 'bin.discovered'
        pairs.discovered = merge(seed.decays, bin.to.add, by='cluster')
        
        pairs.discovered = pairs.discovered[V1 == bin.discovered]

        memb.dt[, numbins := .N, by='cluster']
        ##pairs.discovered[, numbins := .N, by='cluster']
        unique.cluster.numbins = memb.dt %>% unique(by='cluster')
        pairs.discovered = merge(unique.cluster.numbins[, c('cluster','numbins')], pairs.discovered, by='cluster')
        
        check = pairs.discovered[V1 == bin.discovered]

        print('check')
        print(check[cluster==4282])

        prop.check = check[, (choose(numbins, 2) - sum(is.na(bin.to.add))) / choose(numbins, 2), by='cluster']


        cluster.labeled = merge(seeds.loop, dist.decay.test, by='pair.hashes', allow.cartesian=TRUE)
        all.cluster.check = merge(cluster.labeled, bin.to.add, by.x=c('cluster','binterrogate'), by.y=c('cluster','bin.discovered'))

        paircount = seeds.loop[, .(all.count = .N), by='cluster']
        all.cluster.check = merge(all.cluster.check, paircount, by='cluster', all.x=TRUE)

        ##all.cluster.check[pval < pval.thresh & relative.risk>1, enrich.count := .N, by='cluster'] ##I don't think this will be toooooo consequential...
        all.cluster.check[fdr < fdr.expansion.thresh, enrich.count := .N, by='cluster']
        all.cluster.check[is.na(enrich.count), enrich.count := 0]
        all.cluster.check = all.cluster.check[, .(prop = enrich.count / all.count), by='cluster']
        all.cluster.check[, prop := max(prop), by='cluster']
        all.cluster.check[, kill := NA]
        all.cluster.check[prop < expansion.cutoff, kill := TRUE]

        ##prop.check[V1 < 0.5, kill := TRUE]
        print(prop.check)
        
        prop.check = all.cluster.check
        print(prop.check)
        
        pairs.discovered = pairs.discovered[V1 == bin.discovered]


        add.bins.1 = merge(pairs.discovered[, c('i','j','binterrogate','cluster')], pairwise.sym, by.x=c('i','binterrogate'), by.y=c('i','j'), all.x=TRUE)
        add.bins.2 = merge(pairs.discovered[, c('i','j','binterrogate','cluster')], pairwise.sym, by.x=c('j','binterrogate'), by.y=c('i','j'), all.x=TRUE)
        add.bins = rbind(add.bins.1, add.bins.2)


        seeds.loop.chunk = add.bins[, c('pair.hashes','cluster')] %>% unique(by=c('pair.hashes','cluster'))
        print(seeds.loop.chunk)
        
        seeds.loop = rbind(seeds.loop, seeds.loop.chunk) %>% unique(by=c('pair.hashes','cluster'))

        prop.check = unique(prop.check, by='cluster')
        seeds.loop = merge(seeds.loop, prop.check[, c('cluster','kill')], by='cluster')
        seeds.loop = seeds.loop[is.na(kill)] ##stop if you have less than 50% cheunks

        if(dim(seeds.loop)[[1]] == 0){ ##stopping condition
            break
        }
        
        seeds.loop[, kill := NULL]
        
        memb.dt.chunk = merge(seeds.loop, pairwise.sym[, c('pair.hashes','i')], by='pair.hashes')
        ##memb.dt.chunk = memb.dt.chunk[, .(unlist(agg)), by='cluster']
        memb.dt.chunk = unique(memb.dt.chunk, by=c('i','cluster'))
        colnames(memb.dt.chunk)[colnames(memb.dt.chunk) == 'i'] = 'V1'
        memb.dt.chunk$memb = TRUE

        memb.dt.chunk[, numbins.chunk := .N, by='cluster']


        ###kill any bin where the discovered bin is not added to memb dt
        kill.list = merge(memb.dt, unique(memb.dt.chunk[, c('cluster','numbins.chunk')], by='cluster'), by='cluster')
        kill.list.cluster = kill.list[numbins.chunk == numbins]$cluster %>% unique
        seeds.loop = seeds.loop[!(cluster %in% kill.list.cluster)]
        
        memb.dt[, numbins := NULL]
        memb.dt.chunk[, numbins.chunk := NULL]
        memb.dt.chunk[, pair.hashes := NULL]
        memb.dt = rbind(memb.dt, memb.dt.chunk)
        memb.dt = unique(memb.dt, by=c('cluster','V1'))

        seed.decays$iteration = i
        memb.dt.aggregate = rbind(memb.dt.aggregate, seed.decays, fill=TRUE)
        print(memb.dt[cluster==4282])
        print(memb.dt[cluster==30])
    }

    binsets.gr = bins[memb.dt$V1]
    binsets.gr$bid = memb.dt$cluster

    print(binsets.gr)
    print(memb.dt)
    binsets.gr$binid = memb.dt$V1
    if(debug_mode==FALSE){
        return(binsets.gr)
    } else {
        return(list(binsets.gr, memb.dt.aggregate))
    }
    ##merge(dist.decay.test, asdf[, c('pair.hashes','cluster')], by='pair.hashes')
}    






####Some version of this which just outputs gTracks might not be useless...
plot_model_output = function(dt.small.model, pair.hash, bins, plotname='plot.pdf', max.y=NULL, cluster=NULL, view_range=NULL, fdr.thresh=0.25, pairwise=NULL) {
    grange.out = grange_model_prediction(dt.small.model, pair.hash, bins, pairwise, cluster=cluster, fdr.thresh=fdr.thresh)

    model.gr = grange.out[[1]]
    ##view_range = grange.out[[2]]
   
    viewpoint = grange.out[[3]]

    if(is.null(view_range)) {
        ##view_range = viewpoint[[1]] + 1e6
        view_range = (model.gr %Q% (fdr<fdr.thresh) + 1e6) %>% gr.reduce
        dt.range = rbind(gr2dt(view_range), gr2dt(viewpoint[[1]]), gr2dt(viewpoint[[2]]), fill=TRUE)
        view_range = dt2gr(dt.range)
    }
    if(is.null(max.y)){
        max.y = model.gr$num.concats %>% max
    }

    model.gr$log.pval = -log10(model.gr$pval)

    if(is.null(cluster)){
        ppdf(plot(c(gTrack(model.gr, y.field='relative.risk', bars=T, name='log2(O/E)', y0=0), gTrack(model.gr, y.field='log.pval', bars=T, name='-log10(pval)'), gTrack(model.gr, y.field='num.concats.pred', bars=T, name='prediction', y1=max.y), gTrack(model.gr, y.field='num.concats', bars=T, name='actual', y1=max.y), gTrack(viewpoint, name='viewpoint')), view_range+1e6), plotname)
    } else {
        ppdf(plot(c(gTrack(model.gr, y.field='relative.risk', bars=T, name='log2(O/E)', y0=0), gTrack(viewpoint, name='viewpoint')), view_range+5e6), plotname)
    }
    return(list(view_range, max.y))
}



aggregate_synergy_results = function(dir, toplevel=TRUE, strict.check=FALSE) {

    dirs = list.dirs(dir, recursive=FALSE)
    print(dirs)
#####aggregate results function here basically

    ##toplevel=TRUE
    agg.synergy = pbmclapply(dirs, mc.cores = 2, mc.preschedule=FALSE, function(dir) {        
        ##if(dir=="GM12878_dist_decay_all_chr/dist_decay_chr21_knn25_kmin5_resolution50000"){ hardcoded trash
        ##    return(NULL)
        ##}
        if(toplevel==TRUE){
            synergy.chunk = readRDS(paste0(dir, '/synergy_outputs/synergy_results.rds'))
            binsets = readRDS(paste0(dir, '/synergy_outputs/binsets.rds'))
        } else {
            synergy.chunk = readRDS(paste0(dir, '/synergy_results.rds'))
            binsets = readRDS(paste0(dir, '/binsets.rds'))
        }
        if(strict.check==TRUE & toplevel==TRUE) {
            annotate = readRDS(paste0(dir, '/synergy_outputs/chrom_annotate.rds'))
            synergy.chunk = symmetry_check(annotate, synergy.chunk, binsets)
        } else if (strict.check == TRUE & toplevel==FALSE) {
            annotate = readRDS(paste0(dir, '/chrom_annotate.rds'))
            synergy.chunk = symmetry_check(annotate, synergy.chunk, binsets)
        }
        return(synergy.chunk)
    })
    agg.synergy = rbindlist(agg.synergy)
    agg.synergy[, .N, by='annotation']


    agg.synsets = pbmclapply(dirs, mc.cores = 2, function(dir) {
        if(toplevel==TRUE){
            synergy.chunk = readRDS(paste0(dir, '/synergy_outputs/synergy_results.rds'))
            binsets = readRDS(paste0(dir, '/synergy_outputs/binsets.rds'))
        } else {
            synergy.chunk = readRDS(paste0(dir, '/synergy_results.rds'))
            binsets = readRDS(paste0(dir, '/binsets.rds'))
        }
        ##synergy.chunk = readRDS(paste0(dir, '/synergy_results.rds'))
        ##binsets = readRDS(paste0(dir, '/binsets.rds'))
        binsets$bid = factor(binsets$bid)
        synsets = merge(binsets, synergy.chunk[annotation=='chromunity', c('bid','fdr')], by='bid')
        return(synsets)
    })


    
    agg.synsets = rbindlist(agg.synsets)
    synergies = dt2gr(agg.synsets[fdr<.1])
    return(list(agg.synergy, synergies))
}


###CANONICAL FUNCTION


#' @name evaluate_synergy_experimental
#' @description
#'
#' This function bundles together other functions
#' 
#' 
#' @param res Chromunity object
#' @param leave_out_concatemers Concatemers GRange object
#' @param bins GRanges of bins which specify the 
#' @param interchromosomal.dist numeric scalar of "effective" distance for inter chromosomal bins [1e8]
#' @param training.chr Chromosome to use for training distance decay model
#' @param pair.thresh Pairwise contact value used as threshold for considering pairs in analysis
#' @param numchunks Number of chunks to create in annotating higher order distance decay
#' @param mask.bad.regions Will load human telomeres/centromeric regions and remove them from annotated distance deca
#' @param fdr.thresh FDR threshold used in bin-pair network construction
#' @param fdr.expansion.thresh FDR threshold used in performing internally disconnected check to stop bin-set expansion.
#' @param expansion.cutoff Parameter which influences the size of bin-sets produced. Higher values will favor more bins being added.
#' @param num.members Minimum number of nodes in a community for a group to be nominated into a bin-set
#' @param folder Folder to save files
#' @param chromosome Chromosome to nominate bin-sets for. Will run on the entire genome by default.
#' @param model A trained model to score higher order interactions 
#' @pval.thresh
#' @param verbose logical flag
#' @param mc.cores integer how many cores to parallelize
#' @param threads used to set number of data table threads to use with setDTthreads function, segfaults may occur if >1
#' @author Jameson Orvis
#' @export
#' @return 

evaluate_synergy_experimental = function(res, leave_out_concatemers, chid.to.test, chromosome = NULL, filter_binsets = TRUE, folder = NULL, rebin_thresh = 0.85, mc.cores = 20, numchunks = mc.cores*200 + 1, dont_rebin=FALSE, sub.binset.order=3, resolution=1e4, remove.bad=TRUE, genome.to.use =  "BSgenome.Hsapiens.UCSC.hg38::Hsapiens") {

    if (!dir.exists(folder)) {
        stop("output folder does not exist")
    }
    if(is.null(chromosome)){
        chromosome = c(paste0("chr", c(as.character(1:22), "X", "Y"))) ###we should add the Y CHROMOSOME!!!
    }
    
    if(length(chromosome) > 1) {
        all.bad = pbmclapply(chromosome, mc.preschedule=FALSE, function(chr) {
            bad.gr = muffle(load_bad_regions(chr, genome.to.use = genome.to.use))
            bad.dt = gr2dt(bad.gr)
            return(bad.dt)
        })
        this.bad = rbindlist(all.bad) %>% dt2gr
    } else {
        this.bad = load_bad_regions(chromosome, genome.to.use = genome.to.use)
    }

    tiles = gr.tile(hg_seqlengths(genome = genome.to.use), resolution)
    gc5b = readRDS("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/gc.38.rds")
    ## create a list of covariates
    cov_list = list(gc5b)
    ## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
    names(cov_list) <- c("score:gc.cov")
    ## Make the covariate object
    gc_cov = covariate(name = c("gc"), type = c("numeric"), field = c("score"), data = cov_list)
    gc.cov = gc_cov

    if(!(class(res) == 'data.table'))
        this.chrom.dt = gr2dt(res$binsets)
    else
        this.chrom.dt = res
    
    this.chrom.dt[, cardinality := .N, by = chid]
    this.chrom.dt = na.omit(this.chrom.dt)
    ##
    this.all.dat = copy(this.chrom.dt)
    this.all.dat = this.all.dat[cardinality > 2]
    this.all.dat[, bid := chid]
    this.all.dat = this.all.dat[seqnames %in% chromosome]
    this.all.dat[, overall.cardinality := cardinality, by = bid]
    this.chrom.card = unique(this.all.dat[, .(overall.cardinality, bid)])
    ######HELLL NAH
    this.all.dat = this.all.dat[cardinality < 100]
    ####
    this.sub.parq = leave_out_concatemers
    this.sub.parq$cid = this.sub.parq$read_idx
    this.all.dat[, binid := .I]
    ##
    ###filtering out the bad regions
    ###this makes sense why it wouldn't be here for RE chromunity cause you're only looking at annotated regions anyway

###remove only bins which intersect bad region, not all

    this.all.dat$bid = this.all.dat$chid
    this.bad.chrom = unique(as.character((dt2gr(this.all.dat) %&% (this.bad))$bid))                                                 
    print(this.bad.chrom)
    if (remove.bad)
        this.all.dat = this.all.dat[!bid %in% this.bad.chrom]  
                                        #browser()
    #debug(annotate)
    chrom.annot.output = annotate(binsets = dt2gr(this.all.dat[, bid := chid]),
                              k = sub.binset.order,
                              concatemers = this.sub.parq,
                              covariates = gc.cov, resolution = resolution,
                              mc.cores = mc.cores, numchunks = numchunks)

    this.chrom.dat = chrom.annot.output[[1]]
    this.chrom.dat.2 = this.chrom.dat
    this.chrom.dat = merge(this.chrom.dat, this.chrom.card[, bid := as.factor(bid)], by = "bid")
    this.chrom.dat[, annotation := "chromunity"]
    set.seed(198)
    message("generating background binsets for the model")
    ##
                                        #browser()

    n.chrom = length(unique(this.chrom.dat$bid))
    if(length(chromosome) > 1){
        back.dt = help_bins_2(binsets = dt2gr(this.all.dat), n = (n.chrom*3), resolution = resolution, genome.to.use = genome.to.use)#, mc.cores=mc.cores)
    } else {

        back.dt = sliding_window_background(chromosome = chromosome, binsets = dt2gr(this.all.dat), n = (n.chrom*3), resolution = resolution, genome.to.use = genome.to.use)#, mc.cores=mc.cores)
    }
    
    upper.bound = as.data.table(hg_seqlengths(genome = genome.to.use), keep.rownames = T)
    setkeyv(back.dt, c("seqnames", "start"))
    back.dt[, V1 := NULL]
    back.dt = na.omit(back.dt)
    back.dt = back.dt[!bid %in% back.dt[width < (resolution-1)]$bid]
    back.dt = gr2dt(gr.reduce(dt2gr(back.dt), by = "bid"))
    back.dt$bid <- as.factor(back.dt$bid)
    back.dt = merge(back.dt, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.dt = back.dt[end < V2][start < V2]
    back.dt[, overall.cardinality := .N, by = bid]
    back.dt = back.dt[overall.cardinality > 1]
    ##
    this.card = unique(back.dt[, .(bid, overall.cardinality)])
    back_gr = dt2gr(back.dt)
    message("extracting background binsets distances")



    back.train.output = annotate(binsets = dt2gr(back.dt),
                                   interchromosomal.table = NULL, #all.hr.dt.mean,
                                   gg = NULL,
                                   concatemers = this.sub.parq, k = sub.binset.order,
                                   covariates = gc.cov, resolution = resolution,
                                   mc.cores = mc.cores, numchunks = numchunks)
    this.back.train.dat = back.train.output[[1]]
    this.back.train.dat = merge(this.back.train.dat, this.card[, bid := as.factor(bid)], by = "bid")
    this.back.train.dat[, annotation := "random"]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][width <= resolution]$bid)]
    this.back.train.dat = this.back.train.dat[!bid %in% unique(this.back.train.dat[count > 1][min.dist < resolution]$bid)]
    back.dt[, binid := .I]
    this.bad.train = unique(as.character((back_gr %&% (this.bad))$bid))
    this.back.train.dat = this.back.train.dat[!bid %in% this.bad.train]
    this.back.train.dat = this.back.train.dat[!bid %in% this.back.train.dat[, .(sum(count)), by = bid][V1 == 0]$bid]
####    

    ##browser()
    
    back.model = fit(na.omit(this.back.train.dat)[sum.counts > 0][, setdiff(names(this.back.train.dat), c('overall.cardinality', 'chr', 'annotation')), with = F])


    
    this.chrom.dat = sscore(this.chrom.dat, model = back.model)
    message("generating random binsets for testing")
    n.chrom = length(unique(this.chrom.dat$bid))
    this.all.dat = this.all.dat[seqnames %in% chromosome]

    if(length(chromosome) > 1){
        back.test = help_bins_2(binsets = dt2gr(this.all.dat),
                                          n = n.chrom,
                                          resolution = resolution, genome.to.use = genome.to.use)
     } else {
        back.test = sliding_window_background(binsets = dt2gr(this.all.dat), chromosome = chromosome,
                                            n = n.chrom,
                                            resolution = resolution, genome.to.use = genome.to.use)
    }                      

    back.test[, V1 := NULL]
    back.test = na.omit(back.test)
    setkeyv(back.test, c("seqnames", "start"))
    back.test[, start := ifelse(start < 0, 0, start)]
    back.test = back.test[!bid %in% back.test[width < (resolution-1)]$bid]
    back.test = gr2dt(gr.reduce(dt2gr(back.test), by = "bid"))
    back.test = merge(back.test, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
    back.test = back.test[end < V2][start < V2]
    back.test[, overall.cardinality := .N, by = bid]
    back.test = back.test[overall.cardinality > 1]
    back_test = dt2gr(back.test)
###
    back.test$bid <- as.factor(back.test$bid)
    this.back.card = unique(back.test[, .(bid, overall.cardinality)])
    bid.back = unique(back.test$bid)
    message("extracting random binsets distances") 
    back.test.output = annotate(binsets = dt2gr(back.test), k=sub.binset.order,
                                  concatemers = this.sub.parq,
                                  covariates = gc.cov, resolution = resolution,
                                  mc.cores = mc.cores, numchunks=numchunks)
    this.back.test.dat = back.test.output[[1]]
    this.back.test.dat = merge(this.back.test.dat, this.back.card, by = "bid")
    this.back.test.dat[, annotation := "random"]
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][width <= resolution]$bid)]
###
    this.back.test.dat = this.back.test.dat[!bid %in% unique(this.back.test.dat[count > 1][min.dist < resolution]$bid)]
####
####
    this.bad.dat = unique(as.character((back_test %&% (this.bad))$bid))
    this.back.test.dat = this.back.test.dat[!bid %in% this.bad.dat]
    back_test = gr2dt(back_test)[bid %in% unique(this.back.test.dat$bid)]
    #this.back.test.dat = sscore(this.back.test.dat, model = back.model)
####
    back.test[, binid := .I]
#####
    set.seed(178)
    all.bid = unique(this.all.dat$bid)
###
    ##browser()
    message("starting random walks")
    message("extracting shuffled binsets distances")

    ##browser()
    ##nr = 1
    
    ##resolution=1e4

    ##all.bid
    ##
    ##browser()
    s.chrom = synergy(binsets = dt2gr(this.all.dat), #theta = back.model$model$theta,
                      annotated.binsets = this.chrom.dat, model = back.model)
    s.chrom$fdr = signif(p.adjust(s.chrom$p, "BH"), 2)
    print('the hit rate!!!')
    print(s.chrom[, .N])
    print(s.chrom[fdr<0.1, .N])

    saveRDS(s.chrom, paste0(folder,'synergy_results.rds'))
    saveRDS(this.chrom.dat, paste0(folder,'chrom_annotate.rds'))
    
    this.all.dat.shuff5 = pbmclapply(1:length(all.bid), function(nr){
         this.clust = dt2gr(this.all.dat[bid == all.bid[nr]])
         this.chrs = .chr2str(as.character(unique(seqnames(this.clust))))
         this.clust.wind = gr.reduce(this.clust+2e6)
         upper.bound = as.data.table(hg_seqlengths(genome = genome.to.use), keep.rownames = T)
         this.clust.wind = (gr2dt(this.clust.wind)[, start := ifelse(start < 0, 1, start)])
         this.clust.wind = merge(this.clust.wind, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
         this.clust.wind[, end := ifelse(end > V2, V2, end)]
         this.clust.wind = dt2gr(this.clust.wind)

         ##this.clust.wind
         this.sub.win = gr2dt(this.sub.parq %&% this.clust.wind)
         ##this.sub.win[, new.count := .N, by = read_idx]

         ###throw this in there
         

         
         this.tiles.orig = gr.tile(this.clust.wind, resolution)
         
         this.tiles.orig$binid = 1:length(this.tiles.orig)
         this.sub.win = bin_concatemers(dt2gr(this.sub.win), this.tiles.orig)
         this.sub.win = unique(this.sub.win, by=c('read_idx','binid'))
         this.sub.win[, new.count := .N, by = read_idx]
         ##
         card = unique(this.sub.win[, .(read_idx, new.count)])
         this.steps = sum(card$new.count)

         if (nrow(this.sub.win) > 0){
             this.tgm = cocount(dt2gr(this.sub.win), bins = (this.tiles.orig), by = 'read_idx', full = T)
             ##debug(shuffle_concatemers)
             ##shuff.concats = shuffle_concatemers(this.sub.win, this.tgm)
             ##shuff.concats$cid = shuff.concats$read_idx

             ##ppdf(plot(c(shuff.contacts$gtrack(name='shuff', clim=c(0,100)), this.tgm$gtrack(name='normal',clim=c(0,100))), gr.reduce(this.tiles.orig) + 3e5))
             ##ppdf(plot(this.tgm$gtrack(clim=c(0,100)), gr.reduce(this.tiles.orig) + 3e5))
             ##shuff.contacts = cocount(dt2gr(shuff.concats), bins = this.tiles.orig, by='read_idx', full=T)
             A = this.tgm$mat %>% as.matrix
             rownames(A) <- NULL
             colnames(A) <- NULL
             A[cbind(1:nrow(A), 1:nrow(A))] = 0
             A = A + t(A)
             An = A 
             An = round(1+10*An/min(An[An>0]))
             edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
             ## ##
             G = graph.edgelist(edges[, cbind(row, col)])
             RW = random_walk(G, start = 1, steps = sum(card$new.count)) %>% as.numeric
             ## ##
             rm(G)
             rm(edges)
             gc()
             out = this.tgm$gr[RW]%>% gr2dt()
             out = out[1:sum(card$new.count)] ##weird bug
             out$read_idx = card[, rep(read_idx, new.count)] 
             out[, bid := all.bid[nr]]
             out[, cid := read_idx]

             ##debug(annotate)
             sh.output = annotate(binsets = this.clust, verbose = F,
                                          k = sub.binset.order,
                                          concatemers = dt2gr(out),
                                          covariates = gc.cov, resolution = resolution, 
                                          mc.cores = 1, numchunks = numchunks)
             ## shuff.concats$cid = shuff.concats$read_idx
             ## sh.output.2 = annotate(binsets = this.clust, verbose = F,
             ##                              k = 5,
             ##                              concatemers = dt2gr(shuff.concats),
             ##                              covariates = gc.cov, resolution = resolution, 
             ##                              mc.cores = 10, numchunks = numchunks)
             
             this.chrom.sh.dat = sh.output[[1]]
         } else {
             this.chrom.sh.dat = data.table(NA)
         }
         return(this.chrom.sh.dat)
    }, mc.cores = 5, mc.preschedule = T)
    #this.shuff.chrom = tryCatch(rbindlist(the.shuff.list, fill = T), error = function(e) NULL)


    ##browser()

    this.all.dat.shuff5.ne = this.all.dat.shuff5[sapply(this.all.dat.shuff5, function(x) !inherits(x, "try-error"))]
    this.all.sh = rbindlist(this.all.dat.shuff5.ne, fill =T)
###
    this.all.sh = sscore(this.all.sh, model = back.model)
    sh.chrom = tryCatch(synergy(binsets = dt2gr(this.all.dat), annotated.binsets = this.all.sh, model = back.model), error = function(e) NULL)
    sh.chrom$fdr = signif(p.adjust(sh.chrom$p, "BH"), 2)
    ##
#####
    this.back.test.dat = sscore(this.back.test.dat, model = back.model)
    theta = back.model$model$theta
    s.chrom = synergy(binsets = dt2gr(this.all.dat), #theta = back.model$model$theta,
                      annotated.binsets = this.chrom.dat, model = back.model)
    s.chrom$fdr = signif(p.adjust(s.chrom$p, "BH"), 2)
    b.chrom = synergy(binsets = dt2gr(back.test), #theta = back.model$model$theta,
                      annotated.binsets = na.omit(this.back.test.dat), model = back.model)
    b.chrom$fdr = signif(p.adjust(b.chrom$p, "BH"), 2)
    ##

    synergy.inter.EP = rbind(s.chrom[, annotation := "chromunity"],
                             b.chrom[, annotation := "random"], 
                             sh.chrom[, annotation := "shuffled"], fill = T)
    if (!is.null(folder)) {
        saveRDS(this.chrom.dat, paste0(folder,'chrom_annotate.rds'))
        saveRDS(this.back.test.dat, paste0(folder,'back_annotate.rds'))
        saveRDS(this.all.sh, paste0(folder,'shuffled_annotate.rds'))
        saveRDS(this.all.dat, paste0(folder,'binsets.rds'))
        saveRDS(back.model, paste0(folder,'back_model.rds'))
        saveRDS(synergy.inter.EP, paste0(folder,'synergy_results.rds'))
    }
    return(synergy.inter.EP)
}




##subsetted to training chr



train_dist_decay_model = function(iterative.registry, iterative.registry.unlist.filterable, contact_matrix_unique, pairwise, concatemers, bins, pair.thresh=50) {

     dist.decay.train = annotate_distance_decay(concatemers,
                                                   bins,
                                                   min.value = pair.thresh,
                                                   window.size=50,
                                                   simplicial.complex=iterative.registry,
                                                   iterative.registry.unlist.filterable = iterative.registry.unlist.filterable,
                                                   pairwise=pairwise,
                                                   return.training=TRUE)



    ####total contacts within window
    
     total.concats = dist.decay.train[(dist.a < 50 | dist.b < 50), .(total.concats = sum(num.concats)), by=c('pair.hashes')]

     dist.decay.train = merge(dist.decay.train, total.concats, by='pair.hashes')
        ##total.concats = dist.decay.test[(dist.a < 50 | dist.b < 50), .(total.concats = sum(num.concats)), by=c('pair.hashes')]
        
###GLM town

     print('training model')
     covariates = c('value.a.ratio','value.b.ratio')
     fmstring = paste('num.concats ~', paste(paste0('log(', covariates, ')', collapse = ' + ')))
        ##fmstring = paste0(fmstring, " + ", "offset(log(total.concats))") this sometimes does 
     fm = formula(fmstring)


        ##browser()
     num.to.sample = 500000
     if(num.to.sample > dim(dist.decay.train)[[1]]) {
         num.to.sample = dim(dist.decay.train)[[1]]
     }
            
     model = glm.nb(formula = fm, data=dist.decay.train[num.concats>0][sample(.N, num.to.sample)], control=glm.control(maxit=500))
     return(model)
}


count_3way_track = function(pairwise.trimmed, dt.concats.sort, all.pairwise, bins) {
    asdf = merge(pairwise.trimmed[i!=j], dt.concats.sort, by.x='i', by.y='binid', allow.cartesian=TRUE)
    asdf.2 = merge(asdf, dt.concats.sort, by.x=c('j','cidi'), by.y=c('binid','cidi'))
    mergerious = merge(asdf.2[, c('id','cidi','i','j')], dt.concats.sort, by='cidi', allow.cartesian=TRUE)
    mergerious = mergerious[binid != i & binid != j]
    dunka = mergerious[, .(num.concats = .N), by=c('id','binid')]
    return(dunka)
}











create_bin_pair_network_efficient = function(genome.wide.dist.decay, all.pairwise=NULL, rr.thresh=0) {

    if(is.null(all.pairwise)){
        all.pairwise = genome.wide.dist.decay[, c('i','j','pair.hashes')] %>% unique(by='pair.hashes')
    }
    all.pairwise.trimmed.2 = all.pairwise %>% copy
    all.pairwise.trimmed.2$i = all.pairwise$j
    all.pairwise.trimmed.2$j = all.pairwise$i
    pairwise.sym = rbind(all.pairwise, all.pairwise.trimmed.2)



    neighbor.1 = merge(genome.wide.dist.decay[relative.risk >= rr.thresh & pval < 0.05, c('i','j','binterrogate','relative.risk','pair.hashes')], pairwise.sym, by.x=c('i','binterrogate'), by.y=c('i','j'))
    neighbor.2 = merge(genome.wide.dist.decay[relative.risk >= rr.thresh & pval < 0.05, c('i','j','binterrogate','relative.risk','pair.hashes')], pairwise.sym, by.x=c('j','binterrogate'), by.y=c('i','j'))




    all.neighbors = rbind(neighbor.1, neighbor.2)


    all.neighbors$pair.hashes.x = all.neighbors$pair.hashes.x %>% factor
    all.neighbors$pair.hashes.y = all.neighbors$pair.hashes.y %>% factor
    G.kant.10 = graph_from_edgelist(as.matrix(all.neighbors[,c('pair.hashes.x','pair.hashes.y')]), directed=TRUE)
    E(G.kant.10)$weight = all.neighbors$relative.risk
    return(G.kant.10)
}


run_interchr_dist_decay = function(concatemers, folder, bins=NULL, pair.thresh=50, chromosome=NULL, training.chr='chr8', resolution=50000, synergy.resolution = 10000, fdr.thresh=0.1, fdr.expansion.thresh=0.25, expansion.cutoff=0.2, num.members=10, genome.to.use="BSgenome.Hsapiens.UCSC.hg38::Hsapiens", seed=198) {
    set.seed(seed)
    
    if(is.null(chromosome)){
        chromosome = c(paste0("chr", c(as.character(1:22), "X")))
    }

    concatemers = concatemers %Q% (seqnames %in% chromosome)
    
    unique_read_idx = concatemers$read_idx %>% unique
    training.idx = sample(unique_read_idx, length(unique_read_idx)/2)
    training.gr = concatemers %Q% (read_idx %in% training.idx)
    test.gr = concatemers %Q% !(read_idx %in% training.idx)

    chrom = interchr_dist_decay_binsets(training.gr, bins=bins, pair.thresh=pair.thresh, folder=folder, training.chr=training.chr, resolution=resolution, rebin.resolution=synergy.resolution, fdr.thresh=fdr.thresh, fdr.expansion.thresh=fdr.expansion.thresh, expansion.cutoff=expansion.cutoff, num.members=num.members, genome.to.use=genome.to.use)
    
    dir.create(folder)
    saveRDS(chrom, paste0(folder, '/chromunity_results.rds'))
    saveRDS(training.gr, paste0(folder, '/training_concatemers.rds'))

    chid.to.test = (chrom$chid %>% unique)
    syn.interchr = evaluate_synergy_experimental(chrom, test.gr, chid.to.test, folder=folder, chromosome=chromosome, genome.to.use=genome.to.use, resolution=synergy.resolution)
    syns = aggregate_synergy_results_singledir(folder, strict.check=FALSE)
    return(syns)
}


aggregate_synergy_results_singledir = function(dir, strict.check=FALSE, fdr.thresh=0.1) {

    synergy.chunk = readRDS(paste0(dir, '/synergy_results.rds'))
    synergy.chunk = synergy.chunk[!is.na(fdr)]
    binsets = readRDS(paste0(dir, '/binsets.rds'))

    annotate = readRDS(paste0(dir, '/chrom_annotate.rds'))
    if(strict.check)
        synergy.chunk = symmetry_check(annotate, synergy.chunk, binsets)

    binsets$bid = factor(binsets$bid)
    synsets = merge(binsets, synergy.chunk[annotation=='chromunity', c('bid','fdr')], by='bid')
    synsets = synsets[fdr<fdr.thresh]
    return(list(synergy.chunk, dt2gr(synsets)))
}


shuffle_concatemers_spiked = function(concatemers, contact_matrix, spike.set, strength=0.1, noise.parameter = 10000, resolution=10000, trim=TRUE, interval.width=1000) {
    A = contact_matrix$mat %>% as.matrix
    rownames(A) <- NULL
    colnames(A) <- NULL
    A[cbind(1:nrow(A), 1:nrow(A))] = 0
    A = A + t(A)
    An = A 

    An = round(1+10*An/min(An[An>0]))  
    edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
     ##
    G = graph.edgelist(edges[, cbind(row, col)])
    
    concats.binned = bin_concatemers(concatemers, contact_matrix$gr)
    concats.dt = concats.binned
    concats.dt = unique(concats.dt, by=c('binid','cidi')) ###dedupe!!

    concat.counts = concats.dt[, new.count := .N, by='read_idx']


    card = unique(concat.counts[, .(read_idx, new.count)])

    this.steps = sum(card$new.count)

    row.labels = 1:dim(An)[1]
    start.index = row.labels[rowSums(An)>1]
    
    RW = random_walk(G, start = start.index, steps = sum(card$new.count)) %>% as.numeric

    out = contact_matrix$gr[RW] %>% gr2dt()
    ##browser()
    
    out = out[1:sum(card$new.count)]
    out$read_idx = card[, rep(read_idx, new.count)] 

    out[, cardinality := .N, by='read_idx']
##    out = gr2dt(dt2gr(out) - ()) ###Setting width of each monomer to be arbitrarily 1000
    ##browser()
    
    ##spike.concatemers = spike.set
    ##spikable.concats = out[cardinality > 1 & binid %in% spike.set$binid]$read_idx %>% unique  ##no concatemer which only hit a bin-set
    ##concats.to.spike = sample(spikable.concats, number.to.add)

    ##browser()



    ###Split by whatever random binsets you want to add here

    dt = gr2dt(spike.set %&% contact_matrix$gr)[, count := .N, by='bid']
    spike.set = dt2gr(dt[count > 2])
    
    spike.set$binid = gr.match(spike.set, contact_matrix$gr)
    spike.set$chid = spike.set$chid %>% as.character
    binsets.to.add = split(spike.set, spike.set$chid)
    out[, c('query.id','tile.id') := NULL]

###We are going to choose random position within bin to begin monomer from
###And then create fixed width monomer
    
    out$rowid = 1:dim(out)[[1]]
    out[, new.start.interval := sample(resolution - interval.width, 1), by='rowid']
    out[, start := start + new.start.interval]
    out[, end := start + interval.width]
    

    ##spike.set.gr = binsets.to.add$rand_15

    spikes = pbmclapply(binsets.to.add, mc.cores = 5, function(spike.set.gr) {

        ##number.to.add = choose(spike.set.gr[1]$cardinality, 3) * strength
##        spike.set.gr = binsets.to.add$rand_10

        binset.rand = spike.set.gr %>% gr2dt
        binset.rand$bid=1
        binset.rand$index = 1:dim(binset.rand)[[1]]
        power = binset.rand[, powerset(index, 3, 3), by=bid]

        
        num.concats = merge(binset.rand, out, by='binid')$read_idx %>% unique %>% length
        number.to.add = num.concats * strength
        print(number.to.add)
        
        power[, size := .N, by=c('setid','bid')]
        size3 = unique(power[size==3], by=c('bid','setid'))
        size3.samp = sample(size3$setid, number.to.add, replace=TRUE)
        samp = data.table(size3.samp)
        colnames(samp) = 'setid'
        samp.binned = merge(samp, power, by='setid', allow.cartesian=TRUE)
        spike.concats = binset.rand[samp.binned$item]


        

        ###choose spike location randomly from bin-set given
        spike.concats$rowid = 1:dim(spike.concats)[[1]]
        spike.concats[, new.start.interval := sample(resolution - interval.width, 1), by='rowid']
        spike.concats[, start := start + new.start.interval]
        spike.concats[, end := start + interval.width]

        spike.concats.gr = dt2gr(spike.concats)

        dt.spike = gr2dt(spike.concats.gr)
        dt.spike$spike.id = rep(1:dim(samp)[[1]], each=3)
        dt.spike.2 = merge(out, dt.spike, by='binid', allow.cartesian=TRUE)
        sample.spike = dt.spike.2[, .(read_idx = sample(read_idx, 5)), by='spike.id'] ###almost guarantees we don't double select a concatemer
        sample.spike = sample.spike %>% unique(by='read_idx') %>% unique(by='spike.id')
        sample.spike = merge(sample.spike, dt.spike.2[, c('spike.id','read_idx','binid')], by=c('spike.id','read_idx'))
        sample.spike = sample.spike %>% unique(by='read_idx') %>% unique(by='spike.id')
        

        dt.spike[, binid := NULL]
        merge.spike = merge(dt.spike, sample.spike, by='spike.id', allow.cartesian=TRUE)
        merge.spike$cardinality = 3
        
        dt = merge.spike[, c('seqnames','start','end','strand','width','binid','read_idx','cardinality')]


        jitter = rnorm(dim(dt)[[1]], sd=noise.parameter)
        dt$jitter = jitter
        dt[, start := start + jitter]
        dt[, end := end + jitter]
        dt[, jitter := NULL]
        return(dt)
    })
    spikes = spikes %>% rbindlist

####Remove the monomer which allowed concatemer to become spikable in the first place...

    
    out$rowid = 1:dim(out)[[1]]
    out.to.drop = merge(unique(spikes, by='read_idx'), out, by=c('read_idx','binid'))
    out = out[!(rowid %in% out.to.drop$rowid)]
    
    
    ##browser()
    
    spikes$binid = gr.match(dt2gr(spikes), contact_matrix$gr)
    
    spike.out = rbind(out, spikes, fill=TRUE)

    spike.out$cid = spike.out$read_idx
    spike.out[, chid := NULL]

    if(trim == TRUE){
        spiked.concats = spikes[read_idx %in% spikes$read_idx]
        added.pairwise.contacts = cocount(dt2gr(spiked.concats), bins=contact_matrix$gr, by='read_idx')
        dt = added.pairwise.contacts$dat
        to.remove = dt[(i %in% spikes$binid) & (j %in% spikes$binid)][i != j]

        dt.i = merge(out[, c('read_idx','binid')], to.remove, by.x='binid', by.y='i', allow.cartesian=TRUE)
        dt.j = merge(out[, c('read_idx','binid')], dt.i, by.x = c('binid','read_idx'), by.y=c('j','read_idx'))

        dt.j[, count := 1:.N, by='id']
        cid.bins = dt.j[binid != binid.y]
        cid.bins = cid.bins[count <= value]
        
        cid.bins.2 = cid.bins %>% copy

        cid.bins.2$binid = cid.bins$binid.y
        cid.bins.2$binid.y = cid.bins$binid

        dt.bins = rbind(cid.bins, cid.bins.2)[, c('read_idx','binid')] %>% unique(by=c('read_idx','binid'))

        ##browser()
        out$rowid = 1:dim(out)[[1]]
        remove.out = merge(out, dt.bins, by=c('read_idx','binid'))
        remove.out = unique(remove.out, by='read_idx')
        out.trimmed = out[!(rowid %in% remove.out$rowid)]
        out.trimmed[, rowid := NULL]
        
        spike.out.trimmed = rbind(out.trimmed, spikes, fill=TRUE)
        spike.out = spike.out.trimmed
    }
    spike.out[, end := start + interval.width] 
    return(spike.out)
}






analyze_concats_final = function(concatemers, basecalled.bases, num.sequences) {
    stats.table = data.table()

    numchunks = 200
    unique.cids = concatemers$read_idx %>% unique



###Analyze concatemers


    group.cids = data.table(unique.cids, group=ceiling(runif(length(unique.cids))*numchunks))
    ##dt.concats = training.concatemers %>% gr2dt

    cid.split = split(group.cids, by='group')
    cmessage('Reducing concatemers')
    
    reduced.concats.list = pbmclapply(cid.split, mc.cores=10, function(chunk) {
        concat.sub = concatemers %Q% (read_idx %in% chunk$unique.cids)
        reduced.concats = grl.reduce(split(concat.sub, concat.sub$read_idx), pad=100)
        reduced.concats = unlist(reduced.concats)
        reduced.concats$read_idx = names(reduced.concats) %>% as.numeric
        reduced.concats = gr2dt(reduced.concats)
        return(reduced.concats)
    })


    concats.dt = rbindlist(reduced.concats.list)
    concats.dt[, cid := read_idx]

    ##bins = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), width=10000)
    ##concatemers$binid = gr.match(concatemers, bins, max.slice = 1e6, mc.cores = 5, verbose = FALSE)
    ##concatemers = concatemers %Q% (!is.na(binid))

    stats.table = data.table()


    ##concats.dt = gr2dt(concatemers)
    concats.dt[, cidi := as.integer(cid)]

    concats.dt[, bincount := .N, by='read_idx']
    concats.dt[, contacts.count := choose(bincount, 2), by='read_idx']
    undeduped.contact.count = sum(unique(concats.dt, by='cidi')$contacts.count)
    stats.table$undeduped.contact.count = undeduped.contact.count

#####Percent undigested


    ##deduped bin hits
    ##concats.dt = unique(concats.dt, by=c('cidi','binid'))
    concats.dt[, bincount := .N, by='read_idx']
    concats.dt[, contacts.count := choose(bincount, 2), by='read_idx']
    deduped.contact.count = sum(unique(concats.dt, by='cidi')$contacts.count)
    stats.table$contact.count = deduped.contact.count


    new.concats.dt = concats.dt
    new.concats.dt[, bincount := .N, by='read_idx']
    unique.concats = unique(new.concats.dt, by='read_idx')

###percent undigested

    percent.undigested = sum(unique.concats$bincount == 1) / dim(unique.concats)[[1]]
    stats.table$percent.undigested = percent.undigested

    ##number of chromosomes a concatemer hits
    new.concats.dt[, num.chr := length(unique(seqnames)), by='read_idx']

    ##number of times a concatemer hits the given chromosome
    new.concats.dt[, chr.count := .N, by=c('read_idx','seqnames')]

    unique.chr.concats = unique(new.concats.dt, by=c('read_idx','seqnames'))
    unique.chr.concats[, intra.chr.contacts := choose(chr.count, 2)]

    unique.chr.concats[, contacts.count := choose(bincount, 2), by='read_idx']
    unique.concats = unique(unique.chr.concats, by=c('read_idx'))


    intra.vs.inter = unique.chr.concats[, .(intra.contacts = sum(intra.chr.contacts), contacts.count), by='read_idx']
    intra.vs.inter = unique(intra.vs.inter, by='read_idx')

    total.inter.contacts = sum(intra.vs.inter$contacts.count) - sum(intra.vs.inter$intra.contacts)
    intra.vs.inter = unique(intra.vs.inter, by='read_idx')

    ##percent cis
    percent.cis = sum(intra.vs.inter$intra.contacts) / sum(intra.vs.inter$contacts.count)
    stats.table$percent.cis = percent.cis

    contacts.per.gb = deduped.contact.count / ((basecalled.bases) / (10^9))
    stats.table$contacts.per.gb = contacts.per.gb

    ##contact distribution
    contact.order.dist = data.table()
    contact.order.dist$lessone = sum(unique.concats$bincount == 1)

    total.reads = num.sequences


    contact.order.dist$multiway = sum(unique.concats$bincount > 1)
    contact.order.dist$two = sum(unique.concats$bincount == 2) / total.reads
    contact.order.dist$three = sum(unique.concats$bincount == 3)  / total.reads
    contact.order.dist$four = sum(unique.concats$bincount == 4) / total.reads
    contact.order.dist$fivetosix = sum((unique.concats$bincount == 5) | (unique.concats$bincount == 6)) / total.reads
    contact.order.dist$seventoeleven = sum((unique.concats$bincount >= 7) & (unique.concats$bincount <= 11)) / total.reads
    contact.order.dist$twelvetotwentyone = sum((unique.concats$bincount >= 12) & (unique.concats$bincount <= 21)) / total.reads
    contact.order.dist$twentytwotofortynine = sum((unique.concats$bincount >= 22) & (unique.concats$bincount <= 49)) / total.reads
    contact.order.dist$greaterthan50 = sum(unique.concats$bincount >= 50) / total.reads

    stats.table = cbind(stats.table, contact.order.dist)

    distribution = unique.concats[, .N, by=bincount]
    return(list(stats.table, distribution))
    

}




plot_binsets = function(binsets.gr, name='binsets', height=NULL) {
##    binsets.gr$chid = binsets.gr$bid
    binsets.gr$label = NULL
    if(is.null(height)) height=10
    return(gTrack(split(binsets.gr, binsets.gr$chid) %>% unname, name=name, height=height))
}

##Projects/testing/


