## This file will contain functions for running the distance-decay model on Pore-C data

#' @name interchr_dist_decay_binsets
#' @description
#'
#' This function performs the bin-set nomination procedure nominates bin-sets with the distance decay higher order contact modeling procedure
#' 
#'  
#' @param concatemers GRanges of monomers with fields seqnames, start, end, and $cid specifying concatemer id, which will be counted across each binset
#' @param resolution integer specifying the bin width to use for the distance decay model
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
#' @param pairwise.trimmed data table specifying pairs to subset higher order analysis to
#' @param num.to.sample number of pairs to sample for training distance decay model.
#' @param rebin.resolution Resolution to rebin chromunities to for Synergy analysis
#' @param pairs.per.chunk pairs i & j per chunk in parallelization to analyze higher order contacts. overridden by numchunks
#' @param mc.cores cores to use parallelization process. Very mem
#' @param mc.cores integer how many cores to parallelize for distance decay calculation. This is very memory intensive so defaults to 2.
#' @param compressed.representation saves distance decay model output with fewer fields if this is enabled
#' @param numchunks number of computational chunks to create for higher order contact analysis
#' @param genome.to.use Genome model to use for analysis, this argument is passed to hg_seqlengths in the genome argument
#' @param monomer.merge.distance bp distance between which to merge monomers into a single monomer.
#' @author Jameson Orvis
#' @export
#' @return 

interchr_dist_decay_binsets = function(concatemers, resolution=50000, bins=NULL, interchromosomal.distance = 1e8, training.chr = 'chr8', pair.thresh=50, mask.bad.regions = TRUE, fdr.thresh=0.1, fdr.expansion.thresh=0.25, expansion.cutoff=0.2, num.members=10, folder=NULL, chromosome=NULL, model=NULL, pairwise.trimmed=NULL, num.to.sample=250000, rebin.resolution=10000, pairs.per.chunk=1000, mc.cores=2, compressed.representation=FALSE, numchunks=200, genome.to.use = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens", monomer.merge.distance = 100) {

    ##all chromosomes by default
    if(is.null(chromosome)){
        chromosome = c(paste0("chr", c(as.character(1:22), "X", "Y")))
    }
    
    if(is.null(bins)){
        bins = gr.tile(hg_seqlengths(genome = genome.to.use), width=resolution) %Q% (seqnames %in% chromosome)
    }
    bins$binid = 1:length(bins)

    ##Creates virtual pairwise contacts using cocount
    ##This version of cocount removes duplicate contacts from monomers overlapping genomic bins more than once
    contact_matrix_unique = dedupe_cocount(concatemers, bins = bins, by = 'read_idx')
    all.pairwise = contact_matrix_unique$dat
    all.pairwise$id = 1:dim(all.pairwise)[[1]]
    colnames(all.pairwise)[3] = 'pair.value'
    
    concatemers$cid = concatemers$read_idx
    unique.cids = concatemers$read_idx %>% unique

    group.cids = data.table(unique.cids, group=ceiling(runif(length(unique.cids))*numchunks))

    cid.split = split(group.cids, by='group')
    cmessage('Reducing concatemers')

    ##Pre-processing step to merge monomers which are within 100 bp of each other 
    reduced.concats.list = pbmclapply(cid.split, mc.cores=5, function(chunk) {
        concat.sub = concatemers %Q% (read_idx %in% chunk$unique.cids)
        reduced.concats = grl.reduce(split(concat.sub, concat.sub$read_idx), pad=monomer.merge.distance)
        reduced.concats = unlist(reduced.concats)
        reduced.concats$read_idx = names(reduced.concats) %>% as.numeric
        reduced.concats = gr2dt(reduced.concats)
        return(reduced.concats)
    })

    dt.concats.reduced = rbindlist(reduced.concats.list)

    reduced.concatemers = dt2gr(dt.concats.reduced)

    binned.concats = bin_concatemers(reduced.concatemers, bins, max.slice=1e6, mc.cores=5)
    
    if(!is.null(folder)){
        saveRDS(binned.concats, 'binned_concatemers.rds')
    }
    

    dt.concats = unique(binned.concats[, c('cidi','binid')], by=c('cidi','binid'))
    dt.concats.sort = dt.concats[order(binid, cidi)]
    dt.concats.sort[, count := .N, by='cidi']
    
    ##choose subset of bin-pairs S by thresholding
    colnames(all.pairwise)[4] = 'pair.hashes'

    ##allow option of passing in pairwise trimmed manually 
    if(is.null(pairwise.trimmed)) {
        pairwise.trimmed = all.pairwise[pair.value >= pair.thresh]
        pairwise.trimmed[, dist := j-i]
        pairwise.trimmed = pairwise.trimmed[dist > 1] ##ignore pairs directly adjacent
    }

    ##mask genomic regions overlapping centromeres and telomeres
    if(mask.bad.regions==TRUE) {
        ##chromosome = c(paste0("chr", c(as.character(1:22), "X")))
        if(length(chromosome) > 1) {
            all.bad = pbmclapply(chromosome, function(chr) {
                bad.gr = muffle(load_bad_regions(chr, genome.to.use = genome.to.use))
                bad.dt = gr2dt(bad.gr)
                return(bad.dt)
            })
        } else {
            all.bad = muffle(load_bad_region(chromosome, genome.to.use = genome.to.use))
        }
        this.bad = rbindlist(all.bad) %>% dt2gr
        bad.bins = this.bad %*% bins
        pairwise.trimmed = pairwise.trimmed[!(i %in% bad.bins$binid)]
        pairwise.trimmed = pairwise.trimmed[!(j %in% bad.bins$binid)]
    }

    unique.pairs = pairwise.trimmed$pair.hashes %>% unique

    #pre-processing done, now train the model
    ##trains distance decay model using subset of higher order contacts in one chromosome. 
    if(is.null(model)){
        model = train_dist_decay_model_nozero(dt.concats.sort, pairwise.trimmed, bins %Q% (seqnames==training.chr), all.pairwise=all.pairwise, num.to.sample=num.to.sample)
    }

    if(is.null(numchunks))
        numchunks = ceiling(length(unique.pairs) / pairs.per.chunk)  ###Will attempt to process 100 pairs per chunk

    pair.splitting = data.table(unique.pairs, group=ceiling(runif(length(unique.pairs))*numchunks))
    pairwise.trimmed$group = pair.splitting$group

    pairwise.trimmed$id = pairwise.trimmed$pair.hashes
    all.pairwise$id = all.pairwise$pair.hashes

    pairwise.chunks = split(pairwise.trimmed, by='group')


###The most computationally expensive step of this process, analyzes higher order contacts in parallel.
    ##browser()
    scored.chunks = pbmclapply(pairwise.chunks, mc.cores = mc.cores, function(pairwise.chunk) {
        annot.chunk = count_3way_contacts(pairwise.chunk, dt.concats.sort, all.pairwise, bins, interchromosomal.distance = interchromosomal.distance)
        if(mask.bad.regions == TRUE){
            annot.chunk = annot.chunk[!(binterrogate %in% bad.bins$binid)]
        }
        dt.small.scored = score_distance_decay(annot.chunk, model, mode='poisson')

        if(compressed.representation==TRUE){
            dt.small.scored = dt.small.scored[, c('pair.hashes','binterrogate','pval','i','j','num.concats','num.concats.pred')]
        }
        return(dt.small.scored)
    })

    ##browser()
    genome.wide.dist.decay = rbindlist(scored.chunks)

    genome.wide.dist.decay[, relative.risk := log2(num.concats/num.concats.pred)] ###relative.risk is a misnomer, this should probably be relabeled as observed/expected.
    genome.wide.dist.decay$fdr = signif(p.adjust(genome.wide.dist.decay$pval, "BH"), 2)
    trimmed.dist.decay = genome.wide.dist.decay[fdr<fdr.thresh]
    
    
    if(dim(trimmed.dist.decay)[[1]] == 0) {
        stop('Error: No significant three-way contacts discovered with benjamini-hochberg multiple corrections. Try a lower resolution.')
    }
        
    print(trimmed.dist.decay)
    print(all.pairwise)
    
    bin.pair.network = create_bin_pair_network_efficient(trimmed.dist.decay, all.pairwise, rr.thresh=0) ###creates bin-pair network 

    if(!is.null(folder)) {
        if(compressed.representation==TRUE) { 
            saveRDS(genome.wide.dist.decay[pval < fdr.expansion.thresh], paste0(folder,'dist_decay_archive.rds'))
        } else {
            saveRDS(list(genome.wide.dist.decay, bins), paste0(folder,'dist_decay_archive.rds'))
        }
     }

     ###creates chromunity object from bin-pair network
    chrom = derive_binsets_from_network(bin.pair.network, pairwise=all.pairwise, binned.concats=binned.concats, bins=bins, rr.thresh=0, dist.decay.test=genome.wide.dist.decay[fdr < fdr.expansion.thresh], all.pairs.tested=NULL, pairwise.trimmed=pairwise.trimmed, expansion.cutoff = expansion.cutoff, fdr.thresh = fdr.thresh, num.members=num.members, rebin.resolution=rebin.resolution) 
    
    return(chrom)
}

train_dist_decay_model_nozero = function(dt.concats.sort, pairwise.trimmed, bins, all.pairwise, numchunks=NULL, num.to.sample=250000, pairs.to.sample = 10000, mode='poisson', pairs.per.chunk=100){

    unique.pairs = pairwise.trimmed$pair.hashes %>% unique

    #subsetting to make this more efficient
    if(length(unique.pairs) > pairs.to.sample){
        unique.pairs = unique.pairs[sample(pairs.to.sample)]
        pairwise.trimmed = pairwise.trimmed[pair.hashes %in% unique.pairs]
    }

    if(is.null(numchunks))
        numchunks = ceiling(length(unique.pairs) / pairs.per.chunk)  ###Will attempt to process 100 pairs per chunk

    pair.splitting = data.table(unique.pairs, group=ceiling(runif(length(unique.pairs))*numchunks))
    pairwise.trimmed$group = pair.splitting$group

    pairwise.trimmed$id = pairwise.trimmed$pair.hashes
    all.pairwise$id = all.pairwise$pair.hashes

    pairwise.chunks = split(pairwise.trimmed, by='group')

    ##browser()

    scored.chunks = pbmclapply(pairwise.chunks, mc.cores = 5, function(pairwise.chunk) {
        annot.chunk = count_3way_contacts(pairwise.chunk, dt.concats.sort, all.pairwise, bins)
        return(annot.chunk)
    })
    dist.decay.train = rbindlist(scored.chunks)
    
    print('training model')
    covariates = c('value.a.ratio','value.b.ratio')
    fmstring = paste('num.concats ~', paste(paste0('log(', covariates, ')', collapse = ' + ')))
    ##fmstring = paste0(fmstring, " + ", "offset(log(total.concats))") ##this sometimes does 

    fm = formula(fmstring)
    
    if(num.to.sample > dim(dist.decay.train[dist.a <= 50 & dist.b <= 50])[[1]]) {
        train.subset = dist.decay.train[dist.a <= 50 & dist.b <= 50]
    } else { 
        close.subset = dist.decay.train[dist.a <= 50 & dist.b <= 50][sample(.N, num.to.sample)]
        far.subset = dist.decay.train[dist.a > 50 & dist.b > 50][sample(.N, num.to.sample/4)]  ##more readable/straightforward
        train.subset = rbind(close.subset, far.subset)
    }

    if(mode=='poisson'){
        model = glm(formula = fm, data=train.subset[, c('num.concats','value.a.ratio','value.b.ratio')], control=glm.control(maxit=500), family='poisson')
    } else {
        model = glm.nb(formula = fm, data=train.subset[, c('num.concats','value.a.ratio','value.b.ratio')], control=glm.control(maxit=500))
    }
    return(model)
}

bin_concatemers = function(concatemers, bins, max.slice = 1e6, mc.cores=5, verbose=TRUE, hyperedge.thresh=NULL) {
    concatemers$binid = gr.match(concatemers, bins, max.slice = max.slice, mc.cores =  mc.cores, verbose = verbose)

    concatemers$cid = concatemers$read_idx
    ## maybe NA need to be removed

    concatemers = concatemers %Q% (!is.na(binid))
    reads = as.data.table(concatemers)[, `:=`(count, .N), by = cid]    
    reads[, cidi := as.integer(cid)]
    return(reads)
}

score_distance_decay = function(dt.small.model, model, mode='poisson'){
    dt.small.model$num.concats.pred = (predict(model, type = "response", newdata = dt.small.model))

    # TO DO: check that "lower.tail=F" makes sense in this case... 
    if (mode=='poisson'){
        pval = ppois(dt.small.model$num.concats -1, lambda = dt.small.model$num.concats.pred, lower.tail = F)
        pval.right = ppois(dt.small.model$num.concats, lambda = dt.small.model$num.concats.pred, lower.tail = F)
    }else if (mode=='nbinom'){
        pval = pnbinom(dt.small.model$num.concats -1, mu = dt.small.model$num.concats.pred, size=model$theta, lower.tail = F)
        pval.right = pnbinom(dt.small.model$num.concats, mu = dt.small.model$num.concats.pred, size=model$theta, lower.tail = F)
    }

    pval.right = ifelse(is.na(pval.right), 1, pval.right)
    pval = ifelse(is.na(pval), 1, pval)
    dt.small.model$pval = runif(nrow(dt.small.model), min = pval.right, max = pval)
    
    dt.small.model[, enrichment := num.concats / num.concats.pred]
    return(dt.small.model)
}


count_3way_contacts = function(pairwise.trimmed, dt.concats.sort, all.pairwise, bins, interchromosomal.distance = 1e8) {

    ###makes unlisting convenient
    pairwise.trimmed$agg = do.call(Map, c(f = c, pairwise.trimmed[, c('i','j')]))
    
    ###Performing these two joins will give you set of all concatemers which overlap both of i & j
    concats.hitting.i = merge.data.table(pairwise.trimmed[i!=j], dt.concats.sort, by.x='i', by.y='binid', allow.cartesian=TRUE)
    concats.hitting.ij = merge.data.table(concats.hitting.i, dt.concats.sort, by.x=c('j','cidi'), by.y=c('binid','cidi'))


    ###We do not know what else those concatemers are overlapping, join back to dt.concats sort to get the rest of concatemers
    bins.hit.by.ij.concats = merge.data.table(concats.hitting.ij[, c('id','cidi','i','j','agg')], dt.concats.sort, by='cidi', allow.cartesian=TRUE)[binid!=i & binid != j]

    ###Counting three way contacts becomes a matter of counting the number of times binid appears with respect to each i & j
    threeway.contact.counts= bins.hit.by.ij.concats[, .(num.concats = .N, agg), by=c('id','i','j','binid')]
    colnames(threeway.contact.counts)[4]='V1'
    threeway.contact.counts$agg = do.call(Map, c(f = c, threeway.contact.counts[, c('i','j')]))

    ####Now the problem becomes: can we calculate the pairwise contacts between i & the third bin given in V1.
    ####To make this easy with a join we create a new line for i & k and j & k
    dt.sub = threeway.contact.counts[, .(sub.bin = unlist(agg)), by=c('i','j','V1','id','num.concats')]    
    dt.sub[sub.bin < V1, c('sub.bin','V1') := .(V1, sub.bin)] ##Swap these to join with pairwise contact matrix
    dt.sub = merge.data.table(dt.sub, all.pairwise[, c('i','j','pair.value')], by.x=c('V1','sub.bin'), by.y=c('i','j'), all.x=TRUE)
    dt.sub[is.na(pair.value), pair.value := 0]

    ##Find which of sub.bin and binid are not equal to i & j, call this "binterrogate", this is the third bin for whose higher order contacts
    ##we are "interrogating".

    #maybe some refactorring is in order down here
    dt.sub = dt.sub[V1 != sub.bin]
    dt.sub$V1.isinter = !(dt.sub[, V1 == i] | dt.sub[, V1 == j])
    dt.sub$sub.isinter = !(dt.sub[, sub.bin == i] | dt.sub[, sub.bin == j])
    dt.sub$binterrogate = 0
    dt.sub[V1.isinter==TRUE, binterrogate := V1]
    dt.sub[sub.isinter==TRUE, binterrogate := sub.bin]
    dt.sub = dt.sub[binterrogate!=0]

    
    dt.sub$pair.hashes = dt.sub$id
    dt.sub[, sum.pairwise.contacts := sum(pair.value), by=c('pair.hashes','binterrogate')]
    dt.sub[, dist.i := abs(binterrogate - i)]
    dt.sub[, dist.j := abs(binterrogate - j)]
    bins.gr = bins
    bins.gr$binid = 1:length(bins.gr)
    bins.dt = gr2dt(bins.gr)
    setkey(bins.dt, 'binid')
    dt.sub$chr.i = bins.dt[dt.sub$i]$seqnames
    dt.sub$chr.j = bins.dt[dt.sub$j]$seqnames
    dt.sub$chr.binterrogate = bins.dt[dt.sub$binterrogate]$seqnames

    ####interchromosomal distance
    ####we assign interchromosomal distances a value manually passed in

    resolution = median(width(bins))
    inter.dist = interchromosomal.distance / resolution
    dt.sub[chr.i != chr.binterrogate, dist.i := inter.dist]
    dt.sub[chr.j != chr.binterrogate, dist.j := inter.dist]
    dt.sub[, chr.i := NULL]
    dt.sub[, chr.j := NULL]
    dt.sub[, chr.V1 := NULL]
    dt.sub[, diff := j-i]
    dt.sub = dt.sub[diff > 1]
    dt.sub[, pair.value := pair.value + 1] ###for log covariate purposes
    print('calculating close and far pairwise contact values')
    ###pair.value corresponds to pairwise contacts value between V1 and sub.bin
    ###So here the problem is not to determine which of these pairwise contacts should be considered "close" and "far"
    ###For the case where binterrogate is to the left of both i & j (meaning binterrogate < i) this is fairly easy
    dt.sub[binterrogate < i & j == sub.bin, value.b := pair.value] ##CLOSER BIN
    dt.sub[binterrogate < i & i == sub.bin, value.a := pair.value] ##FARTHER BIN
    ###For the case where binterrogate is to the right of both i & j (meaning binterrogate > j) this is also pretty straightforward
    dt.sub[binterrogate > j & i == V1, value.a := pair.value]
    dt.sub[binterrogate > j & j == V1, value.b := pair.value]
    ###The tricky edge case is when binterrogate is between i & j 
    
    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) < (j - binterrogate)) & i==V1, value.a := pair.value]
    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) < (j - binterrogate)) & j==sub.bin, value.b := pair.value]
    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) >= (j - binterrogate)) & j==sub.bin, value.a := pair.value]
    dt.sub[binterrogate > i & binterrogate < j & ((binterrogate - i) >= (j - binterrogate)) & i==V1, value.b := pair.value]

    dt.sub.unique = unique(dt.sub, by=c('pair.hashes','binterrogate'))

    dt.sub.unique[, value.a := NULL]
    dt.sub.unique[, value.b := NULL]

    a.value = dt.sub[!is.na(value.a), c('pair.hashes','value.a','binterrogate')] %>% unique(by=c('pair.hashes','binterrogate'))
    b.value = dt.sub[!is.na(value.b), c('pair.hashes','value.b','binterrogate')] %>% unique(by=c('pair.hashes','binterrogate'))

    dt.sub.unique = merge.data.table(dt.sub.unique, a.value, by=c('pair.hashes','binterrogate'))
    dt.sub.unique = merge.data.table(dt.sub.unique, b.value, by=c('pair.hashes','binterrogate'))

    dt.sub.unique[, dist.a := min(dist.i, dist.j), by=c('pair.hashes','binterrogate')]
    dt.sub.unique[, dist.b := max(dist.i, dist.j), by=c('pair.hashes','binterrogate')]
    dt.small = dt.sub.unique[, c('pair.hashes','num.concats','pair.value','binterrogate','sum.pairwise.contacts','dist.a','dist.b','value.a','value.b')]
    dt.small = unique(dt.small, by=c('pair.hashes','binterrogate'))

    dt.small[, value.a.ratio := value.a / dist.a]
    dt.small[, value.b.ratio := value.b / dist.b]
    dt.small = dt.small[dist.a > 1 & dist.b > 1]

    dt.small = merge.data.table(dt.small, pairwise.trimmed[, c('i','j','id')], by.x='pair.hashes', by.y='id')

    return(dt.small)
}

sparse_n_tensor = function(dt.concats, numchunks=50, cardinality=3, cores=5) {
    unique_cidi = dt.concats$cidi %>% unique
    ucidl = split(unique_cidi, ceiling(runif(length(unique_cidi))*numchunks))

    dt.concats.sort = dt.concats[order(binid, cidi)]
    dt.concats.sort[, count := .N, by='cidi']
    dt.concats.sort = dt.concats.sort[count >= cardinality]

    combinatorics = mclapply(ucidl, mc.cores = cores, function(cidis) {
        combn.chunk = dt.concats.sort[cidi %in% cidis, .(combn(binid,cardinality,simplify=FALSE)), by='cidi']
        combn.chunk$hashes = combn.chunk$V1 %>% paste0
        combn.chunk.count = combn.chunk[, .(count = .N, V1), by='hashes']
        combns = unique(combn.chunk.count, by='hashes')
        return(combns)
    })

    dt = rbindlist(combinatorics)
    return(dt[,.(coords=V1,count)])
}

