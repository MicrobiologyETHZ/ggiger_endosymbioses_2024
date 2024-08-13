count_snps <- function(seq){
    to_count <- c("(?<![\\+-][0-9])A", "(?<![\\+-][0-9])C", "(?<![\\+-][0-9])G", "(?<![\\+-][0-9])T", "\\*", "\\+[0-9][ACGTN]")
    counts <- list()
    for(n in to_count){
        if(grepl(n, toupper(seq), perl=T)){
            counts[[n]] <- length(gregexpr(n, toupper(seq), perl=T)[[1]])
        }else{
            counts[[n]] <- 0
        }
    }
    names(counts) <- c("A", "C", "G", "T", "del", "ins")
    return(unlist(counts))
}

summarise_pile <- function(contig, start, end, samples){
    l = 1+end-start
    paths <- paste("scratch/", samples, "/snippy/", samples, "_ref.bam", sep="")
    write(paths, "samples.txt", ncolumns=1)
    results <- system(paste("samtools mpileup -b samples.txt -r ", contig, ":", start, "-", end, " -aa", sep=""), intern=T)
    results <- strsplit(results, "\t")
    for(x in 1:length(results)){
       names(results[[x]]) <- c("contig", "pos", "ref", paste(rep(samples, each=3), c("cov", "bases", "qual"), sep="_"))
    }
    
    freq_tables <- list()
    for(pos in 1:length(results)){
        freqs <- lapply(results[[pos]][2+1:length(samples)*3], count_snps)
        freq_table <- do.call(rbind, freqs)
        freq_table <- as.data.frame(freq_table)
        freq_table$cov <- results[[pos]][1+1:length(samples)*3]
        rownames(freq_table) <- samples
        freq_tables[[paste(contig, start+pos-1)]] <- freq_table
    }
    return(freq_tables)
}

detail_pile <- function(contig, start, end, samples){
    l = 1+end-start
    paths <- paste("scratch/", samples, "/snippy/", samples, "_ref.bam", sep="")
    write(paths, "samples.txt", ncolumns=1)
    results <- system(paste("samtools mpileup -b samples.txt -r ", contig, ":", start, "-", end, " -aa", sep=""), intern=T)
    results <- strsplit(results, "\t")
    for(x in 1:length(results)){
       names(results[[x]]) <- c("contig", "pos", "ref", paste(rep(samples, each=3), c("cov", "bases", "qual"), sep="_"))
    }
    return(results)
}
