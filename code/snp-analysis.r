library(dplyr)

row_equals <- function(row1, row2){
   return(all(row1==row2, na.rm=T))
}

row_in <- function(row, df){
    for(i in 1:nrow(df)){
        if(row_equals(row, df[i,])){
            return(TRUE)
        }
    }
    return(FALSE)
}

row_count <- function(row, df){
    count = 0
    for(i in 1:nrow(df)){
        if(row_equals(row, df[i,])){
            count = count+1
        }
    }
    return(count)
}

unique_minusc <- function(df){
    udf <- unique(df)
    counts = c()
    for(i in 1:nrow(udf)){
        counts[i] = row_count(udf[i,], df)
    }
    udf <- cbind(udf, count=counts)
    return(udf)
}

sample_meta <- read.table("data/sample_meta.txt", header=T)
bac_meta <- sample_meta[sample_meta$Type=="Inj4",]

bac_snps <- list()
for(bac_sample in bac_meta[,1]){
    bac_snps[[bac_sample]] <- read.table(paste("scratch/", bac_sample, "/snippy/", bac_sample, "_ref.csv", sep=""), sep=",", quote="\"", header=T, fill=T, na.strings="")
    # Remove ones involving N
    bac_snps[[bac_sample]] <- bac_snps[[bac_sample]][!grepl("N", bac_snps[[bac_sample]]$REF) & !grepl("N", bac_snps[[bac_sample]]$ALT),]
}
bac_snps <- lapply(names(bac_snps), function(x) cbind(bac_meta[bac_meta$Sample==x,], bac_snps[[x]]))
names(bac_snps) <- bac_meta[,1]

# Compare two supposedly identical samples G4_01 and GG1
#bac_snps_gained <- anti_join(bac_snps[["G4_01"]][,-c(1:5,11)], bac_snps[["GG1"]][,-c(1:5,11)])
#bac_snps_lost <- anti_join(bac_snps[["GG1"]][,-c(1:5,11)], bac_snps[["G4_01"]][,-c(1:5,11)])
#bac_snps_G401_vs_GG1 <- list(Gain=bac_snps_gained, Loss=bac_snps_lost)
#bac_snps_G401_vs_GG1 <- cbind(Type=c(rep("Gain", nrow(bac_snps_gained)), rep("Loss", nrow(bac_snps_lost))), do.call(rbind, bac_snps_G401_vs_GG1))

# Compare samples to R0 sample GG1
bac_snps_sample_vs_r0 <- list()
bac_samples <- bac_meta[!is.na(bac_meta[,5]),1]
bac_samples_meta <- bac_meta[bac_meta$Sample%in%bac_samples,]
bac_samples <- bac_samples[order(bac_samples_meta$Round, bac_samples_meta$Line)]
for(bac_sample in bac_samples){
    bac_snps_gained <- bac_snps[[bac_sample]][!apply(bac_snps[[bac_sample]][,-c(1:5,11)], 1, function(x) any(apply(bac_snps[["GG1"]][,-c(1:5,11)], 1, function(y) all(x==y, na.rm=T)))),]
    bac_snps_lost <- bac_snps[["GG1"]][!apply(bac_snps[["GG1"]][,-c(1:5,11)], 1, function(x) any(apply(bac_snps[[bac_sample]][,-c(1:5,11)], 1, function(y) all(x==y, na.rm=T)))),]
    bac_snps_combined <- list(Gain=bac_snps_gained, Loss=bac_snps_lost)
    bac_snps_sample_vs_r0[[bac_sample]] <- cbind(SampleCompared=rep(paste(bac_sample, "_vs_GG1", sep=""), nrow(bac_snps_gained)+nrow(bac_snps_lost)), Type=c(rep("Gain", nrow(bac_snps_gained)), rep("Loss", nrow(bac_snps_lost))), do.call(rbind, bac_snps_combined))
}

# Compare samples with their appropriate previous sample
bac_snps_sample_vs_previous <- list()
for(bac_sample in bac_samples){
    previous = bac_meta[bac_meta[,1]==bac_sample,5]
#    bac_snps_gained <- anti_join(bac_snps[[bac_sample]][,-c(1:5,11)], bac_snps[[previous]][,-c(1:5,11)])
#    bac_snps_lost <- anti_join(bac_snps[[previous]][,-c(1:5,11)], bac_snps[[bac_sample]][,-c(1:5,11)])
    bac_snps_gained <- bac_snps[[bac_sample]][!apply(bac_snps[[bac_sample]][,-c(1:5,11)], 1, function(x) any(apply(bac_snps[[previous]][,-c(1:5,11)], 1, function(y) all(x==y, na.rm=T)))),]
    bac_snps_lost <- bac_snps[[previous]][!apply(bac_snps[[previous]][,-c(1:5,11)], 1, function(x) any(apply(bac_snps[[bac_sample]][,-c(1:5,11)], 1, function(y) all(x==y, na.rm=T)))),]
    bac_snps_combined <- list(Gain=bac_snps_gained, Loss=bac_snps_lost)
    bac_snps_sample_vs_previous[[bac_sample]] <- cbind(SampleCompared=rep(paste(bac_sample, "_vs_", previous, sep=""), nrow(bac_snps_gained)+nrow(bac_snps_lost)), Type=c(rep("Gain", nrow(bac_snps_gained)), rep("Loss", nrow(bac_snps_lost))), do.call(rbind, bac_snps_combined))
}
bac_snps_sample_vs_previous[["G4_01"]] <- data.frame()
bac_snps_sample_vs_previous[["G4_05"]] <- data.frame()
bac_snps_sample_vs_previous[["GG1"]] <- data.frame()

# Visual output to make sense of it
bac_sample_coords <- bac_meta[,c("Round", "Line")]
bac_sample_names <- bac_meta[,1]
bac_sample_coords[is.na(bac_sample_coords)] <- 0
bac_sample_coords[bac_sample_coords=="P"]
bac_sample_coords$Line[bac_sample_coords$Line==0] <- 0
bac_sample_coords <- apply(bac_sample_coords, 2, as.numeric)
bac_sample_coords[9,] <- c(7, 4.2) #nudge

pdf("scratch/bac_snp_analysis.pdf", width=24, height=12)
plot(bac_sample_coords, pch=20, xlim=c(0,16), ylim=c(0,12))#, main=titles[i])
seg_coords <- cbind(bac_sample_coords, bac_sample_coords[match(bac_meta[,5], bac_meta[,1]),])
#seg_coords <- seg_coords[apply(seg_coords, 1, function(x) all(!is.na(x))),]
arrows(seg_coords[,3], seg_coords[,4], seg_coords[,1], seg_coords[,2], col=2, length=0.1)
labels <- lapply(bac_snps_sample_vs_previous[bac_meta[,1]], function(x) {ifelse(nrow(x)==0, "", paste(apply(x, 1, paste, collapse=" "), collapse="\n"))})
text(seg_coords[,1]+0.01, seg_coords[,2]+0.01, labels, cex=0.5, adj=c(0,0), srt=20)
text(seg_coords[,1]-0.01, seg_coords[,2]-0.01, bac_sample_names, cex=0.6, adj=c(1,0), srt=20)
dev.off()

# Output tables
all_bac_snps <- do.call(rbind, bac_snps)
#all_bac_snps <- all_bac_snps[order(all_bac_snps$CHROM, all_bac_snps$POS),]
write.table(all_bac_snps, "scratch/all_bac_snps.txt", row.names=F)

all_bac_changes_vs_prev <- do.call(rbind, bac_snps_sample_vs_previous)
#all_bac_changes_vs_prev <- all_bac_changes[order(all_bac_changes_vs_prev$CHROM, all_bac_changes_vs_prev$POS),]
all_bac_changes_vs_prev <- rbind(cbind(SampleCompared="GG1_vs_Ref", Type="Initial", bac_snps[['GG1']]), all_bac_changes_vs_prev)
write.table(all_bac_changes_vs_prev, "scratch/all_bac_changes_vs_prev.txt", row.names=F)

all_bac_changes_vs_r0 <- do.call(rbind, bac_snps_sample_vs_r0)
all_bac_changes_vs_r0 <- rbind(cbind(SampleCompared="GG1_vs_Ref", Type="Initial", bac_snps[['GG1']]), all_bac_changes_vs_r0)
write.table(all_bac_changes_vs_r0, "scratch/all_bac_changes_vs_r0.txt", row.names=F)

# Detailed snp analysis
source("code/summarise_pile.r")

bac_freq_tables <- apply(all_bac_changes_vs_r0, 1, function(x) summarise_pile(x[8], as.numeric(x[9]), as.numeric(x[9]), bac_samples))
sink("scratch/bac_snps_details.txt")
print(bac_freq_tables)
sink()

########################

# Now the fungus samples
fun_meta <- sample_meta[sample_meta$Type=="Fungus",]

fun_snps <- list()
for(fun_sample in fun_meta[,1]){
    fun_snps[[fun_sample]] <- read.table(paste("scratch/", fun_sample, "/snippy/", fun_sample, "_ref.csv", sep=""), sep=",", quote="\"", header=T, fill=T, na.strings="")

    # Remove ones involving N
    fun_snps[[fun_sample]] <- fun_snps[[fun_sample]][!grepl("N", fun_snps[[fun_sample]]$REF) & !grepl("N", fun_snps[[fun_sample]]$ALT),]
}
fun_snps <- lapply(names(fun_snps), function(x) cbind(fun_meta[fun_meta$Sample==x,], fun_snps[[x]]))
names(fun_snps) <- fun_meta[,1]

# BLAST due to lack of annotation
#get_seq_region <- function(snp){
#    cmd <- paste("python ~/scripts/extract_sequence_region.py --contig", snp[1], "--start", as.numeric(snp[2])-500, "--end", as.numeric(snp[2])+500, "scratch/fungus_2/assembly/assembly.fasta")
#    seq <- system(cmd, intern=T)[2]
#    return(seq)
#}
#
#blast_seq <- function(seq){
#    sink("query.fa")
#    cat(">seq\n")
#    cat(seq)
#    cat("\n")
#    sink()
#    cmd <- paste('blastn -query query.fa -db lib/GCF_002708625.1_Rhimi1_1_cds_from_genomic.fna -outfmt "6 stitle" ')
#    cat(paste("Running", cmd, "\n"))
#    hits <- system(cmd, intern=T)
#    return(hits)
#}
#
#seqs <- lapply(fun_snps, function(x) apply(x, 1, function(snp) get_seq_region(snp)))
#blast_hits <- lapply(seqs, function(set) sapply(set, function(seq) blast_seq(seq)))
#reg_hits <- lapply(blast_hits, function(set) sapply(set, function(hit) unique(unlist(regmatches(hit, regexec("protein=(.*?)]", hit, perl=T))))))
#prot_names <- lapply(reg_hits, function(set) sapply(set, function(x) paste(x[1:length(x)%%2==0], collapse="/")))
#
#for(i in 1:length(fun_snps)){
#    fun_snps[[i]]$PRODUCT <- prot_names[[i]]
#}

# Compare samples to R1 sample R0002
fun_snps_sample_vs_r0 <- list()
fun_samples <- fun_meta[!is.na(fun_meta[,5]),1]
fun_samples_meta <- fun_meta[fun_meta$Sample%in%fun_samples,]
fun_samples <- fun_samples[order(fun_samples_meta$Round, fun_samples_meta$Line)]
for(fun_sample in fun_samples){
    fun_snps_gained <- fun_snps[[fun_sample]][!apply(fun_snps[[fun_sample]][,-c(1:5,11)], 1, function(x) any(apply(fun_snps[["R0002"]][,-c(1:5,11)], 1, function(y) all(x==y, na.rm=T)))),]
    fun_snps_lost <- fun_snps[["R0002"]][!apply(fun_snps[["R0002"]][,-c(1:5,11)], 1, function(x) any(apply(fun_snps[[fun_sample]][,-c(1:5,11)], 1, function(y) all(x==y, na.rm=T)))),]
    fun_snps_combined <- list(Gain=fun_snps_gained, Loss=fun_snps_lost)
    fun_snps_sample_vs_r0[[fun_sample]] <- cbind(SampleCompared=rep(paste(fun_sample, "_vs_R0002", sep=""), nrow(fun_snps_gained)+nrow(fun_snps_lost)), Type=c(rep("Gain", nrow(fun_snps_gained)), rep("Loss", nrow(fun_snps_lost))), do.call(rbind, fun_snps_combined))
}

# Compare samples with their appropriate previous sample
fun_snps_sample_vs_previous <- list()
fun_samples <- fun_meta[!is.na(fun_meta[,5]),1]
for(fun_sample in fun_samples){
    previous = fun_meta[fun_meta[,1]==fun_sample,5]
#    fun_snps_gained <- anti_join(fun_snps[[fun_sample]][,-c(1:5,11)], fun_snps[[previous]][,-c(1:5,11)])
#    fun_snps_lost <- anti_join(fun_snps[[previous]][,-c(1:5,11)], fun_snps[[fun_sample]][,-c(1:5,11)])
    fun_snps_gained <- fun_snps[[fun_sample]][!apply(fun_snps[[fun_sample]][,-c(1:5,11)], 1, function(x) any(apply(fun_snps[[previous]][,-c(1:5,11)], 1, function(y) all(x==y, na.rm=T)))),]
    fun_snps_lost <- fun_snps[[previous]][!apply(fun_snps[[previous]][,-c(1:5,11)], 1, function(x) any(apply(fun_snps[[fun_sample]][,-c(1:5,11)], 1, function(y) all(x==y, na.rm=T)))),]
    fun_snps_combined <- list(Gain=fun_snps_gained, Loss=fun_snps_lost)
    fun_snps_sample_vs_previous[[fun_sample]] <- cbind(SampleCompared=rep(paste(fun_sample, "_vs_", previous, sep=""), nrow(fun_snps_gained)+nrow(fun_snps_lost)), Type=c(rep("Gain", nrow(fun_snps_gained)), rep("Loss", nrow(fun_snps_lost))), do.call(rbind, fun_snps_combined))
}

# Visual output to make sense of it
fun_sample_coords <- fun_meta[,c("Round", "Line")]
fun_sample_names <- fun_meta[,1]
fun_sample_coords[is.na(fun_sample_coords)] <- 0
fun_sample_coords[fun_sample_coords=="P"]
fun_sample_coords$Line[fun_sample_coords$Line==0] <- 0
fun_sample_coords <- apply(fun_sample_coords, 2, as.numeric)
fun_sample_coords[9,] <- c(7, 4.2) #nudge

pdf("scratch/fun_snp_analysis.pdf", width=24, height=12)
plot(fun_sample_coords, pch=20, xlim=c(0,16), ylim=c(0,12))#, main=titles[i])
seg_coords <- cbind(fun_sample_coords, fun_sample_coords[match(fun_meta[,5], fun_meta[,1]),])
#seg_coords <- seg_coords[apply(seg_coords, 1, function(x) all(!is.na(x))),]
arrows(seg_coords[,3], seg_coords[,4], seg_coords[,1], seg_coords[,2], col=2, length=0.1)
labels <- lapply(fun_snps_sample_vs_previous[fun_meta[,1]], function(x) {ifelse(nrow(x)==0, "", paste(apply(x, 1, paste, collapse=" "), collapse="\n"))})
text(seg_coords[,1]+0.01, seg_coords[,2]+0.01, labels, cex=0.5, adj=c(0,0), srt=20)
text(seg_coords[,1]-0.01, seg_coords[,2]-0.01, fun_sample_names, cex=0.6, adj=c(1,0), srt=20)
dev.off()

# Output tables
all_fun_snps <- do.call(rbind, fun_snps)
#all_fun_snps <- all_fun_snps[order(all_fun_snps$CHROM, all_fun_snps$POS),]
write.table(all_fun_snps, "scratch/all_fun_snps.txt", row.names=F)

all_fun_changes_vs_prev <- do.call(rbind, fun_snps_sample_vs_previous)
#all_fun_changes <- all_fun_changes[order(all_fun_changes$CHROM, all_fun_changes$POS),]
all_fun_changes_vs_prev <- rbind(cbind(SampleCompared="R0002_vs_Asm", Type="Initial", fun_snps[['R0002']]), all_fun_changes_vs_prev)
write.table(all_fun_changes_vs_prev, "scratch/all_fun_changes_vs_prev.txt", row.names=F)

all_fun_changes_vs_r0 <- do.call(rbind, fun_snps_sample_vs_r0)
all_fun_changes_vs_r0 <- rbind(cbind(SampleCompared="R0002_vs_Asm", Type="Initial", fun_snps[['R0002']]), all_fun_changes_vs_r0)
write.table(all_fun_changes_vs_r0, "scratch/all_fun_changes_vs_r0.txt", row.names=F)

# Detailed snp analysis
source("code/summarise_pile.r")

fun_freq_tables <- apply(all_fun_changes_vs_r0, 1, function(x) summarise_pile(x[8], as.numeric(x[9]), as.numeric(x[9]), fun_samples))
sink("scratch/fun_snps_details.txt")
print(fun_freq_tables)
sink()
