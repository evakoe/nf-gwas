#################################
### annotate top loci with GWAS catalog hits
################################

args <- commandArgs(trailingOnly=T)

# param 1: tsv file of clumped loci, required columns:
# CHROM, POS, 
loci.file <- args[1]

# param 2: bg file of the GWAS catalog
#TODO automatically download latest version from
#GWAS catalog
cat.file <- args[2]

#param 3: reference genome. one of hg19 or hg38
#TODO implement refrence genome matching, lifting
ref <- args[3]

# param 4: search radius for GWAS cat lookup in bp
dist <- args[4]

# param 5: output file
out <- args[5]

# read input files
loci <- read.table(loci.file, header=T, sep="\t")
cat <- read.table(cat.file, header=T, sep="\t", fill=T)


#################################
### functions
################################

parse.pos <- function(string){
    return(unlist(strsplit(string, " x ")))
}

# clean catalog for fast lookup
parse.cat <- function(cat){
    # there can be multiple chr and pos in a field due to SNP interactions
    # read into list
    chr.list <- sapply(cat$CHR_ID, parse.pos)
    pos.list <- sapply(cat$CHR_POS, parse.pos)
    # make an index for the duplicate entries to match
    len <- unlist(lapply(chr.list, length), use.names=F)
    idx <- c()
    for (i in 1:length(len)){
        idx <- c(idx, rep(i, len[i]))
    }
    # save the information if this was a multi SNP entry
    cat.new <- cat[idx,]
    cat.new$CHR_ID <- unlist(chr.list,use.names=F)
    cat.new$CHR_POS <- unlist(lapply(pos.list, as.numeric),use.names=F)
    return(cat.new)
}

quick.clean <- function(cat){
    cat$CHR_POS <- as.numeric(cat$CHR_POS)
    cat$CHR_ID <- as.numeric(cat$CHR_ID)
    return(cat)
}

#TODO check if loci start with "chr" and if yes remove

query.cat <- function(locus, cat){
    chrpos <- unlist(strsplit(locus, "\\:"))
    chr <- chrpos[1]
    pos <- as.numeric(chrpos[2])
    #identify rows from the cat that are +- dist of locus
    lower.bp <- if(pos - dist > 0) pos - dist else 1
    upper.bp <- pos + dist
    idx <- which(cat$CHR_ID == chr & cat$CHR_POS >= lower.bp & cat$CHR_POS <= upper.bp)
    return(idx)
}


#################################
### main
################################

# clean GWAS catalog
cat.clean <- quick.clean(cat)
# get indices to append
idx.list <- sapply(loci$CHRPOS, query.cat, cat.clean)
# add catalog to top loci
df <- c()
for (i in 1:length(idx.list)){
    this.locus <- loci[rep(i, length(idx.list[[i]])),]
    df <- rbind(df, cbind(this.locus, cat.clean[idx.list[[i]], ]))
}
write.table(df, file=out, row.names=F, col.names=T, sep="\t", quote=F)

