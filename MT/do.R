# This script loads pre-generated QDNAseq rda objects (with hg38 bin annotations for various binsizes)
# and then generates and adds the appropriate bins for MT/chrM.
#
# The pre-generated QDNAseq objects are located in MT/rda_orig
# The updated versions with MT/chrM data added are saved to data/
#
# Be sure to download binary version of bigWigAverageOverBed at https://hgdownload.soe.ucsc.edu/admin/exe/ and add it to your $PATH

library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(QDNAseq)
library(future)
library(data.table)

# set working directory
setwd('<PATH TO REPO CLONE>')

## load chrM sequence
mt <- read.table('MT/MT.fa',sep='\n',header=T)[[1]]
mt <- paste(mt, collapse='')
mt_chr <- strsplit(mt,'')[[1]]
mt_len <- length(mt_chr)
mt_dat <- data.table(nt=mt_chr)
mt_dat$pos <- 1:nrow(mt_dat)
mt_dat[,pos2:=pos+1]
mt_dat$chr <- 'MT'
setkey(mt_dat,'chr','pos','pos2')

for (binsize in c(1000, 500, 100, 50, 30, 15, 10, 5, 1)) {
    message(binsize)
    starts <- seq(1, mt_len, by=binsize*1000)
    ends <- c(tail(starts, -1) - 1, mt_len)
    Mbins <- data.table(chr='MT', start=starts, end=ends)
    Mbins[,region:=paste0(chr,':',start,'-',end)]
    setkey(Mbins,'chr','start','end')
    mt_dat_binned <- foverlaps(mt_dat, Mbins, type='any')
    mt_dat_binned <- mt_dat_binned[!duplicated(pos),]

    summarize_bin <- function(mt_dat_binned) { 
        len <- nrow(mt_dat_binned)
        nonN_bases <- sum(toupper(mt_dat_binned$nt) %in% c('A','C','G','T'))
        bases <- round((nonN_bases / len)*100, 4)
        n_gc <- sum(toupper(mt_dat_binned$nt) %in% c('C','G'))
        gc <- round((n_gc / nonN_bases)*100,5)
        data.table(bases=bases, gc=gc, blacklist=0, residual=NaN, use=T, len=len)
    }
    res <- mt_dat_binned[,summarize_bin(.SD), by=c('chr','start','end','region')]
    setnames(res,'chr','chromosome')
    res[,region:=paste0(chromosome,':',start,'-',end)] 
    res <- as.data.frame(res)
    rownames(res) <- res$region
    res$region <- NULL

    ## calculate average mappability for MT
    mappability <- calculateMappability(res,
                                        bigWigFile="MT/MT_mappability.genmap.50mer.bigwig",
                                        bigWigAverageOverBed="bigWigAverageOverBed")
    res$mappability <- mappability

    ## load the rdata object for the given binsize
    rda_file <- paste0('MT/rda_orig/hg38.',binsize,'kbp.SR50.rda')
    load(rda_file)
    obj <- eval(parse(text=paste0('hg38.',binsize,'kbp.SR50')))
    bins <- obj@data
    res <- res[,names(bins)]
    bins <- rbind(bins, res)
    obj@data <- bins  
    save(obj, file=paste0("data/hg38.", binsize, "kbp.SR50.rda"), compress='xz')
}





