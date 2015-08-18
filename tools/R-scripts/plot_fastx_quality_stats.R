library(ggplot2)
library(plyr)
library(reshape2)
library(stringr)

fastx.quality.stats.dir <- c("example-configurations/example-out/ChIPseq/fastx_quality_stats-202e/")
fastx.quality.stats.files <- list.files(fastx.quality.stats.dir)
fastx.quality.stats.path <- paste(fastx.quality.stats.dir, fastx.quality.stats.files, sep="/")
fastx.quality.stats <- llply(fastx.quality.stats.path, read.csv, sep="\t", fill=FALSE, header=TRUE)

sample.names <- strsplit(fastx.quality.stats.files, ".fastq.quality.tsv")
sample.names <- ldply(sample.names, function(x) x[1])
names(fastx.quality.stats) <- sample.names[[1]]

raw.read.counts <- ldply(fastx.quality.stats, function(list.name) list.name$ALL_count[1])
names(raw.read.counts) <- c("Sample", "Raw.reads")

i <- ggplot(raw.read.counts, aes_string(x="Sample", y="Raw.reads", fill="Sample"))
i <- i + geom_bar(stat="identity", position="dodge")
i <- i + theme_bw(base_size=10)
i <- i + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#               legend.position=c(0,1),
#               legend.justification=c(0,1))
i <- i + geom_hline(aes(yintercept=c(10000000, 20000000)))
#i <- i + geom_text(aes_string(label="NGS_QC_SCORE"), size=3)
i <- i + scale_y_log10("Reads")
i
ggsave(plot=i, file="raw_reads_per_sample.pdf", width=50, height=16, units="cm")


sample.info.file <- c("../statistics/sample_info.txt")
sample.info.df <- read.csv(sample.info.file, sep=",", fill=FALSE, header=TRUE)
raw.reads.stat <- c("../statistics/count_fastq.raw_read_statistics.txt")
raw.reads.df <- read.csv(raw.reads.stat, sep=",", fill=FALSE, header=TRUE)
mapped.reads.stat <- c("../statistics/count_tophat2.mapped_reads_statistics.txt")
mapped.reads.df <- read.csv(mapped.reads.stat, sep=",", fill=FALSE, header=TRUE)
duplicates.stat <- c("../statistics/MarkDuplicates_statistics.txt")
duplicates.df <- read.csv(duplicates.stat, sep=",", fill=FALSE, header=TRUE)
ngsqc.rmdup.stat <- c("../statistics/ngs-qc_scores-rm-dup.txt")
ngsqc.rmdup.df <- read.csv(ngsqc.rmdup.stat, sep=",", fill=FALSE, header=TRUE)
experiment.stat <- c("../statistics/experiment_statistics.txt")
experiment.df <- read.csv(experiment.stat, sep=",", fill=FALSE, header=TRUE)
ratios.stat <- c("../statistics/fragment_ratio_100-500_vs_50-9000.csv")
ratios.df <- read.csv(ratios.stat, sep=",", fill=FALSE, header=TRUE)

dfs2merge <- list(sample.info.df, raw.reads.df, mapped.reads.df, duplicates.df, ngsqc.rmdup.df, experiment.df)

merge.dfs <- function(x,y){
    return(merge(x, y, by="SAMPLE_NAME"))
}

# prepare data.frame "all" by adding data
all <- Reduce(merge, dfs2merge)
all$FRAGMENTS_EXAMINED <- all$UNPAIRED_READS_EXAMINED + all$READ_PAIRS_EXAMINED
all$DUPLICATE_FRAGMENTS <- all$UNPAIRED_READ_DUPLICATES + all$READ_PAIR_DUPLICATES
all$UNIQ_FRAGMENTS <- all$FRAGMENTS_EXAMINED - all$DUPLICATE_FRAGMENTS
q <- as.character(c("A", "B", "C", "D"))
p <- expand.grid(q, q, q)
o <- paste(p[[3]], p[[2]], p[[1]], sep="")
all$NGS_QC_SCORE <- factor(all$NGS_QC_SCORE, levels=rev(o), ordered=TRUE)

source("~/develop/bioinf-scripts/R/ggplot2/ggplot2_convenience_functions.R")

xfp <- all[, c("SAMPLE_NAME", "BC", "TIMEPOINT", "ANTIBODY", "UNIQ_FRAGMENTS", "DUPLICATE_FRAGMENTS", "NGS_QC_SCORE")]
xfp <- melt(xfp, id.vars=c("SAMPLE_NAME", "BC", "TIMEPOINT", "ANTIBODY", "NGS_QC_SCORE"))

i <- ggplot(xfp, aes_string(x="SAMPLE_NAME", y="value", fill="variable"))
i <- i + geom_bar(stat="identity", position="dodge")
i <- i + theme_bw(base_size=10)
i <- i + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#               legend.position=c(0,1),
#               legend.justification=c(0,1))
i <- i + geom_hline(aes(yintercept=c(10000000, 20000000)))
#i <- i + geom_text(aes_string(label="NGS_QC_SCORE"), size=3)
i <- i + scale_y_log10("Reads")
i
ggsave(plot=i, file="unique_vs_duplicate_reads_per_sample.pdf", width=25, height=16, units="cm")

yfp <- all[, c("SAMPLE_NAME", "BC", "TIMEPOINT", "ANTIBODY", "RAW_READS", "MAPPED_READS", "NGS_QC_SCORE")]
yfp <- melt(yfp, id.vars=c("SAMPLE_NAME", "BC", "TIMEPOINT", "ANTIBODY", "NGS_QC_SCORE"))
y <- ggplot(yfp, aes_string(x="SAMPLE_NAME", y="value", fill="variable"))
y <- y + geom_bar(stat="identity", position="dodge")
y <- y + theme_bw(base_size=10)
y <- y + theme(axis.text.x = element_text(angle = 90, hjust = 1),
               legend.position=c(0,1),
               legend.justification=c(0,1))
y <- y + scale_y_continuous("Reads")
y <- y + scale_x_discrete("Sample Name")
y
ggsave(plot=y, file="raw_vs_mapped_reads_per_sample.pdf", width=25, height=16, units="cm")


jfp <- all[, c("SAMPLE_NAME", "BC", "TIMEPOINT", "ANTIBODY", "FRAGMENTS_EXAMINED", "MAPPED_READS")]
jfp <- melt(jfp, id.vars=c("SAMPLE_NAME", "BC", "TIMEPOINT", "ANTIBODY"))

j <- ggplot(jfp, aes_string(x="SAMPLE_NAME", y="value", fill="variable"))
j <- j + geom_bar(stat="identity", position="dodge")
j <- j + theme_bw(base_size=10)
j <- j + theme(axis.text.x = element_text(angle = 90, hjust = 1))
j
ggsave(plot=j, file="mapped_reads_vs_examined_fragments_per_sample.pdf", width=25, height=10, units="cm")



#################################################################################
# Percent Duplication plots
#################################################################################

m <- ggplot(all, aes_string(x="PERCENT_DUPLICATION", y="PRECIPITATED_DNA_ng.µl", color="TIMEPOINT", group=1))
m <- m + geom_point(stat="identity")
m <- m + theme_bw(base_size=10)
m <- m + scale_y_continuous( "Immunoprecipitated DNA [ng/micro l]" )
m <- m + xlim(0, 1)
m <- m + scale_x_continuous("Fraction Duplicated Fragments")
#m <- m + geom_smooth(method="lm")
m
ggsave(plot=m, file="percent_duplicates_vs_IP_DNA.pdf", width=25, height=10, units="cm")

summary( lm(all$PERCENT_DUPLICATION ~ all$PRECIPITATED_DNA_ng.µl) )

a <- ggplot(all, aes_string(x="PERCENT_DUPLICATION", y="RAW_READS", color="TIMEPOINT"))
a <- a + geom_point()
a <- a + scale_x_log10()
a <- a + scale_y_continuous()
a <- a + theme_bw(base_size=10)
a <- a + ylab("log10(Nr. of Raw Reads)")
a <- a + xlim(0, 1)
a <- a + xlab("Fraction Duplicated Fragments")
a
ggsave(plot=a, file=c("percent_duplicates_vs_raw_reads.pdf"), width=25, height=10, units="cm")

duplicates.info <- merge(sample.info.df, duplicates.df, by="SAMPLE_NAME")
duplicates.ratios <- merge(duplicates.info, ratios.df, by=c("BC", "TIMEPOINT"))
b <- ggplot(duplicates.ratios, aes_string(x="PERCENT_DUPLICATION", y="pg_µl_ratio..100.500nt..50.9000nt..", color="TIMEPOINT", group=1))
b <- b + geom_point()
b <- b + scale_x_continuous()
b <- b + scale_y_continuous()
b <- b + theme_bw(base_size=10)
b <- b + ylim(0, 1)
b <- b + ylab("Ratio Fragment Conc. [c(100nt to 500nt)/c(50nt to 9000nt)]")
b <- b + xlim(0, 1)
b <- b + xlab("Fraction Duplicated Fragments")
# b <- b + geom_smooth(method="lm")
b
ggsave(plot=b, file=c("percent_duplicates_vs_ratio_fragment_conc.pdf"), width=25, height=10, units="cm")

c <- ggplot(duplicates.ratios, aes_string(x="PERCENT_DUPLICATION", y="pmol_l.ratio..100.500nt..50.9000nt..", color="TIMEPOINT", group=1))
c <- c + geom_point()
c <- c + scale_x_continuous()
c <- c + scale_y_continuous()
c <- c + theme_bw(base_size=10)
c <- c + ylim(0, 1)
c <- c + ylab("Ratio Fragment Mol. [M(100nt to 500nt)/M(50nt to 9000nt)]")
c <- c + xlim(0, 1)
c <- c + xlab("Fraction Duplicated Fragments")
c
ggsave(plot=c, file=c("percent_duplicates_vs_ratio_fragment_mol.pdf"), width=25, height=10, units="cm")



# summary( lm(duplicates.ratios$PERCENT_DUPLICATION ~ duplicates.ratios$pg_µl_ratio..100.500nt..50.9000nt..) )

#################################################################################
# Antibody Quality
#################################################################################

g <- ggplot(all, aes_string(x="ANTIBODY", y="PERCENT_DUPLICATION", fill="TIMEPOINT"))
g <- g + geom_boxplot()
g <- g + geom_jitter(aes_string(color="TIMEPOINT"))
g <- g + scale_x_discrete()
g <- g + scale_y_continuous()
g <- g + theme_bw(base_size=10)
g <- g + xlab("Antibody")
g <- g + ylim(0, 1)
g <- g + ylab("Fraction Duplicated Fragments")
g
ggsave(plot=g, file="percent_duplicates_per_antibody_and_timepoint.pdf", width=32, height=16, units="cm")

k <- ggplot(all, aes_string(x="ANTIBODY", y="PRECIPITATED_DNA_ng.µl", fill="TIMEPOINT"))
k <- k + geom_boxplot()
k <- k + geom_jitter(aes_string(color="TIMEPOINT"))
k <- k + theme_bw(base_size=10)
k <- k + scale_y_continuous("Immunoprecipitated DNA [ng/µl]")
k <- k + xlim(0, 1)
k <- k + scale_x_discrete("Antibody")
k
ggsave(plot=k, file="antibody_vs_IP_DNA.pdf", width=25, height=10, units="cm")

l <- ggplot(all, aes_string(x="ANTIBODY", y="NGS_QC_SCORE", color="TIMEPOINT"))
l <- l + geom_jitter(position=position_jitter(w = 0.25, h = 0.25))
l <- l + theme_bw(base_size=10)
l <- l + scale_y_discrete("NGS-QC Scores", drop=FALSE)
l <- l + scale_x_discrete("Antibody")
l
ggsave(plot=l, file="antibody_vs_ngs_qc_score.pdf", width=25, height=16, units="cm")

h <- ggplot(all, aes_string(x="ANTIBODY", y="RAW_READS", color="TIMEPOINT"))
h <- h + geom_boxplot()
h <- h + geom_jitter()
h <- h + scale_x_discrete("Antibody")
h <- h + scale_y_log10("log10(Nr. of Raw Reads)")
h <- h + theme_bw(base_size=10)
h
ggsave(plot=h, file="antibody_vs_raw_reads.pdf", width=25, height=10, units="cm")
