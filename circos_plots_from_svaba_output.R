library("VariantAnnotation")
library("StructuralVariantAnnotation")
library("circlize")
library('tools')
library('dplyr')


# R script to take the vcf files from the output of svaba and creates 
# ciros plots for the structural variants
# genome mapped against is hg38 (this can be changed to hg19)
# Works in the current working directory and puts plots in dir called plots
# plots dir will have to be created 
# mkdir plots


circos_from_svaba_output <- function(x) {
  file_name <- strsplit(x, "\\.")[[1]]
  new_name <- paste(file_name[[1]], file_name[[4]], sep="_")
  out_name <- paste(new_name, '.jpeg', sep='')
  out_path <- paste('./plots/', out_name, sep='')
  vcf <- VariantAnnotation::readVcf(x, "hg38")
  g <- breakpointRanges(vcf)
  
  colo829_bpgr_with_chr_prefix <- g
  seqlevelsStyle(colo829_bpgr_with_chr_prefix) <- "UCSC"
  pairs <- breakpointgr2pairs(colo829_bpgr_with_chr_prefix)

  seq1 <- (as.data.frame(S4Vectors::first(pairs)))
  seq2 <- (as.data.frame(S4Vectors::second(pairs)))

  seq1 <- cbind(index = rownames(seq1), seq1)
  seq2 <- cbind(index = rownames(seq2), seq2)

  seq1 <- seq1 %>% 
    select('index','seqnames','start', 'end')  %>% 
    mutate(id = row_number())

  seq2 <- seq2 %>% 
    select('index', 'seqnames','start', 'end') %>% 
    mutate(id = row_number())

  df <- merge(x = seq1, y = seq2, by = "id", all = TRUE)
  df <- df[!grepl("chrUn", df$seqnames.x),]
  df <- df[!grepl("random", df$seqnames.x),]
  df <- df[!grepl("alt", df$seqnames.x),]
  df <- df[!grepl("chrUn", df$seqnames.y),]
  df <- df[!grepl("random", df$seqnames.y),]
  df <- df[!grepl("alt", df$seqnames.y),]

  seq1 <- df[c('index.x', 'seqnames.x', 'start.x', 'end.x')]

  seq1 <- seq1 %>% 
    rename(
      index = index.x,
      seqnames = seqnames.x,
      start = start.x,
      end = end.x
    )

  seq2 <- df[c('index.y', 'seqnames.y', 'start.y', 'end.y')]
  seq2 <- seq2 %>% 
    rename(
      index = index.y,
      seqnames = seqnames.y,
      start = start.y,
      end = end.y
    )

  rownames(seq1) <- seq1$index
  seq1$index <- NULL
  rownames(seq2) <- seq2$index
  seq2$index <- NULL  


  jpeg(out_path, width = 6, height = 6, units = 'in', res = 200)
  plot <-circos.initializeWithIdeogram(species = "hg38")
  plot <- circos.genomicLink(seq1, seq2)
  dev.off()
}

file.list <- dir(pattern = "\\.vcf$")

sapply(file.list, circos_from_svaba_output)