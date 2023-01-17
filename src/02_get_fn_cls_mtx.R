pacman::p_load(rtracklayer,ggplot2, reshape2, Rsamtools, ShortRead)
gff3File = "hBak_seg_genome_14122022.gff3"
ann = import(gff3File)
files = list.files(pattern = ".*.mapped.Rds$")
libs = lapply(files, readRDS)
names(libs) = files

get_fn_mtx = function(x){
  lib  = subset(x, seqnames(x) != "MtDNA")
  fnuc = substr(lib$seq,1,1)
  lens = split(fnuc, width(lib))
  lis = lapply(lens,table)
  lis = lapply(lis, function(x){x[c("C","T","A","G")]})
  counts = as.data.frame(do.call(rbind, lis))
  return(counts)
}


get_cls_mtx = function(lib) {
  rs_lib = resize(lib, 1, "center")
  cls_lib =  mergeByOverlaps(rs_lib, ann)
  cls_lib$n_type = "Other"
  cls_lib$n_type = ifelse(grepl("exons_As|introns_As", cls_lib$ID), "Gene_As", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("DNA.*_AS|RC.*_As", cls_lib$ID), "Transposon_As", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("LINE.*_As|SINE.*_As|LTR.*_As", cls_lib$ID), "Retrotransposon_As", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("Unk", cls_lib$ID), "Unknown", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("Sat|Low|Simple", cls_lib$ID), "Other_repeats", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("miRNA_S", cls_lib$ID), "miRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("rRNA_S", cls_lib$ID), "rRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("tRNA_S", cls_lib$ID), "tRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("yRNA_S", cls_lib$ID), "yRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("piRNA_S|snoRNA_S|snRNA_S|other_ncRNA_S", cls_lib$ID), "Other_ncRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("inter", cls_lib$ID), "Intergenic", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("exon_S|intron_S|DNA.*_S|RC.*_S|LINE.*_S|SINE.*_S|LTR.*_S|miRNA_As|rRNA_As|tRNA_As|yRNA_As|piRNA_As|snoRNA_As|snRNA_As|other_ncRNA_As", cls_lib$ID), "Other", cls_lib$n_type)
  cls_lib$n_type = factor(cls_lib$n_type)
  len_lib = split(cls_lib, lib$len)
  cls_len = lapply(len_lib, function(x){table(x$n_type)})
  mtx = do.call(rbind, cls_len)
  mtx_ed = apply(mtx, 2,as.numeric)
  rownames(mtx_ed) = rownames(mtx)
  return(mtx_ed)
}


fil_cls_mtx = function(x){
  cats = c("Gene_As", "Transposon_As", "Retrotransposon_As", "Unknown", "Other_repeats",
           "miRNA", "rRNA", "tRNA", "yRNA", "Other_ncRNA", "Intergenic", "Other")
  l = x
  ns = setdiff(cats, colnames(l))
  if (identical(ns, character(0))) {
      return(l)
  }else{
      vs =  t(matrix(0, length(ns), nrow(l)))
      colnames(vs) =  ns
      rownames(vs) = rownames(l)
      l = cbind(l, vs)
      l =  l[,cats]
      return(l)  
  }
  }

cls_mtx =  lapply(libs, get_cls_mtx)
cls_mtx = lapply(cls_mtx, fil_cls_mtx)
fn_mtx = lapply(libs, get_fn_mtx)
lst = list(fn_mtx,cls_mtx)
saveRDS(lst, "Hb_libs_fn_cls_mtx.R")


