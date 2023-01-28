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
  cls_lib$n_type = ifelse(grepl("exons|introns", cls_lib$ID), "Gene", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("DNA|RC", cls_lib$ID), "Transposon", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("LINE|SINE.*_As|LTR", cls_lib$ID), "Retrotransposon", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("Unk", cls_lib$ID), "Unknown", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("Sat|Low|Simple", cls_lib$ID), "Other_repeats", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("miRNA", cls_lib$ID), "miRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("rRNA", cls_lib$ID), "rRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("tRNA", cls_lib$ID), "tRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("yRNA", cls_lib$ID), "yRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("piRNA|snoRNA|snRNA|other_ncRNA", cls_lib$ID), "Other_ncRNA", cls_lib$n_type)
  cls_lib$n_type = ifelse(grepl("inter", cls_lib$ID), "Intergenic", cls_lib$n_type)
  cls_lib$n_type = factor(cls_lib$n_type)
  len_lib = split(cls_lib, lib$len)
  cls_len = lapply(len_lib, function(x){table(x$n_type)})
  mtx = do.call(rbind, cls_len)
  mtx_ed = apply(mtx, 2,as.numeric)
  rownames(mtx_ed) = rownames(mtx)
  return(mtx_ed)
}


fil_cls_mtx = function(x){
  cats = c("Gene", "Transposon", "Retrotransposon", "Unknown", "Other_repeats",
           "miRNA", "rRNA", "tRNA", "yRNA", "Other_ncRNA", "Intergenic")
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
saveRDS(lst, "Hb_libs_fn_cls_mtx.Rds")


