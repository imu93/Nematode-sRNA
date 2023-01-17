pacman::p_load(rtracklayer, Rsamtools)

system('samtools split -f %!.%. -@ 20 -v merged_alignments.bam')
system('for file in *trim.bam ; do outFile=$(echo $file| sed "s|.trim.bam|.trim.mapped.bam|"); samtools view -F 4  -b -@ 20 $file > $outFile; done')
system('rm -rf *.trim.bam')

infiles = list.files(pattern = ".*trim.bam$")
outfiles = sub("WIP", "IP", infiles)
outfiles = sub("unbound", "Ub", outfiles)
outfiles = sub("_mono|_mp", "_monoP", outfiles)
outfiles = sub("_poly|_pp", "_polyP", outfiles)
outfiles = sub("Hb_adult", "Adult_total", outfiles)
outfiles = sub("Hb_vesicle", "EV_total", outfiles)
file.rename(from_files,outfiles)

what_to_extract = c("rname", "pos", "qwidth", "strand", "flag", "seq", "qname")
bams = list.files(pattern = ".*trim.bam$")
for (file in bams) {
  bamf = file
  print(bamf)
  example = as.data.frame(scanBam(bamf, param = ScanBamParam(what = what_to_extract, reverseComplement = TRUE, flag = scanBamFlag(isUnmappedQuery = FALSE))))
  toGR = GRanges(seqnames = example$rname,IRanges(start = example$pos, width = example$qwidth), strand = example$strand)
  toGR$len = width(toGR)
  toGR$Id = as.character(example$qname)
  toGR$flag = as.numeric(example$flag)
  toGR$seq = as.character(example$seq)
  #toGR <- toGR[toGR$len >= 20 & toGR$len <= 25,]
  saveRDS(toGR, file = gsub(".bam", ".Rds", bamf))
}
