library(data.table)
library(rtracklayer)

tir=fread('~/Downloads/tir_B73_2018-07-26_extra.txt')

tir=tir[,c('mtec', 'chr', 'start.adj', 'end.adj', 'tirseqSingle', 'tsdadjacentup', 'fam', 'sup')]
tir=tir[tir$start.adj-tir$end.adj<0,]

#names(tir)=c('mtec', 'chr', 'start', 'end', 'tirscore', 'V6', 'V7', 'TIR1', 'TIR2')
tir.gr=GRanges(seqnames=tir$chr, ranges=IRanges(start=tir$start.adj, end=tir$end.adj))
#mcols(tir.gr)$score=tir$tirscore

## okay to overlap if entirely within others - but not adjacent overlaps
selfOver=findOverlaps(tir.gr, drop.self=T, ignore.strand=T, type='equal')
rmRows=sapply(1:length(selfOver), function(x){ 
#  scoreA=mcols(tir.gr)$score[queryHits(selfOver)[x]]
#  scoreB=mcols(tir.gr)$score[subjectHits(selfOver)[x]]
#  if(scoreA>=scoreB){   ## arbitrarily prioritize the first entry if the scores are equal.
  return(max(subjectHits(selfOver)[x],queryHits(selfOver)[x])) ## return the larger number to be removed.
#  }else if (scoreB>scoreA){
#    return(subjectHits(selfOver)[x])
#  }
}
)



## need to assign family based on REDUCED MTEC FAMS!
#### there are way too many similar seqs, so same region called to different families.


