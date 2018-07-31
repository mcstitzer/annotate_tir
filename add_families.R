library(data.table)
library(rtracklayer)
library(stringr)
library(dplyr)
library(stringr)

#tir=fread('~/Downloads/tir_B73_2018-07-26_extra.txt')

#tir=tir[,c('mtec', 'chr', 'start.adj', 'end.adj', 'tirseqSingle', 'tsdadjacentup', 'fam', 'sup')]
#tir=tir[tir$start.adj-tir$end.adj<0,]

#names(tir)=c('mtec', 'chr', 'start', 'end', 'tirscore', 'V6', 'V7', 'TIR1', 'TIR2')
tir.gr=GRanges(seqnames=tir$chr, ranges=IRanges(start=tir$start.adj, end=tir$end.adj))[tir$tirsmatch & tir$tsdadjacentequal]
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


tir.gr=tir.gr[-rmRows,]
tir.gr$mtec=tir$mtec[tir$tirsmatch & tir$tsdadjacentequal][-rmRows]
tir.gr$mtecfamnum=substr(tir.gr$mtec, 7,11)
tir.gr$sup=substr(tir.gr$mtec, 1, 3)
tir.gr$famname=paste(tir.gr$sup, tir.gr$mtecfamnum, sep='')
tir.gr$TSD=tir$tsdadjacentup[tir$tirsmatch & tir$tsdadjacentequal][-rmRows]
tir.gr$TIR=tir$tirseqSingle[tir$tirsmatch & tir$tsdadjacentequal][-rmRows]

strand(tir.gr)=tir$strand[tir$tirsmatch & tir$tsdadjacentequal][-rmRows]


## okay, now get the ones that overlap partially
selfOver=findOverlaps(tir.gr, drop.self=T, ignore.strand=T, type='start')
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
tir.gr=tir.gr[-rmRows,]

selfOver=findOverlaps(tir.gr, drop.self=T, ignore.strand=T, type='end')
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
tir.gr=tir.gr[-rmRows,]

tirNO=anti_join(rbind(data.frame(findOverlaps(tir.gr, ignore.strand=T, type='start')), 
                data.frame(findOverlaps(tir.gr, ignore.strand=T, type='end')), 
                data.frame(findOverlaps(tir.gr, ignore.strand=T, type='any'))), 
              rbind(data.frame(findOverlaps(tir.gr, ignore.strand=T, type='within')),
                setNames(data.frame(findOverlaps(tir.gr, ignore.strand=T, type='within')), c('subjectHits', 'queryHits'))))$subjectHits

tir.gr=tir.gr[-unique(tirNO),]

SHORTNAME='Zm00001d' ## for b73
#SHORTNAME='Zm00004b'  ## for w22

## assign 8digit family name (RL.XXXXX) and 10 digit copy name (B73v4XXXXX)
mcols(tir.gr)$Name=NULL
for (x in names(table(mcols(tir.gr)$famname))){
  mcols(tir.gr)$Name[mcols(tir.gr)$famname==x & !is.na(mcols(tir.gr)$famname)]=paste(SHORTNAME, str_pad(1:sum(mcols(tir.gr)$famname[!is.na(mcols(tir.gr)$famname)]==x), 5, pad='0'),sep='')
}

tir.gr$sup=substr(tir.gr$mtec, 1,3)
tir.gr$ID=paste0(tir.gr$famname, tir.gr$Name)


GENOMENAME='B73'
#GENOMENAME='W22'
### end -1 for gff3 format!
d=data.frame(chr=seqnames(tir.gr), 'TARGeT', 'terminal_inverted_repeat_element', start(tir.gr), end(tir.gr), '.', strand(tir.gr), '.', Name=paste0('ID=', tir.gr$ID, ';Name=', tir.gr$ID, '_', tir.gr$TSD, '_', tir.gr$TIR))
#d=d[tir$tsdadjacentequal & tir$tirsmatch,]
#write.table(d[!is.na(tir$whichrule) & d[,4]<d[,5],], file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
#write.table(d[d[,4]<d[,5],], file=paste0(GENOMENAME, '_unfiltered_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
write.table(d, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.', SHORTNAME, '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)

## for maizegdb with Chr
dd=d						
levels(dd$chr)[1:10]=paste0('Chr', levels(dd$chr)[1:10])	
write.table(dd, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.', SHORTNAME, '.Chr.gff3'), col.names=F, row.names=F, sep='\t', quote=F)


## need to assign family based on REDUCED MTEC FAMS!
#### there are way too many similar seqs, so same region called to different families.


