library(data.table)
library(rtracklayer)
library(stringr)
library(dplyr)
library(stringr)


ti=import.gff3('B73.2018-08-02.allTE.gff3')


#tir=fread('~/Downloads/tir_B73_2018-07-26_extra.txt')

#tir=tir[,c('mtec', 'chr', 'start.adj', 'end.adj', 'tirseqSingle', 'tsdadjacentup', 'fam', 'sup')]
#tir=tir[tir$start.adj-tir$end.adj<0,]

## filters
filtersm=!is.na(tirm$closestTSDoffset) & tirm$seqdist<(nchar(tirm$adjustedTIRup)*0.2) & nchar(tirm$tirseqSingle)>=5 & tirm$end.adj-tirm$start.adj>=40

tirmf=tirm[filtersm,]

#names(tir)=c('mtec', 'chr', 'start', 'end', 'tirscore', 'V6', 'V7', 'TIR1', 'TIR2')
#tir.gr=GRanges(seqnames=tir$chr, ranges=IRanges(start=tir$start.adj, end=tir$end.adj))[filters]## add 5 bp minimum for TIR length, that is Bergamo's
tir.gr=GRanges(seqnames=tirmf$chr, ranges=IRanges(start=tirmf$start.adj-1, end=tirmf$end.adj-1))
#mcols(tir.gr)$score=tir$tirscore

### length filters - incorporated into filters above!
#tir.gr=tir.gr[width(tir.gr)>=40,]## stowaway mites can be small, using 40 bp as size cutoff. This removes 154 copies from B73.


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
tir.gr$mtec=tirmf$mtec[-rmRows]
tir.gr$mtecfamnum=substr(tir.gr$mtec, 7,11)
tir.gr$sup=substr(tir.gr$mtec, 1, 3)
tir.gr$famname=paste(tir.gr$sup, tir.gr$mtecfamnum, sep='')
tir.gr$TSD=tirmf$tsdadjacentup[-rmRows]
tir.gr$TIR=tirmf$tirseqSingle[-rmRows]

strand(tir.gr)=tirmf$strand[-rmRows]


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
#mcols(tir.gr)$Name=NULL
#for (x in names(table(mcols(tir.gr)$famname))){
#  mcols(tir.gr)$Name[mcols(tir.gr)$famname==x & !is.na(mcols(tir.gr)$famname)]=paste(SHORTNAME, str_pad(1:sum(mcols(tir.gr)$famname[!is.na(mcols(tir.gr)$famname)]==x), 5, pad='0'),sep='')
#}
#
#tir.gr$sup=substr(tir.gr$mtec, 1,3)
#tir.gr$ID=paste0(tir.gr$famname, tir.gr$Name)
#


## need to assign family based on REDUCED MTEC FAMS!
#### there are way too many similar seqs, so same region called to different families.

mtec=read.table('~/projects/tir_modified_mcs/mtec/MTEC_TIR.8080.fnodes')
mtec$famname=paste0(substr(mtec$V2,1,3), substr(mtec$V2,7,11))
mtec$collapsedfamname=sapply(1:nrow(mtec), function(x) sort(mtec$famname[mtec$V1==mtec$V1[x]])[1])

                             
## this is concerning - some families collapse across superfamilies.
mtec[substr(mtec$famname,1,3)!=substr(mtec$collapsedfamname,1,3),]
#### looking carefully, these sequences are nearly identical :(
## and more look like mutator (have a 9 bp TSD next door)                             
sum(tir$fam=='DTA00166' & tir$tsdadjacentequal)
sum(tir$fam=='DTA00273' & tir$tsdadjacentequal)
sum(tir$fam=='DTM00555' & tir$tsdadjacentequal)                            
### this fam doesn't get pulled in so maybe we don't care??
sum(tir$fam=='DTA00050' & tir$tsdadjacentequal)                           
sum(tir$fam=='DTM06396' & tir$tsdadjacentequal)                           

## manually change the first one to Mutator, as it appears it must be that (9 bp tsd, more G-rich terminal inverted repeats)
mtec$collapsedfamname[mtec$V1=='MTEC000019']='DTM00555'
                             
                             
### output info on collapsed families
collapsed=mtec[mtec$famname!=mtec$collapsedfamname,2:4]
names(collapsed)=c('mtec_fasta_name', 'mtec_name', 'family_added_to')
write.table(collapsed, '~/projects/tir_modified_mcs/mtec/MTEC_families_collapsed_into_other_families.txt', row.names=F, col.names=T, quote=F, sep='\t')


## change families
tir.gr$updatedfamname=tir.gr$famname
tir.gr$updatedfamname[tir.gr$famname %in% collapsed$mtec_name]=mapvalues(tir.gr$famname[tir.gr$famname %in% collapsed$mtec_name],
                                                                  from=collapsed$mtec_name,
                                                                  to=collapsed$family_added_to, warn_missing=F)
  
## assign 8digit family name (RL.XXXXX) and 10 digit copy name (B73v4XXXXX)
mcols(tir.gr)$Name=NULL
for (x in names(table(mcols(tir.gr)$updatedfamname))){
  mcols(tir.gr)$Name[mcols(tir.gr)$updatedfamname==x & !is.na(mcols(tir.gr)$updatedfamname)]=paste(SHORTNAME, str_pad(1:sum(mcols(tir.gr)$updatedfamname[!is.na(mcols(tir.gr)$updatedfamname)]==x), 5, pad='0'),sep='')
}

## confirm this is 0 to make sure switched correctly                             
sum(tir.gr$updatedfamname!=substr(tir.gr$ID,1,8))
                             
tir.gr$sup=substr(tir.gr$updatedfamname, 1,3)
tir.gr$ID=paste0(tir.gr$updatedfamname, tir.gr$Name)
  
  
  
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
