library(data.table)
library(rtracklayer)
library(stringr)
library(dplyr)
library(stringr)
library(plyr)

#tir=fread('~/Downloads/tir_B73_2018-07-26_extra.txt')
#tir=fread('all_tir_B73_2018-08-25_extra.txt')
#tir=tir[,c('mtec', 'chr', 'start.adj', 'end.adj', 'tirseqSingle', 'tsdadjacentup', 'fam', 'sup')]
#tir=tir[tir$start.adj-tir$end.adj<0,]

#tir=import.gff3('B73_tir_2018-08-31.gff3')
tir=import.gff3('W22_tir_2018-08-31.gff3')
end(tir)=end(tir)+1  ## for W22 2018-08-31.gff3 where the end needs to be incremented
filt=nchar(str_split_fixed(tir$ID, '_', 7)[,5])>=5 & width(tir)>=40
tir=tir[filt,]
tircols=str_split_fixed(tir$ID, '_', 7)
tir$mtec=sapply(1:length(tir), function(x) paste(unlist(tircols[x,1:3]), collapse='_'))
tir$sup=substr(tir$mtec, 1, 3)
tir$mtecfamnum=substr(tir$mtec,7,11)
tir$famname=paste0(tir$sup,tir$mtecfamnum)
tir$TSD=tircols[,4]
tir$TIRup=tircols[,5]
tir$TIRmm=tircols[,6]
tir$TIRdown=tircols[,7]
tir.gr=tir
#mcols(tir.gr)=NULL





#tirmm=import.gff3('B73_tir_2018-08-31.mismatchAll.gff3')
tirmm=import.gff3('W22_tir_2018-09-01.mismatchAll.gff3')
mmfilt=nchar(str_split_fixed(tirmm$ID, '_', 7)[,5])>=5 & width(tirmm)>=40
tirmm=tirmm[mmfilt,]
tirmmcols=str_split_fixed(tirmm$ID, '_', 7)
tirmm$mtec=sapply(1:length(tirmm), function(x) paste(unlist(tirmmcols[x,1:3]), collapse='_'))
tirmm$sup=substr(tirmm$mtec, 1, 3)
tirmm$mtecfamnum=substr(tirmm$mtec,7,11)
tirmm$famname=paste0(tirmm$sup,tirmm$mtecfamnum)
tirmm$TSD=tirmmcols[,4]
tirmm$TIRup=tirmmcols[,5]
tirmm$TIRmm=tirmmcols[,6]
tirmm$TIRdown=tirmmcols[,7]
tirmm.gr=tirmm
#mcols(tirmm.gr)=NULL
                 
tir.gr=c(tir.gr,tirmm.gr)
           
                  
### yikes! I know this says B73 but trust these are W22 - it's hard coded in the original script I used to run. But these are DEFINITELY w22, and in the w22 directory!                  
#for(dmSplit in c('B73_detectMITE_2018-09-04.1.gff3', 'B73_detectMITE_2018-09-04.1.gff3', 'B73_detectMITE_2018-09-04.1.gff3')){                 
for(dmSplit in c(paste0('B73_detectMITE_2018-09-04.',1:8,'.gff3'), paste0('B73_detectMITE_2018-09-05.',9:20,'.gff3'))){
#tirdm=import.gff3('B73_tir_2018-08-31.mismatchAll.gff3')
tirdm=import.gff3(dmSplit)
end(tirdm)=end(tirdm)+1
dmfilt=nchar(str_split_fixed(tirdm$ID, '_', 7)[,5])>=5 & width(tirdm)>=40
tirdm=tirdm[dmfilt,]
tirdmcols=str_split_fixed(tirdm$ID, '_', 7)
tirdm$mtec=sapply(1:length(tirdm), function(x) paste(unlist(tirdmcols[x,1:3]), collapse='_'))
tirdm$sup=substr(tirdm$mtec, 1, 3)
tirdm$mtecfamnum=substr(tirdm$mtec,4,8)
tirdm$famname=paste0(tirdm$sup,tirdm$mtecfamnum)
tirdm$TSD=tirdmcols[,4]
tirdm$TIRup=tirdmcols[,5]
tirdm$TIRmm=tirdmcols[,6]
tirdm$TIRdown=tirdmcols[,7]
tirdm.gr=tirdm
#mcols(tirmm.gr)=NULL

tir.gr=c(tir.gr,tirdm.gr)

} 
                  
## filters
#filters=tir$tirsmatch & tir$tsdadjacentequal & nchar(tir$tirseqSingle)>=5 & tir$end.adj-tir$start.adj>=40

#tir=tir[filters,]
  
#names(tir)=c('mtec', 'chr', 'start', 'end', 'tirscore', 'V6', 'V7', 'TIR1', 'TIR2')
#tir.gr=GRanges(seqnames=tir$chr, ranges=IRanges(start=tir$start.adj, end=tir$end.adj))## add 5 bp minimum for TIR length, that is Bergamo's
#mcols(tir.gr)$score=tir$tirscore

#strand(tir.gr)=tir$strand
                  
                  
origlen=length(tir.gr)
mmlen=length(tirmm.gr)

### length filters - incorporated into filters above!
#tir.gr=tir.gr[width(tir.gr)>=40,]## stowaway mites can be small, using 40 bp as size cutoff. This removes 154 copies from B73.

#tir.gr$mtec=c(tir$mtec,tirmm$mtec)
#tir.gr$mtecfamnum=substr(tir.gr$mtec, 7,11)
#tir.gr$sup=substr(tir.gr$mtec, 1, 3)
#tir.gr$famname=paste(tir.gr$sup, tir.gr$mtecfamnum, sep='')
#tir.gr$TSD=c(tir$tsdadjacentup,tirmm$TSD)
#tir.gr$TIRup=c(tir$tirseqSingle, tirmm$TIRup)
#tir.gr$TIRdown=c(rep('NA',nrow(tir)), tirmm$TIRdown)
#tir.gr$TIRmm=c(rep('mismatch=0',nrow(tir)), tirmm$TIRmm)



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

#SHORTNAME='Zm00001d' ## for b73
SHORTNAME='Zm00004b'  ## for w22

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
#write.table(collapsed, '~/projects/tir_modified_mcs/mtec/MTEC_families_collapsed_into_other_families.txt', row.names=F, col.names=T, quote=F, sep='\t')


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
  
  
  
#GENOMENAME='B73'
GENOMENAME='W22'
### end -1 for gff3 format!
d=data.frame(chr=seqnames(tir.gr), 'TARGeT', 'terminal_inverted_repeat_element', start(tir.gr), end(tir.gr), '.', strand(tir.gr), '.', Name=paste0('ID=', tir.gr$ID, ';Name=', tir.gr$ID, '_', tir.gr$TSD, '_', tir.gr$TIRup, '_', tir.gr$TIRmm, '_', tir.gr$TIRdown))
#d=d[tir$tsdadjacentequal & tir$tirsmatch,]
#write.table(d[!is.na(tir$whichrule) & d[,4]<d[,5],], file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
#write.table(d[d[,4]<d[,5],], file=paste0(GENOMENAME, '_unfiltered_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
write.table(d, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.', SHORTNAME, '.withDetectMITE.gff3'), col.names=F, row.names=F, sep='\t', quote=F)

## for maizegdb with Chr
dd=d
levels(dd$chr)[1:10]=paste0('Chr', levels(dd$chr)[1:10])	
write.table(dd, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.', SHORTNAME, '.withDetectMITE.Chr.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
