## mcs find TIR/TSD in sequence from TARGeT

library(Rlibstree)
library(plyr)
library(Biostrings)     # Provides DNAString, DNAStringSet, etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges, etc
library(rtracklayer)    # Provides import() and export()
library(data.table)


## read in refgen
seqs=readDNAStringSet('../agpv4_te_annotation/ncbi_pseudomolecule/B73V4.both_pseudo_AND_unplaced.fa')
#seqs=readDNAStringSet('../../../W22__Ver12.fasta')

## read in all output 
## ex line >1_B73V4.both_pseudo_AND_unplaced.fa Query:DTC_ZM00060_consensus Sbjct:10 Length:4429 Location:(93697655 - 93710698) Direction:minus
## could do with fread
## or outside: grep -h ">" */*/*.genomic > all_b73_tir_target_matches.txt
### or grep -h ">" */*/*.genomic > all_w22_tir_target_matches.txt

a=fread('all_b73_tir_target_matches.txt', header=F)
#a=fread('all_w22_tir_target_matches.txt', header=F)
#a=fread('all_w22_mtec_tir_target_matches.txt', header=F)
a$mtec=gsub('Query:', '', a$V2)
a$chr=gsub('Sbjct:', '', a$V3)
a$start=as.numeric(gsub('Location:\\(', '', a$V5))
a$end=as.numeric(gsub('\\)', '', a$V7))
a$direction=gsub('Direction:', '', a$V8)

## turn this into tir!
tir=a[,c('mtec', 'chr', 'start', 'end', 'direction')]

## make bioconductor like chromosome names
#tir$chrnew=as.character(tir$chr)
#tir$chrnew[is.na(as.numeric(as.character(tir$chr)))]=paste0('B73V4_', tir$chrnew[is.na(as.numeric(as.character(tir$chr)))])

tir$chrnew=tir$chr ## not needed anymore, but so as not to change the names later, changing here.

## just changed DTM to 9 bp TSD ala Mobile DNA II
tsdlens=data.frame(sup=c('DTA', 'DTC', 'DTH', 'DTM', 'DTT'), tsdlen=c(8, 3, 3, 9, 2))



#tir$tirlen=nchar(as.character(tir$TIR1))

tir$tsdlen=mapvalues(substr(tir$mtec,1,3), from=tsdlens$sup, to=tsdlens$tsdlen)
tir$tsdlen=as.numeric(tir$tsdlen)

tir$tirlen=mapvalues(substr(tir$mtec,1,3), from=tsdlens$sup, to=tsdlens$tsdlen)
tir$tirlen=as.numeric(tir$tirlen)


#tir$upstream=as.character(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$start-tir$tsdlen, end=tir$start+tir$tirlen))))
#tir$downstream=as.character(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$end-tir$tirlen, end=tir$end+tir$tsdlen))))

tir$downstreamRC=as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$end-tir$tirlen, end=tir$end+tir$tsdlen)))))

######
########## Introduce way to not look for overlaps within the TE

tir$origlen=tir$end-tir$start

offset=200
tir$upstreamExtra=as.character(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$start-offset, end=tir$start+offset))))
tir$downstreamExtra=as.character(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$end-offset, end=tir$end+offset))))
tir$downstreamExtraRC=as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$end-offset, end=tir$end+offset)))))
tir$upstreamExtraRC=as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$start-offset, end=tir$start+offset)))))

tir$upstreamExtra[tir$origlen <= 2*offset ] = as.character(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$start[tir$origlen <= 2*offset ] + (tir$origlen[tir$origlen <= 2*offset ]/2)))))
tir$downstreamExtra[tir$origlen <= 2*offset ] = as.character(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$end[tir$origlen <= 2*offset ]- ( tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$end[tir$origlen <= 2*offset ] + (2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)))))
tir$upstreamExtraRC[tir$origlen <= 2*offset ] = as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$start[tir$origlen <= 2*offset ] + (tir$origlen[tir$origlen <= 2*offset ]/2))))))
tir$downstreamExtraRC[tir$origlen <= 2*offset ] = as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$end[tir$origlen <= 2*offset ]- ( tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$end[tir$origlen <= 2*offset ] + (2*offset-tir$origlen[tir$origlen <= 2*offset ]/2))))))


## now do for everything!!!

## get sequence of TIRs
## note that this can get two different, comma separated TIRs. Be careful, or deal with both? For now, this just drops those.
## probably best to use the closest one to the coordinates? or test both and see if there are TSDs?
tir$tirseq=sapply(1:nrow(tir), function(x) getLongestCommonSubstring(c(tir$upstreamExtra[x], tir$downstreamExtraRC[x])))
tir$tirseqRC=sapply(1:nrow(tir), function(x) getLongestCommonSubstring(c(tir$upstreamExtraRC[x], tir$downstreamExtra[x])))

## deal with pesky double counts by removing them here. May want to revisit, especially if some of these don't generate a TSD!!!!
tir$tirseqSingle=unlist(lapply(tir$tirseq, function(l) l[[1]]))
tir$tirseqRCSingle=unlist(lapply(tir$tirseqRC, function(l) l[[length(l)]]))


## get position of TIR in forward orientation of upstream extract.
tir$tirstartup=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqSingle[x], tir$upstreamExtra[x])))

## get position of TIR in forward orientation of downstream extract.
tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x])))

## having an issue with weird characters, removing those that don't have a tir present here.
## removes about 150 copies where no match is found. 
tirorig=tir
tir=tir[tir$tirstartup != -1,]

## get TSD sequence - need to adjust offset!!!! right now is ten which is not enough for mutator
tir$tsdseq=sapply(1:nrow(tir), function(x) getLongestCommonSubstring(c(substr(tir$upstreamExtra[x], tir$tirstartup[x]-10-tir$tsdlen[x], tir$tirstartup[x]-1),  ## need the minus one to exclude the first base of the TIR
						  substr(tir$downstreamExtra[x], tir$tirstartdown[x]+ nchar(tir$tirseqSingle[x]), tir$tirstartdown[x]+nchar(tir$tirseqSingle[x])+9+tir$tsdlen[x])))
				  )
tir$tsdseq[sapply(tir$tsdseq, is.null)]=NA
tir$tsdseqSingle=unlist(lapply(tir$tsdseq, function(l) l[[1]]))

## try with a larger TSD window
tir$tsdseq20=sapply(1:nrow(tir), function(x) getLongestCommonSubstring(c(substr(tir$upstreamExtra[x], tir$tirstartup[x]-20-tir$tsdlen[x], tir$tirstartup[x]-1),  ## need the minus one to exclude the first base of the TIR
						  substr(tir$downstreamExtra[x], tir$tirstartdown[x]+ nchar(tir$tirseqSingle[x]), tir$tirstartdown[x]+nchar(tir$tirseqSingle[x])+19+tir$tsdlen[x])))
				  )
tir$tsdseq20[sapply(tir$tsdseq20, is.null)]=NA
tir$tsdseqSingle20=unlist(lapply(tir$tsdseq20, function(l) l[[1]]))


tir$tsdadjacentup=sapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tirstartup[x] - tir$tsdlen[x], tir$tirstartup[x]-1))
tir$tsdadjacentdown=sapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]), tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) + tir$tsdlen[x] -1 ))
tir$tsdadjacentequal=tir$tsdadjacentup == tir$tsdadjacentdown


### is inferred TSD the expected TSD length?

tir$tsdrightlength=nchar(tir$tsdseqSingle)==tir$tsdlen
tir$tirrightlength=nchar(tir$tirseqSingle)==tir$tirlen

## check if TSD abuts TIR
tir$tsdtirjunctionpresent=sapply(1:nrow(tir), function(x)  ## useBytes=T in grepl so there aren't locale errors when introducing weird characters with negative ranges
								grepl(paste0(tir$tsdseqSingle[x], tir$tirseqSingle[x]), tir$upstreamExtra[x], useBytes=T) & 
								grepl(paste0(tir$tirseqRCSingle[x], tir$tsdseqSingle[x]), tir$downstreamExtra[x], useBytes=T))
							
tir$tsdtirjunctionpresent20=sapply(1:nrow(tir), function(x)  ## useBytes=T in grepl so there aren't locale errors when introducing weird characters with negative ranges
								grepl(paste0(tir$tsdseqSingle20[x], tir$tirseqSingle[x]), tir$upstreamExtra[x], useBytes=T) & 
								grepl(paste0(tir$tirseqRCSingle[x], tir$tsdseqSingle20[x]), tir$downstreamExtra[x], useBytes=T))
							
										
## and that it is doing so in the expected position
##### running out of memory- abandon??
#tir$tsdtirjunctioninrightspot=sapply(1:nrow(tir), function(x)
#				regexpr(paste0(tir$tsdseq[x], tir$tirseqSingle[x]), tir$upstreamExtra[x])==tir$tirstartup[x]-nchar(tir$tsdseq) &
#				regexpr(paste0(tir$tirseqRCSingle[x], tir$tsdseq[x]), tir$downstreamExtra[x])==tir$tirstartdown)

##adjust positions
tir$start.adj=sapply(1:nrow(tir), function(x) (tir$start[x]-200) + tir$tirstartup[x] - 1 )## parentheses puts on same scale as upstreamExtra
tir$end.adj= sapply(1:nrow(tir), function(x) (tir$end[x]-200) + tir$tirstartdown[x] + nchar(tir$tirseqSingle[x]) - 2)

## add strand (relative to mtec)
tir$strand=ifelse(tir$direction=='plus', '+', '-')



tir$sup=substr(tir$mtec,1,3)
tir$fam=paste0(tir$sup, substr(tir$mtec, 7,11))



library(stringdist)
library(stringr)

tir$tsdseqSingle[is.na(tir$tsdseqSingle)]='N'

## get position of TSD in forward orientation of upstream extract.
### in a nongreedy way .*? gets minimal???
tir$tsdstartup=sapply(1:nrow(tir), function(x) as.numeric(regexpr(paste0(tir$tsdseqSingle[x], '[^', tir$tsdseqSingle[x], ']*', tir$tirseqSingle[x]), tir$upstreamExtra[x])))

## get position of TSD in forward orientation of downstream extract.
#tir$tsdstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(paste0(tir$tirseqRCSingle[x], '.*?', tir$tsdseqSingle[x]), tir$downstreamExtra[x])))

tir$tsdstartdown=sapply(1:nrow(tir), function(x) tir$tirstartdown[x] + ## this is the start position of the match
												 nchar(str_extract(tir$downstreamExtra[x], paste0(tir$tirseqRCSingle[x], '.*?', tir$tsdseqSingle[x]) ))- nchar(tir$tsdseqSingle[x])) # adjusting to get to end of new TIR (right inside of TSD)


#Check candidate TIRs
tir$adjustedTIRup=sapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tsdstartup[x]+nchar(tir$tsdseqSingle[x]), tir$tirstartup[x]+nchar(tir$tirseqSingle[x])-1))  ## need the minus one to exclude the first base of the TIR
						  
tir$adjustedTIRdown=sapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x], tir$tsdstartdown[x]-1))

tir$adjustedTIRdownRC=NA
tir$adjustedTIRdownRC[which(!is.na(tir$adjustedTIRdown))]=sapply(which(!is.na(tir$adjustedTIRdown)), function(x) as.character(reverseComplement(DNAString(tir$adjustedTIRdown[x]))))

tir$seqdist=sapply(1:nrow(tir), function(x) stringdist(tir$adjustedTIRup[x], tir$adjustedTIRdownRC[x], method='h'))

####################################
## repeat with TSD search in 20 bp up and downstream

tir$tsdseqSingle20[is.na(tir$tsdseqSingle20)]='N'
## get position of TSD in forward orientation of upstream extract.
### in a nongreedy way .*? gets minimal???
tir$tsdstartup20=sapply(1:nrow(tir), function(x) as.numeric(regexpr(paste0(tir$tsdseqSingle20[x], '[^', tir$tsdseqSingle20[x], ']*', tir$tirseqSingle[x]), tir$upstreamExtra[x])))

## get position of TSD in forward orientation of downstream extract.
tir$tsdstartdown20=sapply(1:nrow(tir), function(x) tir$tirstartdown[x] + ## this is the start position of the match
												 nchar(str_extract(tir$downstreamExtra[x], paste0(tir$tirseqRCSingle[x], '.*?', tir$tsdseqSingle20[x]) ))- nchar(tir$tsdseqSingle20[x])) # adjusting to get to end of new TIR (right inside of TSD)


#Check candidate TIRs
tir$adjustedTIRup20=sapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tsdstartup20[x]+nchar(tir$tsdseqSingle20[x]), tir$tirstartup[x]+nchar(tir$tirseqSingle[x])-1))  ## need the minus one to exclude the first base of the TIR
						  
tir$adjustedTIRdown20=sapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x], tir$tsdstartdown20[x]-1))

tir$adjustedTIRdownRC20=NA
tir$adjustedTIRdownRC20[which(!is.na(tir$adjustedTIRdown20))]=sapply(which(!is.na(tir$adjustedTIRdown20)), function(x) as.character(reverseComplement(DNAString(tir$adjustedTIRdown20[x]))))

tir$seqdist20=sapply(1:nrow(tir), function(x) stringdist(tir$adjustedTIRup20[x], tir$adjustedTIRdownRC20[x], method='h'))


## so now apply some rules

## 1. If TSD of expected length and adjacent to longest identical TIR, keep.
tir$whichrule=NA
tir$whichrule[tir$tsdadjacentequal]=1
## 2. If TSD of unexpected length and adjacent to longest identical TIR, keep.
###### amendment - but TSD must be >80% of expected length?? problem of two bp tsds for Mutator, which is wrong.
tir$whichrule[tir$tsdtirjunctionpresent & is.na(tir$whichrule) & nchar(tir$tsdseqSingle)>=0.8*tir$tsdlen]=2
## 2.5 same as 2 except with 20 bp to search for TSD. 
tir$whichrule[tir$tsdtirjunctionpresent20 & is.na(tir$whichrule) & nchar(tir$tsdseqSingle20)>=0.8*tir$tsdlen]=2.5
## 3. If extended TIR has seqdist between extended TIRs of < 0.2, keep.
tir$whichrule[tir$seqdist/nchar(tir$adjustedTIRdownRC) <= 0.2 & is.na(tir$whichrule) & nchar(tir$adjustedTIRdownRC)==nchar(tir$adjustedTIRup)& nchar(tir$tsdseqSingle)>=0.8*tir$tsdlen]=3
## 3.5 same as 3 except with 20 bp to search for TSD. 
tir$whichrule[tir$seqdist20/nchar(tir$adjustedTIRdownRC20) <= 0.2 & is.na(tir$whichrule) & nchar(tir$adjustedTIRdownRC20)==nchar(tir$adjustedTIRup20) & nchar(tir$tsdseqSingle20)>=0.8*tir$tsdlen]=3.5



## TODO --- DONE!!!!!!!!!!!!!!!!!!!
# remove negative widths, or account for TE matches shorter than 400 bp internal (such that going 200 bp in on either side gets the same sequence)



## output gffs
GENOMENAME='B73'
#GENOMENAME='W22'
d=data.frame(tir$chrnew, 'TARGeT', 'terminal_inverted_repeat_element', tir$start.adj, tir$end.adj, '.', tir$strand, '.', paste0('ID=', tir$mtec, '_', tir$tsdseqSingle, '_', tir$tirseqSingle, '_rule=', tir$whichrule))
write.table(d[!is.na(tir$whichrule),], file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
write.table(d, file=paste0(GENOMENAME, '_unfiltered_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)



 write.table(tir[,-c('tsdseq20', 'tsdseq', 'tirseqRC', 'tirseq')], paste0('all_tir_', GENOMENAME, '_extra.txt'), quote=F, sep='\t', col.names=T, row.names=F)

pdf('b73_tir_compare.pdf')
#pdf('w22_tir_compare.pdf')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj)) + geom_density()
ggplot(tir, aes(x=tir$end.adj-tir$start.adj)) + geom_density() + scale_x_log10()
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdadjacentequal)) + geom_density(alpha=0.3)
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdadjacentequal)) + geom_density(alpha=0.3) + scale_x_log10()
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdtirjunctionpresent)) + geom_density(alpha=0.3)
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdtirjunctionpresent)) + geom_density(alpha=0.3) + scale_x_log10()
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=paste0(tir$tsdtirjunctionpresent, tir$tsdadjacentequal))) + geom_density(alpha=0.3) + labs(fill='junction, adjacent')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=paste0(tir$tsdtirjunctionpresent, tir$tsdadjacentequal))) + geom_density(alpha=0.3) + labs(fill='junction, adjacent') + scale_x_log10()
## facet on superfam
ggplot(tir, aes(x=tir$end.adj-tir$start.adj)) + geom_density() + facet_wrap(~substr(mtec,1,3), ncol=1, scales='free_y')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj)) + geom_density() + scale_x_log10()+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free_y')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdadjacentequal)) + geom_density(alpha=0.3)+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free_y')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdadjacentequal)) + geom_density(alpha=0.3) + scale_x_log10()+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free_y')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdtirjunctionpresent)) + geom_density(alpha=0.3)+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free_y')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdtirjunctionpresent)) + geom_density(alpha=0.3) + scale_x_log10()+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free_y')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=paste0(tir$tsdtirjunctionpresent, tir$tsdadjacentequal))) + geom_density(alpha=0.3) + labs(fill='junction, adjacent')+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free_y')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=paste0(tir$tsdtirjunctionpresent, tir$tsdadjacentequal))) + geom_density(alpha=0.3) + labs(fill='junction, adjacent') + scale_x_log10()+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free_y')

## facet on superfam, freex
ggplot(tir, aes(x=tir$end.adj-tir$start.adj)) + geom_density() + facet_wrap(~substr(mtec,1,3), ncol=1, scales='free')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj)) + geom_density() + scale_x_log10()+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdadjacentequal)) + geom_density(alpha=0.3)+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdadjacentequal)) + geom_density(alpha=0.3) + scale_x_log10()+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdtirjunctionpresent)) + geom_density(alpha=0.3)+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=tir$tsdtirjunctionpresent)) + geom_density(alpha=0.3) + scale_x_log10()+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=paste0(tir$tsdtirjunctionpresent, tir$tsdadjacentequal))) + geom_density(alpha=0.3) + labs(fill='junction, adjacent')+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free')
ggplot(tir, aes(x=tir$end.adj-tir$start.adj, fill=paste0(tir$tsdtirjunctionpresent, tir$tsdadjacentequal))) + geom_density(alpha=0.3) + labs(fill='junction, adjacent') + scale_x_log10()+ facet_wrap(~substr(mtec,1,3), ncol=1, scales='free')


dev.off()

tira=tir %>% group_by(fam, tirseqSingle, sup) %>% summarize(n=n())

pdf('b73_tir_summary.pdf',16,8)
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle))), 100),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity') +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle))), 100),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTA']))), 100),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTA')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTA']))), 100),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTA')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTC']))), 100),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTC')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTC']))), 100),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTC')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTH']))), 100),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTH')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTH']))), 100),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTH')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTM']))), 100),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTM')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTM']))), 100),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTM')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTT']))), 100),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTT')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTT']))), 100),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTT')



dev.off()





pdf('b73_tir_summary.30.pdf',12,8)
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle))), 30),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity') +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle))), 30),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTA']))), 30),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTA')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTA']))), 30),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTA')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTC']))), 30),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTC')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTC']))), 30),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTC')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTH']))), 30),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTH')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTH']))), 30),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTH')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTM']))), 30),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTM')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTM']))), 30),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTM')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTT']))), 30),], aes(x=tirseqSingle, y=n, fill=sup)) + geom_bar(stat='identity')  +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTT')
ggplot(tira[tira$tirseqSingle %in% tail(names(sort(table(tir$tirseqSingle[tir$sup=='DTT']))), 30),], aes(x=tirseqSingle, y=n, fill=fam)) + geom_bar(stat='identity') + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('DTT')



dev.off()



