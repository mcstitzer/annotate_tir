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

## read in gff output
# 1       TARGeT  terminal_inverted_repeat_element        37697   42585   .       .       .       ID=DTA00167Zm00001d00006;Name=DTA00167Zm00001d00006_TSDlen8_TIRlen11
# 1       detectMITE      terminal_inverted_repeat_element        975732  975820  .       .       .       ID=DTT10880Zm00001d00001;Name=DTT10880Zm00001d00001_TSDlen3_TIRlen13

##  need to be able to get out superfamily, so I get predicted TSD
## This expects ID field in the 9th column to be:
# ID=DTA00167Zm00001d00006;Name=DTA00167Zm00001d00006_TSDlen8_TIRlen11


a=fread('~/projects/b73_ecology/RawData/tes/B73v4.Zm00001d.allTE.gff3', header=F)
#a=fread('all_w22_tir_target_matches.txt', header=F)
#a=fread('all_w22_mtec_tir_target_matches.txt', header=F)
a$fam=substr(a$V9, 4,11)
a$sup=substr(a$V9, 4,6)
a$chr=as.character(a$V1)
a$start=as.numeric(as.character(a$V4))
a$end=as.numeric(as.character(a$V5))
#a$direction=gsub('Direction:', '', a$V8) # could get this, but absent in this gff3
a$direction=NA

## turn this into tir!
tir=a[a$V2=='TARGeT' & a$V3=='terminal_inverted_repeat_element' ,c('sup', 'fam', 'chr', 'start', 'end', 'direction')]  ### for now, remove all but TARGeT - should be able to parse detectMITE too. 

tir$chrnew=tir$chr ## not needed anymore, but so as not to change the names later, changing here.

## which TSD lengths do we expect for each superfamily?
## just changed DTM to 9 bp TSD ala Mobile DNA II
tsdlens=data.frame(sup=c('DTA', 'DTC', 'DTH', 'DTM', 'DTT'), tsdlen=c(8, 3, 3, 9, 2))
tir$tsdlen=mapvalues(tir$sup, from=tsdlens$sup, to=tsdlens$tsdlen)
tir$tsdlen=as.numeric(tir$tsdlen)

########## 
########## Define a region surrounding each end of the TE to search for TIRs
#### this is 400 bp at each TE boundary, centered on the edge of the TARGeT call
### note that when the TE is less than 400 bp long, I still search in 400 bp, just assign half of the TE to upstream and half to downstream.
########## 


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

########## 
## get sequence of TIRs
## note that this can get multiple different, comma separated TIRs.
## this tests all candidates to see if there are TSDs
########## 
tir$tirseq=sapply(1:nrow(tir), function(x) getLongestCommonSubstring(c(tir$upstreamExtra[x], tir$downstreamExtraRC[x])))
tir$tirseqRC=sapply(1:nrow(tir), function(x) getLongestCommonSubstring(c(tir$upstreamExtraRC[x], tir$downstreamExtra[x])))

## deal with pesky multiple TIRs by removing them here. 
tir$tirseqSingle=unlist(lapply(tir$tirseq, function(l) l[[1]]))
tir$tirseqRCSingle=unlist(lapply(tir$tirseqRC, function(l) l[[length(l)]]))  ## not good - this can be either first or second! AAAHHHHHHHHHH
tir=tir[!grepl('N', tir$tirseqSingle),]

## get position of TIR in forward orientation of upstream extract.
tir$tirstartup=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqSingle[x], tir$upstreamExtra[x])))

## get position of TIR in forward orientation of downstream extract.
tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x])))
				 
## having an issue with weird characters, removing those that don't have a tir present here.
## removes about 150 copies where no match is found. 
#tirorig=tir
tir=tir[tir$tirstartup != -1,]
				 
				 
				 
tir$tsdadjacentup=sapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tirstartup[x] - tir$tsdlen[x], tir$tirstartup[x]-1))
tir$tsdadjacentdown=sapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]), tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) + tir$tsdlen[x] -1 ))
tir$tsdadjacentequal=tir$tsdadjacentup == tir$tsdadjacentdown

## shorten TIR by 1 bp in case the TSD is pallindromic and getting incorporated!
### start with 1 bp to see how much it changes - in theory, could reduce more.
tir$tsdadjacentup1=sapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tirstartup[x] + 1 - tir$tsdlen[x], tir$tirstartup[x]-1 +1 ))
tir$tsdadjacentdown1=sapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) -1 , tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) + tir$tsdlen[x] -1 -1))
tir$tsdadjacentequal1=tir$tsdadjacentup1 == tir$tsdadjacentdown1
			   
tir$tsdadjacentup2=sapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tirstartup[x] + 2 - tir$tsdlen[x], tir$tirstartup[x]-1 +2 ))
tir$tsdadjacentdown2=sapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) -2 , tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) + tir$tsdlen[x] -1 -2))
tir$tsdadjacentequal2=tir$tsdadjacentup2 == tir$tsdadjacentdown2


#### replace with these new values of semi-shifted TIR coordinates - if TSD is pallindromic, it will be incorporated to the TIR. We don't want this!		    
tir$tsdadjacentup[tir$tsdadjacentequal1]=tir$tsdadjacentup1[tir$tsdadjacentequal1]			    			    
tir$tsdadjacentdown[tir$tsdadjacentequal1]=tir$tsdadjacentdown1[tir$tsdadjacentequal1]			    			    
tir$tsdadjacentequal[tir$tsdadjacentequal1]=T
tir$tirseqSingle[tir$tsdadjacentequal1]=sapply(tir$tirseqSingle[tir$tsdadjacentequal1], function(tirtoshorten) substr(tirtoshorten, 2, nchar(tirtoshorten)))    
tir$tirseqRCSingle[tir$tsdadjacentequal1]=sapply(tir$tirseqRCSingle[tir$tsdadjacentequal1], function(tirtoshorten) substr(tirtoshorten, 1, nchar(tirtoshorten)-1))
## and 2 offset				       
tir$tsdadjacentup[tir$tsdadjacentequal2]=tir$tsdadjacentup2[tir$tsdadjacentequal2]			    			    
tir$tsdadjacentdown[tir$tsdadjacentequal2]=tir$tsdadjacentdown2[tir$tsdadjacentequal2]			    			    
tir$tsdadjacentequal[tir$tsdadjacentequal2]=T
tir$tirseqSingle[tir$tsdadjacentequal2]=sapply(tir$tirseqSingle[tir$tsdadjacentequal2], function(tirtoshorten) substr(tirtoshorten, 3, nchar(tirtoshorten)))    
tir$tirseqRCSingle[tir$tsdadjacentequal2]=sapply(tir$tirseqRCSingle[tir$tsdadjacentequal2], function(tirtoshorten) substr(tirtoshorten, 1, nchar(tirtoshorten)-2))






### actually deal with pesky multiple TIRs by searching each for a TSD here.


temp=lapply(which(sapply(tir$tirseq, length)>1), function(x) {
#		     print(x)
		    tirseq=as.character(tir$tirseq[[x]])
		    tirseqRC=sapply(tirseq, function(tirF) tryCatch({as.character(reverseComplement(DNAString(tirF)))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}))
                    upposns = as.numeric(sapply(tirseq, function(tirF) regexpr(tirF, tir$upstreamExtra[x])))
#                    downposns=as.numeric(sapply(tir$tirseqRC[[x]], function(te) regexpr(te, tir$downstreamExtra[x])))
                    downposns=as.numeric(sapply(tirseqRC, function(tirR) regexpr(tirR,tir$downstreamExtra[x]))) + sapply(tirseq, function(te) nchar(te)) ## this should be the END of the TIR/TE!!!
                    uptsds=sapply(upposns, function(uppos) substr(tir$upstreamExtra[x], uppos - tir$tsdlen[x], uppos-1))
                    downtsds=sapply(1:length(downposns), function(downpos) substr(tir$downstreamExtra[x], downposns[downpos] , downposns[downpos]  + tir$tsdlen[x] -1)) # was double counting tir length??
		    tsdsequal= uptsds==downtsds & uptsds!='' & downtsds!='' ## account for empty seqs
                    tsdtirjunctionpresent=sapply(1:length(upposns), function(index)  ## useBytes=T in grepl so there aren't locale errors when introducing weird characters with negative ranges
								grepl(paste0(uptsds[[index]], tirseq[index]), tir$upstreamExtra[x], useBytes=T) & 
								grepl(paste0(tryCatch({as.character(reverseComplement(DNAString(tirseq[index])))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}), downtsds[index]), tir$downstreamExtra[x], useBytes=T))
                    uptsds1=sapply(upposns, function(uppos) substr(tir$upstreamExtra[x], uppos - tir$tsdlen[x] +1, uppos-1 +1))
                    downtsds1=sapply(1:length(downposns), function(downpos) substr(tir$downstreamExtra[x], downposns[downpos] -1, downposns[downpos]  + tir$tsdlen[x] -1 -1)) # was double counting tir length??
		    tsdsequal1= uptsds1==downtsds1 & uptsds1!='' & downtsds1!='' ## account for empty seqs
                    uptsds2=sapply(upposns, function(uppos) substr(tir$upstreamExtra[x], uppos - tir$tsdlen[x] +2, uppos-1 +2))
                    downtsds2=sapply(1:length(downposns), function(downpos) substr(tir$downstreamExtra[x], downposns[downpos] -2, downposns[downpos]  + tir$tsdlen[x] -1 -2)) # was double counting tir length??
		    tsdsequal2= uptsds2==downtsds2 & uptsds2!='' & downtsds2!='' ## account for empty seqs
		     
				     
		    if(sum(tsdsequal1)==1){uptsds=uptsds1; downtsds=downtsds1; tsdsequal=tsdsequal1
					  tirseq=sapply(tirseq, function(tirtoshorten) substr(tirtoshorten, 2, nchar(tirtoshorten)))
					  tirseqRC=sapply(tirseq, function(tirtoshorten) substr(tirtoshorten, 1, nchar(tirtoshorten)-1))
					  upposns=upposns[tsdsequal1 & sum(tsdsequal1)==1] +1
					  downposns=as.numeric(sapply(tirseqRC, function(tirR) regexpr(tirR,tir$downstreamExtra[x]))) + sapply(tirseq, function(te) nchar(te)) -1
					}else if(sum(tsdsequal2)==1){uptsds=uptsds2; downtsds=downtsds2; tsdsequal=tsdsequal2
					  tirseq=sapply(tirseq, function(tirtoshorten) substr(tirtoshorten, 3, nchar(tirtoshorten)))
					  tirseqRC=sapply(tirseq, function(tirtoshorten) substr(tirtoshorten, 1, nchar(tirtoshorten)-2))
					  upposns=upposns[tsdsequal1 & sum(tsdsequal1)==1] +2
					  downposns=as.numeric(sapply(tirseqRC, function(tirR) regexpr(tirR,tir$downstreamExtra[x]))) + sapply(tirseq, function(te) nchar(te)) -2
						}
		    return(list(tirseq, upposns, downposns, uptsds, downtsds, tsdsequal, tsdtirjunctionpresent))
                    })

## okay, so not inconsequential number of copies with multiple possible adjacent TIRs
table(sapply(1:length(temp), function(x) sum(temp[[x]][[7]]))) ## this is always at least 2, because there are at least 2 TIR candidates here!
table(sapply(1:length(temp), function(x) any(temp[[x]][[7]]))) # the only one left is the exclamation point disaster - ignore it??
## oh wait, this needs to be both adjacent AND equal to each other
table(sapply(1:length(temp), function(x) sum(temp[[x]][[7]] & temp[[x]][[6]])))


######## PICK ONE from temp to add to tir's columns!
## replace multis with the best TIR, abandoning those with >1 best TIR/TSD pair
tir$tirseqSingle[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[1]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
tir$tirseqRCSingle[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[2]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
#tir$tsdseq[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[5]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
tir$tsdadjacentup[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[5]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
tir$tsdadjacentdown[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[5]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
tir$tsdadjacentequal[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, T, F))
#tir$tsdtirjunctionpresent[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, T, F))
#tir$tsdstartup[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[3]][temp[[x]][[7]] & temp[[x]][[6]]], -1))
#tir$tsdstartdown[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[4]][temp[[x]][[7]] & temp[[x]][[6]]]-nchar(temp[[x]][[1]][1]), -1))


###########
### OKAY! Now we have real TIRs, and can do stuff with them!!!
###########


## get position of TIR in forward orientation of upstream extract.
tir$tirstartup=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqSingle[x], tir$upstreamExtra[x])))

## get position of TIR in forward orientation of downstream extract.
tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x])))

## having an issue with weird characters, removing those that don't have a tir present here.
## removes about 150 copies where no match is found. 
tirorig=tir
tir=tir[tir$tirstartup != -1,]



############
## Now, adjust positions to put them back on the same scale as the genome!!!
############
##adjust positions
tir$start.adj=sapply(1:nrow(tir), function(x) (tir$start[x]-offset) + tir$tirstartup[x] - 1 )## parentheses puts on same scale as upstreamExtra
tir$end.adj= sapply(1:nrow(tir), function(x) (tir$end[x]-offset) + tir$tirstartdown[x] + nchar(tir$tirseqSingle[x]) - 2)

## NEED TO CONFIRM THIS FLOOR IS WORKING RIGHT!!!!
tir$start.adj[tir$origlen <= 2*offset ]= sapply(which(tir$origlen <= 2*offset), function(x) (tir$start[x]-  offset - floor(tir$origlen[x]/2) + tir$tirstartup[x] - 1 ))
tir$end.adj[tir$origlen <= 2*offset ]=   sapply(which(tir$origlen <= 2*offset), function(x) (tir$end[x] - floor(tir$origlen[x]/2) + tir$tirstartdown[x] + nchar(tir$tirseqSingle[x]) - 2))
#as.character(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$start[tir$origlen <= 2*offset ] + (tir$origlen[tir$origlen <= 2*offset ]/2)))))
##tir$downstreamExtra[tir$origlen <= 2*offset ] = as.character(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$end[tir$origlen <= 2*offset ]- ( tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$end[tir$origlen <= 2*offset ] + (2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)))))



## add strand (relative to mtec)
tir$strand=ifelse(tir$direction=='plus', '+', ifelse(tir$direction=='minus', '-', NA))




#tir$sup=substr(tir$mtec,1,3)
#tir$fam=paste0(tir$sup, substr(tir$mtec, 7,11))

################
## output gffs
################
## only keep those where TSDs are equal, adjacent to TIR!!!
tir=tir[tir$tsdadjacentequal,]

GENOMENAME='B73'
#GENOMENAME='W22'
d=data.frame(tir$chrnew, 'TARGeT', 'terminal_inverted_repeat_element', tir$start.adj, tir$end.adj, '.', tir$strand, '.', paste0('ID=', tir$mtec, '_', tir$tsdadjacentup, '_', tir$tirseqSingle))
write.table(d[!is.na(tir$whichrule) & d[,4]<d[,5],], file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
write.table(d[d[,4]<d[,5],], file=paste0(GENOMENAME, '_unfiltered_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
write.table(d[d[,4]<d[,5] & tir$tsdadjacentequal,], file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)


 write.table(tir[,-c( 'tirseqRC', 'tirseq')], paste0('all_tir_', GENOMENAME, '_', Sys.Date(), '_extra.txt'), quote=F, sep='\t', col.names=T, row.names=F)


### 
