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

tir$chrnew=tir$chr ## not needed anymore, but so as not to change the names later, changing here.

## which TSD lengths do we expect for each superfamily?
## just changed DTM to 9 bp TSD ala Mobile DNA II
tsdlens=data.frame(sup=c('DTA', 'DTC', 'DTH', 'DTM', 'DTT'), tsdlen=c(8, 3, 3, 9, 2))
tir$tsdlen=mapvalues(substr(tir$mtec,1,3), from=tsdlens$sup, to=tsdlens$tsdlen)
tir$tsdlen=as.numeric(tir$tsdlen)

########## 
########## Define a region surrounding each end of the TE to search for TIRs
#### this is 400 bp at each TE boundary, centered on the edge of the TARGeT call
### note that when the TE is less than 400 bp long, I still search in 400 bp, just assign half of the TE to upstream and half to downstream.
########## 


tir$origlen=tir$end-tir$start ## note that none of these are negative!

offset=200
tir$upstreamExtra=as.character(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$start-offset, end=tir$start+offset))))
tir$downstreamExtra=as.character(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$end-offset, end=tir$end+offset))))
tir$downstreamExtraRC=as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$end-offset, end=tir$end+offset)))))
tir$upstreamExtraRC=as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$start-offset, end=tir$start+offset)))))

tir$upstreamExtra[tir$origlen <= 2*offset ] = as.character(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$start[tir$origlen <= 2*offset ] + (tir$origlen[tir$origlen <= 2*offset ]/2)))))
tir$downstreamExtra[tir$origlen <= 2*offset ] = as.character(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$end[tir$origlen <= 2*offset ]- ( tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$end[tir$origlen <= 2*offset ] + (2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)))))
tir$upstreamExtraRC[tir$origlen <= 2*offset ] = as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$start[tir$origlen <= 2*offset ] + (tir$origlen[tir$origlen <= 2*offset ]/2))))))
tir$downstreamExtraRC[tir$origlen <= 2*offset ] = as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$end[tir$origlen <= 2*offset ]- ( tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$end[tir$origlen <= 2*offset ] + (2*offset-tir$origlen[tir$origlen <= 2*offset ]/2))))))

### the previous all look great, don't think I use the RC's???

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
## removes about 150 copies where no match is found. ## now 871 for full b73
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

#tir$tsdadjacentup3=sapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tirstartup[x] + 3 - tir$tsdlen[x], tir$tirstartup[x]-1 +3 ))
#tir$tsdadjacentdown3=sapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) -3 , tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) + tir$tsdlen[x] -1 -3))
#tir$tsdadjacentequal3=tir$tsdadjacentup3 == tir$tsdadjacentdown3

			    
			    
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
		    return(list(tirseq, upposns, downposns, uptsds, downtsds, tsdsequal, tsdtirjunctionpresent, tirseqRC))
                    })

## okay, so not inconsequential number of copies with multiple possible adjacent TIRs
table(sapply(1:length(temp), function(x) sum(temp[[x]][[7]]))) ## this is always at least 2, because there are at least 2 TIR candidates here!
table(sapply(1:length(temp), function(x) any(temp[[x]][[7]]))) # the only one left is the exclamation point disaster - ignore it??
## oh wait, this needs to be both adjacent AND equal to each other
table(sapply(1:length(temp), function(x) sum(temp[[x]][[7]] & temp[[x]][[6]])))


######## PICK ONE from temp to add to tir's columns!
## replace multis with the best TIR, abandoning those with >1 best TIR/TSD pair
tir$tirseqSingle[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[1]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
tir$tirseqRCSingle[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[8]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
#tir$tsdseq[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[5]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
tir$tsdadjacentup[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[5]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
tir$tsdadjacentdown[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[5]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
tir$tsdadjacentequal[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, T, F))
#tir$tsdtirjunctionpresent[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, T, F))
#tir$tsdstartup[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[3]][temp[[x]][[7]] & temp[[x]][[6]]], -1))
#tir$tsdstartdown[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[4]][temp[[x]][[7]] & temp[[x]][[6]]]-nchar(temp[[x]][[1]][1]), -1))

							   
###########
### OKAY! Now we have real TIRs, and can do stuff with them!!!
###    use TIR's RC to find in reverse sequence??
###########


## get position of TIR in forward orientation of upstream extract.
tir$tirstartup=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqSingle[x], tir$upstreamExtra[x])))

## get position of TIR in forward orientation of downstream extract.
# to do RC: as.character(reverseComplement(DNAString(tirF)))
## and weird char introduced?: sapply(tirseq, function(tirF) tryCatch({as.character(reverseComplement(DNAString(tirF)))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}))
#tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x]))-1)
#tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tryCatch({as.character(reverseComplement(DNAString(tir$tirseqRCSingle[x])))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}), tir$downstreamExtra[x]))-1)
tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x]))-1)


## having an issue with weird characters, removing those that don't have a tir present here.
## removes about 150 copies where no match is found. 
#tirorig=tir
#tir=tir[tir$tirstartup != -1,]
#tir=tir[tir$tirstartdown != -2,]


############
## Now, adjust positions to put them back on the same scale as the genome!!!
############
##adjust positions
tir$start.adj=sapply(1:nrow(tir), function(x) (tir$start[x]-offset) + tir$tirstartup[x] - 1 )## parentheses puts on same scale as upstreamExtra
tir$end.adj= sapply(1:nrow(tir), function(x) (tir$end[x]-offset) + tir$tirstartdown[x] + nchar(tir$tirseqSingle[x]) - 1)

## floor is what GRanges does to decimal values, so replicate this here (e.g. 0.5 becomes 0, 3.5 becomes 3)
#  tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)
tir$start.adj[tir$origlen <= 2*offset ]= sapply(which(tir$origlen <= 2*offset), function(x) (tir$start[x] - floor(2*offset - tir$origlen[x]/2) + tir$tirstartup[x] - 1 ))
tir$end.adj[tir$origlen <= 2*offset ]=   sapply(which(tir$origlen <= 2*offset), function(x) (tir$end[x] - floor(tir$origlen[x]/2) + tir$tirstartdown[x] + nchar(tir$tirseqSingle[x]) - 1))
#as.character(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$start[tir$origlen <= 2*offset ] + (tir$origlen[tir$origlen <= 2*offset ]/2)))))
##tir$downstreamExtra[tir$origlen <= 2*offset ] = as.character(getSeq(seqs, GRanges(tir$chrnew[tir$origlen <= 2*offset ], IRanges(start=tir$end[tir$origlen <= 2*offset ]- ( tir$origlen[tir$origlen <= 2*offset ]/2), end=tir$end[tir$origlen <= 2*offset ] + (2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)))))



## add strand (relative to mtec)
tir$strand=ifelse(tir$direction=='plus', '+', '-')

tir$sup=substr(tir$mtec,1,3)
tir$fam=paste0(tir$sup, substr(tir$mtec, 7,11))

						
						
################
## check that i have what i think i have
################						

tir$tirupinseq=	as.character(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$start.adj, end=tir$start.adj+nchar(tir$tirseqSingle)-1))))				
tir$tirdowninseqRC=as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$end.adj-nchar(tir$tirseqSingle)+1, end=tir$end.adj)))))
																		      
tir$tirsmatch=tir$tirseqSingle==tir$tirupinseq & tir$tirseqSingle==tir$tirdowninseq & tir$tirupinseq!=''	


################
## output gffs
################
## only keep those where TSDs are equal, adjacent to TIR!!!
#tir=tir[tir$tsdadjacentequal,]

GENOMENAME='B73'
#GENOMENAME='W22'
### end -1 for gff3 format!
d=data.frame(tir$chrnew, 'TARGeT', 'terminal_inverted_repeat_element', tir$start.adj, tir$end.adj-1, '.', tir$strand, '.', paste0('ID=', tir$mtec, '_', tir$tsdadjacentup, '_', tir$tirseqSingle))
d=d[tir$tsdadjacentequal & tir$tirsmatch,]
#write.table(d[!is.na(tir$whichrule) & d[,4]<d[,5],], file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
#write.table(d[d[,4]<d[,5],], file=paste0(GENOMENAME, '_unfiltered_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
write.table(d, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)

## for maizegdb with Chr
dd=d						
levels(dd$tir.chrnew)[1:10]=paste0('Chr', levels(dd$tir.chrnew)[1:10])	
write.table(dd, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.Chr.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
						
 write.table(tir[,-c( 'tirseqRC', 'tirseq')], paste0('all_tir_', GENOMENAME, '_', Sys.Date(), '_extra.txt'), quote=F, sep='\t', col.names=T, row.names=F)


### 
						
						
#########################						
#### TIR mismatches  ####
#########################
tirm=tir[!tir$tsdadjacentequal & !tir$tirsmatch,]

## deal with pesky multiple TIRs by removing them here. 
tirm$tirseqSingle=unlist(lapply(tirm$tirseq, function(l) l[[1]]))
tirm$tirseqRCSingle=unlist(lapply(tirm$tirseqRC, function(l) l[[length(l)]]))  ## not good - this can be either first or second! AAAHHHHHHHHHH
tirm=tirm[!grepl('N', tirm$tirseqSingle),]

## get position of TIR in forward orientation of upstream extract.
tirm$tirstartup=sapply(1:nrow(tirm), function(x) as.numeric(regexpr(tirm$tirseqSingle[x], tirm$upstreamExtra[x])))

## get position of TIR in forward orientation of downstream extract.
tirm$tirstartdown=sapply(1:nrow(tirm), function(x) as.numeric(regexpr(tirm$tirseqRCSingle[x], tirm$downstreamExtra[x])))
				 
## having an issue with weird characters, removing those that don't have a tir present here.
## removes about 150 copies where no match is found. ## now 871 for full b73
#tirorig=tir
tirm=tirm[tirm$tirstartup != -1,]
			 
			 
#### whereas above, I searched inwards to find TSD accidentally called TIR, here			 
#### search outside to find a matching TSD in a loop!	
### so essentially, find closest TSD possible by stepping out one bp at a time.
			 
tirm$closestTSDseq=NA    ## this is the sequence of the matching TSD
tirm$closestTSDoffset=NA ## this is the offset from the perfect TIR edge
for(tirposoffset in 0:20){
	tsdadjacentup=sapply(1:nrow(tirm), function(x) substr(tirm$upstreamExtra[x], tirm$tirstartup[x] - tirm$tsdlen[x]-tirposoffset , tirm$tirstartup[x]-1 -tirposoffset ))
	tsdadjacentdown=sapply(1:nrow(tirm), function(x) substr(tirm$downstreamExtra[x], tirm$tirstartdown[x]  + nchar(tirm$tirseqSingle[x])+tirposoffset , tirm$tirstartdown[x]  + nchar(tirm$tirseqSingle[x]) + tirm$tsdlen[x] -1 +tirposoffset ))
	tsdadjacentequal=tsdadjacentup == tsdadjacentdown & tsdadjacentequal !=''
	print(sum(tsdadjacentequal))
	tirm$closestTSDseq[tsdadjacentequal & is.na(tirm$closestTSDseq)]=tsdadjacentup[tsdadjacentequal & is.na(tirm$closestTSDseq)]
	tirm$closestTSDoffset[tsdadjacentequal & is.na(tirm$closestTSDoffset)]=tirposoffset
	}
			 
### next, have to put in these offsets into the actual positions! and include the TSD sequence
### then, check the TIRs they suggest, and see how many mismatches			       
			       
## start to check TIR proposed by these new TSDs
#### replace with these new values of semi-shifted TIR coordinates - if TSD is pallindromic, it will be incorporated to the TIR. We don't want this!		    
tir$tsdadjacentup[tir$tsdadjacentequal1]=tir$tsdadjacentup1[tir$tsdadjacentequal1]			    			    
tir$tsdadjacentdown[tir$tsdadjacentequal1]=tir$tsdadjacentdown1[tir$tsdadjacentequal1]			    			    
tir$tsdadjacentequal[tir$tsdadjacentequal1]=T
tir$tirseqSingle[tir$tsdadjacentequal1]=sapply(tir$tirseqSingle[tir$tsdadjacentequal1], function(tirtoshorten) substr(tirtoshorten, 2, nchar(tirtoshorten)))    
tir$tirseqRCSingle[tir$tsdadjacentequal1]=sapply(tir$tirseqRCSingle[tir$tsdadjacentequal1], function(tirtoshorten) substr(tirtoshorten, 1, nchar(tirtoshorten)-1))
			       
			       
			       
			       
