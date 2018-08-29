## mcs find TIR/TSD in sequence from TARGeT

library(Rlibstree)
library(plyr)
library(Biostrings)     # Provides DNAString, DNAStringSet, etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges, etc
library(rtracklayer)    # Provides import() and export()
library(data.table)
library(stringdist)
library(stringr)
library(parallel)

ncores=64


## read in refgen
#seqs=readDNAStringSet('../agpv4_te_annotation/ncbi_pseudomolecule/B73V4.both_pseudo_AND_unplaced.fa')
seqs=readDNAStringSet('~/projects/maize_te_counts/mapping/Zea_mays.AGPv4.dna.toplevel.fa')
names(seqs)=str_split_fixed(names(seqs), ' ', 2)[,1]
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
#tir$tirstartup=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqSingle[x], tir$upstreamExtra[x])))
tir$tirstartup=unlist(mclapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqSingle[x], tir$upstreamExtra[x])), mc.cores=ncores))

				 
## get position of TIR in forward orientation of downstream extract.
#tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x])))
tir$tirstartdown=unlist(mclapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x])), mc.cores=ncores))
				 
## having an issue with weird characters, removing those that don't have a tir present here.
## removes about 150 copies where no match is found. ## now 871 for full b73
#tirorig=tir
tir=tir[tir$tirstartup != -1,]
				 
				 
				 
tir$tsdadjacentup=unlist(mclapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tirstartup[x] - tir$tsdlen[x], tir$tirstartup[x]-1), mc.cores=ncores))
tir$tsdadjacentdown=unlist(mclapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]), tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) + tir$tsdlen[x] -1 ), mc.cores=ncores))
tir$tsdadjacentequal=tir$tsdadjacentup == tir$tsdadjacentdown

## shorten TIR by 1 bp in case the TSD is pallindromic and getting incorporated!
### start with 1 bp to see how much it changes - in theory, could reduce more.
tir$tsdadjacentup1=unlist(mclapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tirstartup[x] + 1 - tir$tsdlen[x], tir$tirstartup[x]-1 +1 ), mc.cores=ncores))
tir$tsdadjacentdown1=unlist(mclapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) -1 , tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) + tir$tsdlen[x] -1 -1), mc.cores=ncores))
tir$tsdadjacentequal1=tir$tsdadjacentup1 == tir$tsdadjacentdown1
			   
tir$tsdadjacentup2=unlist(mclapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tirstartup[x] + 2 - tir$tsdlen[x], tir$tirstartup[x]-1 +2 ), mc.cores=ncores))
tir$tsdadjacentdown2=unlist(mclapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) -2 , tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) + tir$tsdlen[x] -1 -2), mc.cores=ncores))
tir$tsdadjacentequal2=tir$tsdadjacentup2 == tir$tsdadjacentdown2

#tir$tsdadjacentup3=sapply(1:nrow(tir), function(x) substr(tir$upstreamExtra[x], tir$tirstartup[x] + 3 - tir$tsdlen[x], tir$tirstartup[x]-1 +3 ))
#tir$tsdadjacentdown3=sapply(1:nrow(tir), function(x) substr(tir$downstreamExtra[x], tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) -3 , tir$tirstartdown[x]  + nchar(tir$tirseqSingle[x]) + tir$tsdlen[x] -1 -3))
#tir$tsdadjacentequal3=tir$tsdadjacentup3 == tir$tsdadjacentdown3

			    
			    
#### replace with these new values of semi-shifted TIR coordinates - if TSD is pallindromic, it will be incorporated to the TIR. We don't want this!		    
tir$tsdadjacentup[tir$tsdadjacentequal1]=tir$tsdadjacentup1[tir$tsdadjacentequal1]			    			    
tir$tsdadjacentdown[tir$tsdadjacentequal1]=tir$tsdadjacentdown1[tir$tsdadjacentequal1]			    			    
tir$tsdadjacentequal[tir$tsdadjacentequal1]=T
tir$tirseqSingle[tir$tsdadjacentequal1]=unlist(mclapply(tir$tirseqSingle[tir$tsdadjacentequal1], function(tirtoshorten) substr(tirtoshorten, 2, nchar(tirtoshorten)), mc.cores=ncores))    
tir$tirseqRCSingle[tir$tsdadjacentequal1]=unlist(mclapply(tir$tirseqRCSingle[tir$tsdadjacentequal1], function(tirtoshorten) substr(tirtoshorten, 1, nchar(tirtoshorten)-1), mc.cores=ncores))
## and 2 offset				       
tir$tsdadjacentup[tir$tsdadjacentequal2]=tir$tsdadjacentup2[tir$tsdadjacentequal2]			    			    
tir$tsdadjacentdown[tir$tsdadjacentequal2]=tir$tsdadjacentdown2[tir$tsdadjacentequal2]			    			    
tir$tsdadjacentequal[tir$tsdadjacentequal2]=T
tir$tirseqSingle[tir$tsdadjacentequal2]=unlist(mclapply(tir$tirseqSingle[tir$tsdadjacentequal2], function(tirtoshorten) substr(tirtoshorten, 3, nchar(tirtoshorten)), mc.cores=ncores))    
tir$tirseqRCSingle[tir$tsdadjacentequal2]=unlist(mclapply(tir$tirseqRCSingle[tir$tsdadjacentequal2], function(tirtoshorten) substr(tirtoshorten, 1, nchar(tirtoshorten)-2), mc.cores=ncores))






### actually deal with pesky multiple TIRs by searching each for a TSD here.


temp=mclapply(which(sapply(tir$tirseq, length)>1), function(x) {
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
                    }, mc.cores=ncores)

## okay, so not inconsequential number of copies with multiple possible adjacent TIRs
table(sapply(1:length(temp), function(x) sum(temp[[x]][[7]]))) ## this is always at least 2, because there are at least 2 TIR candidates here!
table(sapply(1:length(temp), function(x) any(temp[[x]][[7]]))) # the only one left is the exclamation point disaster - ignore it??
## oh wait, this needs to be both adjacent AND equal to each other
table(sapply(1:length(temp), function(x) sum(temp[[x]][[7]] & temp[[x]][[6]])))


######## PICK ONE from temp to add to tir's columns!
## replace multis with the best TIR, abandoning those with >1 best TIR/TSD pair
tir$tirseqSingle[which(sapply(tir$tirseq, length)>1)] = unlist(mclapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[1]][temp[[x]][[7]] & temp[[x]][[6]]], ''), mc.cores=ncores))
tir$tirseqRCSingle[which(sapply(tir$tirseq, length)>1)] = unlist(mclapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[8]][temp[[x]][[7]] & temp[[x]][[6]]], ''), mc.cores=ncores))
#tir$tsdseq[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[5]][temp[[x]][[7]] & temp[[x]][[6]]], ''))
tir$tsdadjacentup[which(sapply(tir$tirseq, length)>1)] = unlist(mclapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[5]][temp[[x]][[7]] & temp[[x]][[6]]], ''), mc.cores=ncores))
tir$tsdadjacentdown[which(sapply(tir$tirseq, length)>1)] = unlist(mclapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[5]][temp[[x]][[7]] & temp[[x]][[6]]], ''), mc.cores=ncores))
tir$tsdadjacentequal[which(sapply(tir$tirseq, length)>1)] = unlist(mclapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, T, F), mc.cores=ncores))
#tir$tsdtirjunctionpresent[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, T, F))
#tir$tsdstartup[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[3]][temp[[x]][[7]] & temp[[x]][[6]]], -1))
#tir$tsdstartdown[which(sapply(tir$tirseq, length)>1)] = sapply(1:length(temp), function(x) ifelse(sum(temp[[x]][[7]] & temp[[x]][[6]])==1, temp[[x]][[4]][temp[[x]][[7]] & temp[[x]][[6]]]-nchar(temp[[x]][[1]][1]), -1))

							   
###########
### OKAY! Now we have real TIRs, and can do stuff with them!!!
###    use TIR's RC to find in reverse sequence??
###########

### ADDING IN TSD TO FIND THIS, AS SHORTENED TIR CAN BE IN MULTIPLE PLACES!!!!
## get position of TIR in forward orientation of upstream extract.
#tir$tirstartup=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqSingle[x], tir$upstreamExtra[x])))
tir$tirstartup=unlist(mclapply(1:nrow(tir), function(x) as.numeric(regexpr(paste0(tir$tsdadjacentup, tir$tirseqSingle[x]), tir$upstreamExtra[x])), mc.cores=ncores))
		      
		      
## get position of TIR in forward orientation of downstream extract.
# to do RC: as.character(reverseComplement(DNAString(tirF)))
## and weird char introduced?: sapply(tirseq, function(tirF) tryCatch({as.character(reverseComplement(DNAString(tirF)))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}))
#tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x]))-1)
#tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tryCatch({as.character(reverseComplement(DNAString(tir$tirseqRCSingle[x])))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}), tir$downstreamExtra[x]))-1)
tir$tirstartdown=unlist(mclapply(1:nrow(tir), function(x) as.numeric(regexpr(paste0(tir$tirseqRCSingle[x], tir$tsdadjacentdown), tir$downstreamExtra[x]))-1, mc.cores=ncores))
			
			
			

## having an issue with weird characters, removing those that don't have a tir present here.
## removes about 150 copies where no match is found. 
#tirorig=tir
#tir=tir[tir$tirstartup != -1,]
#tir=tir[tir$tirstartdown != -2,]


############
## Now, adjust positions to put them back on the same scale as the genome!!!
############
##adjust positions
tir$start.adj=unlist(mclapply(1:nrow(tir), function(x) (tir$start[x]-offset) + tir$tirstartup[x] - 1 , mc.cores=ncores))## parentheses puts on same scale as upstreamExtra
tir$end.adj= unlist(mclapply(1:nrow(tir), function(x) (tir$end[x]-offset) + tir$tirstartdown[x] + nchar(tir$tirseqSingle[x]) - 1, mc.cores=ncores))

## floor is what GRanges does to decimal values, so replicate this here (e.g. 0.5 becomes 0, 3.5 becomes 3)
#  tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)
tir$start.adj[tir$origlen <= 2*offset ]= unlist(mclapply(which(tir$origlen <= 2*offset), function(x) (tir$start[x] - floor(2*offset - tir$origlen[x]/2) + tir$tirstartup[x] - 1 ), mc.cores=ncores))
tir$end.adj[tir$origlen <= 2*offset ]=   unlist(mclapply(which(tir$origlen <= 2*offset), function(x) (tir$end[x] - floor(tir$origlen[x]/2) + tir$tirstartdown[x] + nchar(tir$tirseqSingle[x]) - 1), mc.cores=ncores))
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
						
#############################################################################################################################						
#############################################################################################################################					
#### TIR mismatches  ########################################################################################################
#############################################################################################################################
#############################################################################################################################
tirm=tir[!tir$tsdadjacentequal & !tir$tirsmatch,]

## deal with pesky multiple TIRs by removing them here. 
tirm$tirseqSingle=unlist(lapply(tirm$tirseq, function(l) l[[1]]))
tirm$tirseqRCSingle=unlist(lapply(tirm$tirseqRC, function(l) l[[length(l)]]))  ## not good - this can be either first or second! AAAHHHHHHHHHH
tirm=tirm[!grepl('N', tirm$tirseqSingle),]

## get position of TIR in forward orientation of upstream extract.
tirm$tirstartup=unlist(mclapply(1:nrow(tirm), function(x) as.numeric(regexpr(tirm$tirseqSingle[x], tirm$upstreamExtra[x])), mc.cores=ncores))

## get position of TIR in forward orientation of downstream extract.
tirm$tirstartdown=unlist(mclapply(1:nrow(tirm), function(x) as.numeric(regexpr(tirm$tirseqRCSingle[x], tirm$downstreamExtra[x])), mc.cores=ncores)) ## NO DON'T!  -1 ## be consistent!
				 
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
	tsdadjacentup=unlist(mclapply(1:nrow(tirm), function(x) substr(tirm$upstreamExtra[x], tirm$tirstartup[x] - tirm$tsdlen[x]-tirposoffset , tirm$tirstartup[x]-1 -tirposoffset ), mc.cores=ncores))
	tsdadjacentdown=unlist(mclapply(1:nrow(tirm), function(x) substr(tirm$downstreamExtra[x], tirm$tirstartdown[x]  + nchar(tirm$tirseqSingle[x])+tirposoffset , tirm$tirstartdown[x]  + nchar(tirm$tirseqSingle[x]) + tirm$tsdlen[x] -1 +tirposoffset ), mc.cores=ncores))
	tsdadjacentequal=tsdadjacentup == tsdadjacentdown & tirm$tsdadjacentequal !=''
	print(sum(tsdadjacentequal))
	tirm$closestTSDseq[tsdadjacentequal & is.na(tirm$closestTSDseq)]=tsdadjacentup[tsdadjacentequal & is.na(tirm$closestTSDseq)]
	tirm$closestTSDoffset[tsdadjacentequal & is.na(tirm$closestTSDoffset)]=tirposoffset
	}
			 
### next, have to put in these offsets into the actual positions! and include the TSD sequence
### then, check the TIRs they suggest, and see how many mismatches			       

			       
tirm$tirstartup.adj=tirm$tirstartup
tirm$tirstartup.adj[!is.na(tirm$closestTSDoffset)]=(tirm$tirstartup-tirm$closestTSDoffset)[!is.na(tirm$closestTSDoffset)]
#### Don't do this - dumb! the TIR length changes in the positive direction, so adjusting downstream doesn't matter.
#tirm$tirstartdown.adj=tirm$tirstartdown
#tirm$tirstartdown.adj[!is.na(tirm$closestTSDoffset)]=(tirm$tirstartdown-tirm$closestTSDoffset)[!is.na(tirm$closestTSDoffset)]
### if less than 1, remove
tirm$tirstartup.adj[tirm$tirstartup.adj<1]=NA
tirm$tirstartup.adj[tirm$tirstartup.adj>2*offset]=NA

			       
#Check candidate TIRs
tirm$adjustedTIRup=unlist(mclapply(1:nrow(tirm), function(x) substr(tirm$upstreamExtra[x], tirm$tirstartup.adj[x], tirm$tirstartup[x]+nchar(tirm$tirseqSingle[x])-1), mc.cores=ncores))  ## need the minus one to exclude the first base of the TIR
						  
tirm$adjustedTIRdown=unlist(mclapply(1:nrow(tirm), function(x) substr(tirm$downstreamExtra[x], tirm$tirstartdown[x], tirm$tirstartdown[x]+nchar(tirm$adjustedTIRup[x])-1), mc.cores=ncores))
#tirm$adjustedTIRdown=unlist(mclapply(1:nrow(tirm), function(x) substr(tirm$downstreamExtra[x], tirm$tirstartdown[x]+1, tirm$tirstartdown[x]+nchar(tirm$adjustedTIRup[x])), mc.cores=ncores)) ## this minus one is already subracted in the regex! see diffs between tir adjustedtirdown before and after regex - i subtract AFTER getting TIR sequence.

tirm$adjustedTIRdownRC=NA
tirm$adjustedTIRdownRC[which(!is.na(tirm$adjustedTIRdown))]=unlist(mclapply(which(!is.na(tirm$adjustedTIRdown)), function(x) as.character(reverseComplement(DNAString(tirm$adjustedTIRdown[x]))), mc.cores=ncores))

tirm$seqdist=unlist(mclapply(1:nrow(tirm), function(x) stringdist(tirm$adjustedTIRup[x], tirm$adjustedTIRdownRC[x], method='h'), mc.cores=ncores))

		    
## either filter by those with TIRs <0.2 distant from each other		    
tirm[!is.na(tirm$closestTSDoffset) & tirm$seqdist<(nchar(tirm$adjustedTIRup)-nchar(tirm$tirseqSingle))*0.2,]
## or those with new extension of TIR <0.2 distant from each other
tirm[!is.na(tirm$closestTSDoffset) & tirm$seqdist<(nchar(tirm$adjustedTIRup))*0.2,]

###############################################
#### REDEFINE JUST AS FOR PERFECT TIR #########
###############################################
tirm$closestTSDseq[tirm$closestTSDseq=='']=NA	## don't want an empty TSD seq to match!!!		     
			     
tirm$tirstartup.regex=unlist(mclapply(1:nrow(tirm), function(x) as.numeric(regexpr(paste0(tirm$closestTSDseq[x], tirm$adjustedTIRup[x]), tirm$upstreamExtra[x]))+nchar(tirm$closestTSDseq[x]), mc.cores=ncores))
		      
		      
## get position of TIR in forward orientation of downstream extract.
# to do RC: as.character(reverseComplement(DNAString(tirF)))
## and weird char introduced?: sapply(tirseq, function(tirF) tryCatch({as.character(reverseComplement(DNAString(tirF)))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}))
#tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x]))-1)
#tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tryCatch({as.character(reverseComplement(DNAString(tir$tirseqRCSingle[x])))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}), tir$downstreamExtra[x]))-1)
tirm$tirstartdown.regex=unlist(mclapply(1:nrow(tirm), function(x) as.numeric(regexpr(paste0(tirm$adjustedTIRdown[x], tirm$closestTSDseq[x]), tirm$downstreamExtra[x]))-1, mc.cores=ncores))
tirm$tirstartup.regex[tirm$tirstartdown.regex==0]=NA ## i use startup.regex as a filter below
tirm$tirstartdown.regex[tirm$tirstartdown.regex==0]=NA ## you can't subset in R with a 0!     
			     
############
## Now, adjust positions to put them back on the same scale as the genome!!!
############
## also have to redo as above!				      
tirm$start.adj=unlist(mclapply(1:nrow(tirm), function(x) (tirm$start[x]-offset) + tirm$tirstartup.regex[x] - 1 , mc.cores=ncores))## parentheses puts on same scale as upstreamExtra
tirm$end.adj= unlist(mclapply(1:nrow(tirm), function(x) (tirm$end[x]-offset) + tirm$tirstartdown.regex[x] + nchar(tirm$adjustedTIRup[x]) - 1, mc.cores=ncores))

## floor is what GRanges does to decimal values, so replicate this here (e.g. 0.5 becomes 0, 3.5 becomes 3)
#  tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)
tirm$start.adj[tirm$origlen <= 2*offset ]= unlist(mclapply(which(tirm$origlen <= 2*offset), function(x) (tirm$start[x] - floor(2*offset - tirm$origlen[x]/2) + tirm$tirstartup.regex[x] - 1 ), mc.cores=ncores))
tirm$end.adj[tirm$origlen <= 2*offset ]=   unlist(mclapply(which(tirm$origlen <= 2*offset), function(x) (tirm$end[x] - floor(tirm$origlen[x]/2) + tirm$tirstartdown.regex[x] + nchar(tirm$adjustedTIRup[x]) - 1), mc.cores=ncores))

## fix short odds
is.odd <- function(x) x %% 2 != 0				   
tirm$start.adj[tirm$origlen <= 2*offset & is.odd(tirm$origlen)]= unlist(mclapply(which(tirm$origlen <= 2*offset & is.odd(tirm$origlen)), function(x) (tirm$start[x] - floor(2*offset - tirm$origlen[x]/2) + tirm$tirstartup.regex[x] - 2 ), mc.cores=ncores))
tirm$end.adj[tirm$origlen <= 2*offset & is.odd(tirm$origlen)]=   unlist(mclapply(which(tirm$origlen <= 2*offset & is.odd(tirm$origlen)), function(x) (tirm$end[x] - floor(tirm$origlen[x]/2) + tirm$tirstartdown.regex[x] + nchar(tirm$adjustedTIRup[x]) - 2), mc.cores=ncores))
							   
							   
##adjust positions
#tirm$start.adj=unlist(mclapply(1:nrow(tirm), function(x) (tirm$start[x]-offset) + tirm$tirstartup.adj[x] - 1 , mc.cores=ncores))## parentheses puts on same scale as upstreamExtra
######### ADDED A MINUS ONE HERE - A BIT WORRIED WHY THIS MATTERS FOR THIS ADJ BUT NOT THE FIRST (maybe because I used the adjusted one?)
#tirm$end.adj= unlist(mclapply(1:nrow(tirm), function(x) (tirm$end[x]-offset) + tirm$tirstartdown.adj[x] + nchar(tirm$adjustedTIRup[x]) - 1, mc.cores=ncores))
###### no longer needed with new search above
##is.odd <- function(x) x %% 2 != 0
### floor is what GRanges does to decimal values, so replicate this here (e.g. 0.5 becomes 0, 3.5 becomes 3)
##  tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)
########### EVENS WORK FINE AS IS!!!!
#tirm$start.adj[tirm$origlen <= 2*offset & !is.odd(tirm$origlen)]= unlist(mclapply(which(tirm$origlen <= 2*offset& !is.odd(tirm$origlen)), function(x) (tirm$start[x] - floor(2*offset - tirm$origlen[x]/2) + tirm$tirstartup.adj[x] - 1 ), mc.cores=ncores))
########## ADDED A MINUS ONE HERE - A BIT WORRIED WHY THIS MATTERS FOR THIS ADJ BUT NOT THE FIRST (maybe because I used the adjusted one?)
#tirm$end.adj[tirm$origlen <= 2*offset & !is.odd(tirm$origlen)]=   unlist(mclapply(which(tirm$origlen <= 2*offset& !is.odd(tirm$origlen)), function(x) (tirm$end[x] - floor(tirm$origlen[x]/2) + tirm$tirstartdown[x] + nchar(tirm$adjustedTIRup[x]) - 1), mc.cores=ncores))
########## ODDS NEED A MODIFICATION
#tirm$start.adj[tirm$origlen <= 2*offset & is.odd(tirm$origlen)]= unlist(mclapply(which(tirm$origlen <= 2*offset & is.odd(tirm$origlen)), function(x) (tirm$start[x] - floor(2*offset - tirm$origlen[x]/2) + tirm$tirstartup.adj[x] - 2 ), mc.cores=ncores))
########## ADDED A MINUS ONE HERE - A BIT WORRIED WHY THIS MATTERS FOR THIS ADJ BUT NOT THE FIRST (maybe because I used the adjusted one?)
#tirm$end.adj[tirm$origlen <= 2*offset & is.odd(tirm$origlen)]=   unlist(mclapply(which(tirm$origlen <= 2*offset & is.odd(tirm$origlen)), function(x) (tirm$end[x] - floor(tirm$origlen[x]/2) + tirm$tirstartdown[x] + nchar(tirm$adjustedTIRup[x]) - 3), mc.cores=ncores))

						  
################
## check that i have what i think i have
################						
#tir$tirupinseq=	as.character(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$start.adj, end=tir$start.adj+nchar(tir$tirseqSingle)-1))))				
#tir$tirdowninseqRC=as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$end.adj-nchar(tir$tirseqSingle)+1, end=tir$end.adj)))))
tirm$tiradjustedupinseq=NA
tirm$tiradjusteddowninseqRC=NA
tirm$tiradjustedupinseq[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex) & tirm$tirstartup.regex!=-1]=as.character(getSeq(seqs, GRanges(tirm$chrnew[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex)& tirm$tirstartup.regex!=-1], IRanges(start=tirm$start.adj[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex)& tirm$tirstartup.regex!=-1], end=tirm$start.adj[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex)& tirm$tirstartup.regex!=-1]+nchar(tirm$adjustedTIRup[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex)& tirm$tirstartup.regex!=-1])-1))))				
tirm$tiradjusteddowninseqRC[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex) & tirm$tirstartup.regex!=-1]=as.character(reverseComplement(getSeq(seqs, GRanges(tirm$chrnew[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex)& tirm$tirstartup.regex!=-1], IRanges(start=tirm$end.adj[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex)& tirm$tirstartup.regex!=-1]-nchar(tirm$adjustedTIRup[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex)& tirm$tirstartup.regex!=-1])+1, end=tirm$end.adj[!is.na(tirm$closestTSDseq) & !is.na(tirm$tirstartup.regex)& tirm$tirstartup.regex!=-1])))))
																		      
tirm$tirsadjustedmatch=tirm$adjustedTIRup==tirm$tiradjustedupinseq & tirm$adjustedTIRdownRC==tirm$tiradjusteddowninseqRC & tirm$tiradjustedupinseq!=''	

table(tirm$tirsadjustedmatch)						  
						  
dm=data.frame(tirm$chrnew, 'TARGeT', 'terminal_inverted_repeat_element', tirm$start.adj, tirm$end.adj-1, '.', tirm$strand, '.', paste0('ID=', tirm$mtec, '_', tirm$closestTSDseq, '_', tirm$adjustedTIRup, '_mismatch=', tirm$seqdist, '_', tirm$adjustedTIRdown))
### this concerns me!!!
dm[,4:5]=dm[,4:5]-1
#dm=dm[!is.na(dm[,1]),]
dm=dm[!is.na(tirm$closestTSDoffset) & tirm$seqdist<(nchar(tirm$adjustedTIRup)*0.2),]
dm=dm[complete.cases(dm),]
#write.table(d[!is.na(tir$whichrule) & d[,4]<d[,5],], file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
#write.table(d[d[,4]<d[,5],], file=paste0(GENOMENAME, '_unfiltered_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
write.table(dm, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.mismatch.gff3'), col.names=F, row.names=F, sep='\t', quote=F)

## for maizegdb with Chr
ddm=dm				
levels(ddm$tirm.chrnew)[1:10]=paste0('Chr', levels(ddm$tirm.chrnew)[1:10])	
write.table(ddm, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.mismatch.Chr.gff3'), col.names=F, row.names=F, sep='\t', quote=F)

### need to redo with the multiple possible TIRs			       

						  
						  
						  
#######################################################################################
### actually deal with pesky multiple TIRs by searching each for a TSD here.     ######
#######################################################################################

						  
checkTIRcandidateforTSD=function(tirseq, upstreamExtra, downstreamExtra, tsdlen, offset=0){
	tirseqRC=tryCatch({as.character(reverseComplement(DNAString(tirseq)))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')})
	upposn=as.numeric(regexpr(tirseq, upstreamExtra))-offset
	downposn=as.numeric(regexpr(tirseqRC, downstreamExtra)) + nchar(tirseq)+offset
	uptsd=substr(upstreamExtra, upposn-tsdlen, upposn-1)
	downtsd=substr(downstreamExtra, downposn, downposn + tsdlen -1)
	tsdsequal=uptsd==downtsd & uptsd != ''
#	tsdtirjunctionpresent=grepl(paste0(uptsd, tirseq), upstreamExtra, useBytes=T) & 
#			      grepl(paste0(tirseqRC, downtsd), downstreamExtra, useBytes=T)
	if(tsdsequal){
		return(list(tirseq, tirseqRC, upposn, downposn, uptsd, offset))
		}else{
		return(list('', '', -1, -1, '', 0))
		}
	}

checkTIRcandidateforOffsetTSD=function(tirseqSingle, upstreamExtra, downstreamExtra, offset=0, tsdlen){
	tirseqRCSingle=tryCatch({as.character(reverseComplement(DNAString(tirseqSingle)))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')})
	tirstartup=as.numeric(regexpr(tirseqSingle, upstreamExtra))
	tirstartdown=as.numeric(regexpr(tirseqRCSingle, downstreamExtra))
	upposn=tirstartup-offset
	downposn=tirstartdown + nchar(tirseqSingle) + offset
	uptsd=substr(upstreamExtra, upposn-tsdlen, upposn-1)
	downtsd=substr(downstreamExtra, downposn, downposn + tsdlen -1)
	tsdsequal= uptsd==downtsd & uptsd!=''
	if(tsdsequal){
		return(list(closestTSDoffset=offset, closestTSDseq=uptsd, tirstartup=upposn, tirstartdown=tirstartdown))
		}else{
		return(list(closestTSDoffset=NA, closestTSDseq=NA, tirstartup=NA, tirstartdown=NA))
		}
	}

tirm$closestTSDseq[which(sapply(tirm$tirseq, length)>1)]=NA    ## this is the sequence of the matching TSD
tirm$closestTSDoffset[which(sapply(tirm$tirseq, length)>1)]=NA ## this is the offset from the perfect TIR edge

						  
						  
## this works, but is running out of memory?? So loop through and do in chunks.
for(i in split(which(sapply(tirm$tirseq, length)>1), ceiling(seq_along(which(sapply(tirm$tirseq, length)>1))/1000))){
tempm=mclapply(i, function(x) {
	allcombos=expand.grid(as.character(tirm$tirseq[[x]]), 0:20)
	colnames(allcombos)=c('tirseq', 'offset')
	allcombos$tirseq=as.character(allcombos$tirseq)
	allcombos$closestTSDseq=NA
	tircands=data.frame(t(data.frame(sapply(1:nrow(allcombos), function(tirseqcand) checkTIRcandidateforOffsetTSD(tirseqSingle=allcombos$tirseq[tirseqcand], tirm$upstreamExtra[x], tirm$downstreamExtra[x], offset=allcombos$offset[tirseqcand], tsdlen=tirm$tsdlen[x])))))
	if(sum(!is.na(tircands$closestTSDseq))==1){
		tempstore=cbind(tircands[!is.na(tircands$closestTSDseq),], tirseq=allcombos$tirseq[!is.na(tircands$closestTSDseq)])
	}else{tempstore=tircands[1,]
	      tempstore[1,]=c(NA,NA,NA,NA,NA)
	     }
	}, mc.cores=ncores
	)
## this is a data frame of closestTSDoffset, closestTSDseq, tirstartup, and tirstartdown.
##  there's one entry for each tirm that has more than one potential TIR, but many are NA
tirmworking=as.data.frame(do.call(rbind, tempm))
tirmworking$tirseq=as.character(tirmworking$tirseq)
					 
		
tirm$closestTSDseq[i]=tirmworking[,'closestTSDseq']
tirm$closestTSDoffset[i]=tirmworking[,'closestTSDoffset']
tirm$tirstartup[i]=tirmworking$tirstartup
tirm$tirstartdown[i]=tirmworking$tirstartdown
tirm$tirseqSingle[i]=tirmworking$tirseq
					 
}					 

					 
tirm$closestTSDseq=unlist(tirm$closestTSDseq)
tirm$closestTSDoffset=unlist(tirm$closestTSDoffset)					 
tirm$tirstartup=unlist(tirm$tirstartup)
tirm$tirstartdown=unlist(tirm$tirstartdown)
tirm$tirseqSingle=unlist(tirm$tirseqSingle)		 
					 
					 
					 
tirm$tirstartup.adj[which(sapply(tirm$tirseq, length)>1)]=tirm$tirstartup[which(sapply(tirm$tirseq, length)>1)]
## the actual adjustment was already done in the loop!!!!
#tirm$tirstartup.adj[!is.na(tirm$closestTSDoffset)]=(tirm$tirstartup-tirm$closestTSDoffset)[!is.na(tirm$closestTSDoffset)]
tirm$tirstartup.adj[tirm$tirstartup.adj<1]=NA
tirm$tirstartup.adj[tirm$tirstartup.adj>2*offset]=NA

## essentially redoing this now that we have better candidates
#tirm$tirstartup.adj=(tirm$tirstartup-as.numeric(tirm$closestTSDoffset))
#tirm$tirstartup.adj[tirm$tirstartup.adj<1]=NA
### HERE, I can redo everything! And not do it above. BUT because I sutracted from tirstartup.adj above, don't redo those 1 tir copies here?
#Check candidate TIRs
tirm$adjustedTIRup[which(sapply(tirm$tirseq, length)>1)]=unlist(mclapply(which(sapply(tirm$tirseq, length)>1), function(x) substr(tirm$upstreamExtra[x], tirm$tirstartup.adj[x], tirm$tirstartup[x]+nchar(tirm$tirseqSingle[x]) + tirm$closestTSDoffset[x] -1), mc.cores=ncores))  ## need the minus one to exclude the first base of the TIR
tirm$adjustedTIRdown[which(sapply(tirm$tirseq, length)>1)]=unlist(mclapply(which(sapply(tirm$tirseq, length)>1), function(x) substr(tirm$downstreamExtra[x], tirm$tirstartdown[x], tirm$tirstartdown[x]+nchar(tirm$tirseqSingle[x]) + tirm$closestTSDoffset[x] -1), mc.cores=ncores))
#tirm$adjustedTIRdown=unlist(mclapply(1:nrow(tirm), function(x) substr(tirm$downstreamExtra[x], tirm$tirstartdown[x]+1, tirm$tirstartdown[x]+nchar(tirm$adjustedTIRup[x])), mc.cores=ncores)) ## this minus 1 already subtracted in the regex!!
tirm$adjustedTIRdownRC=NA
tirm$adjustedTIRdownRC[which(!is.na(tirm$adjustedTIRdown))]=unlist(mclapply(which(!is.na(tirm$adjustedTIRdown)), function(x) as.character(reverseComplement(DNAString(tirm$adjustedTIRdown[x]))), mc.cores=ncores))
tirm$seqdist[which(sapply(tirm$tirseq, length)>1)]=unlist(mclapply(which(sapply(tirm$tirseq, length)>1), function(x) stringdist(tirm$adjustedTIRup[x], tirm$adjustedTIRdownRC[x], method='h'), mc.cores=ncores))
		
#### SO NOW, I THINK I HAVE IT ALL!!!!!!!!!!!!!
		    
		    
#####################
		    
## either filter by those with TIRs <0.2 distant from each other		    
#tirm[!is.na(tirm$closestTSDoffset) & tirm$seqdist<(nchar(tirm$adjustedTIRup)-nchar(tirm$tirseqSingle))*0.2,]
## or those with new extension of TIR <0.2 distant from each other
head(tirm[!is.na(tirm$closestTSDoffset) & tirm$seqdist<(nchar(tirm$adjustedTIRup))*0.2,])
		    
		    
############
## Now, adjust positions to put them back on the same scale as the genome!!!
############
##adjust positions
#tirm$start.adj=unlist(mclapply(1:nrow(tirm), function(x) (tirm$start[x]-offset) + tirm$tirstartup.adj[x] - 1 , mc.cores=ncores))## parentheses puts on same scale as upstreamExtra
######### ADDED A MINUS ONE HERE - A BIT WORRIED WHY THIS MATTERS FOR THIS ADJ BUT NOT THE FIRST (maybe because I used the adjusted one?)
#tirm$end.adj= unlist(mclapply(1:nrow(tirm), function(x) (tirm$end[x]-offset) + tirm$tirstartdown[x] + nchar(tirm$adjustedTIRup[x]) - 1, mc.cores=ncores))-1 ## worried there is another -1 here - is this because I'm defining differently??

## floor is what GRanges does to decimal values, so replicate this here (e.g. 0.5 becomes 0, 3.5 becomes 3)
#  tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)
#tirm$start.adj[tirm$origlen <= 2*offset ]= unlist(mclapply(which(tirm$origlen <= 2*offset), function(x) (tirm$start[x] - floor(2*offset - tirm$origlen[x]/2) + tirm$tirstartup.adj[x] - 1 ), mc.cores=ncores))
######### ADDED A MINUS ONE HERE - A BIT WORRIED WHY THIS MATTERS FOR THIS ADJ BUT NOT THE FIRST (maybe because I used the adjusted one?)
#tirm$end.adj[tirm$origlen <= 2*offset ]=   unlist(mclapply(which(tirm$origlen <= 2*offset), function(x) (tirm$end[x] - floor(tirm$origlen[x]/2) + tirm$tirstartdown[x] + nchar(tirm$adjustedTIRup[x]) - 1), mc.cores=ncores))-1

###############################################
#### REDEFINE JUST AS FOR PERFECT TIR #########
###############################################
tirm$closestTSDseq[tirm$closestTSDseq=='']=NA	## don't want an empty TSD seq to match!!!		     
			     
tirm$tirstartup.regex[which(sapply(tirm$tirseq, length)>1)]=unlist(mclapply(which(sapply(tirm$tirseq, length)>1), function(x) as.numeric(regexpr(paste0(tirm$closestTSDseq[x], tirm$adjustedTIRup[x]), tirm$upstreamExtra[x]))+nchar(tirm$closestTSDseq[x]), mc.cores=ncores))
		      
		      
## get position of TIR in forward orientation of downstream extract.
# to do RC: as.character(reverseComplement(DNAString(tirF)))
## and weird char introduced?: sapply(tirseq, function(tirF) tryCatch({as.character(reverseComplement(DNAString(tirF)))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}))
#tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tir$tirseqRCSingle[x], tir$downstreamExtra[x]))-1)
#tir$tirstartdown=sapply(1:nrow(tir), function(x) as.numeric(regexpr(tryCatch({as.character(reverseComplement(DNAString(tir$tirseqRCSingle[x])))}, error=function(e){print(paste('line not working', x, 'error is', e)); return('NNNNN')}), tir$downstreamExtra[x]))-1)
tirm$tirstartdown.regex[which(sapply(tirm$tirseq, length)>1)]=unlist(mclapply(which(sapply(tirm$tirseq, length)>1), function(x) as.numeric(regexpr(paste0(tirm$adjustedTIRdown[x], tirm$closestTSDseq[x]), tirm$downstreamExtra[x]))-1, mc.cores=ncores))
			     
			     
############
## Now, adjust positions to put them back on the same scale as the genome!!!
############
## also have to redo as above!				      
tirm$start.adj[which(sapply(tirm$tirseq, length)>1)]=unlist(mclapply(which(sapply(tirm$tirseq, length)>1), function(x) (tirm$start[x]-offset) + tirm$tirstartup.regex[x] - 1 , mc.cores=ncores))## parentheses puts on same scale as upstreamExtra
tirm$end.adj[which(sapply(tirm$tirseq, length)>1)]= unlist(mclapply(which(sapply(tirm$tirseq, length)>1), function(x) (tirm$end[x]-offset) + tirm$tirstartdown.regex[x] + nchar(tirm$adjustedTIRup[x]) - 1, mc.cores=ncores))

## floor is what GRanges does to decimal values, so replicate this here (e.g. 0.5 becomes 0, 3.5 becomes 3)
#  tir$start[tir$origlen <= 2*offset ]- ( 2*offset-tir$origlen[tir$origlen <= 2*offset ]/2)
tirm$start.adj[tirm$origlen <= 2*offset & sapply(tirm$tirseq, length)>1 ]= unlist(mclapply(which(sapply(tirm$tirseq, length)>1 & tirm$origlen <= 2*offset), function(x) (tirm$start[x] - floor(2*offset - tirm$origlen[x]/2) + tirm$tirstartup.regex[x] - 1 ), mc.cores=ncores))
tirm$end.adj[tirm$origlen <= 2*offset & sapply(tirm$tirseq, length)>1]=   unlist(mclapply(which(sapply(tirm$tirseq, length)>1 & tirm$origlen <= 2*offset), function(x) (tirm$end[x] - floor(tirm$origlen[x]/2) + tirm$tirstartdown.regex[x] + nchar(tirm$adjustedTIRup[x]) - 1), mc.cores=ncores))
				      				      
					  
################
## check that i have what i think i have
################						
#tir$tirupinseq=	as.character(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$start.adj, end=tir$start.adj+nchar(tir$tirseqSingle)-1))))				
#tir$tirdowninseqRC=as.character(reverseComplement(getSeq(seqs, GRanges(tir$chrnew, IRanges(start=tir$end.adj-nchar(tir$tirseqSingle)+1, end=tir$end.adj)))))
tirm$tiradjustedupinseq[which(sapply(tirm$tirseq, length)>1)]=NA
tirm$tiradjusteddowninseqRC[which(sapply(tirm$tirseq, length)>1)]=NA
tirm$tiradjustedupinseq[!is.na(tirm$closestTSDoffset) & sapply(tirm$tirseq, length)>1]=as.character(getSeq(seqs, GRanges(tirm$chrnew[!is.na(tirm$closestTSDseq) & sapply(tirm$tirseq, length)>1], IRanges(start=tirm$start.adj[!is.na(tirm$closestTSDseq) & sapply(tirm$tirseq, length)>1], end=tirm$start.adj[!is.na(tirm$closestTSDseq) & sapply(tirm$tirseq, length)>1]+nchar(tirm$adjustedTIRup[!is.na(tirm$closestTSDseq) & sapply(tirm$tirseq, length)>1])-1))))				
tirm$tiradjusteddowninseqRC[!is.na(tirm$closestTSDseq) & sapply(tirm$tirseq, length)>1]=as.character(reverseComplement(getSeq(seqs, GRanges(tirm$chrnew[!is.na(tirm$closestTSDseq) & sapply(tirm$tirseq, length)>1], IRanges(start=tirm$end.adj[!is.na(tirm$closestTSDseq) & sapply(tirm$tirseq, length)>1]-nchar(tirm$adjustedTIRup[!is.na(tirm$closestTSDseq) & sapply(tirm$tirseq, length)>1])+1, end=tirm$end.adj[!is.na(tirm$closestTSDseq) & sapply(tirm$tirseq, length)>1])))))
																		      
tirm$tirsadjustedmatch=tirm$adjustedTIRup==tirm$tiradjustedupinseq & tirm$adjustedTIRdownRC==tirm$tiradjusteddowninseqRC & tirm$tiradjustedupinseq!=''	
						  
################
## check that i have what i think i have
################						

#tirm$tiradjustedupinseq[!is.na(tirm$adjustedTIRup)]=as.character(getSeq(seqs, GRanges(tirm$chrnew[!is.na(tirm$adjustedTIRup)], IRanges(start=tirm$start.adj[!is.na(tirm$adjustedTIRup)], end=tirm$start.adj[!is.na(tirm$adjustedTIRup)]+nchar(tirm$adjustedTIRup[!is.na(tirm$adjustedTIRup)])-1))))				
#tirm$tiradjusteddowninseqRC[!is.na(tirm$adjustedTIRup)]=as.character(reverseComplement(getSeq(seqs, GRanges(tirm$chrnew[!is.na(tirm$adjustedTIRup)], IRanges(start=tirm$end.adj[!is.na(tirm$adjustedTIRup)]-nchar(tirm$adjustedTIRup[!is.na(tirm$adjustedTIRup)])+1, end=tirm$end.adj[!is.na(tirm$adjustedTIRup)])))))																		      
#tirm$tirsadjustedmatch=tirm$adjustedTIRup==tirm$tiradjustedupinseq & tirm$adjustedTIRdownRC==tirm$tiradjusteddowninseqRC & tirm$tiradjustedupinseq!=''	

tirm[!is.na(tirm$closestTSDoffset) & !tirm$tirsadjustedmatch,]

### this is bad, but I can't come up with a solution. 
#tirm$start.adj[!tirm$tirsadjustedmatch & !is.na(tirm$adjustedTIRup)]=tirm$start.adj[!tirm$tirsadjustedmatch & !is.na(tirm$adjustedTIRup)] +1
#tirm$end.adj[!tirm$tirsadjustedmatch & !is.na(tirm$adjustedTIRup)]=tirm$end.adj[!tirm$tirsadjustedmatch & !is.na(tirm$adjustedTIRup)] +1
	
## now recalculate
#tirm$tiradjustedupinseq[!is.na(tirm$adjustedTIRup)]=as.character(getSeq(seqs, GRanges(tirm$chrnew[!is.na(tirm$adjustedTIRup)], IRanges(start=tirm$start.adj[!is.na(tirm$adjustedTIRup)], end=tirm$start.adj[!is.na(tirm$adjustedTIRup)]+nchar(tirm$adjustedTIRup[!is.na(tirm$adjustedTIRup)])-1))))				
#tirm$tiradjusteddowninseqRC[!is.na(tirm$adjustedTIRup)]=as.character(reverseComplement(getSeq(seqs, GRanges(tirm$chrnew[!is.na(tirm$adjustedTIRup)], IRanges(start=tirm$end.adj[!is.na(tirm$adjustedTIRup)]-nchar(tirm$adjustedTIRup[!is.na(tirm$adjustedTIRup)])+1, end=tirm$end.adj[!is.na(tirm$adjustedTIRup)])))))																		      
#tirm$tirsadjustedmatch=tirm$adjustedTIRup==tirm$tiradjustedupinseq & tirm$adjustedTIRdownRC==tirm$tiradjusteddowninseqRC & tirm$tiradjustedupinseq!=''	
						  
						  
						  
dm=data.frame(tirm$chrnew, 'TARGeT', 'terminal_inverted_repeat_element', tirm$start.adj, tirm$end.adj-1, '.', tirm$strand, '.', paste0('ID=', tirm$mtec, '_', tirm$closestTSDseq, '_', tirm$adjustedTIRup, '_mismatch=', tirm$seqdist, '_', tirm$adjustedTIRdown))
### this concerns me!!!
dm[,4:5]=dm[,4:5]-1
dm=dm[!is.na(tirm$closestTSDoffset) & tirm$seqdist<(nchar(tirm$adjustedTIRup)*0.2),]
#write.table(d[!is.na(tir$whichrule) & d[,4]<d[,5],], file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
#write.table(d[d[,4]<d[,5],], file=paste0(GENOMENAME, '_unfiltered_tir_', Sys.Date(), '.gff3'), col.names=F, row.names=F, sep='\t', quote=F)
write.table(dm, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.mismatchAll.gff3'), col.names=F, row.names=F, sep='\t', quote=F)

## for maizegdb with Chr
ddm=dm				
levels(ddm$tirm.chrnew)[1:10]=paste0('Chr', levels(ddm$tirm.chrnew)[1:10])	
write.table(ddm, file=paste0(GENOMENAME, '_tir_', Sys.Date(), '.mismatchAll.Chr.gff3'), col.names=F, row.names=F, sep='\t', quote=F)

