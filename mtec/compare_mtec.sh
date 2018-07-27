
cat ~/te_reference_fasta/individual_fasta/DT* > MTEC_TIR.fa

SILIX=~/software/bin/silix
USEARCH=/home/mstitzer/software/vsearch/bin/vsearch


srun -p bigmemh --time=1-00:00 ${USEARCH} -allpairs_global MTEC_TIR.fa -blast6out MTEC_TIR.allvall.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads 1


srun -p bigmemh --time=1-00:00 ${SILIX} MTEC_TIR.fa MTEC_TIR.allvall.out -f MTEC -i 0.8 -r 0.8 --net > MTEC_TIR.8080.fnodes
