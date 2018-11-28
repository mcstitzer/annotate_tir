# annotate_tir
find terminal inverted repeats (TIRs)


#### not worth removing terminal seq difference
Ac/Ds is terminated by imperfect TIRs, at the terminus.

cAGGGATGAAA
tAGGGATGAAA

[Kunze and Weil (2002), Mobile DNA II](http://www.asmscience.org/content/book/10.1128/9781555817954.chap24)



#### Some TIRs can be short

Bergamo has 5 bp TIRs

CAGGG

[Kunze and Weil (2002), Mobile DNA II](http://www.asmscience.org/content/book/10.1128/9781555817954.chap24)

Rlibstree can be difficult to install, best way:
`install.packages("devtools")`
`library(devtools)`
`install_github("omegahat/Rlibstree")`

### path right now is: (need to clean up)

1. `tir_search_simple.mismatchTIR.R` - NOTE that this has to be on UNIX (even using R) because I'm having locale issues I haven't been able to resolve on my Mac.
2. `add_families.R` uses output from #1, but needs to be adapted to read in only gff's not the crazy large text file I output

This generates a gff3 with all TIR copies, not two (identical TIR and mismatch TIR) as I was doing before.

### old path was (before Aug 29)
1. `tir_search_simple.R`
2. `add_families.R` not great documentation of how to get source
3. `tir_search_simple.mismatchTIR.R` - NOTE that this has to be on UNIX (even using R) because I'm having locale issues I haven't been able to resolve on my Mac.
4. `add_families.mismatch.R` uses output from #2


