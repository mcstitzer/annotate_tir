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

### path right now is: (need to clean up)

1. `tir_search_simple.R`
2. `add_families.R` not great documentation of how to get source
3. `tir_search_simple.mismatchTIR.R` 
4. `add_families.mismatch.R` uses output from #2
