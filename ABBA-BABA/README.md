# ABBA-BABA analyses

- `calcD.R`: helper functions to calculate the number of ABBA and BABA patterns
  from each gene, bootstrap genes, calculate the D statistic and its bootstrap
  standard deviation.
- `reorder_fasta.py`: python script to read in a large fasta alignment,
  and extract a sub-aligment of a desired set of 4, in a specific order.
- `ABBA-BABA.Rmd`: scripts to perform the ABBA-BABA analyses to test
  specific questions about in the baobabs. Uses `reorder_fasta.py` and `calcD.R`
  repeatedly.
- `ABBA-BABA-Dstat.csv`: result file created by `ABBA-BABA.Rmd` with the D statistic,
  its standard deviation and associated z score, for many 4-taxon sets.
