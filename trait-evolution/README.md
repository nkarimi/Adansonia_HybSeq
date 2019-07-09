# network calibration and trait evolution analyses

## network calibration

julia script in `calibration.jl`  
results:
- pairwise genetic distances in `genetreePairwiseDistance_*.csv`
- calibrated network in `*_calibrated.tre`, with branch lengths
  scaled to a height of 1 between the crown node of the `Adansonia`
  clade and the tips.

## trait evolution

- trait data in `Adansonia_floralTraits_accession.csv`
- julia script in `traitanalysis.jl`
