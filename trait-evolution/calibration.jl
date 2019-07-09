# software versions: julia v1.1.1, see details at the end
using Statistics
using StatsBase
using DataFrames
using CSV
using PhyloNetworks
using RCall
using PhyloPlots

# naming: net1, net2 = networks from primary / haplotype data
#         o = with outgroups, i = ingroup only
#         f = full taxon names
#         tre1, tre2 = primary / haplotype gene trees

fulltaxnamedict = Dict(
  "Pcr070"=>"Pseudobombax crozatii", "Smi165"=>"Scleronema micrantha", # outgroups
  "Bce020"=>"Bombax ceiba",
  "Adi003"=>"A. digitata", "Adi002"=>"A. digitata", "Adi001"=>"A. digitata", # Africa
  "Age001"=>"A. gregorii", # Australia
  "Aga001"=>"A. grandidieri", "Aga002"=>"A. grandidieri", # brevi
  "Asu001"=>"A. suarezensis",
  "Aru001"=>"A. rubrostipa", "Aru127"=>"A. rubrostipa",   # longi
  "Ama018"=>"A. madagascariensis", "Ama006"=>"A. madagascariensis",
  "Aza037"=>"A. za",          "Aza135"=>"A. za",
  "Ape009"=>"A. perrieri",    "Ape001"=>"A. perrieri",
)
outgroup = Set(["Pcr070","Smi165","Bce020"])
ingroup  = setdiff(keys(fulltaxnamedict), outgroup)
taxa = ["Pcr070","Smi165","Bce020", # outgroups
        "Adi001","Adi002","Adi003", "Age001", # Africa, Australia
        "Aga001","Aga002","Asu001",  # brevi
        "Aru001","Aru127", "Ama006","Ama018", "Aza037","Aza135", "Ape001","Ape009"]
Set(taxa) == keys(fulltaxnamedict) || error("incorrect taxa or fulltaxnamedict")
ntips = length(taxa) # 18

#--------------------------------------------------
#     pairwise distances: helper functions
#--------------------------------------------------

"""
    getPairwiseDistances(genetrees)

Use `pairwiseTaxonDistanceMatrix(tree)` from `PhyloNetworks`, which outputs
a matrix in which the rows correspond to taxa in the order in which they
come in `tipLabels(tree)`

Tricky part: the taxa may not come in the same order in `tipLabels` across
all trees. first build a taxon list: to get a matrix of taxa
in the same order for each gene. Some genes don't even have all taxa.
"""
function getpairwisedistances(genetrees)
  # ntips & taxa defined outside: see top of file
  D = Array{Float64,2}[]; # empty matrix. will contain distances
  ngenes = zeros(Int, ntips, ntips) # number of genes that have data for each pair
  geneind = Int[];        # indices of genes with missing taxa
  istaxonmissing = Vector{Bool}(undef, ntips) # to be modified in place for each gene
  for i in 1:length(genetrees)
    #println("getting distances from gene $i...")
    g = genetrees[i]
    M = zeros(ntips,ntips) # initialized at 0.0: for missing pairs
    taxnames = tipLabels(g)
    tipind = Int[]
    for k in 1:ntips
      j = findfirst(isequal(taxa[k]), taxnames)
      istaxonmissing[k] = j==nothing # modified in place
      if j!=nothing
          push!(tipind, j) # tipind[k] = j if not missing data
      end
    end
    M[.!istaxonmissing, .!istaxonmissing] = pairwiseTaxonDistanceMatrix(g)[tipind,tipind]
    ngenes[.!istaxonmissing, .!istaxonmissing] .+= 1
    if any(istaxonmissing)
      # @warn "taxa not found in gene $i: $(taxa[findall(istaxonmissing)])"
      push!(geneind,i)
    end
    push!(D, M)
  end
  return D, ngenes, geneind
end

"""
averagepairwisedistances(D, ngenes)

Average pairwise distance matrices, weighted by number of genes
with data on the pair.
For each pair of taxa `i,j`, only `ngenes[i,j]` contributed data
to `D[i,j]`. The other genes contributed 0.
"""
function averagepairwisedistances(D, ngenes)
  return sum(D) ./ ngenes
end

"""
    normalizedistances_outgroup2ingroup!(D;
        taxa=taxa, ingroup=ingroup, outgroup=outgroup)

Rescale each input distance matrix `D[k]`, such that all have the same
median patristic distance between outgroup taxa and ingroup taxa.
Input: `D` should be a vector of pairwise distances matrices,
one per gene tree.

*median*: such one taxon or one small clade with an unusually
large (or low) substitution rate does have an undue influence on the
scaling factor.

Assumptions:
- all trees have at least 1 outgroup and 1 ingroup.
- row & column `i` in Ds correspond to `taxa[i]` (`taxa` is external)
- `D[k][i,j]` = 0 if gene `k` doesn't have both taxa `i` and `j`,
- `ingroup` and `outgroup` are sets (external variables as well)
"""
function normalizedistances_outgroup2ingroup!(D;
            taxa=taxa, ingroup=ingroup, outgroup=outgroup)

  ntax = length(taxa) # should also be the nrow & ncol of each D[k]!!
  inding = findall(in(ingroup),  taxa) # indices of ingroup taxa. 4:18
  indout = findall(in(outgroup), taxa) # 1:3
  medianingroup2outgroup = Float64[]
  for dm in D # dm = distance matrix
    absent = findall([all(dm[:,i] .== 0.0) for i in 1:ntax])
    push!(medianingroup2outgroup,
          median(dm[setdiff(inding, absent), setdiff(indout, absent)]) )
  end
  mi2o = mean(medianingroup2outgroup)
  for i in 1:length(D)
    D[i] .*= mi2o/medianingroup2outgroup[i]
  end
  return medianingroup2outgroup
end

"""
    normalizetrees_totaltreelength!(genetrees)

UNUSED! Rescale each input gene tree to have a total length that is
the original mean of all gene trees' total tree length.
Problem: smaller tree length for trees that are missing taxa.
"""
function normalizetrees_totaltreelength!(genetrees)
  tl = [sum(e.length for e in n.edge) for n in genetrees] # tree lengths
  describe(tl)
  # haphunt, boxplot: 0.07175---[0.138622-0.18593-0.248549]---0.9167
  #          mean=0.205243
  ntax = [length(tipLabels(n)) for n in genetrees]
  # R"plot($ntax, $tl)"; R"cor($ntax, $tl)"; # -0.001099542
  tlmean = mean(tl)
  for i in 1:length(genetrees)
    for e in genetrees[i].edge
      e.length *= tlmean/tl[i]
    end
  end
end

#--------------------------------------------------
#   primary: pairwise distances from 372 gene trees
#--------------------------------------------------

tre1 = readMultiTopology("372genes_raxml/primary_372.tre");
length(tre1) == 372 # trees
countmap([n.numNodes for n in tre1]) # in StatsBase
countmap([length(tipLabels(n)) for n in tre1])
# 16 taxa (30 nodes) => 1 tree; 17 (32) => 14; 18 (34) => 357
findall(n -> isempty(outgroup ∩ tipLabels(n)), tre1)
# none: all genes have at least 1 outgroup (and 1 ingroup) taxon

# calculate the pairwise distances from each tree
D1, ngenes1, geneind1 = getpairwisedistances(tre1);
# consistency checks
length(geneind1) # 15 = 372-357 genes are missing 1 or more taxa
minimum([minimum(ngenes1[:,i]) for i in 1:ntips])
# 360: all pairs are informed by a minimum of 360 trees that have both taxa
median([d[i,j] for d in D1 for i in 2:ntips for j in 1:(i-1)]) # 0.020074336240907806
Adi12 = [d[4,5] for d in D1];
Aga12 = [d[8,9] for d in D1];
Ama6Adi3 = [d[13,6] for d in D1];
median(Adi12) # 0.0037664: good, sisters
median(Aga12) # 0.0029353: good, sisters
median(Ama6Adi3) # 0.01546609:  good, distant
# normalize distances, to same median ingroup-outgroup dist across all genes
med_in2out1 = normalizedistances_outgroup2ingroup!(D1)
describe(med_in2out1) # 0.026617---[0.044601-0.055564-0.074690]---0.367189 mean:0.067741
median([d[4,5] for d in D1]) # Adi12: 0.00460401
median([d[8,9] for d in D1]) # Aga12: 0.00348860
median([d[13,6] for d in D1])# Ama6Adi3: 0.0194310
# average across genes
avD1 = averagepairwisedistances(D1, ngenes1)
avD1[[4,6],[4,5,9]]
avD1[8,9]  # 0.00973443
avD1[13,6] # 0.02425564
# save average pairwise distances to file
CSV.write("traitanalysis/genetreePairwiseDistance_primary.csv",
          DataFrame([avD1[:,j] for j in 1:ntips], [Symbol(t) for t in taxa]))

#--------------------------------------------------
# haplotype: pairwise distances from 344 gene trees
#--------------------------------------------------

tre2 = readMultiTopology("344genes_haphunt/haphunt_subset1_344.tre");
length(tre2) == 344 # trees
countmap([n.numNodes for n in tre2]) # in StatsBase
countmap([length(tipLabels(n)) for n in tre2])
# 13 taxa (24 nodes) => 1 tree; 14 (26) => 3; 15 (28) => 3
# 16 (30) => 22; 17 (32) => 73; 18 (34) => 242

# pattern of taxon missingness / taxon sharing among gene trees
sharedtaxa = tipLabels(tre2[1])
for n in tre2 intersect!(sharedtaxa, tipLabels(n)); end
sharedtaxa # "Aru001", "Ape001" only!!
findall(n -> !("Aga001" in tipLabels(n)), tre2) # 278
length(tipLabels(tre2[278])) # 17 taxa: all but 1!
[count(n -> !(tax in tipLabels(n)), tre2) for tax in taxa]
# [11, 5, 68, 6, 5, 5, 16, 1, 3, 1, 0, 4, 4, 2, 3, 7, 0, 2]
taxa[3] # "Bce020" missing from 68 gene trees
findall(n -> isempty(outgroup ∩ tipLabels(n)), tre2)
# none: all genes have at least 1 outgroup (and 1 ingroup) taxon

# calculate the pairwise distances from each tree
D2, ngenes2, geneind2 = getpairwisedistances(tre2);
# consistency checks
length(geneind2) # 102 = 344-242 genes are missing 1 or more taxa
minimum([minimum(ngenes2[:,i]) for i in 1:ntips]) # 269
# some pairs are together in as few as 269 genes
maximum([minimum(ngenes2[:,i]) for i in 1:ntips]) # 276
# one taxon is together with every other taxon in at most 276 genes
minimum([maximum(ngenes2[:,i]) for i in 1:ntips]) # 276
maximum([maximum(ngenes2[:,i]) for i in 1:ntips]) # 344
# some pairs are together in all 344 genes
maximum([d[i,i] for d in D2 for i in 1:ntips]) # 0.0
minimum([d[i,j] for d in D2 for i in 2:ntips for j in 1:i]) # 0.0
median([d[i,j] for d in D2 for i in 2:ntips for j in 1:(i-1)]) # 0.024737784161420247
Adi12 = [d[4,5] for d in D2];
Aga12 = [d[8,9] for d in D2];
Ama6Adi3 = [d[13,6] for d in D2];
median(Adi12) # 0.007630487881: good, sisters
median(Aga12) # 0.002770531561: good, sisters
median(Ama6Adi3) # 0.031402587:  good, distant
# normalize distances, to same median ingroup-outgroup dist across all genes
med_in2out = normalizedistances_outgroup2ingroup!(D2)
describe(med_in2out) # 0.020649---[0.043737-0.053208-0.066023]---0.188097 mean:0.058690
median([d[4,5] for d in D2]) # Adi12: 0.008494699
median([d[8,9] for d in D2]) # Aga12: 0.002953929
median([d[13,6] for d in D2])# Ama6Adi3: 0.032502604
# average across genes
avD2 = averagepairwisedistances(D2, ngenes2)
avD2[[4,6],[4,5,9]]
avD2[8,9]  # 0.0064783
avD2[13,6] # 0.0441279
# save average pairwise distances to file
CSV.write("traitanalysis/genetreePairwiseDistance_haphunt.csv",
          DataFrame([avD2[:,j] for j in 1:ntips], [Symbol(t) for t in taxa]))

#------------------------------------------------
#    primary data set (372 genes), with outgroups
#------------------------------------------------

net1o = readTopology("snaq_results/BestH1_372g.out")
rootonedge!(net1o, 7); # 2 outgroups here only:
# Scleronema micrantha (Smi165) Pseudobombax crozatii (Pcr070)
# Bombax ceiba (Bce020) missing.
# plot(net1o, :R, useEdgeLength=true, showEdgeNumber=true);
# plot(net1o, :R, useEdgeLength=true, showNodeNumber=true);
for n in [-7,-6,-5,-8,-10,-17,-11,-13] rotate!(net1o, n); end
plot(net1o, :R, useEdgeLength=true);

net1of = deepcopy(net1o);
for n in net1of.leaf n.name = fulltaxnamedict[n.name]; end
# plot(net1of, :R, useEdgeLength=true);

# nice plot, showing topology only (not branch lengths)
hi = findfirst([!e.isMajor for e in net1of.edge]) # hi for hybrid index: 32
R"pdf"("traitanalysis/primary372_withOutgroups.pdf", width=7, height=5);
R"par"(mar=[0,0,0,0]);
res = plot(net1of, :R, xlim=[1.0,17.0], tipOffset=0.1);
hx1, hx2, hy1, hy2 = (res[i][hi] for i in 9:12); # coordinates for minor hybrid edge
R"arrows"(hx1, hy1, hx2, hy2, col="deepskyblue", length=0.08, angle=20);
#R"text"(x=1, y=17.5, labels="no branch lengths", adj=[0,0], cex=0.8);
x0=14.7; x1=15.2
R"text"(res[14][hi,:x]+0.6, res[14][hi,:y], res[14][hi,:gam], col="deepskyblue", cex=0.75);
R"lines"(x=[x0,x0], y=[16,17], lwd=2); R"text"(x=x1, y=16.5,"outgroups", adj=0);
R"lines"(x=[x0,x0], y=[13,15], lwd=2); R"text"(x=x1, y=14,"Africa", adj=0);
R"lines"(x=[x0,x0], y=[2, 12], lwd=2); R"text"(x=x1, y=6, "Madagascar", adj=0);
R"lines"(x=[x0,x0], y=[0.5,1.2], lwd=2); R"text"(x=x1, y=1,"Australia", adj=0);
R"dev.off()";

# load pairwise distance, to prepare for network calibration
avD1 = disallowmissing(convert(Matrix, CSV.read("traitanalysis/genetreePairwiseDistance_primary.csv")))
# check that the column names are the same (and in the same order) as in taxa:
all(taxa .== String.(names(CSV.read("traitanalysis/genetreePairwiseDistance_primary.csv"))))
# modify branch lengths for good starting values
maximum(avD1) # maximum distance: 0.10849
# plot(net1o, :R, useEdgeLength=true, showEdgeNumber=true); # to see which edges to add to get max distance
sum(net1o.edge[i].length for i in [7,35,8,31,32,27,25]) # 13.498. 0.10849/13.498=0.00804
calnet1o = deepcopy(net1o)
for e in calnet1o.edge
  if e.length < 0.0 e.length  = 0.001; else e.length *= 0.008; end
end
# calibration: remove Bce020 from taxa and distance matrix, because not in network
# taxa1o = deepcopy(taxa); deleteat!(taxa1o, 3)
ind1o = deleteat!(collect(1:18),3)
taxa1o = taxa[ind1o]
avD1o = avD1[ind1o,ind1o]
# run the calibration 3 times to make sure we optimized branch lengths thoroughly:
calibrateFromPairwiseDistances!(calnet1o, avD1o, taxa1o, forceMinorLength0=false, NLoptMethod=:LD_MMA)
calibrateFromPairwiseDistances!(calnet1o, avD1o, taxa1o, forceMinorLength0=false, NLoptMethod=:LN_COBYLA) # score got worse
calibrateFromPairwiseDistances!(calnet1o, avD1o, taxa1o, forceMinorLength0=false, NLoptMethod=:LD_MMA)
# scores: 0.0017388694714770016, 0.0017388695273561135, 0.0017388694990967278

# almost 2 polytomies: 1 edge of estimated length 2.7e-11, 1 other quite small
sort([e.length for e in calnet1o.edge]) # 2.7e-11,3.8e-5,0.00033,0.00040,...
[e.length for e in calnet1o.edge if isapprox(e.length, 0., atol=1e-5)]
[e.length for e in calnet1o.edge if e.hybrid && !e.isMajor]
# the minor hybrid edge is NOT small: ~ 0.0046848
R"par(mar=c(0,0,0,0))";
# plot(calnet1o, :R, useEdgeLength=true, showNodeNumber=true); # ingroup crown: -6
# normalize calibrated network to unit height of Adansonia crown
crownind = indexin([-6], [n.number for n in calnet1o.nodes_changed])[1] # 2
crownage = getNodeAges(calnet1o)[crownind] # 0.01974636374183969
for e in calnet1o.edge e.length /= crownage; end
#= check the normalization: plot with scale bar
plot(calnet1o, :R, useEdgeLength=true)
R"lines"(x=[1.66,2.66], y=[0.5,0.5], lwd=2);
R"text"(x=2.16, y=0.7, labels="1", adj=[0,0], cex=0.8);
=#
# save calibrated network
writeTopology(calnet1o, "snaq_results/BestH1_372g_calibrated.tre")
calnet1of = deepcopy(calnet1o)
for n in calnet1of.leaf
  n.name = replace(replace(fulltaxnamedict[n.name], ". "=>"_"), " "=>"_")
end
writeTopology(calnet1of, "snaq_results/BestH1_372g_calibrated_fullnames.tre")

#------------------------------------------------
#    primary data set (372 genes), ingroup only
#------------------------------------------------

net1i = readTopology("snaq_results/BestH1_noOut_372g.out")
rootonedge!(net1i, 9); # assuming (Africa,(Australia,Madagascar)), as when outgroups are included
# plot(net1i, :R, showEdgeLength=true, showEdgeNumber=true);
# plot(net1i, :R, useEdgeLength=true, showNodeNumber=true);
for n in [-6,-4,-2,-5,-12,-13] rotate!(net1i, n); end
# plot(net1i, :R, useEdgeLength=true);

net1if = deepcopy(net1i);
for n in net1if.leaf n.name = fulltaxnamedict[n.name]; end
hi = findfirst([!e.isMajor for e in net1if.edge]) # hi for hybrid index: 32
R"pdf"("traitanalysis/primary372_onlyIngroup.pdf", width=7, height=4.5);
R"par"(mar=[0,0,0,0]);
res = plot(net1if, :R, xlim=[1.0,17.0], tipOffset=0.1);
hx1, hx2, hy1, hy2 = (res[i][hi] for i in 9:12); # coordinates for minor hybrid edge
R"arrows"(hx1, hy1, hx2, hy2, col="deepskyblue", length=0.08, angle=20);
#R"text"(x=1, y=17.5, labels="no branch lengths", adj=[0,0], cex=0.8);
x0=14.7; x1=15.2
R"text"(res[14][hi,:x]+0.6, res[14][hi,:y], res[14][hi,:gam], col="deepskyblue", cex=0.75);
R"lines"(x=[x0,x0], y=[13,15], lwd=2); R"text"(x=x1, y=14,"Africa", adj=0);
R"lines"(x=[x0,x0], y=[2, 12], lwd=2); R"text"(x=x1, y=6, "Madagascar", adj=0);
R"lines"(x=[x0,x0], y=[0.5,1.2], lwd=2); R"text"(x=x1, y=1,"Australia", adj=0);
R"dev.off()";

# calibration prep: modify branch lengths for good starting values
calnet1i = deepcopy(net1i)
for e in calnet1i.edge
  if e.length < 0.0 e.length  = 0.001; else e.length *= 0.008; end
end
# calibration: all 3 outgroups from taxa and distance matrix, because not in network
# taxa1o = deepcopy(taxa); deleteat!(taxa1o, 3)
ind1i = sort(indexin(ingroup, taxa)) # same as collect(4:18), but more general
taxa1i = taxa[ind1i]
avD1i = avD1[ind1i,ind1i]
# run the calibration 3 times to make sure we optimized branch lengths thoroughly:
calibrateFromPairwiseDistances!(calnet1i, avD1i, taxa1i, forceMinorLength0=false, NLoptMethod=:LD_MMA, verbose=true)
calibrateFromPairwiseDistances!(calnet1i, avD1i, taxa1i, forceMinorLength0=false, NLoptMethod=:LN_COBYLA) # score got worse
calibrateFromPairwiseDistances!(calnet1i, avD1i, taxa1i, forceMinorLength0=false, NLoptMethod=:LD_MMA)
# scores: 0.002369945913157443, 0.002370139803006534, 0.002369945952360207

# 2 polytomies & horizontal minor hybrid edge:
sort([e.length for e in calnet1i.edge]) # -1.3e-9,-3.8e-11,-2.9e-11,.00040,.00059,...
# set to 0 the 3 negative branch lengths
for e in calnet1i.edge if e.length <0.0 e.length=0.0; end; end
[e.length for e in calnet1i.edge if e.hybrid && !e.isMajor] # minor hybrid edge ~ 0.0
R"par(mar=c(0,0,0,0))";
# plot(calnet1i, :R, useEdgeLength=true, showNodeNumber=true); # ingroup=root: 17
# polytomies: at crown divergence, and right after gene flow
# normalize calibrated network to unit height of Adansonia crown
crownind = indexin([17], [n.number for n in calnet1i.nodes_changed])[1] # 1: root!
crownage = getNodeAges(calnet1i)[crownind] # 0.014384170739935843
for e in calnet1i.edge e.length /= crownage; end
#= check the normalization: plot with scale bar
plot(calnet1i, :R, useEdgeLength=true);
R"lines"(x=[1,2], y=[0.5,0.5], lwd=2);
R"text"(x=1.5, y=0.7, labels="1", adj=[0,0], cex=0.8);
=#
# save calibrated network
writeTopology(calnet1i, "snaq_results/BestH1_noOut_372g_calibrated.tre")
calnet1if = deepcopy(calnet1i)
for n in calnet1if.leaf
  n.name = replace(replace(fulltaxnamedict[n.name], ". "=>"_"), " "=>"_")
end
writeTopology(calnet1if, "snaq_results/BestH1_noOut_372g_calibrated_fullnames.tre")

#------------------------------------------------
#  haplotype data set (344 genes), with outgroups
#------------------------------------------------

net2o = readTopology("snaq_results/BestH1_HHcombined.out")
rootonedge!(net2o, 26); # all 3 outgroups present
# plot(net2o, :R, useEdgeLength=true, showEdgeNumber=true);
# plot(net2o, :R, useEdgeLength=false, showNodeNumber=true);
for n in [-14,-15,-18,-4,-11,-5,-7] rotate!(net2o, n); end
net2of = deepcopy(net2o);
for n in net2of.leaf n.name = fulltaxnamedict[n.name]; end
# plot(net2of, :R, useEdgeLength=true);

# nice plot, showing topology only (not branch lengths)
hi = findfirst([!e.isMajor for e in net2of.edge]) # hi for hybrid index: 32
R"pdf"("traitanalysis/haplotype344_withOutgroups.pdf", width=7, height=5);
R"par"(mar=[0,0,0,0]);
res = plot(net2of, :R, xlim=[1.0,17.0], ylim=[0.8,18.2], tipOffset=0.1);
hx1, hx2, hy1, hy2 = (res[i][hi] for i in 9:12); # coordinates for minor hybrid edge
R"arrows"(hx1, hy1, hx2, hy2, col="deepskyblue", length=0.08, angle=20);
#R"text"(x=1, y=17.5, labels="no branch lengths", adj=[0,0], cex=0.8);
x0=14.7; x1=15.2
R"text"(res[14][hi,:x]+0.6, res[14][hi,:y], res[14][hi,:gam], col="deepskyblue", cex=0.75);
R"lines"(x=[x0,x0], y=[16,18], lwd=2); R"text"(x=x1, y=17,"outgroups", adj=0);
R"lines"(x=[x0,x0], y=[13,15], lwd=2); R"text"(x=x1, y=14,"Africa", adj=0);
R"lines"(x=[x0,x0], y=[1, 11], lwd=2); R"text"(x=x1, y=5, "Madagascar", adj=0);
R"lines"(x=[x0,x0], y=[11.8,12.2], lwd=2); R"text"(x=x1, y=12,"Australia", adj=0);
R"dev.off()";

# load pairwise distance, to prepare for network calibration
avD2 = disallowmissing(convert(Matrix, CSV.read("traitanalysis/genetreePairwiseDistance_haphunt.csv")))
# check that the column names are the same (and in the same order) as in taxa:
all(taxa .== String.(names(CSV.read("traitanalysis/genetreePairwiseDistance_haphunt.csv"))))
# modify branch lengths for good starting values
maximum(avD2) # maximum distance: 0.08294
# plot(net2o, :R, useEdgeLength=true, showEdgeNumber=true); # to see which edges to add to get max distance
sum(net2o.edge[i].length for i in [24,26,37,35,21,20,18,16]) # 7.345. 0.08294/7.345=0.011
calnet2o = deepcopy(net2o)
for e in calnet2o.edge
  if e.length < 0.0 e.length  = 0.008; else e.length *= 0.01; end
end
# run the calibration 3 times to make sure we optimized branch lengths thoroughly:
calibrateFromPairwiseDistances!(calnet2o, avD2, taxa, forceMinorLength0=false, NLoptMethod=:LD_MMA)
calibrateFromPairwiseDistances!(calnet2o, avD2, taxa, forceMinorLength0=false, NLoptMethod=:LN_COBYLA) # score got worse
calibrateFromPairwiseDistances!(calnet2o, avD2, taxa, forceMinorLength0=false, NLoptMethod=:LD_MMA)
# scores: 0.006979654701045381, 0.0069797352144146665, 0.006979654640760858

# 7 edges of estimated length < 3e-9, 1.2e-5,.00036,.001,...
sort([e.length for e in calnet2o.edge]) # 1 is negative
for e in calnet2o.edge if e.length < 0.0 e.length=0.0; end; end
[e.length for e in calnet2o.edge if isapprox(e.length, 0., atol=1e-5)]
[e.length for e in calnet2o.edge if e.hybrid && !e.isMajor]
# the minor hybrid edge is ~ 0
R"par(mar=c(0,0,0,0))";
# plot(calnet2o, :R, useEdgeLength=true, showNodeNumber=true); # ingroup crown: -13
# normalize calibrated network to unit height of Adansonia crown
crownind = indexin([-13], [n.number for n in calnet2o.nodes_changed])[1]
crownage = getNodeAges(calnet2o)[crownind] # 0.022535188873004656
for e in calnet2o.edge e.length /= crownage; end
#= check the normalization: plot with scale bar
plot(calnet2o, :R, useEdgeLength=true)
R"lines"(x=[1.44,2.44], y=[0.5,0.5], lwd=2);
R"text"(x=1.94, y=0.7, labels="1", adj=[0,0], cex=0.8);
=#
# save calibrated network
writeTopology(calnet2o, "snaq_results/BestH1_HHcombined_calibrated.tre")
calnet2of = deepcopy(calnet2o)
for n in calnet2of.leaf
  n.name = replace(replace(fulltaxnamedict[n.name], ". "=>"_"), " "=>"_")
end
writeTopology(calnet2of, "snaq_results/BestH1_HHcombined_calibrated_fullnames.tre")

#------------------------------------------------
#  haplotype data set (344 genes), ingroup only
#------------------------------------------------

net2i = readTopology("snaq_results/BestH1_noOutgroups_HHcombined.out")
rootonedge!(net2i, 27);
# assuming ((Africa,Australia),Madagascar), as when outgroups are included
# root on edge 25 for (Africa,(Australia,Madagascar))
# plot(net2i, :R, showEdgeLength=true, showEdgeNumber=true);
# plot(net2i, :R, useEdgeLength=true, showNodeNumber=true);
for n in [-15,-6,-9,-10,-12] rotate!(net2i, n); end
net2if = deepcopy(net2i);
for n in net2if.leaf n.name = fulltaxnamedict[n.name]; end
# plot(net2i, :R, useEdgeLength=false);

# nice plot, showing topology only (not branch lengths)
hi = findfirst([!e.isMajor for e in net2if.edge]) # hi for hybrid index: 32
R"pdf"("traitanalysis/haplotype344_onlyIngroup.pdf", width=7, height=5);
R"par"(mar=[0,0,0,0]);
res = plot(net2if, :R, xlim=[1.0,17.0], ylim=[0.8,15.2], tipOffset=0.1);
hx1, hx2, hy1, hy2 = (res[i][hi] for i in 9:12); # coordinates for minor hybrid edge
R"arrows"(hx1, hy1, hx2, hy2, col="deepskyblue", length=0.08, angle=20);
#R"text"(x=1, y=17.5, labels="no branch lengths", adj=[0,0], cex=0.8);
x0=14.7; x1=15.2
R"text"(res[14][hi,:x]+0.6, res[14][hi,:y], res[14][hi,:gam], col="deepskyblue", cex=0.75);
R"lines"(x=[x0,x0], y=[13,15], lwd=2); R"text"(x=x1, y=14,"Africa", adj=0);
R"lines"(x=[x0,x0], y=[1, 11], lwd=2); R"text"(x=x1, y=5, "Madagascar", adj=0);
R"lines"(x=[x0,x0], y=[11.8,12.2], lwd=2); R"text"(x=x1, y=12,"Australia", adj=0);
R"dev.off()";

# calibration prep: modify branch lengths for good starting values
calnet2i = deepcopy(net2i)
for e in calnet2i.edge
  if e.length < 0.0 e.length  = 0.008; else e.length *= 0.01; end
end
# prune all 3 outgroups from taxa and distance matrix, because not in network
ind2i = sort(indexin(ingroup, taxa)) # same as collect(4:18), but more general
taxa2i = taxa[ind2i]
avD2i = avD2[ind2i,ind2i]
# run the calibration 3 times to make sure we optimized branch lengths thoroughly:
calibrateFromPairwiseDistances!(calnet2i, avD2i, taxa2i, forceMinorLength0=false, NLoptMethod=:LD_MMA)
calibrateFromPairwiseDistances!(calnet2i, avD2i, taxa2i, forceMinorLength0=false, NLoptMethod=:LN_COBYLA) # score got worse
calibrateFromPairwiseDistances!(calnet2i, avD2i, taxa2i, forceMinorLength0=false, NLoptMethod=:LD_MMA)
# scores: 7.24570473e-5, 7.26570834e-5, 7.2457021984e-5

# 4 polytomies: branch lengths of -1.6e-11,1.8e-11,3.4e-11,1.2e-5,0.00035,...
sort([e.length for e in calnet2i.edge]) # 1 is negative
for e in calnet2i.edge if e.length < 0.0 e.length=0.0; end; end
[e.length for e in calnet2i.edge if isapprox(e.length, 0., atol=1e-5)]
[e.length for e in calnet2i.edge if e.hybrid && !e.isMajor]
# the minor hybrid edge is ~ 0
R"par(mar=c(0,0,0,0))";
# plot(calnet2i, :R, useEdgeLength=true, showNodeNumber=true); # ingroup=root=17
# normalize calibrated network to unit height of Adansonia crown
crownind = indexin([17], [n.number for n in calnet2i.nodes_changed])[1] # 1: root!
crownage = getNodeAges(calnet2i)[crownind] # 0.0216794187731583
for e in calnet2i.edge e.length /= crownage; end
#= check the normalization: plot with scale bar
plot(calnet2i, :R, useEdgeLength=true)
R"lines"(x=[1,2], y=[0.5,0.5], lwd=2);
R"text"(x=1.5, y=0.7, labels="1", adj=[0,0], cex=0.8);
=#
# save calibrated network
writeTopology(calnet2i, "snaq_results/BestH1_noOutgroups_HHcombined_calibrated.tre")
calnet2if = deepcopy(calnet2i)
for n in calnet2if.leaf
  n.name = replace(replace(fulltaxnamedict[n.name], ". "=>"_"), " "=>"_")
end
writeTopology(calnet2if, "snaq_results/BestH1_noOutgroups_HHcombined_calibrated_fullnames.tre")

#---------------------------------
#      package version details
#---------------------------------

"""
(v1.1) pkg> status
    Status `~/.julia/environments/v1.1/Project.toml`
  [7e6ae17a] BioSequences v1.1.0
  [3c28c6f8] BioSymbols v3.1.0
  [336ed68f] CSV v0.5.5
  [a93c6f00] DataFrames v0.18.3
  [31c24e10] Distributions v0.20.0
  [33ad39ac] PhyloNetworks v0.9.1
  [c0d5b6db] PhyloPlots v0.2.0
  [6f49c342] RCall v0.13.2
  [2913bbd2] StatsBase v0.30.0
  [3eaba693] StatsModels v0.5.0
"""
