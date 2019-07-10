# software versions: julia v1.1.1, see details at the end
using Statistics
using StatsBase
using DataFrames
using CSV
using PhyloNetworks
using RCall
using PhyloPlots
# @rlibrary ape
# ERROR: invalid redefinition of constant node_depth_edgelength

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

#--------------------------------------------------
#     load the various networks and data
#--------------------------------------------------

# load the various phylogenies
calnet1o = readTopology("snaq_results/BestH1_372g_calibrated.tre")
calnet1i = readTopology("snaq_results/BestH1_noOut_372g_calibrated.tre")
calnet2o = readTopology("snaq_results/BestH1_HHcombined_calibrated.tre")
calnet2i = readTopology("snaq_results/BestH1_noOutgroups_HHcombined_calibrated.tre")

# load the trait data: flower color and pollinator
dat = CSV.read("traitanalysis/Adansonia_floralTraits_accession.csv",
               categorical=false, copycols=true)
dropmissing!(dat, :flowercolor) # disallowmissing!

# plot the data on each network
"""
    plotnet(network, pdf file name)

assume:
- a single hybrid node in the network
- trait data in `dat` (external object)
"""
function plotnet(net, file; xlim=[1.0,4.0], tipoffset=0.15, xleg=1.2, onfile=true)
  hi = findfirst([!e.isMajor for e in net.edge]) # "h"ybrid "i"ndex
  o = [findfirst(isequal(tax), dat[:accession]) for tax in tipLabels(net)]
  msgn = o.==nothing
  pollab = dat[:pollination][o[.!msgn]];
  polpch = map(x -> (ismissing(x) ? 24 : (x=="mammal" ? 21 : 22)), pollab);
  flolab = dat[:flowercolor][o[.!msgn]];
  flocol = map(x -> (ismissing(x) ? "grey" : x), flolab);
  if onfile R"pdf"(file, width=6, height=5); end
  R"par"(mar=[0,0,0,0]);
  res = plot(net, :R, useEdgeLength=true, xlim=xlim, tipOffset=tipoffset);
  hx1, hx2, hy1, hy2 = (res[i][hi] for i in 9:12); # coordinates for minor hybrid edge
  R"arrows"(hx1, hy1, hx2, hy2, col="deepskyblue", length=0.08, angle=20);
  R"text"(res[14][hi,:x]+0.03*(xlim[2]-xlim[1]), res[14][hi,:y], res[14][hi,:gam], col="deepskyblue", cex=0.75);
  R"points"(x=res[13][:x][.!msgn].+0.05, y=res[13][:y][.!msgn], pch=polpch, bg=flocol, cex=1.5);
  R"legend"(x=xleg, y=5, legend=["mammal","hawkmoth"], pch=21:22, bty="n", title="pollinator", var"title.adj"=0);
  colnames = ["white","yellow","red"];
  R"legend"(x=xleg, y=9, legend=colnames, pch=22, var"pt.bg"=colnames, bty="n",
          title="flower color", var"title.adj"=0);
  if onfile R"dev.off"(); end
end

#=
plotnet(calnet1o, "traitanalysis/figures/primary372_withOutgroups_traits.pdf"; xlim=[1,3]);
plotnet(calnet1i, "traitanalysis/figures/primary372_onlyIngroup_traits.pdf"; xlim=[1,2.3], xleg=1.01);
plotnet(calnet2o, "traitanalysis/figures/haplotype344_withOutgroups_traits.pdf"; xlim=[1,2.8]);
plotnet(calnet2i, "traitanalysis/figures/haplotype344_onlyIngroup_traits.pdf"; xlim=[1,2.3], xleg=1.01);
=#

"""
    plotasr(network, ancestral state predictions, pdf file name)

assume:
- a single hybrid node in the network
- full taxon names taken from external dictionary `fulltaxnamedict`
"""
function plotasr(fit, file; xlim=[1,3], ylim=[0.8,15.2], tipoffset=0.15,
                 onfile=true, coltrait=["brown","grey"], radius=0.016,
                 legpos="topright", legendtitle="", legend=true)
    asr = ancestralStateReconstruction(fit)
    colnames = names(asr)[3:end] # column names
    asr[:fake] = "";
    net = deepcopy(fit.net)
    for n in net.leaf
        n.name = fulltaxnamedict[n.name]
    end
    R"library(ape)"; # floating.pie.asp not visible from within R otherwise
    if onfile R"pdf"(file, width=6.5, height=5); end
    R"par"(mar=[0,0,0,0]);
    res = plot(net, :R, nodeLabel = asr[[:nodenumber, :fake]], ylim=ylim,
            useEdgeLength=true, xlim=xlim, tipOffset=tipoffset);
    hi = findfirst([!e.isMajor for e in net.edge]) # "h"ybrid "i"ndex
    hx1, hx2, hy1, hy2 = (res[i][hi] for i in 9:12); # coordinates for minor hybrid edge
    R"arrows"(hx1, hy1, hx1+0.95*(hx2-hx1), hy1+0.95*(hy2-hy1), col="deepskyblue", length=0.08, angle=20);
    R"text"(res[14][hi,:x]+0.03*(xlim[2]-xlim[1]), res[14][hi,:y], res[14][hi,:gam], col="deepskyblue", cex=0.75);
    for i in 1:nrow(asr)
        ii = findfirst(isequal(string(asr[:nodenumber][i])), res[13][:num]);
        colpp = convert(Vector{Float64}, asr[i,colnames]);
        R"floating.pie.asp"(res[13][ii,:x], res[13][ii,:y], colpp, radius=radius, col=coltrait);
    end
    if legend
    R"legend"(legpos, legend=colnames, pch=21, var"pt.bg"=coltrait,
       bty="n", title=legendtitle, var"title.adj"=0, var"pt.cex"=1.5);
    end
    if onfile R"dev.off()"; end
    return nothing
end

"""
Bayes factor for inheritance via the minor parent.

Assuming h=1: only 2 displayed trees, tree 2 and tree 1 corresponding to the
major & minor edges at one reticulation, so we want the
BF for trait inherited via tree 2 versus tree 1.

Bayes factor:

P(data|parent1)/P(data|parent2) =
posterior(parent1|data)/posterior(parent2|data) * prior(parent2)/prior(parent1)
"""
function geneflowBF(fit)
    exp(fit.postltw[2] - fit.postltw[1] + fit.priorltw[1] - fit.priorltw[2])
end
#--------------------------------------------------
#     fit trait models: preparation
#--------------------------------------------------


# preparation for pollinator
ipol = completecases(dat, :pollination);
pollevels = disallowmissing(unique(dat[:pollination][ipol]))
mpol = BinaryTraitSubstitutionModel([0.1, 0.1], pollevels);
mpol_equal = EqualRatesSubstitutionModel(2, 0.1, pollevels);
datpol = dropmissing!(deepcopy(dat[[:accession,:pollination]]))
# Set(tipLabels(calnet1o_pruned)) == Set(datpol[:accession]) # true

# preparation for flower color
dat[:pigmentation] = map(x -> x=="white" ? "white" : "pigmented", dat[:flowercolor])
collevels = unique(dat[:pigmentation])
mcol = BinaryTraitSubstitutionModel([0.1, 0.1], collevels);
mcol_equal = EqualRatesSubstitutionModel(2, 0.1, collevels);
datcol = dat[[:accession,:pigmentation]]

"""
    deleteoutgroups!(net)

Delete outgroup with custom names "Pcr070","Smi165","Bce020"
from a network, keeping the network rooted at the ingroup crown node.

This is a fix for a bug with the pruned network. The bug is that
this does not work: the loglikelihood score is always log(1/2) regardless
of the rates, so the rates don't change during the optimization.

```julia
fitDiscrete(calnet1o, mpol, datpol, verbose=true) # fails
calnet1o_pruned = deepcopy(calnet1o);
for species in ["Pcr070","Smi165"] deleteleaf!(net, species); end
fitDiscrete(calnet1o_pruned, mpol, datpol, verbose=true) # fails
# the plot below shows that Age001 is the root
plot(calnet1o_pruned, :R, showEdgeLength=true, showEdgeNumber=true, useEdgeLength=true)
calnet1o_pruned.node[calnet1o_pruned.root].name # root = leaf "Age001"!
```

`deleteoutgroups!(net)` deletes the 2 or 3 outgroups
without destroying the branch lengths at the base of the ingroup.
`deleteleaf!` currently (v0.9.1) messes up the branch lengths
at the base of the ingroup & unroots the network, that is:
- one ingroup edge disappears,
- the other ingroup edge (sister originally) is too long.
One edge is restored with `rootonedge!`, then both branch lengths
are manually restored.
"""
function deleteoutgroups!(net)
    directEdges!(net)
    sminum = net.node[findfirst(n -> n.name=="Smi165", net.node)].number
    rootedge = net.node[net.root].edge
    stemingroup = findfirst(e -> sminum ∉ PhyloNetworks.descendants(e), rootedge)
    crown = PhyloNetworks.getChild(rootedge[stemingroup]) # crown node
    basaledge = [e for e in crown.edge if PhyloNetworks.getParent(e)==crown]
    basaledgenum = [e.number for e in basaledge]
    basaledgeBL = [e.length for e in basaledge]
    for species in ["Pcr070","Smi165","Bce020"]
        if species in tipLabels(net) deleteleaf!(net, species); end
    end
    edge1_ind = findfirst(in(basaledgenum), [e.number for e in net.edge])
    edge1 = net.edge[edge1_ind]
    rootonedge!(net, edge1)
    edge1_o = findfirst(isequal(edge1.number), basaledgenum)
    edge1.length = basaledgeBL[edge1_o]
    edge2_o   = (edge1_o  ==1 ? 2 : 1)
    rootedge = net.node[net.root].edge
    edge2_ind = findfirst(e -> e.number != edge1.number, rootedge)
    edge2 = rootedge[edge2_ind]
    edge2.length = basaledgeBL[edge2_o]
    return net
end

#--------------------------------------------------
#     fit pollinator models
#--------------------------------------------------

# pollinator: primary with outgroups
#-----------------------------------

calnet1o_pruned = deepcopy(calnet1o);
deleteoutgroups!(calnet1o_pruned)
#plot(calnet1o, :R, showEdgeLength=true, showEdgeNumber=true)
#plot(calnet1o_pruned, :R, showEdgeLength=true, showEdgeNumber=true, useEdgeLength=true)

# unconstrained rates:
fitp1o = fitDiscrete(calnet1o_pruned, mpol, datpol)
aic(fitp1o)
#=
rate mammal→hawkmoth α=0.8462
rate hawkmoth→mammal β=0.0
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -6.8396
aic: 17.679203122276146
=#

# equal rates: fits better
fitp1o_equal = fitDiscrete(calnet1o_pruned, mpol_equal, datpol)
aic(fitp1o_equal)
show(ancestralStateReconstruction(fitp1o_equal), allrows=true)
plotasr(fitp1o_equal, "traitanalysis/figures/primary372_withOutgroups_ase_pol.pdf")
#=
            mammal hawkmoth
   mammal        *   0.4804
 hawkmoth   0.4804        *
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -7.10813
aic: 16.21625857933925
31×4 DataFrame
│ Row │ nodenumber │ nodelabel │ mammal      │ hawkmoth   │
│     │ Int64      │ String    │ Float64     │ Float64    │
├─────┼────────────┼───────────┼─────────────┼────────────┤
│ 1   │ 1          │ Adi001    │ 1.0         │ 0.0        │
│ 2   │ 2          │ Adi003    │ 1.0         │ 0.0        │
│ 3   │ 3          │ Adi002    │ 1.0         │ 0.0        │
│ 4   │ 4          │ Asu001    │ 1.0         │ 0.0        │
│ 5   │ 5          │ Aga001    │ 1.0         │ 0.0        │
│ 6   │ 6          │ Aga002    │ 1.0         │ 0.0        │
│ 7   │ 7          │ Aza135    │ 0.0         │ 1.0        │
│ 8   │ 8          │ Aza037    │ 0.0         │ 1.0        │
│ 9   │ 9          │ Ama018    │ 0.0         │ 1.0        │
│ 10  │ 10         │ Ama006    │ 0.0         │ 1.0        │
│ 11  │ 11         │ Ape009    │ 0.0         │ 1.0        │
│ 12  │ 12         │ Ape001    │ 0.0         │ 1.0        │
│ 13  │ 13         │ Aru001    │ 0.0         │ 1.0        │
│ 14  │ 14         │ Aru127    │ 0.0         │ 1.0        │
│ 15  │ 15         │ Age001    │ 0.0         │ 1.0        │
│ 16  │ 16         │ 16        │ 0.369518    │ 0.630482   │
│ 17  │ 17         │ 17        │ 0.464133    │ 0.535867   │
│ 18  │ 18         │ 18        │ 0.230936    │ 0.769064   │
│ 19  │ 19         │ 19        │ 0.0210102   │ 0.97899    │
│ 20  │ 20         │ 20        │ 0.230255    │ 0.769745   │
│ 21  │ 21         │ 21        │ 0.0138985   │ 0.986102   │
│ 22  │ 22         │ 22        │ 0.000645436 │ 0.999355   │
│ 23  │ 23         │ 23        │ 0.000601462 │ 0.999399   │
│ 24  │ 24         │ 24        │ 0.000239543 │ 0.99976    │
│ 25  │ 25         │ 25        │ 0.00037092  │ 0.999629   │
│ 26  │ 26         │ 26        │ 0.464133    │ 0.535867   │
│ 27  │ 27         │ 27        │ 0.991841    │ 0.00815857 │
│ 28  │ 28         │ 28        │ 0.99538     │ 0.00461983 │
│ 29  │ 29         │ #H18      │ 0.877681    │ 0.122319   │
│ 30  │ 30         │ 30        │ 0.885554    │ 0.114446   │
│ 31  │ 31         │ 31        │ 0.978093    │ 0.0219073  │
=#

# contribution of gene flow to pollinator
exp.(fitp1o_equal.postltw)
exp.(fitp1o_equal.priorltw)
geneflowBF(fitp1o_equal) # 3.0077799559628824
#=
marginal (posterior) probability that the trait evolved on each displayed tree:
0.7006004378846591,0.2993995621153409 vs prior: 0.8756,0.1244
=#

# pollinator: primary ingroup only
#-----------------------------------

# unconstrained rates:
fitp1i = fitDiscrete(calnet1i, mpol, dat[[:accession,:pollination]])
aic(fitp1i)
#=
rate mammal→hawkmoth α=0.0
rate hawkmoth→mammal β=0.4161
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -6.62297
aic: 17.245933822077305
=#

# equal rates: fits better
fitp1i_equal = fitDiscrete(calnet1i, mpol_equal, dat[[:accession,:pollination]])
aic(fitp1i_equal)
show(ancestralStateReconstruction(fitp1i_equal), allrows=true)
plotasr(fitp1i_equal, "traitanalysis/figures/primary372_onlyIngroup_ase_pol.pdf")
exp.(fitp1i_equal.postltw)  # 0.6684740218348195, 0.3315259781651807
exp.(fitp1i_equal.priorltw) # 0.8304256934254225 0.16957430657457748
geneflowBF(fitp1i_equal)    # 2.428699361757745
#=
            mammal hawkmoth
   mammal        *   0.3439
 hawkmoth   0.3439        *
log-likelihood: -6.88785
aic: 15.775704888161254
31×4 DataFrame
│ Row │ nodenumber │ nodelabel │ mammal      │ hawkmoth   │
│     │ Int64      │ String    │ Float64     │ Float64    │
├─────┼────────────┼───────────┼─────────────┼────────────┤
│ 1   │ 1          │ Adi002    │ 1.0         │ 0.0        │
│ 2   │ 2          │ Adi003    │ 1.0         │ 0.0        │
│ 3   │ 3          │ Adi001    │ 1.0         │ 0.0        │
│ 4   │ 4          │ Asu001    │ 1.0         │ 0.0        │
│ 5   │ 5          │ Aga001    │ 1.0         │ 0.0        │
│ 6   │ 6          │ Aga002    │ 1.0         │ 0.0        │
│ 7   │ 7          │ Aza135    │ 0.0         │ 1.0        │
│ 8   │ 8          │ Aza037    │ 0.0         │ 1.0        │
│ 9   │ 9          │ Ama018    │ 0.0         │ 1.0        │
│ 10  │ 10         │ Ama006    │ 0.0         │ 1.0        │
│ 11  │ 11         │ Ape009    │ 0.0         │ 1.0        │
│ 12  │ 12         │ Ape001    │ 0.0         │ 1.0        │
│ 13  │ 13         │ Aru127    │ 0.0         │ 1.0        │
│ 14  │ 14         │ Aru001    │ 0.0         │ 1.0        │
│ 15  │ 15         │ Age001    │ 0.0         │ 1.0        │
│ 16  │ 16         │ 16        │ 0.372088    │ 0.627912   │
│ 17  │ 17         │ 17        │ 0.372088    │ 0.627912   │
│ 18  │ 18         │ 18        │ 0.303884    │ 0.696116   │
│ 19  │ 19         │ 19        │ 0.0572633   │ 0.942737   │
│ 20  │ 20         │ 20        │ 0.0101631   │ 0.989837   │
│ 21  │ 21         │ 21        │ 0.334003    │ 0.665997   │
│ 22  │ 22         │ #H16      │ 0.00725974  │ 0.99274    │
│ 23  │ 23         │ 23        │ 0.00725974  │ 0.99274    │
│ 24  │ 24         │ 24        │ 0.000339542 │ 0.99966    │
│ 25  │ 25         │ 25        │ 0.000511302 │ 0.999489   │
│ 26  │ 26         │ 26        │ 0.000155894 │ 0.999844   │
│ 27  │ 27         │ 27        │ 0.000330931 │ 0.999669   │
│ 28  │ 28         │ 28        │ 0.875103    │ 0.124897   │
│ 29  │ 29         │ 29        │ 0.97655     │ 0.0234502  │
│ 30  │ 30         │ 30        │ 0.992811    │ 0.00718883 │
│ 31  │ 31         │ 31        │ 0.995962    │ 0.00403829 │
=#


# pollinator: haplotype with outgroups
#-------------------------------------

calnet2o_pruned = deepcopy(calnet2o);
deleteoutgroups!(calnet2o_pruned)

# unconstrained rates:
fitp2o = fitDiscrete(calnet2o_pruned, mpol, datpol)
aic(fitp2o)
#=
rate mammal→hawkmoth α=0.80376
rate hawkmoth→mammal β=0.6995
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -7.35788
aic: 18.71575574701439
=#

# equal rates: fits better
fitp2o_equal = fitDiscrete(calnet2o_pruned, mpol_equal, datpol)
aic(fitp2o_equal)
show(ancestralStateReconstruction(fitp2o_equal), allrows=true)
plotasr(fitp2o_equal, "traitanalysis/figures/haplotype344_withOutgroups_ase_pol.pdf")
exp.(fitp2o_equal.postltw)  # 0.9649790114890178, 0.03502098851098212
exp.(fitp2o_equal.priorltw) # 0.928069201055275, 0.07193079894472504
geneflowBF(fitp2o_equal)    # BF = 0.468248084321128
#=
            mammal hawkmoth
   mammal        *   0.7565
 hawkmoth   0.7565        *
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -7.36448
aic: 16.728951397896637
31×4 DataFrame
│ Row │ nodenumber │ nodelabel │ mammal      │ hawkmoth  │
│     │ Int64      │ String    │ Float64     │ Float64   │
├─────┼────────────┼───────────┼─────────────┼───────────┤
│ 1   │ 1          │ Adi001    │ 1.0         │ 0.0       │
│ 2   │ 2          │ Adi002    │ 1.0         │ 0.0       │
│ 3   │ 3          │ Adi003    │ 1.0         │ 0.0       │
│ 4   │ 4          │ Age001    │ 0.0         │ 1.0       │
│ 5   │ 5          │ Asu001    │ 1.0         │ 0.0       │
│ 6   │ 6          │ Aga001    │ 1.0         │ 0.0       │
│ 7   │ 7          │ Aga002    │ 1.0         │ 0.0       │
│ 8   │ 8          │ Aza135    │ 0.0         │ 1.0       │
│ 9   │ 9          │ Aza037    │ 0.0         │ 1.0       │
│ 10  │ 10         │ Ama006    │ 0.0         │ 1.0       │
│ 11  │ 11         │ Ama018    │ 0.0         │ 1.0       │
│ 12  │ 12         │ Ape001    │ 0.0         │ 1.0       │
│ 13  │ 13         │ Ape009    │ 0.0         │ 1.0       │
│ 14  │ 14         │ Aru001    │ 0.0         │ 1.0       │
│ 15  │ 15         │ Aru127    │ 0.0         │ 1.0       │
│ 16  │ 16         │ 16        │ 0.53424     │ 0.46576   │
│ 17  │ 17         │ 17        │ 0.212056    │ 0.787944  │
│ 18  │ 18         │ 18        │ 0.0230316   │ 0.976968  │
│ 19  │ 19         │ 19        │ 0.212056    │ 0.787944  │
│ 20  │ 20         │ 20        │ 0.00336996  │ 0.99663   │
│ 21  │ 21         │ 21        │ 0.000142936 │ 0.999857  │
│ 22  │ 22         │ 22        │ 0.000578102 │ 0.999422  │
│ 23  │ 23         │ 23        │ 0.000140004 │ 0.99986   │
│ 24  │ 24         │ 24        │ 0.000573917 │ 0.999426  │
│ 25  │ 25         │ 25        │ 0.212056    │ 0.787944  │
│ 26  │ 26         │ #H19      │ 0.890839    │ 0.109161  │
│ 27  │ 27         │ 27        │ 0.890839    │ 0.109161  │
│ 28  │ 28         │ 28        │ 0.890839    │ 0.109161  │
│ 29  │ 29         │ 29        │ 0.890839    │ 0.109161  │
│ 30  │ 30         │ 30        │ 0.925283    │ 0.0747174 │
│ 31  │ 31         │ 31        │ 0.988154    │ 0.011846  │
=#

# pollinator: haplotype ingroup only
#-------------------------------------

# unconstrained rates:
fitp2i = fitDiscrete(calnet2i, mpol, datpol)
aic(fitp2i)
#=
rate mammal→hawkmoth α=0.0
rate hawkmoth→mammal β=0.55952
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -6.62245
aic: 17.244893933108358
=#

# equal rates: fits better
fitp2i_equal = fitDiscrete(calnet2i, mpol_equal, datpol)
aic(fitp2i_equal)
show(ancestralStateReconstruction(fitp2i_equal), allrows=true)
plotasr(fitp2i_equal, "traitanalysis/figures/haplotype344_onlyIngroup_ase_pol.pdf")
exp.(fitp2i_equal.postltw)  # 0.7215360254195731, 0.2784639745804271
exp.(fitp2i_equal.priorltw) # 0.8723075169958534, 0.12769248300414657
geneflowBF(fitp2i_equal)    # 2.6364241479560104
#=
            mammal hawkmoth
   mammal        *   0.4982
 hawkmoth   0.4982        *
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -6.75842
aic: 15.516840448250345
31×4 DataFrame
│ Row │ nodenumber │ nodelabel │ mammal      │ hawkmoth   │
│     │ Int64      │ String    │ Float64     │ Float64    │
├─────┼────────────┼───────────┼─────────────┼────────────┤
│ 1   │ 1          │ Adi001    │ 1.0         │ 0.0        │
│ 2   │ 2          │ Adi002    │ 1.0         │ 0.0        │
│ 3   │ 3          │ Adi003    │ 1.0         │ 0.0        │
│ 4   │ 4          │ Age001    │ 0.0         │ 1.0        │
│ 5   │ 5          │ Asu001    │ 1.0         │ 0.0        │
│ 6   │ 6          │ Aga001    │ 1.0         │ 0.0        │
│ 7   │ 7          │ Aga002    │ 1.0         │ 0.0        │
│ 8   │ 8          │ Aza135    │ 0.0         │ 1.0        │
│ 9   │ 9          │ Aza037    │ 0.0         │ 1.0        │
│ 10  │ 10         │ Ama018    │ 0.0         │ 1.0        │
│ 11  │ 11         │ Ama006    │ 0.0         │ 1.0        │
│ 12  │ 12         │ Ape001    │ 0.0         │ 1.0        │
│ 13  │ 13         │ Ape009    │ 0.0         │ 1.0        │
│ 14  │ 14         │ Aru127    │ 0.0         │ 1.0        │
│ 15  │ 15         │ Aru001    │ 0.0         │ 1.0        │
│ 16  │ 16         │ 16        │ 0.462012    │ 0.537988   │
│ 17  │ 17         │ 17        │ 0.248735    │ 0.751265   │
│ 18  │ 18         │ 18        │ 0.0269252   │ 0.973075   │
│ 19  │ 19         │ 19        │ 0.00785307  │ 0.992147   │
│ 20  │ 20         │ 20        │ 0.281056    │ 0.718944   │
│ 21  │ 21         │ #H16      │ 0.00157456  │ 0.998425   │
│ 22  │ 22         │ 22        │ 0.00157456  │ 0.998425   │
│ 23  │ 23         │ 23        │ 2.90668e-5  │ 0.999971   │
│ 24  │ 24         │ 24        │ 0.00018089  │ 0.999819   │
│ 25  │ 25         │ 25        │ 2.81216e-5  │ 0.999972   │
│ 26  │ 26         │ 26        │ 0.000179662 │ 0.99982    │
│ 27  │ 27         │ 27        │ 0.929591    │ 0.0704091  │
│ 28  │ 28         │ 28        │ 0.992217    │ 0.00778312 │
│ 29  │ 29         │ 29        │ 0.462012    │ 0.537988   │
│ 30  │ 30         │ 30        │ 0.97042     │ 0.0295801  │
│ 31  │ 31         │ 31        │ 0.977831    │ 0.0221687  │
=#

#--------------------------------------------------
#     flower color models: 2 states
#--------------------------------------------------

# pigmentation: primary with outgroups
#-------------------------------------

calnet1o_pruned = deepcopy(calnet1o);
deleteleaf!(calnet1o_pruned,"Pcr070")

# unconstrained rates:
fitc1o = fitDiscrete(calnet1o_pruned, mcol, datcol)
aic(fitc1o)
#=
rate white→pigmented α=0.22505
rate pigmented→white β=0.0
1 traits, 16 species, on a network with 1 reticulations
log-likelihood: -7.35164
aic: 18.70328778409157
=#

# equal rates: fits better
fitc1o_equal = fitDiscrete(calnet1o_pruned, mcol_equal, datcol)
aic(fitc1o_equal)
show(ancestralStateReconstruction(fitc1o_equal), allrows=true)
plotasr(fitc1o_equal, "traitanalysis/figures/primary372_withOutgroups_ase_col.pdf",
        coltrait=["white","orange"], legpos="topleft", ylim=[0.8,16.2], xlim=[1,3.4])
exp.(fitc1o_equal.postltw)  # 0.4679756876353036, 0.5320243123646963
exp.(fitc1o_equal.priorltw) # 0.8755949666644305, 0.12440503333556951
geneflowBF(fitc1o_equal)    # BF = 8.001538374191426
#=
               white pigmented
     white         *    0.2578
 pigmented    0.2578         *
1 traits, 16 species, on a network with 1 reticulations
log-likelihood: -7.52447
aic: 17.048930650086135
33×4 DataFrame
│ Row │ nodenumber │ nodelabel │ white       │ pigmented   │
│     │ Int64      │ String    │ Float64     │ Float64     │
├─────┼────────────┼───────────┼─────────────┼─────────────┤
│ 1   │ 1          │ Smi165    │ 1.0         │ 0.0         │
│ 2   │ 2          │ Adi001    │ 1.0         │ 0.0         │
│ 3   │ 3          │ Adi003    │ 1.0         │ 0.0         │
│ 4   │ 4          │ Adi002    │ 1.0         │ 0.0         │
│ 5   │ 5          │ Asu001    │ 1.0         │ 0.0         │
│ 6   │ 6          │ Aga001    │ 1.0         │ 0.0         │
│ 7   │ 7          │ Aga002    │ 1.0         │ 0.0         │
│ 8   │ 8          │ Aza135    │ 0.0         │ 1.0         │
│ 9   │ 9          │ Aza037    │ 0.0         │ 1.0         │
│ 10  │ 10         │ Ama018    │ 0.0         │ 1.0         │
│ 11  │ 11         │ Ama006    │ 0.0         │ 1.0         │
│ 12  │ 12         │ Ape009    │ 0.0         │ 1.0         │
│ 13  │ 13         │ Ape001    │ 0.0         │ 1.0         │
│ 14  │ 14         │ Aru001    │ 0.0         │ 1.0         │
│ 15  │ 15         │ Aru127    │ 0.0         │ 1.0         │
│ 16  │ 16         │ Age001    │ 1.0         │ 0.0         │
│ 17  │ 17         │ 17        │ 0.891483    │ 0.108517    │
│ 18  │ 18         │ 18        │ 0.930281    │ 0.0697194   │
│ 19  │ 19         │ 19        │ 0.899755    │ 0.100245    │
│ 20  │ 20         │ 20        │ 0.327796    │ 0.672204    │
│ 21  │ 21         │ 21        │ 0.0159165   │ 0.984084    │
│ 22  │ 22         │ 22        │ 0.322711    │ 0.677289    │
│ 23  │ 23         │ 23        │ 0.0104269   │ 0.989573    │
│ 24  │ 24         │ 24        │ 0.000203307 │ 0.999797    │
│ 25  │ 25         │ 25        │ 9.82146e-5  │ 0.999902    │
│ 26  │ 26         │ 26        │ 3.81599e-5  │ 0.999962    │
│ 27  │ 27         │ 27        │ 5.37961e-5  │ 0.999946    │
│ 28  │ 28         │ 28        │ 0.899755    │ 0.100245    │
│ 29  │ 29         │ 29        │ 0.999354    │ 0.000645909 │
│ 30  │ 30         │ 30        │ 0.999724    │ 0.000276274 │
│ 31  │ 31         │ #H18      │ 0.974545    │ 0.0254554   │
│ 32  │ 32         │ 32        │ 0.97765     │ 0.0223505   │
│ 33  │ 33         │ 33        │ 0.997433    │ 0.00256721  │
=#

# pigmentation: primary ingroup only
#-------------------------------------

# unconstrained rates:
fitc1i = fitDiscrete(calnet1i, mcol, datcol)
aic(fitc1i)
#=
rate white→pigmented α=0.29114
rate pigmented→white β=0.0
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -6.10205
aic: 16.204103200981013
=#

# equal rates: fits better
fitc1i_equal = fitDiscrete(calnet1i, mcol_equal, datcol)
aic(fitc1i_equal)
show(ancestralStateReconstruction(fitc1i_equal), allrows=true)
plotasr(fitc1i_equal, "traitanalysis/figures/primary372_onlyIngroup_ase_col.pdf",
        coltrait=["white","orange"], ylim=[0.8,15.2], xlim=[1,3.2])
exp.(fitc1i_equal.postltw)  # 0.3998050420528541, 0.6001949579471464
exp.(fitc1i_equal.priorltw) # 0.8304256934254225, 0.16957430657457748
geneflowBF(fitc1i_equal)    # BF = 7.351649679664247
#=
               white pigmented
     white         *    0.2191
 pigmented    0.2191         *
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -6.49775
aic: 14.99549554466643
31×4 DataFrame
│ Row │ nodenumber │ nodelabel │ white       │ pigmented   │
│     │ Int64      │ String    │ Float64     │ Float64     │
├─────┼────────────┼───────────┼─────────────┼─────────────┤
│ 1   │ 1          │ Adi002    │ 1.0         │ 0.0         │
│ 2   │ 2          │ Adi003    │ 1.0         │ 0.0         │
│ 3   │ 3          │ Adi001    │ 1.0         │ 0.0         │
│ 4   │ 4          │ Asu001    │ 1.0         │ 0.0         │
│ 5   │ 5          │ Aga001    │ 1.0         │ 0.0         │
│ 6   │ 6          │ Aga002    │ 1.0         │ 0.0         │
│ 7   │ 7          │ Aza135    │ 0.0         │ 1.0         │
│ 8   │ 8          │ Aza037    │ 0.0         │ 1.0         │
│ 9   │ 9          │ Ama018    │ 0.0         │ 1.0         │
│ 10  │ 10         │ Ama006    │ 0.0         │ 1.0         │
│ 11  │ 11         │ Ape009    │ 0.0         │ 1.0         │
│ 12  │ 12         │ Ape001    │ 0.0         │ 1.0         │
│ 13  │ 13         │ Aru127    │ 0.0         │ 1.0         │
│ 14  │ 14         │ Aru001    │ 0.0         │ 1.0         │
│ 15  │ 15         │ Age001    │ 1.0         │ 0.0         │
│ 16  │ 16         │ 16        │ 0.917471    │ 0.0825293   │
│ 17  │ 17         │ 17        │ 0.917471    │ 0.0825293   │
│ 18  │ 18         │ 18        │ 0.757313    │ 0.242687    │
│ 19  │ 19         │ 19        │ 0.0944873   │ 0.905513    │
│ 20  │ 20         │ 20        │ 0.0108531   │ 0.989147    │
│ 21  │ 21         │ 21        │ 0.766186    │ 0.233814    │
│ 22  │ 22         │ #H16      │ 0.00781162  │ 0.992188    │
│ 23  │ 23         │ 23        │ 0.00781162  │ 0.992188    │
│ 24  │ 24         │ 24        │ 0.000188413 │ 0.999812    │
│ 25  │ 25         │ 25        │ 0.00014197  │ 0.999858    │
│ 26  │ 26         │ 26        │ 4.68619e-5  │ 0.999953    │
│ 27  │ 27         │ 27        │ 8.40095e-5  │ 0.999916    │
│ 28  │ 28         │ 28        │ 0.9715      │ 0.0284995   │
│ 29  │ 29         │ 29        │ 0.99624     │ 0.00375979  │
│ 30  │ 30         │ 30        │ 0.999424    │ 0.000576369 │
│ 31  │ 31         │ 31        │ 0.99972     │ 0.000280042 │
=#

# pigmentation: haplotype with outgroups
#-------------------------------------

calnet2o_pruned = deepcopy(calnet2o);
for n in ["Bce020","Pcr070"] deleteleaf!(calnet2o_pruned, n); end

# unconstrained rates:
fitc2o = fitDiscrete(calnet2o_pruned, mcol, datcol)
aic(fitc2o)
#=
rate white→pigmented α=0.17652
rate pigmented→white β=0.75842
1 traits, 16 species, on a network with 1 reticulations
log-likelihood: -6.20307
aic: 16.406144695367196
=#

# equal rates: fits better
fitc2o_equal = fitDiscrete(calnet2o_pruned, mcol_equal, datcol)
aic(fitc2o_equal)
show(ancestralStateReconstruction(fitc2o_equal), allrows=true)
plotasr(fitc2o_equal, "traitanalysis/figures/haplotype344_withOutgroups_ase_col.pdf",
        coltrait=["white","orange"], legpos="topleft", ylim=[0.8,16.2], xlim=[1,3.2])
exp.(fitc2o_equal.postltw)  # 0.9779840550872033, 0.02201594491279644
exp.(fitc2o_equal.priorltw) # 0.928069201055275, 0.07193079894472504
geneflowBF(fitc2o_equal)    # BF = 0.29044976673246325
#=
               white pigmented
     white         *    0.4192
 pigmented    0.4192         *
1 traits, 16 species, on a network with 1 reticulations
log-likelihood: -6.63987
aic: 15.279740501324724
33×4 DataFrame
│ Row │ nodenumber │ nodelabel │ white      │ pigmented  │
│     │ Int64      │ String    │ Float64    │ Float64    │
├─────┼────────────┼───────────┼────────────┼────────────┤
│ 1   │ 1          │ Smi165    │ 1.0        │ 0.0        │
│ 2   │ 2          │ Adi001    │ 1.0        │ 0.0        │
│ 3   │ 3          │ Adi002    │ 1.0        │ 0.0        │
│ 4   │ 4          │ Adi003    │ 1.0        │ 0.0        │
│ 5   │ 5          │ Age001    │ 1.0        │ 0.0        │
│ 6   │ 6          │ Asu001    │ 1.0        │ 0.0        │
│ 7   │ 7          │ Aga001    │ 1.0        │ 0.0        │
│ 8   │ 8          │ Aga002    │ 1.0        │ 0.0        │
│ 9   │ 9          │ Aza135    │ 0.0        │ 1.0        │
│ 10  │ 10         │ Aza037    │ 0.0        │ 1.0        │
│ 11  │ 11         │ Ama006    │ 0.0        │ 1.0        │
│ 12  │ 12         │ Ama018    │ 0.0        │ 1.0        │
│ 13  │ 13         │ Ape001    │ 0.0        │ 1.0        │
│ 14  │ 14         │ Ape009    │ 0.0        │ 1.0        │
│ 15  │ 15         │ Aru001    │ 0.0        │ 1.0        │
│ 16  │ 16         │ Aru127    │ 0.0        │ 1.0        │
│ 17  │ 17         │ 17        │ 0.700961   │ 0.299039   │
│ 18  │ 18         │ 18        │ 0.681549   │ 0.318451   │
│ 19  │ 19         │ 19        │ 0.22532    │ 0.77468    │
│ 20  │ 20         │ 20        │ 0.0131006  │ 0.986899   │
│ 21  │ 21         │ 21        │ 0.22532    │ 0.77468    │
│ 22  │ 22         │ 22        │ 0.00184694 │ 0.998153   │
│ 23  │ 23         │ 23        │ 2.08118e-5 │ 0.999979   │
│ 24  │ 24         │ 24        │ 9.64094e-5 │ 0.999904   │
│ 25  │ 25         │ 25        │ 1.95566e-5 │ 0.99998    │
│ 26  │ 26         │ 26        │ 9.56649e-5 │ 0.999904   │
│ 27  │ 27         │ 27        │ 0.225321   │ 0.774679   │
│ 28  │ 28         │ #H19      │ 0.998671   │ 0.00132889 │
│ 29  │ 29         │ 29        │ 0.998671   │ 0.00132887 │
│ 30  │ 30         │ 30        │ 0.998671   │ 0.00132885 │
│ 31  │ 31         │ 31        │ 0.998671   │ 0.00132885 │
│ 32  │ 32         │ 32        │ 0.960131   │ 0.0398691  │
│ 33  │ 33         │ 33        │ 0.996336   │ 0.00366424 │
=#

# pigmentation: haplotype ingroup only
#-------------------------------------

# unconstrained rates:
fitc2i = fitDiscrete(calnet2i, mcol, datcol)
aic(fitc2i)
#=
rate white→pigmented α=0.21776
rate pigmented→white β=0.79936
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -6.27119
aic: 16.542380322802295
=#

# equal rates: fits better
fitc2i_equal = fitDiscrete(calnet2i, mcol_equal, datcol)
aic(fitc2i_equal)
show(ancestralStateReconstruction(fitc2i_equal), allrows=true)
plotasr(fitc2i_equal, "traitanalysis/figures/haplotype344_onlyIngroup_ase_col.pdf",
        coltrait=["white","orange"], ylim=[0.8,15.2], xlim=[1,3.2])
exp.(fitc2i_equal.postltw)  # 0.5942097074178719, 0.40579029258212806
exp.(fitc2i_equal.priorltw) # 0.8723075169958534, 0.12769248300414657
geneflowBF(fitc2i_equal)    # BF = 4.665156248910101
#=
               white pigmented
     white         *    0.3612
 pigmented    0.3612         *
1 traits, 15 species, on a network with 1 reticulations
log-likelihood: -6.47541
aic: 14.950810664666887
│ Row │ nodenumber │ nodelabel │ white      │ pigmented  │
│     │ Int64      │ String    │ Float64    │ Float64    │
├─────┼────────────┼───────────┼────────────┼────────────┤
│ 1   │ 1          │ Adi001    │ 1.0        │ 0.0        │
│ 2   │ 2          │ Adi002    │ 1.0        │ 0.0        │
│ 3   │ 3          │ Adi003    │ 1.0        │ 0.0        │
│ 4   │ 4          │ Age001    │ 1.0        │ 0.0        │
│ 5   │ 5          │ Asu001    │ 1.0        │ 0.0        │
│ 6   │ 6          │ Aga001    │ 1.0        │ 0.0        │
│ 7   │ 7          │ Aga002    │ 1.0        │ 0.0        │
│ 8   │ 8          │ Aza135    │ 0.0        │ 1.0        │
│ 9   │ 9          │ Aza037    │ 0.0        │ 1.0        │
│ 10  │ 10         │ Ama018    │ 0.0        │ 1.0        │
│ 11  │ 11         │ Ama006    │ 0.0        │ 1.0        │
│ 12  │ 12         │ Ape001    │ 0.0        │ 1.0        │
│ 13  │ 13         │ Ape009    │ 0.0        │ 1.0        │
│ 14  │ 14         │ Aru127    │ 0.0        │ 1.0        │
│ 15  │ 15         │ Aru001    │ 0.0        │ 1.0        │
│ 16  │ 16         │ 16        │ 0.862235   │ 0.137765   │
│ 17  │ 17         │ 17        │ 0.448113   │ 0.551887   │
│ 18  │ 18         │ 18        │ 0.034968   │ 0.965032   │
│ 19  │ 19         │ 19        │ 0.00792135 │ 0.992079   │
│ 20  │ 20         │ 20        │ 0.469911   │ 0.530089   │
│ 21  │ 21         │ #H16      │ 0.00151575 │ 0.998484   │
│ 22  │ 22         │ 22        │ 0.00151575 │ 0.998484   │
│ 23  │ 23         │ 23        │ 1.37262e-5 │ 0.999986   │
│ 24  │ 24         │ 24        │ 6.89614e-5 │ 0.999931   │
│ 25  │ 25         │ 25        │ 1.27793e-5 │ 0.999987   │
│ 26  │ 26         │ 26        │ 6.84331e-5 │ 0.999932   │
│ 27  │ 27         │ 27        │ 0.962292   │ 0.0377079  │
│ 28  │ 28         │ 28        │ 0.996884   │ 0.00311646 │
│ 29  │ 29         │ 29        │ 0.862235   │ 0.137765   │
│ 30  │ 30         │ 30        │ 0.994735   │ 0.00526501 │
│ 31  │ 31         │ 31        │ 0.996235   │ 0.00376515 │
=#

#--------------------------------------------------
#     summary table
#--------------------------------------------------

allfit = [fitp1o, fitp1o_equal, fitp1i, fitp1i_equal,
          fitp2o, fitp2o_equal, fitp2i, fitp2i_equal,
          fitc1o, fitc1o_equal, fitc1i, fitc1i_equal,
          fitc2o, fitc2o_equal, fitc2i, fitc2i_equal]
df = DataFrame(
    trait = repeat(["pollinator","flower_color"], inner=8),
    network_data = repeat(["primary","haplotype"], inner=4, outer=2),
    network_sampling = repeat(["with_outgroups","only_ingroup"], inner=2, outer=4),
    rates = repeat(["unconstrained","equal"], inner=1, outer=8),
    likelihood = [fit.loglik for fit in allfit],
    aic = [aic(fit) for fit in allfit],
    bf = [geneflowBF(fit) for fit in allfit],
    rate1 = [fit.model.rate[1] for fit in allfit],
    rate2 = map(fit -> length(fit.model.rate)>1 ? fit.model.rate[2] : missing , allfit)
)
df |> CSV.write("traitanalysis/modelsummary.csv")


#---------------------------------
#      figure for manuscript
#---------------------------------

# main:
R"pdf"("traitanalysis/figures/Figure7.ase.pdf", width=7, height=7)
R"layout"([1 1 1 1 2 2 2; 3 3 3 3 4 4 4]);
plotasr(fitc1o_equal, ""; onfile=false, legend=false, ylim=[0.8,16.2], xlim=[1,3.5], coltrait=["white","orange"], radius=0.025)
R"legend"(x=1.2, y=16, legend=["white", "pigmented"], pch=21, var"pt.bg"=["white","orange"],
        bty="n", title="flower color", var"title.adj"=0, var"pt.cex"=1.5);
plotasr(fitc1i_equal, ""; onfile=false, legend=false, ylim=[0.8,15.2], xlim=[1,3.0], coltrait=["white","orange"], radius=0.025)
plotasr(fitp1o_equal, ""; onfile=false, legendtitle="pollinator", radius=0.020)
plotasr(fitp1i_equal, ""; onfile=false, legend=false, radius=0.025) # radius 0.016 by default
R"dev.off"()

# supplemental:
R"pdf"("traitanalysis/figures/FigureS7.ase.pdf", width=7, height=7)
R"layout"([1 1 1 1 2 2 2; 3 3 3 3 4 4 4]);
plotasr(fitc2o_equal, ""; onfile=false, legend=false, ylim=[0.8,16.2], xlim=[0.9,3.2], coltrait=["white","orange"], radius=0.024)
R"legend"(x=0.95, y=16, legend=["white", "pigmented"], pch=21, var"pt.bg"=["white","orange"],
        bty="n", title="flower color", var"title.adj"=0, var"pt.cex"=1.5);
plotasr(fitc2i_equal, ""; onfile=false, legend=false, ylim=[0.8,15.2], xlim=[1,3.0], coltrait=["white","orange"], radius=0.026)
plotasr(fitp2o_equal, ""; onfile=false, legendtitle="pollinator", radius=0.021)
plotasr(fitp2i_equal, ""; onfile=false, legend=false, radius=0.026)
R"dev.off"()

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
