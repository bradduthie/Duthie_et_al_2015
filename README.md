Trade-offs and coexistence in fluctuating environments
=======================================================

The files in this repository are also available in the Dryad Digital Repository:

Duthie AB, Abbott KC, Nason JD (2015) Data from: Tradeoffs and coexistence in fluctuating environments: evidence for a key dispersal-fecundity tradeoff in five nonpollinating fig wasps. Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.4dj10

Overview
-----------------------------------------

The subsequent code replicates the analyses used in Duthie et al. (2015) *American Naturalist*. All data were collected in 2010 in the field from populations of Sonoran Desert rock fig (*Ficus petiolaris*) trees. For this analysis, three files are needed:

The first file includes counts of all female non-pollinating fig wasps from the syconia of the Sonoran Desert rock fig, *Ficus petiolaris* ('wasp_data.csv'). The data file includes the site, tree and syconia ('Fruit') labels from each syconia sampled (rows); the data also includes columns indicating tree lat-lon coordinates, pollinating foundresses arriving to the syconia (corpses remain within syconia), fruit volume, counts of all wasp species, and the number of neighbours a tree has within a 1 km radius. The pollinator ('Poll') is an unnamed (as of 2015) species of *Pegoscapus*. Wasps 'LO1', 'SO1', and 'SO2' are all unnamed species of the genus *Idarnes*, and wasps 'Het1' and 'Het2' are both unnamed species of the genus *Heterandrium*. The wasps 'Physo' and 'Aepoc' are wasps of unnamed species of the genera *Physothorax* and *Aepocerus*, respectively; 'NA' denotes missing values.

The second file ('wing_loadings.csv') includes measurements to estimate wing loading for 83 wasps (rows) of the species of interest in Duthie et al. (2015): LO1, SO1, SO2, Het1, and Het2. Each row is a single wasps, and columns show the site, tree, and fruit from which the wasp was sampled. Lengths and widths of wasp heads, thoraxes, and abdomens are included as measured at their widest points (e.g., for the wasp abdomen, the whole length of the segment is reported, along with an estimate of its width at the widest point). Wing areas were calculated using wing images and ImageJ software.

The third file ('egg_loads.csv') includes estimates of wasp egg loads from 54 wasps (rows) of the following species:  LO1, SO1, SO2, Het1, and Het2. Both mature and immature egg counts were estimated, and columns include the site, tree, and fruit from which the wasp was sampled.

The code below will replicate the statistical analysis of Duthie et al. (2015). This analysis includes: 1) The estimation of a colonization index for each species, and the correlation of this index with species mean egg load. 2) The correlation between mean species wing loading and mean species fecundity. 3) The correlation betwen species wing loading and the effect (regression slope) that the number of neighbouring fig trees within a 1 km radius of the tree from which wasps were sampled had on the species abundance; the same correlation is estimated with egg load instead of wing loading. Finally, 4) the differences between species fecundities/wing loadings as correlated with species among-tree density correlations (i.e., individual points are absolute differences in fecundities or wing loadings, and are compared with individual points that are correlations between the densities of species on a fig tree).

After the analyses are replicated, the remaining code re-builds the figures presented in Duthie et al., 2015. 

Any enquiries about these data or the analysis and code that follows can be made to Brad Duthie (aduthie@abdn.ac.uk; brad.duthie@gmail.com).

Duthie, A. B., Abbott, K. C., & Nason, J. D. (2015). Trade-offs and coexistence in fluctuating environments: evidence for a key dispersal-fecundity trade-off in five nonpollinating fig wasps. American Naturalist, 186(1), 151–158. doi:10.1086/681621
