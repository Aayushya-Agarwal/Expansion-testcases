# Expansion-testcases
The following repository contains supplementary material for [1] including scripts using KNITRO [2] and Gurobi [3] to perform expansion as well as data for two cases of expansion. The two cases for expansion that are considered are:
1) Reconstruction of a 14 bus testcase in which each transmission line and generator is considered for expansion. This testcase is labeled "case_14_expansion.m"
2) Expansion of a synthetic Eastern Interconnection modified from [4] that decommissions certain generators as well as increases the loading factor in a certain area. This testcase is labeled, "ACTIVSg70k_expansion.raw"

The repository contains two scripts to consider the expansion of the 14 bus testcase. The first script, labeled "knitro_14_bus_example.m" considers AC constraints using KNITRO while looking at the expansion and operation of the 14 bus testcase. The second script, labeled "dc_expansion_14_bus.m" considers DC constraints using Gurobi. Both of these scripts are run using Matlab.

[1]

[2][3]	R. H. Byrd, J. Nocedal, and R.A. Waltz, “KNITRO: An integrated package for nonlinear optimization”, In G. di Pillo and M. Roma, editors, Large-Scale Nonlinear Optimization, pages 35–59. Springer, 2006.

[3]	Gurobi Optimization, Incorporate. "Gurobi optimizer reference manual." (2018).

[4][23]	 Xu; A. B. Birchfield; K. S. Shetye; T. J. Overbye, “Creation of Synthetic Electric Grid Models for Transient Stability Studies,”  2017 IREP Symposium Bulk Power System Dynamics and Control, Espinho, Portugal, 2017.
