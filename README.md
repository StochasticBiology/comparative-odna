# pgls-odna
Phylogenetic comparative approaches with challenging (oDNA) datasets

This project explores the use of different phylogenetic comparative methods, including phylogenetic generalised least squares (PGLS), to investigate evolutionary correlations between variables with awkward data structures.

Synthetic controls
-----
All these scripts simulate an evolutionary process meant to mimic organelle DNA reduction, and explore the power and precision of comparative methods to detect correlations from awkward resulting data.

`method-comparison.R` compares PGLS, naive correlation, and other non-parametric approaches to characterise correlations.
`pgls-high-obs-noise.R` explores the PGLS effect of a substantial false negative rate in observations of the predictor variable.
`pgls-branches.R` explores the PGLS effect of accounting for, or not accounting for, branch lengths.
`pgls-cor-structs.R` and `pgls-branches-cor-structs.R` do as the previous two but test different models for correlation structures.

Organelle DNA
-----
`MTFull22.txt` and `PTFull22.txt` are tab-separated datafiles of organismal traits by species, compiled from a variety of sources. These datafiles are very sparse -- one of the issues we are investigating. `tree-for-traits-clean-mt.phy` and `tree-for-traits-clean-pt.phy` are Newick trees from NCBI's Common Taxonomy Tree tool linking these observed species.

`mt-corr.R` and `pt-corr.R` use the PGLS pipeline above to explore correlations between different traits and organelle DNA counts across species. This is computationally demanding, given the large correlation structure linking observations. A machine with 4GB RAM (or maybe more) is probably necessary; some of the traits will take perhaps dozens of minutes on a modern machine. As there are around 80 traits we are probably talking a 24-hour run time.

`plot-odna-corr.R` plots the results (which are written to files by the above scripts).
