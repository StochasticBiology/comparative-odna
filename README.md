# comparative-odna
Phylogenetic comparative approaches with challenging (organelle DNA) datasets

Permanent link to post-peer-review release: [![DOI](https://zenodo.org/badge/657160573.svg)](https://zenodo.org/doi/10.5281/zenodo.10693678)

This project explores the use of different phylogenetic comparative methods, including phylogenetic generalised least squares (PGLS) and phylogenetic generalised linear models (PGLM), as well as approaches accounting for relatedness in a more coarse-grained way, to investigate evolutionary correlations between variables with awkward data structures.

![image](https://github.com/StochasticBiology/pgls-odna/assets/50171196/6dca8ca1-e609-4226-ac7b-c8a0acc95451)

Synthetic controls
-----
All these scripts simulate an evolutionary process meant to mimic organelle DNA reduction, and explore the power and precision of comparative methods to detect correlations from awkward resulting data.

`method-comparison.R` compares PGLS, naive correlation, PGLM, and other non-parametric approaches to characterise correlations.
`pgls-high-obs-noise.R` explores the PGLS effect of a substantial false negative rate in observations of the predictor variable.
`pgls-branch-lengths.R` explores the PGLS and PGLM effect of accounting for, or not accounting for, branch lengths.
`pgls-branch-cor-structs.R` does as the previous but tests different models for correlation structures.
`pgls-illustration.R` illustrates some simulations.

Organelle DNA -- cross-eukaryotic data, taxonomy tree
-----
`MTFull22.txt` and `PTFull22.txt` are tab-separated datafiles of organismal traits by species, compiled from a variety of sources (see Sources below). These datafiles are very sparse -- one of the issues we are investigating. `tree-for-traits-clean-mt.phy` and `tree-for-traits-clean-pt.phy` are Newick trees from NCBI's Common Taxonomy Tree tool linking these observed species. `odna-illustration.R` uses these files to plot trees illustrating a particular trait (parasitism) with oDNA gene count.

`[mt/pt]-corr.R` use the PGLS pipeline above to explore correlations between different traits and organelle DNA counts across species. This is computationally demanding, given the large correlation structure linking observations. A machine with 4GB RAM (or maybe more) is probably necessary; some of the traits will take perhaps dozens of minutes on a modern machine. As there are around 80 traits, some of which have multiple factor levels, we are probably talking a 48-hour run time for the mtDNA data* (the ptDNA data is smaller and sparser and will probably run in a handful of hours). There is a final "habitat" feature that has many dozen levels. This will take a possible couple of days by itself! It's all very parallelisable though if you have the memory. 
`[mt/pt]-test.R` do PLM and PGLM for oDNA protein-coding gene counts, also producing plots according to latest significance coding; `[mt/pt]-test-ncbi.R` do the same for NCBI-derived CDS counts; `mt-test-[ncbi-]no-metazoa.R` for the MT data with metazoans removed (many datapoints, little diversity). These approaches are much faster than the PGLS above and should run in seconds on a modern machines. Summary plots are produced.

`[mt/pt]-test-blockonly.R` use Kruskal-Wallis or Scheirer-Ray-Hare tests to block eukaryotic clade and analyse the remaining variance due to a factor of interest. `[mt/pt]-test-mixed.R` use LMM and GLMM approaches to assign random effects associated with eukaryotic clade and assess the remaining link with the feature of interest. `[mt/pt]-test-clademean.R` subtract each clade's mean gene count from all members of that clade before running PLM.

The separate script `plot-odna-corrs.R` runs most of these analyses then structures summary plots for manuscript layout.

Organelle DNA -- plant kingdom data, estimated phylogeny
-----
We provide scripts to process an external (vascular) plant phylogeny, output the necessary txt files and then plot the PGLM and PLM volcano plots as in `(mt/pt)-test(-ncbi).R`, once with the branch lengths and once with uniform branch lengths (neglecting specific branch length information; labelled NOBL).

`plantPhyloMT.R` is the script that fetches the phylogeny data from https://github.com/megatrees and tries to find the matches with our `MTFull22.txt` dataset. Then, it exports two files: `MTspeciesHits.txt` with the subset of `MTFull22.txt` and `mt-pruned-plantPhylo.phy`, which is the the phylogenetic tree for the former list of species. `plantPhyloPT.R` works the same for the plastid case, that is the `PTFull22.txt` input dataset.

`mt-test-Plants.R` reads the output files from plantPhyloMT.R and plots the significance and coefficients for (filtered) traits in a similar fashion to `mt-test.R` for PLM and PGLM. `mt-test-Plants-NOBL.R` does to same, only with  removing the branch length information from the plant phylogeny (i.e., the underlying tree has uniform edge lengths, but the topology remains the same). `pt-test-Plants.R` and `pt-test-Plants-NOBL.R` do the same for plastid gene counts.

Assisted manual parsing of Wikipedia pages
-----

`wiki-html.py` searches an XML dump (for example, from Wikipedia's Special:Export https://en.wikipedia.org/wiki/Special:Export ) for a given regular expression, and produces an HTML file designed to assist manual parsing of all the entries matching that pattern. The HTML page consists of the snippets of text surrounding each pattern match, labelled by page name. Links and other interesting text are accompanied by checkboxes, which when clicked populate a text box on the right-hand side of the page with the text immediately before the checkbox. Of course, text can also be manually entered into these boxes.

At the bottom of the HTML page there is a summary button which compiles all the text from other boxes into a single list.

The idea is that the user can quickly check any examples where a given bit of text (e.g. a species or taxon name) genuinely corresponds to the feature of interest, then compile all these positive cases into a summary list.

Here's an example. `parasitaxus.xml` is the XML output from a Special:Export search for _"parasitaxus"_. We can search for occurrences of the string _"arasit"_ within this XML with

`python3 wiki-html.py arasit parasitaxus.xml form-1.html`

`firefox form-1.html`

In the form we'll see lots of matches. One key one is

_"It is generally mentioned that ''Parasitaxus usta'' is the only known ''parasitic'' ''gymnosperm''. The species remarkably lacks ''root''s and is always found attached to the roots of ''''Falcatifolium taxoides'''' (another member of"_

This seems to provide a compelling suggestion that _Parasitaxus usta_ is a parasite. We thus click the checkbox next to _''Parasitaxus usta''_, which now appears in the right-hand text box.

After scrolling through the rest, we come to _"Produce summary"_ at the bottom, and clicking this gives us a list of all the text snippets we've checked.

The script supports regexs but you'll need to quote them in the command line call e.g.
 
`python3 wiki-html.py "[Uu]nicell|[Uu]ni-cell|[Ss]ingle-cell" mt-wiki.xml form-2.html`

Sources
----

Please see bibliography in https://www.biorxiv.org/content/10.1101/2023.08.11.553003v1.full

Species names and genetic information (such as organellar genome size, GC content, and taxonomical information) were collected via NCBI’s Organelle database (O’Leary et al. 2016) and processed according to the pipeline in (Giannakis et al. 2022). After curation, we retained 9296 mtDNA sequences (8835 of which are metazoans) and 4264 ptDNA sequences. The majority of the ecological trait entries were acquired via Encyclopedia of Life (Parr et al. 2014) through individual queries for traits of interest, and then cross-referenced with some domain-specific databases (like GloBI (Poelen, Simons, and Mungall 2014), but also others, see following references). Encyclopedia of Life is an open encyclopedia that provides a catalogue with plenty of traits and taxonomic information on organisms, which consists of contributions and curation by experts and organisations. Palm data were collected through (Kissling et al. 2019). Organismal habitats and traits related to cross-species interactions (including parasitism) were taken from (Cohen et al. 2020), a large ecological scale survey on species interactions. For the majority of plant data, ecological traits and characters were taken from Encyclopedia of Life and crosschecked with the TRY database (2018 update) (Kattge et al. 2020).

For the cross-eukaryote data, NCBI’s Common Taxonomy Tool was used to estimate phylogenetic topology (Federhen 2012). More specific phylogenetic information was obtained through the U.PhyloMaker R package (Jin and Qian 2023; 2022) which uses the Plant megaphylogeny found in (Jin and Qian 2022) comprised by phylogenies from (S. A. Smith and Brown 2018) and (Zanne et al. 2014). The plant phylogeny assembled in (Jin and Qian 2022) and sourced from (S. A. Smith and Brown 2018; Zanne et al. 2014) gave us a sufficient number of matches (up to the level of genus) to our list of species (120 for the species with mitochondrial gene count and 3851 for plastid) and we were able to test different scenaria.
 
