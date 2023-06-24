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

`mt-corr.R` and `pt-corr.R` use the PGLS pipeline above to explore correlations between different traits and organelle DNA counts across species. This is computationally demanding, given the large correlation structure linking observations. A machine with 4GB RAM (or maybe more) is probably necessary; some of the traits will take perhaps dozens of minutes on a modern machine. As there are around 80 traits, some of which have multiple factor levels, we are probably talking a 48-hour run time for the mtDNA data (the ptDNA data is smaller and sparser and will probably run in a handful of hours). EDIT: there is now a final "habitat" feature that has many dozen levels. This will take a possible couple of days by itself! It's all very parallelisable though if you have the memory.

`plot-odna-corr.R` plots the results (which are written to files by the above scripts).

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
 
