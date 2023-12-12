# ArcPP - Archaeal Proteome Project

Modern proteomics approaches can explore whole proteomes within a single mass spectrometry (MS) run. However, the enormous amount of MS data generated often remains incompletely analyzed due to a lack of sophisticated bioinformatic tools and expertise needed from a diverse array of fields. In particular, in the field of microbiology, efforts to combine large-scale proteomic datasets have so far largely been missing. Thus, despite their relatively small genomes, the proteomes of most archaea remain incompletely characterized. This in turn undermines our ability to gain a greater understanding of archaeal cell biology. 

Therefore, we have initiated the Archaeal Proteome Project (ArcPP), a community effort that works towards a comprehensive analysis of archaeal proteomes. Starting with the model archaeon Haloferax volcanii, using state-of-the-art bioinformatic tools, we have:
* reanalyzed 26 Mio. spectra
* optimized the analysis using parameter sweeps, multiple search engines implemented in Ursgal, and the combination of resultsthrough the combined PEP approach
* thoroughly controlled false discovery rates for high confidence protein identifications using the picked protein FDR approach and limiting FDR to 0.5%
* identified more than 45k peptides, corresponding to 3069 proteins (>75% of the proteome) with a median sequence coverage of 55%.
* analyzed N-terminal protein processing, including N-terminal acetylation and signal peptide cleavage
* performed a detailed glycoproteomic analysis, identifying >230 glycopeptides corresponding to 45 glycoproteins

Benefiting from the established bioinformatic infrastructure, we will follow up on this analysis focusing on H. volcanii proteogenomics as well as the characterization of various post-translational modifications. Furthermore, ArcPP will integrate quantitative results obtained from the individual datasets in order to identify common regulatory mechanisms. These studies on the H. volcanii proteome can serve as a blueprint for comprehensive proteomic analyses performed on a diverse range of archaea and bacteria.


## Citation

***Initial publication:***
Schulze, S.; Adams, Z.; Cerletti, M. et al. (2020). The Archaeal Proteome Project advances knowledge about archaeal cell biology through comprehensive proteomics. Nat Commun 11, 3145. https://doi.org/10.1038/s41467-020-16784-7

**Glycoproteomic analysis using ArcPP:**
Schulze, S.; Pfeiffer, F.; Garcia, B.A.; Pohlschroder, M. (2021). Comprehensive glycoproteomics shines new light on the complexity and extent of glycosylation in archaea. PLOS Biol.  https://doi.org/10.1371/journal.pbio.3001277

**Step-by-step description of sample preparation and bioinformatic analysis:**
Schulze, S.; Pohlschroder, M. (2022). Proteomic sample preparation and data analysis in line with the Archaeal Proteome. Meth Mol Biol. https://doi.org/10.1007/978-1-0716-2445-6_18


## Download

Get the latest version via GitHub: https://github.com/arcPP/ArcPP

as a .zip package: https://github.com/arcPP/ArcPP/archive/master.zip

or via git clone::

    git clone https://github.com/ArcPP/ArcPP.git

Result files can found here: https://doi.org/10.5281/zenodo.3825856.

An interactive, searchable representation of protein and peptide identifications can be found at: https://archaealproteomeproject.org 


## How to use

For looking at the meta data, database and result files, of course no specific software is required.

However, if you want to execute the analysis scripts, Python 3.6 or higher is required as well as the
Python package Ursgal. Further information on the installation and usage of Ursgal can be found here:
https://github.com/ursgal/ursgal

Furthermore, docstrings in the scripts explain the basic functions and steps in the analysis pipeline.


## Questions and Participation

If you encounter any problems you can open up issues at GitHub, or contact us directly by email.
For any contributions, fork us at https://github.com/arcPP/ArcPP and open up pull requests.


## Changelog

### Version 1.4.0 (12-2023)
* Full integration of PXD040781

### Version 1.3.0 (06-2021)
* Full integration of PXD021827

### Version 1.2.0 (05-2021)
* Integration of PXD021827 (meta data)
* Integration of PXD021874 (meta data + results)
* Including glycoproteomic analysis scripts
* Including general_example_workflow.py as a generalized example analysis script

### Version 1.1.0 (05-2020)
* Including analysis of Natrialba magadii as part of PXD009116

### Version 1.0.0 (05-2020)
* Initial ArcPP release

## Copyrights

Copyright 2019-today by authors and contributors in alphabetical order

* Zach Adams
* Micaela Cerletti
* Rosana De Castro
* Sébastien Ferreira-Cerca
* Christian Fufezan
* Benjamin A. Garcia
* María Inés Giménez
* Michael Hippler
* Zivojin Jevtic
* Robert Knüppel
* Georgio Legerme
* Christof Lenz
* Anita Marchfelder
* Julie Maupin-Furlow
* Roberto A. Paggi
* Friedhelm Pfeiffer
* Ansgar Poetsch
* Mechthild Pohlschroder
* Stefan Schulze
* Henning Urlaub

## Contact

Dr. Stefan Schulze
Gosnell School of Life Sciences
Rochester Institute of Technology 
Brown Hall 86-1121, Rochester, NY 14623
email: stefan.schulze@rit.edu

## Disclaimer

While all analyses for ArcPP have been performed with greatest care, errors in the underlying software, raw datasets or database cannot be ruled out completely. Therefore, the source code used for this project is available here and should be inspected if in doubt, or used to verify results. As common practice in science: be critical and never trust results blindly.

## License

The content of this project itself is licensed under the [Creative Commons Attribution 3.0 Unported license](https://creativecommons.org/licenses/by/3.0/), and the underlying source code used to format and display that content is licensed under the [GNU Lesser General Public License v3.0](https://github.com/arcPP/ArcPP/blob/master/LICENSE).
