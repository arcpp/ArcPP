# ArcPP - Archaeal Proteome Project

Modern proteomics approaches can explore whole proteomes within a single mass spectrometry (MS) run. However, the enormous amount of MS data generated often remains incompletely analyzed due to a lack of sophisticated bioinformatic tools and expertise needed from a diverse array of fields. In particular, in the field of microbiology, efforts to aggregate and integrate proteomic datasets have so far been missing. Thus, despite their relatively small genomes, the proteomes of most archaea remain incompletely characterized. This in turn undermines our ability to gain a greater understanding of archaeal cell biology. 

Therefore, we have initiated the Archaeal Proteome Project (ArcPP), a community effort that works towards a comprehensive analysis of archaeal proteomes. Starting with the model archaeon Haloferax volcanii, using state-of-the-art bioinformatic tools, we have:
* reanalyzed more than 23 Mio. spectra
* optimized the analysis using parameter sweeps, multiple search engines implemented in Ursgal, and the combination of resultsthrough the combined PEP approach
* thoroughly controlled false discovery rates for high confidence protein identifications using the picked protein FDR approach and limiting FDR to 0.5%
* identified more than 40k peptides, corresponding to over 3000 proteins (>75% of the proteome) with a median sequence coverage of 50%

Benefiting from the established bioinformatic infrastructure, we will follow up on this analysis focusing on H. volcanii proteogenomics as well as the characterization of various post-translational modifications. Furthermore, ArcPP will integrate quantitative results obtained from the individual datasets in order to identify common regulatory mechanisms. These studies on the H. volcanii proteome can serve as a blueprint for comprehensive proteomic analyses performed on a diverse range of archaea and bacteria.


## Citation

Coming soon.

## Installation

Get the latest version via GitHub: https://github.com/StSchulze/ArcPP

as a .zip package: https://github.com/StSchulze/ArcPP/archive/master.zip

or via git clone::

 	git clone https://github.com/StSchulze/ArcPP.git

For looking at the meta data, database and result files, of course no specific software is required.
However, if you want to use the Dash app to search for proteins,
you will need Python (tested for Python 3.5 and higher) with the packages Ursgal, Pandas and Dash.
Download Python here: https://www.python.org/downloads/
It comes with pip, which you can use to install the dependencies::

	pip install -r requirements.txt

## How to use

To start the Dash app, simply execute the script Arcpp_results_dash_app.py in the main folder::

	python Arcpp_results_dash_app.py

Dash automatically starts a local server and you can view the results opening the following page in you browser:
http://127.0.0.1:8050/

## Documentation

As the project develops, a more extensive documentation via readthedocs will be included.

## Questions and Participation

If you encounter any problems you can open up issues at GitHub, or contact us directly by email.
For any contributions, fork us at https://github.com/StSchulze/ArcPP. and open up pull requests.

## Copyrights

Copyright 2019-today by authors and contributors in alphabetical order

* Zach Adams
* Rosana de Castro
* Micaela Cerletti
* Sébastien Ferreira-Cerca
* Michael Hippler
* Zivojin Jevtic
* Robert Knüppel
* Christof Lenz
* Anita Marchfelder
* Julie Maupin-Furlow
* Friedhelm Pfeiffer
* Mechthild Pohlschroder
* Ansgar Poetsch
* Stefan Schulze
* Henning Urlaub

## Contact

Dr. Stefan Schulze
Department of Biology, University of Pennsylvania
201 Leidy Laboratories, 3740 Hamilton Walk, Philadelphia, PA 19104
email: sschulze@sas.upenn.edu

## Disclaimer

While all analyses for ArcPP have been performed with greatest care, errors in the underlying software, raw datasets or database cannot be ruled out completely. Therefore, the source code used for this project is available here and should be inspected if in doubt or used to verify results. As common practice in science: be critical and never trust results blindly.

## License

The content of this project itself is licensed under the [Creative Commons Attribution 3.0 Unported license](https://creativecommons.org/licenses/by/3.0/), and the underlying source code used to format and display that content is licensed under the [GNU Lesser General Public License v3.0](https://github.com/StSchulze/ArcPP/blob/master/LICENSE).
