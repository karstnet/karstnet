## V1.3.0 (...) - Julien Straubhaar & Celia Trunz
- Extract functions that apply to Kg.graph only to make the functions usable outside of the Kg object. While keeping it in the method
- Add new packages:
	- io: import export and functions (previous import export functions used as methods exist now as functions under different names, but are still callable as the original method)
	- clean: functions used to clean the graph before initating the kg object. 
	- view: function for plotting or other viewing tools, that are applied to the graph only. Previously existing methods are transformed in functions with a different name, but methods is kept by calling the function in the 

## V1.2.5 (30/08/2024) - Philippe Renard

- Modifying package structure for distribution via pypi
- Replacing setup.py by pyproject.toml

## V1.2.4 (23/07/2024) - Philippe Renard

- Updated the automated construction of the documentation
- Minor bug fix in the computation of the cv length

## V1.2.1 (19/10/2023) - Philippe Renard

- Added verbosity option in creator functions
- Corrected import from Therion
- Added changelog.md

## V1.2.0 (18/05/2021) - Pauline Collon

- Modification of entropy computation

## V1.1.0 (25/11/2019) - Philippe Vernant & Pauline Collon

- Adding stereonet features

## V1.0.0 (2019) - Pauline Collon & Philippe Renard

- First complete release of Karsnet

## V0.1.0 (10/12/2018) - Philippe Renard

- Pre-release of Karstnet
