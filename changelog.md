
## V2.0.0 (xx/12/2026) - Julien Straubhaar & Celia Trunz
- PACKAGE REDESIGN 
- New modules / packages:
	- io: import / export functions, various format are handled (previous import export functions still available)
	- clean: functions used to clean the graph before instantiating the KGraph object
	- view: functions for plotting / viewing graphs (in 2D and 3D, based on networkx / matplotlib / pyvista) (main argument is a Networkx graph for more generality - can be used with kg.graph for a KGraph object)
	- tools: various tools and algorithm for graphs
	- utils: basic utilities
- New notebooks: numerous notebooks has been added for illustrating different aspects of karstnet

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
