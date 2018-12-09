te

# Karstnet TODO list

This document intents to register and organize what we still have to do.



## Objectives

- [ ] **version 1.0** : allows computing the statistics of a karst network
  - [ ] reasonably debugged code
    - [ ] to solve : remove warning with histograms for aspl
    - [ ] check entropy with Pauline
    - [ ] test cases with several connex components (is it already done?)
  - [x] basic documentation
    - [x] seting up the sphinx infrastructure
    - [x] an intro page with link to paper and philosophy of the package
    - [x] add example : load data / compute characteristics
    - [x] one or two simple notebooks
  - [ ] unitest
    - [x] tests on simple cases 
    - [ ] test example with real data set
    - [ ] test for import and export functions
  - [x] no data base in that version
- [ ] **version 1.1** : allows in addition to test if a simulated network is plausible
  - [ ] includes  
    - [ ] clean data base
    - [ ] function to test a simulated network
  - [ ] supported by a minimum of additional testing of the methodology
  - [ ] this is the version required for the paper
- [ ] **version 1.2** : nicer version with additional features, such as
  - [ ] functions to add data to database
  - [ ] function to export in additional formats (Shapefile, etc.)
  - [ ] additional metrics (Hendricks, Jouves ?)
  - [ ] stereonets ?
  - [ ] etc.



## CODING


### Priority

- [ ] distribution clean
  - [x] structure pour install avec pip (en cours)
  - [ ] readme / makefile and co
- [ ] unit tests
  - [x] coding simple unitary tests: see test_karstnet.py
  - [ ] coding unit test for a pair of larger open source networks
  - [ ] review the code: list functions of karstnet and check that they are tested
  - [ ] chose public data set to be distributed with the tests / code corresponding tests
- [ ] test graph simplification (not clear to Philippe - should be discussed)

  - [ ] in case of connex components 
  - [ ] loop combined with a cycle !
- [ ] code function to test if a network is plausible
- [ ] database maintenance and extension 

  - [ ] recompute the complete set of indicators and rebuild a table with the new code
  - [ ] function to add a network to the database

### Non priority

- [ ] function to damage the graphs in order to test the robustness of the statistics:
   * random point of degree 2 suppression
   * suppression of one part of the network : suppress one conex component
- [ ] additional characteristics that could be added
   - [ ] fractal dimension (see Martin's code)
   - [ ] other metrics from Jouves ?
- [ ] ensure that plot2 and plot3 are not creating new figures / this is to allow subplot
- [ ] generate the stereonets
- [ ] enlever pos3D et ne garder que z pour eviter les doublons
- [ ] metadata management thanks to YML ?
- [ ] export en shapefile (et dxf?)
- [ ] changer les plots pour mettre par défaut les 2 et choisir l'un et/ou l'autre avec un parametre
- [ ] Pour la pondération par les longueurs, peut-etre que ce serait plus intéressant 
  pour comparer les réseau de diviser par la longueur moyenne d'un edge ?? 
  Je ne sais pas je n'ai pas encore vu de biblio la dessus... ??
  Philippe -> I do not think that it is useful for the moment.



## Documentation



- [x] prepare the setup for a simple documentation including a basic tutorial using Sphinx
- [ ] add reference to original paper Colon et al
- [ ] list the networks, the input format, and initial provider
- [ ] add text / decide if GitHub or documentation to thank the data providers 
  - [ ] P Le Roux pour Foussoubie



## Numerical experiments



- [ ] Evaluate if the statistical tests developed with David are efficient

- [x] Add all networks - complete data base

-

## Paper



- [ ] Write errata on correlation des vertex degree : update table
- [ ] Start writing new paper 
  - [ ] Environmental Modelling & Software ?



## Files

### running

- test_karstnet.ipynb    
- test_networkx.ipynb
- assortativity.ipynb
- assortativity_simple_cases.ipynb
- KGraph tests.ipynb
- edge_attribute.ipynb   
- simplify_graph.ipynb
- testPauline.py
- testPlot.py
- testCharacterize.py

### Files that are not running any more 

(not compatible with current version):

- connected_component_tests.ipynb




## DONE

### Code

- [x] add entropy of lengths (PC : partially done, still differences with previous works)
- [x] add entropy of orientations
- [x] read/import directly Plines (PC) 
- [x] (PC) - import Plines with properties on vertices
- [x] (PC) - export to gocad
- [x] (PC) - changer la structure pour faire une classe kgraph qui possèderait un graphe complet et un graphe simplifi
- [x] (PC) - corriger l'average SPL pour gérer les réseaux à plusieurs composantes connexes
- [x] (PC) - rajouter au averageSPL la possibilité de pondérer par les longueurs