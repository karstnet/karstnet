# Karstnet - statistics of karstic networks

Karstnet is a python3 project providing tools for the statistical analysis of karstic networks.

[![Documentation Status](https://readthedocs.org/projects/karstnet/badge/?version=latest)](https://karstnet.readthedocs.io/en/latest/?badge=latest)



Version 1.2.5 - August 2024 

Please check the file [changelog.md](changelog.md) to track the novel functionalities of karstnet.


## Simple installation from pypi.org

If you just need to use karstnet (e.g. in Jupyter notebooks), the simplest way to install it is to get it from [Pypi.org](https://pypi.org/) by typing the command:

`pip install karstnet`

This command should install directly all the required dependencies for a fully functional version of karstnet and you do not need to download manually anything.

## Installation from github

If you want to access the source code and potentially contribute. You should follow the following steps.

### 1. Download

Download karstnet from the [Github repository](https://github.com/karstnet/karstnet): green button "clone or download". Then, unzip it on your computer. 


### 2. Go in the directory

Once this has been done, open a Python prompt (like the Anaconda prompt), and go to your directory location (ex: 
`cd C:\Users\YourName\Documents\KARSTNET\karstnet`).

### 3. Launch the local installation

You have three options: the minimal one will install the code.

`pip install .`

Do not forget the point "." at the end of the command

If you want to run it directly from the source files (useful for development):

` pip install -e .`

You can also install all the dependencies used in karstnet to have complete development environment:

` pip install -e .[dev]`

All these options can also be run without moving previously in the karstnet folder. 
In that case, just type in your Anaconda prompt :

` pip install -e your\path\to\karstnet` 

### 4. Testing

If you start modifying the code, you should regularly check that you did not break some key features. For that you can run the unit tests.  From source directory, and after instaling **karstnet** run:

`pytest tests/test_karstnet.py`

## In Jupyter notebooks

Example of jupyter notebooks are provided to help you use **karstnet**. 
To use **karstnet** in notebooks, you just have to write

`import karstnet as kn`

A call-test function is available to help you check if the package is ready to use : just type: 
`kn.test_kn()`

## Documentation

The html documentation is available in the sub directory:  [docs/_build/html/index.html](docs/_build/html/index.html)
 and it is available online at: https://karstnet.readthedocs.io/

A french version of a guide for students willing to code karstnet extensions is available:
[FR_GuideDebutant_karstnet_jupyter_github.pdf](docs/pdf/FR_GuideDebutant_karstnet_jupyter_github.pdf)


**Remark on ENTROPIES**:
Note that, by default, Karstnet computes Entropies as described in the paper (mode = "default") : 
on normalized values ranged on 10 bins for branch lengths and on 18 bins of 10° for orientations

**If you want to compute entropies using Sturges'rule use** : 
l_entrop = myKGraph.length_entropy(mode = "sturges")
or_entropy = myKGraph.orientation_entropy(mode = "sturges")


## Reference and Corrigendum

The karstnet package implements some of the statistical metrics that were
investigated and discussed in:
Collon, P., Bernasconi D., Vuilleumier C., and Renard P., 2017, Statistical
metrics for the characterization of karst network geometry and topology.
Geomorphology. 283: 122-142 doi:10.1016/j.geomorph.2017.01.034
<http://dx.doi.org/doi:10.1016/j.geomorph.2017.01.034>

An **updated paper** (see remarks below) is available in the "doc" folder of this github and 
can be downloaded here  <https://hal.univ-lorraine.fr/hal-01468055v3/document>. 
 (complete link : <https://hal.univ-lorraine.fr/hal-01468055>)

**Concerning the paper, important remarks should be made :** 

There was some **errors** in the old Matlab implementation (the one used for the paper) **that have been corrected in Karstnet**. 
A **corrigendum** has been published in Geomorphology journal : Geomorphology 389, 107848. <http://dx.doi.org/doi:10.1016/j.geomorph.2021.107848>.
The **results** obtained on the same 34 networks than the ones used for the paper but 
with the **implementation of Karstnet** are proposed for information in the **doc part** (New_Statistics_results.xls) 
as well as the updated author version of the paper in pdf : **2016Pap_Collon_Geomorphology_Autho_Upd_2021.pdf**.

Here we summarize the main differences : 

**Correlation of Vertex Degrees, rk** : The corrected values of the correlation of vertex degree, rk, are all negative, indicating that the karstic networks in our data set are disassortative as it was reported for other natural networks by Newman (2002). 
	
**Branch lengths entropy, Hlen** : In the previous Matlab code, we computed the branch length entropy on 11 bins instead of the 10 bins described in the paper. Correcting this (computing on 10 bins) slightly decreases the values, ranging now from 0.07 to 0.67 instead of 0.18 to 0.74 (page 9). The Karstnet values are now correct. 

**CVlengths** : just an error of transcription in the table 2 of the paper where CVlen was supposed to be provided in %, but was provided in standard number

**branch lengths** : In the previous code, when computing the mean length of the branches, the looping branches (branches that closes on the starting point) were ignored. If this is meaningful for tortuosity computation, these branches should not be ignored for the mean length computation. This has been corrected and explains the minor differences observable for the mean length values of Agen Allwed, Daren Cilau, Foussoubie Goule, Krubera, Lechuguilla, MammuthHöhle, Ratasse, SaintMarcel, Sakany, Shuanghe, SiebenHengsteFull and SiebenHengsteSP2.
	
**SPL** : In addition, the previous code also ignored independent connected components of two nodes. This is not justified. This correction slightly impacts the values of the (SPL) coefficient of Agen Allwed, Arphidia Robinet, Daren Cilau, Grotte du Roy, Krubera, Lechuguilla, Llangattwg, MammuthHöhle Ratasse, SaintMarcel, Sakany, SiebenHengsteUpPart

We are sorry for any inconvenience.
