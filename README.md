# Karstnet - statistics of karstic networks

Karstnet is a python3 project providing tools for the statistical analysis of karstic networks.

[![Documentation Status](https://readthedocs.org/projects/karstnet/badge/?version=latest)](https://karstnet.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/UniNE-CHYN/karstnet.svg?branch=master)](https://travis-ci.org/UniNE-CHYN/karstnet)


Version 1.2.0 - May 2021 - Updated by Pauline Collon



## Installation

Download karstnet from this Github platform : green button "clone or download". Then, unzip it on your computer. 

Once this have been done, you can open a Python prompt (like the Anaconda prompt) to install it 

Install from source (from project main directory), just to use it (e.g. in Jupyter notebooks): 

1)In your prompt, go to your directory location (ex: 
`cd C:\Users\YourName\Documents\KARSTNET\karstnet`

2) then launch karstnet installation by taping:
`pip install .`
Do not forget the point "." at the end of the command

If you want to run it directly from source (useful for development):
` pip install -e .`

Both these options can also be run without moving previously in the karstnet folder. 
In that case, just type in your Anaconda prompt :

` pip install -e your\path\to\karstnet` 

## Testing

From source directory, and after instaling **karstnet** run:

`pytest tests/test_karstnet.py`

## In Jupyter notebooks

Example of jupyter notebooks are provided to help you use **karstnet**. 
To use **karstnet** in notebooks, you just have to write

`import karstnet as kn`

A call-test function is available to help you check if the package is ready to use : just type: 
`kn.test_kn()`

## Documentation

The html documentation is available in the sub directory:  ``docs/_build/html/index.html``

Also available online at: https://karstnet.readthedocs.io/

## Reference and Corrigendum

The karstnet package implements some of the statistical metrics that were
investigated and discussed in:
Collon, P., Bernasconi D., Vuilleumier C., and Renard P., 2017, Statistical
metrics for the characterization of karst network geometry and topology.
Geomorphology. 283: 122-142 doi:10.1016/j.geomorph.2017.01.034
<http://dx.doi.org/doi:10.1016/j.geomorph.2017.01.034>

An **updated paper** (see remarks below) is available in the "doc" folder of this github and 
can be downloaded here  <https://hal.univ-lorraine.fr/hal-01468055v2/document>. 
 (complete link : <https://hal.univ-lorraine.fr/hal-01468055>)

**Concerning the paper, important remarks should be made :** 

There was some **errors** in the old Matlab implementation (the one used for the paper) that have been corrected in Karstnet. 
A **corrigendum** is currently submitted to the journal but the updated author version of the paper has 
already been downloaded on the HAL platform (link above).
The **results** obtained on the same 34 networks than the ones used for the paper but 
with the **implementation of Karstnet** are proposed for information in the **doc part** (New_Statistics_results.xls). 

Here we summarize the main differences : 

	- **Correlation of Vertex Degrees** : contrary to what was previously computed, all studied karstic networks appear to be disassortative ( rk < 0). 
	
	- **Branch lengths entropy** : in the old matlab code, we computed entropy on 11 bins instead of the 10 wanted bins. 
The Karstnet values are now correct. 

Note that, by default, Karstnet computes Entropies as described in the paper (mode = "default") : 
on normalized values ranged on 10 bins for branch lengths and on 18 bins of 10° for orientations
If you want to compute entropies using Sturges'rule use : 
l_entrop = myKGraph.length_entropy(mode = "sturges")
or_entropy = myKGraph.orientation_entropy(mode = "sturges")

	- **CVlengths** : just an error of transcription in the table 2 of the paper where CVlen was supposed to be provided in %, 
	but was provided in standard number

	- **branch lengths** : in the previous Matlab code, loops were ignore to compute the average branch lengths
	or the CV lengths, this induces very slight differences on some networks computation (ex : Agen Allwed)
	
	- **SPL** : slight differences are also noticeable for some networks as isolated connected components were ignored
	in the previous code, while we decided to keep them in this clean updated Python version (ex :  Agen Allwed is now 11.17 instead 11.72)

We are sorry for any inconvenience.
