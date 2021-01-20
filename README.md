# Karstnet - statistics of karstic networks

Karstnet is a python3 project providing tools for the statistical analysis of karstic networks.

[![Documentation Status](https://readthedocs.org/projects/karstnet/badge/?version=latest)](https://karstnet.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/UniNE-CHYN/karstnet.svg?branch=master)](https://travis-ci.org/UniNE-CHYN/karstnet)


Version 1.1.0 - August 2020



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

The paper can be downloaded here
<http://members.unine.ch/philippe.renard/articles/collon2017.pdf> or here <https://hal.univ-lorraine.fr/hal-01468055v1>. 

**Concerning this publication, important remarks should be made :** 

- There was an **error** in the Matlab implementation of **Correlation of Vertex Degrees** used for the 2017 paper. 
The Karstnet implementation corrects this, and, as a result, all studied karstic networks appears to be disassortative ( rk < 0)
contrary to what was initially found. We wait to finish identify all potential errors before sendin a corrigendum to the Journal.

- The **implementation of the entropy is not the same** in Karstnet that the ones used in the paper. For the moment in Karstnet, 
the number of bins is computed using Sturges'rule and is thus varying between networks. As a result, the values obtained with Karstnet 
are different from the ones presented in the paper. We are studying the possibility of including an option for the user 
to choose between various implementations. 

- Some unitary tests are still under development. The **results** obtained on the same networks than the ones used for the paper but 
with the **implementation of Karstnet** are proposed for information in the **doc part**. 

