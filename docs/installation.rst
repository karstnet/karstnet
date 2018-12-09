Installation guide
==================

Installation
------------

Install Karstnet from source. From the project's main directory, type::

    pip install .

If you want to run it directly from source (useful for development)::

    pip install -e .


Testing
-------

To check if the installation was done properly, from the source directory,
and after instaling **karstnet**, run::

    pytest tests/test_karstnet.py


Dependencies
------------
karstnet requires the following python packages to function properly:
 * networkx
 * numpy
 * scipy-stats
 * matplotlib.pyplot
