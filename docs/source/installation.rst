Installation guide
==================

Karstnet is available on:
    
    - `PyPI <https://pypi.org/project/karstnet>`_ (The Python Package Index)
    - `Github <https://github.com/karstnet/karsnet>`_ (repository)

Installation from `PyPI <https://pypi.org/project/karstnet>`_
-------------------------------------------------------------

In a terminal type::

    pip install karstnet

Or, equivalently: `python -m pip install karstnet`.


Installation from the `Github repository <https://github.com/karstnet/karsnet>`_
--------------------------------------------------------------------------------

In a terminal, change directory where to download Karstnet, and type::

    git clone https://github.com/karstnet/karstnet.git
    cd karstnet
    pip install .

.. tip::
    Use `pip install . --verbose` or `pip install . -v` for printing (more) messages during the installation.

Alternatively:

- Instead of `git clone ...`, you can download Karstnet from the `Github repository <https://github.com/karstnet/karstnet>`_: click on the green button "code" and choose "Download ZIP". 
- Then, unzip the archive on your computer.
- Finally, in a terminal, go into the unzipped directory, and type `pip install .`

.. tip::
    To check if the installation was done properly, from the source directory, and after instaling **karstnet**, run::

    pytest tests/test_karstnet.py


Dependencies
------------
karstnet requires the following python packages to function properly:
 * geopandas
 * matplotlib
 * mplstereonet
 * networkx
 * numpy
 * pandas
 * pyvista
 * scipy
