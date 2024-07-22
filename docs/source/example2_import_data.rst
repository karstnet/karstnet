A second example
==================

In this page, we show how one can import data in Karstnet and compute
the statistical properties of the karstic network.

Loading a network
------------------

First, let us import the karstnet module::

    import karstnet as kn


There are two built in functions to import data from ASCII files.


The **First**
one assumes that the data describing the geometry of the network
are stored in two ASCII files: one containing the position of the nodes,
the other containing the description of the links between the nodes
(the edges).

The convention used by Karstnet is that the two files have the same basename,
and have different extensions:

- `basename_nodes.dat`: the matrix of 3D node coordinates, one line per node,
  one column per coordinates

- `basename_links.dat`: the list of links (edges), one edge per line,
  two nodes ID for each link

To load the karstic network geometry from these files, use this command::

    myobject = kn.from_nodlink_dat('basename')

After reading the two data files,
karstnet creates the internal structure of the graph and
precomputed all the information that it requires, like for example the
distance between the nodes, the orientation of the edges, and stored
everything in the object.

In particular, the simplified graph is created during the import. The
simplified graph is used to quantify the topological characteristics
of the network.

A **second** possibility is to import the network from a GOCAD ASCII
Pline file::

  myobject = kn.from_pline('filename.pl')

More information in the API documentation :py:func:`base.from_pline`

Characterizing the network
--------------------------

To check the import and the network simplification,
it is possible to plot the map of the original
and simplified networks side by side in 2D::

    myobject.plot()

Other functions are available to plot the network in 3D::

    myobject.plot3()

Finally, to compute the statistical
characteristics of this network, just call the characterize function::

    results = myobject.characterize_graph( verbose = True )
