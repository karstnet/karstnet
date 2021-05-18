#    Copyright (C) 2018 by
#    Philippe Renard <philippe.renard@unine.ch>
#    Pauline Collon <pauline.collon@univ-lorraine.fr>
#    All rights reserved.
#    MIT license.
#

"""
Test module for the Karstnet class.
Execute with pytest : `pytest test_karstnet.py`
"""

import karstnet as kn
import numpy as np
import pytest

"""
Utilities to construct the tests
"""


def float_eq(a, b, tolerance=1e-3):
    """
    Returns True if the difference between a and b is lower
    than tolerance.
    """
    return abs(a-b) < tolerance


@pytest.fixture
def periodic():
    """ Returns a simple periodic graph

        The graph contains 9 nodes on regular grid. The extreme edges are
        connected to reproduce a periodic connectivity pattern.
        In this way, all nodes have a degree 4.

        Used to test:
        central point dominance = 0,
        average shortest path = 1.5,
        mean degree = 4,
        cv degree = 0,
        correlation vertex degree = 1,
        tortuosity = 1,

    """
    pos = {1: (-1, 0), 2: (-0.2, -0.2), 3: (1, 0), 4: (-1.5, 1.5), 5: (0, 1),
           6: (1.5, 1.5), 7: (-1.5, -1.5), 8: (0, -1), 9: (1.5, -1.5)}
    edges = [(1, 2), (2, 3), (4, 5), (5, 6), (7, 8), (8, 9), (1, 4), (1, 7),
             (2, 5), (2, 8), (3, 6), (3, 9), (5, 8), (1, 3), (6, 4), (4, 7),
             (6, 9), (9, 7)]
    k = kn.KGraph(edges, pos)
    return k


@pytest.fixture
def disassortative():
    """ Returns a simple disassortative graph

        The graph contains 2 connected central nodes. Each of these nodes
        is surrounded by 10 other nodes. All distances are 1.

        Used to test:
        mean length = 1.000
        cv length = 0.000
        length entropy = 0.000
        tortuosity = 1.000
        correlation vertex degree = -0.909
        central point dominance = 0.703

    """
    npt = 10
    pos = {}
    pos[0] = (0.5, 0)
    pos[npt+1] = (-0.5, 0)
    edges = [(0, npt+1)]

    for i in range(npt):
        angle = 1.2 * np.pi * i / (npt - 1) - np.pi / 2 - 0.1 * np.pi
        pos[i+1] = (0.5 + np.cos(angle), np.sin(angle))
        pos[i+npt+2] = (-0.5 - np.cos(angle), np.sin(angle))
        edges.append((0, i+1))
        edges.append((npt+1, i+npt+2))

    k = kn.KGraph(edges, pos)
    return k


@pytest.fixture
def assortative():
    """
    Returns a simple assortative graph

    The graph contains a core central grid and lateral branches.

    Used to test:
        mean length = 1.118
        length entropy = 0.206
        orientation entropy = 0.698
        aspl = 2.947
        cpd = 0.155
        mean degree = 3.520
        cv degree = 0.662
        correlation vertex degree = 0.353
    """
    pos = {0: [0.0, 0.0], 1: [0.0, 1.0], 2: [1.0, 0.0], 3: [0.0, -1.0],
           4: [-1.0, 0.0], 5: [1.0, 1.0], 6: [1.0, -1.0], 7: [-1.0, -1.0],
           8: [-1.0, 1.0], 9: [-0.5, 2.0], 10: [0.5, 2.0], 11: [2.0, 0.5],
           12: [2.0, -0.5], 13: [0.5, -2.0], 14: [-0.5, -2.0],
           15: [-2.0, -0.5], 16: [-2.0, 0.5], 17: [-0.5, 3.0],
           18: [0.5, 3.0], 19: [3.0, 0.5], 20: [3.0, -0.5],
           21: [0.5, -3.0], 22: [-0.5, -3.0], 23: [-3.0, -0.5],
           24: [-3.0, 0.5]}
    edges = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4),
             (4, 1), (0, 5), (0, 6), (0, 7), (0, 8), (1, 5), (5, 2),
             (2, 6), (6, 3), (3, 7), (7, 4), (4, 8), (8, 1), (8, 9),
             (1, 9), (1, 10), (5, 10), (5, 11), (2, 11), (2, 12),
             (6, 12), (6, 13), (3, 13), (3, 14), (7, 14), (7, 15),
             (4, 15), (4, 16), (8, 16), (9, 17), (10, 18), (11, 19),
             (12, 20), (13, 21), (14, 22), (15, 23), (16, 24))

    k = kn.KGraph(edges, pos)
    return k


@pytest.fixture
def complete():
    """ Returns a complete graph with 10 nodes

        The graph contains a core central grid and lateral branches.

        Used to test:
        length entropy = 0.635
        orientation entropy = 0.841
        aspl = 1.0
        cpd = 0.0
        mean degree = 9.0
        cv degree = 0.0
        correlation vertex degree = 1.0
    """
    npt = 10
    pos = {}
    edges = []

    for i in range(npt):
        angle = 2*np.pi * i / (npt - 1)
        pos[i] = (np.cos(angle), np.sin(angle))

    for i in range(npt):
        for j in range(npt):
            if i != j:
                edges.append((i, j))

    return kn.KGraph(edges, pos)


@pytest.fixture
def semibinary():
    """
    Returns a (semi) binary tree graph (one node has 3 edges)

    Used to test:
        aspl = 3.128
        cpd = 0.463
        mean degree = 1.846
        correlation vertex degree = -0.6
    """
    pos = {1: (0, 0), 2: (0, 1), 3: (-1, 2), 4: (1, 2), 5: (-1.5, 3),
           6: (-0.5, 3), 7: (0.5, 3), 8: (1.5, 3), 9: (-2, 4), 10: (-1, 4),
           11: (0, 4), 12: (1, 4), 13: (0.5, 5)}
    edges = [(1, 2), (2, 3), (2, 4), (3, 5), (3, 6), (4, 7), (4, 8),
             (5, 9), (5, 10), (7, 11), (7, 12), (7, 13)]
    return kn.KGraph(edges, pos)


"""
Main tests
"""


def test_mean_tortuosity(periodic):
    assert periodic.mean_tortuosity() == 1


def test_central_point_dominance(periodic, disassortative,
                                 assortative, complete):
    assert periodic.central_point_dominance() == 0
    assert complete.central_point_dominance() == 0
    assert float_eq(disassortative.central_point_dominance(), 0.703)
    assert float_eq(assortative.central_point_dominance(), 0.155)


def test_average_SPL(periodic, assortative, complete, semibinary):
    assert float_eq(complete.average_SPL(), 1.0)
    assert float_eq(periodic.average_SPL(), 1.5)
    assert float_eq(assortative.average_SPL(), 2.947)
    assert float_eq(semibinary.average_SPL(), 3.128)


def test_mean_degree_and_CV(periodic, assortative, complete, semibinary):
    md, cvde = periodic.mean_degree_and_CV()
    assert md == 4
    assert cvde == 0

    md, cvde = assortative.mean_degree_and_CV()
    assert float_eq(md, 3.520)
    assert float_eq(cvde, 0.662)

    md, cvde = semibinary.mean_degree_and_CV()
    assert float_eq(md, 1.846)  # Lower than 2 for a tree
    assert float_eq(cvde, 0.619)

    md, cvde = complete.mean_degree_and_CV()
    assert md == 9
    assert cvde == 0


def test_correlation_vertex_degree(periodic, disassortative,
                                   assortative, complete, semibinary):
    cvd = periodic.correlation_vertex_degree()
    assert cvd == 1
    cvd = disassortative.correlation_vertex_degree()
    assert float_eq(cvd, -0.909)
    cvd = assortative.correlation_vertex_degree()
    assert float_eq(cvd, 0.353)
    cvd = complete.correlation_vertex_degree()
    assert float_eq(cvd, 1.0)
    cvd = semibinary.correlation_vertex_degree()
    assert float_eq(cvd, -0.6)


def test_mean_length(disassortative, assortative):
    assert disassortative.mean_length() == 1
    assert float_eq(assortative.mean_length(), 1.118)


def test_coef_variation_length(disassortative):
    assert float_eq(disassortative.coef_variation_length(), 0, 1e-8)


def test_length_entropy(disassortative, assortative, complete):
    assert disassortative.length_entropy() == 0
    assert float_eq(assortative.length_entropy(), 0.206)
    assert float_eq(complete.length_entropy(), 0.635)


def test_orientation_entropy(assortative, complete):
    assert float_eq(assortative.orientation_entropy(), 0.698)
    assert float_eq(complete.orientation_entropy(), 0.841)
