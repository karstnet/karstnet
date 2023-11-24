"""
Karstnet
========

Karstnet is a Python package for the analysis of karstic networks.

This file aims to compare karstic networks thank to characteristics compute by
the Karsntet code in KGraph class. The networks are together by groups and
showed in a dendrogram

"""

# ----External libraries importations

import glob
import pandas as pd
import openpyxl as op
from scipy.cluster import hierarchy
from scipy.spatial import distance
import matplotlib.pyplot as plt
import numpy as np
from typing import Optional, List, Tuple, Union, Any

# ----Internal module dependent
import karstnet as kn


# -----------------------------------------------------------------------------
#                                PUBLIC FUNCTIONS (Louise Houbre 2021)
# -----------------------------------------------------------------------------
def clustering_complet(db_dirname: str, nt_dirname: str, extension: str,
                       charact_not_used: Optional[List[str]] = None,
                       figsize: Tuple[float] = (10, 10)):
    """
    Function to compare a network to  a databased and plot the dendrogramm
    associated and the comparaison lines.

    Used the  functions organize_dataset and dendrogram and _bar_comparaison

    Parameters
    ----------
        db_dirname : string
            dirname of the database record

        nt_dirname : string
            dirname of the candidate network to compare to the database

        extension : string
            the extension of the file describing the candidate network (.dat,
            .pl, .sql)

        charact_not_used : list[string]
            a list composed of the caracteristiques not used to compare the
            networks and make the dendrogram
            For example 'l', 'H0', ...

        figsize : tuple
            contains the (x,y) dimension of the figure

    Examples
    --------
    >>> nt_dirname_ = '../data/toadd'
    >>> db_dirname_ = '../network/'
    >>> extension_ = '.dat'
    >>> clustering_complet(db_dirname_, nt_dirname_, extension_)
    >>> charact_not_used_ = ['l']
    >>> clustering_complet(db_dirname_, nt_dirname_, extension_,
    >>> charact_not_used_)
    """
    organize_dataset(db_dirname, nt_dirname, extension)

    dendrogram(db_dirname + 'compute_db_nt.xlsx', charact_not_used, figsize)

    name_charact = ['l', 'CVl', 'Hl', 't', 'H0', 'kmoy', 'CVk', 'rk', 'SPL',
                    'CPD']
    plt.show()
    _bar_comparaison(name_charact, nt_dirname, extension, db_dirname)


def dendrogram(name_file: str, charact_not_used: Optional[List[str]] = None,
               figsize: Tuple[float] = (6, 6)):
    """
    Plot the dendrogram

    Used in the function plot_clustering

    Parameters
    ----------
        name_file : str
            a xlsx file columns are characters and lines are the network, used
            to created the dendrogram

        charact_not_used : list[string]
            a list composed of the caracteristiques not used to compare the
            networks and make the dendrogram
            For example 'l', 'H0', ...

        figsize : tuple
            contains the (x,y) dimension of the figure

    Examples
    --------
    >>> name_file_ = '../network/compute_db.xlsx'
    >>> kn.dendrogram(name_file_)
    >>> name_file_ = '../network/compute_db_nt.xlsx'
    >>> kn.dendrogram(name_file_)
    """
    typ = name_file.split('_')
    cpt = 0
    print('type = ', typ)
    if typ[-1] == 'nt.xlsx':
        cpt = 1
    name_network, name_charact, list_values = _read_file_xls(name_file, cpt)

    list_val = _normalize(list_values)[1]
    print('Values are normalized')
    print('------------------------------------------------------------')
    print('All datas are ready to be compared')

    if charact_not_used is not None:
        for charact in charact_not_used:
            list_val, name_charact = _del_charact(charact, name_charact,
                                                  list_val)

    plt.subplots(figsize=figsize)
    hierarchy.set_link_color_palette(['r', 'b', 'y', 'g', 'p'])
    dist = hierarchy.ward(distance.pdist(list_val))
    print('------------------------------------------------------------')
    print('Dendrogram of ' + name_file + '  file ')
    hierarchy.dendrogram(dist, above_threshold_color='black',
                         labels=name_network,
                         orientation='left')


def organize_dataset(db_dirname: str, nt_dirname: str, extension: str):
    """
    Importation of network database
    Creation of  compute_db.xlsx, compute_db_nt.xlsx and extremum.dat files

    Used in the function clustering_complet

    If the user want to plot only the dendrogram of a database, it's necessary
    to create the files .Xlsx and .dat before run the function dendrogram

    Parameters
    ----------
        db_dirname : string
            dirname of the database record

        nt_dirname : string
            dirname of the candidate network to compare to the database

        extension : string
            the extension of the file describing the candidate network (.dat,
            .pl, .sql)
    Examples
    -------
    >>> nt_dirname_ = '../data/candidate'
    >>> db_dirname_ = '../network/'
    >>> extension_ = '.dat'
    >>> organize_dataset(db_dirname_, nt_dirname_, extension_)
    """
    _create_file_charact(db_dirname)
    print('------------------------------------------------------------')
    _create_file_extremum(db_dirname)
    print('------------------------------------------------------------')
    print('Creation of the KGraph corresponding to the network to compare')
    _network_to_compare(db_dirname, nt_dirname, extension)
    print('------------------------------------------------------------')
    print(' Files are created ')


# -----------------------------------------------------------------------------
#                             NOT PUBLIC FUNCTIONS (Louise Houbre 2021)
# -----------------------------------------------------------------------------

# ***********************************************
# Private function for plot the comparaison lines
# ***********************************************

def _bar_comparaison(name_charact: List[str], nt_dirname: str, extension: str,
                     db_dirname: str):
    """
    NOT PUBLIC

    Plot the barres to compare the values normalized of the network to the
    extremum values of the database

    Parameters
    ----------
        name_charact : list(string)
            list composed of the caracteristiques used to compare the networks
            and make the dendrogram

        nt_dirname : string
            dirname of the candidate network to compare to the database

        extension : string
            the extension of the file describing the candidate network (.dat,
            .pl, .sql)

        db_dirname : string
            dirname of the database record
    """
    if 'name' in name_charact:
        name_charact.remove('name')

    list_values = _read_file_xls(db_dirname + 'compute_db.xlsx')[-1]
    normalize = _normalize(list_values)
    list_norm = normalize[-1]
    list_inv = _inverse_list(list_norm)
    list_min, list_max = list(), list()

    list_val = _network_to_compare(db_dirname, nt_dirname, extension)[1][2:]
    list_val_norm = _normalize_one(list_val, normalize[0])

    for i in range(len(list_inv)):
        min_, max_ = min(list_inv[i]), max(list_inv[i])
        list_min.append(min_)
        list_max.append(max_)

    M = len(name_charact)
    for i in range(M):
        min_, max_, val = list_min[M - i - 1], list_max[M - i - 1], \
                          list_val_norm[M - i - 1]

        # Le caractère charact ne correspond pas aux valeurs utilisées à chaque
        # boucle, mais permet un affichage de la légende bon: c'est à dire que
        # la bonne chaine de caractère est affichée en face de la ligne
        # correspondante
        charact = name_charact[i]
        x = np.linspace(min_, max_, 10)
        y = [i] * 10
        plt.plot(x, y, label=charact, color='navy')
        plt.plot(val, i, 'D', color='crimson')
        plt.legend(bbox_to_anchor=(1, 0., 0.4, 1.), loc='right',
                   labelspacing=1.2, frameon=False)

    # ************************************************
    # Private functions for read and used extern files (Louise Houbre 2021)
    # ************************************************


def _create_file_charact(dirname: str) -> str:
    """
    NOT PUBLIC

    Création du fichier pour la base de bonnées : calcul des caractéristiques
    des réseaux présents dans l abase de données permettant enseuite la
    comparaison avec la réseau suplémentaire

    Parameters
    ----------
        dirname : string
            dirname of a  record composed of networks

    Returns
    -------
    string
        new_dirname : the name and direction of the .xlsx file created

    """
    new_dirname = dirname + 'compute_db.xlsx'

    try:
        open(new_dirname, 'r')

    except:
        # This list is used to replace long names of attributes with short names
        # in the final excel file
        list_metrics = [["l", "mean length"], ["CVl", "cv length"],
                        ["Hl", "length entropy"],
                        ["t", "tortuosity"], ["H0", "orientation entropy"],
                        ["kmoy", "mean degree"],
                        ["CVk", "cv degree"],
                        ["rk", "correlation vertex degree"], ["SPL", "aspl"],
                        ["CPD", "cpd"]]

        # Generates a list of all base names for the nodes and links data file
        listfname = glob.glob(dirname + '*_links.dat')
        listbname = []
        for tmp in listfname:
            listbname.append(tmp[11:-10])

            # Prints the list to check
        print(listbname)
        # Creates an empty dictionnary with the proper keys to store the results
        pdict = dict()
        pdict['name'] = []
        for short, long in list_metrics:
            pdict[short] = []

            # This is the main loop over all the data sets
        for basename in listbname:
            filename = dirname + basename
            print(
                '------------------------------------------------------------')
            print('Reading network:', filename)
            print(
                '------------------------------------------------------------')

            karst = kn.from_nodlink_dat(filename)
            results = karst.characterize_graph()
            pdict['name'].append(basename)
            for short, long in list_metrics:
                pdict[short].append(results[long])

            # Transforms the data in a pandas dataframe to facilitate export
        df = pd.DataFrame(data=pdict)

        # Show the results
        # print('------------------------------------------------------------')
        # print(df)

        # Export them in an excel file
        df.to_excel(new_dirname)
        print('Export results to file compute_db.xlsx')

    return new_dirname


def _create_file_extremum(dirname: str):
    """
    NOT PUBLIC

    Create a .dat file composed of the minimum, maximum and mean values of each
    characteristics

     Parameters
     ---------
        dirname : string
            the name and direction of a .xlsx file where we want to compute the
            maximum, minimum and mean values of each columns.
    """
    # Check if the file compute_db.xlsx exists. If its the cas, we dont write
    # it twice
    _create_file_charact(dirname)

    # Compute the min, max and meand of each charact of the database and write
    # them on the .dat file
    name_line, list_charact, list_values = _read_file_xls(
        dirname + 'compute_db.xlsx')
    list_stat = _normalize(list_values)[0]
    name_stat = ['min', 'max', 'mean']

    name = dirname + 'extremum.dat'
    file = open(name, 'w')
    file.write('name ')

    for charact in list_charact:
        file.write(str(charact) + ' ')

    for i in range(3):
        file.write('\n' + name_stat[i] + ' ')
        for j in range(len(list_charact)):
            file.write(str(list_stat[j][i]) + ' ')

    file.close()
    print('File extremum.dat created')


def _network_to_compare(db_dirname: str, nt_dirname: str, extension: str) -> \
        Tuple[str, List[Union[str, Any]]]:
    """
    NOT PUBLIC

    Parameters
    ----------
        nt_dirname : string
            dirname of the candidate network to compare to the database

        extension : string
            the extension of the file describing the candidate network (.dat,
            .pl, .sql)

        db_dirname : string
            dirname of the database record

    Returns
    -------
    string
        nwfile_dirname : the dirname of the new file created : a xlsx with the
        data of the database and of the candidate network at the end of the file
    """

    file_name = db_dirname + 'compute_db.xlsx'

    file = op.load_workbook(file_name)
    source = file.active
    source = pd.DataFrame(source.values)

    if extension == '.pl' or extension == 'pl':
        graph = kn.from_pline(nt_dirname)
    if extension == '.dat' or extension == 'dat':
        graph = kn.from_nodlink_dat(nt_dirname)
    if extension == '.sql' or extension == 'sql':
        graph = kn.from_therion_sql(nt_dirname)
    else:
        print('The extension ' + extension +
              ' is not considerated by the Karstnet code')

    results = graph.characterize_graph()
    val = list(results.values())
    name = nt_dirname.split('/')
    network = name[-1]
    ntw = ['-1', network]
    values = [val[i] for i in range(4)]
    values.extend([val[8], val[11], val[12], val[13], val[9], val[10]])
    val_ntw = [round(value, 6) for value in values]
    ntw.extend(val_ntw)

    line_network = pd.Series(ntw)
    new_file = source.append(line_network, ignore_index=True)
    nwfile_dirname = db_dirname + 'compute_db_nt.xlsx'
    new_file.to_excel(nwfile_dirname)

    print('File compute_db_nt.xlsx is created')
    return nwfile_dirname, ntw


def _read_file_xls(file_name: str, cpt: float = 0) -> \
        Tuple[List[str], List[str], List[List[float]]]:
    """
    NOT PUBLIC

    Read the xls file and convert the table in a List of List

    Parameters
    ----------
        file_name : string
            name of the file to read

        cpt : float 0 or 1
            1 : when the candidate network is considerated in the file
            0 : when the file contains only the values of the database

    Returns
    -------
    list
        name_network[1:] : name networks list

    list
        row_list[0] : name of characteriqtique (l, CVl,HO,...)

    list
         row_list[1:] : values
    """
    file = op.load_workbook(file_name)
    source = file.active
    wsheet = file.copy_worksheet(source)

    row_list = list()
    name_lines = list()

    if cpt == 0:

        for row in wsheet.values:
            row_list.append(list(row[2:]))
            name_lines.append(row[1])

        return name_lines[1:], row_list[0], row_list[1:]

    if cpt == 1:
        for row in wsheet.values:
            row_list.append(list(row[2:]))
            name_lines.append(row[2])

        row_list_ = [row_list[1]]
        for line in range(2, len(row_list)):
            li_ = list()
            for column in range(1, len(row_list[line])):
                li_.append(float(row_list[line][column]))
            row_list_.append(li_)

        return name_lines[2:], row_list_[0], row_list_[1:]

    # ************************************
    # Private functions to compute values (Louise Houbre 2021)
    # ************************************


def _del_charact(charact: str, name_charact: List[str],
                 list_values: List[List[float]]) -> \
        Tuple[List[List[float]], List[str]]:
    """
    NOT PUBLIC

    Return a new list_values without the values corresponding to the
    charact imput

    Parameters
    ----------
        charact : string
            the character to del of the name_charact list

        name_charact : list
            a list composed of the caracteristiques used to compare the
            networks

        list_values : list
             Each list correspond to a line of the compute_db_nt.xlsx file

    Returns
    -------

    list
        list_values_del : the list_values without the columns correspondig to
        the character inputeach list correspond to a line of the
        compute_db_nt.xlsx file

    list
        name_charact : a list composed of the caracteristiques used to compare
        the networks whitout the charact input
    """
    if charact not in name_charact:
        print("The character " + charact + " can't be delete. \
        It is not a character used to make the dendrogram.")
        return list_values, name_charact

    else:
        list_values_del = list()
        pos = name_charact.index(charact)
        name_charact.remove(charact)

        for row in list_values:
            row.pop(pos)
            list_values_del.append(row)
        return list_values_del, name_charact


def _normalize(list_values: List[List[float]]):
    """
    NOT PUBLIC

    Normalize the values for each characteristics

    Parameters
    ----------

        list_values : list
            list of list corresponding of the file read : each under list is a
            line

    Returns
    -------

    list
        list_stat : a list of list. Each under line is composed of the min,
        the max and the mean value of the character considered

    list
        list_norm : the list_values normalized
        using the formul  : (val - mean) / (max - min)
    """
    l1 = _inverse_list(list_values)
    list_stat = list()
    l2 = list()
    for list_ in l1:
        min_ = min(list_)
        max_ = max(list_)
        mean_ = sum(list_) / len(list_values)
        list_stat.append([min_, max_, mean_])
        val_norm = [(val - mean_) / (max_ - min_) for val in list_]
        l2.append(val_norm)

    list_norm = _inverse_list(l2)
    return list_stat, list_norm


def _normalize_one(value_ntw: List[float], list_stat: List[List[float]]) -> \
        List[float]:
    """
    NOT PUBLIC

    Normalize the values for each characteristics for one network

    Parameters
    ----------

        value_ntw : list
            list of the value of the candidate network in the same order as in
            the initial .xlsx file, with all chaarcteristiscs

        list_stat : list
            a list of list. Each under line is composed of the min, the max
            and the mean value of the character considered, compute only in the
            dataset

    Returns
    -------

    list
        list_norm : list of the value of the candidate network normalized by the
        min, max and mean value of the database
    """
    list_norm = list()
    for i in range(len(value_ntw)):
        list_norm.append((value_ntw[i] - list_stat[i][2]) /
                         (list_stat[i][1] - list_stat[i][0]))
    return list_norm


def _inverse_list(list_: List[List[float]]) -> List[List[float]]:
    """
    NOT PUBLIC

    Inverse de list of lists.

    For example
    Inverse a list as : [v1.1, v1.2, v1.3],[v2.1, v2.2, v2.3]]
    to return [[v1.1, v2.1], [v1.2, v2.2], [v1.3, v2.3]]

    Parameters
    ----------
        list_ : list
            list of list

    Returns
    -------

    list
       list_int : list of list
    """
    nbr = len(list_[0])
    list_int = list()
    for i in range(nbr):
        val = list()
        for li in list_:
            val.append(li[i])
        list_int.append(val)

    return list_int
