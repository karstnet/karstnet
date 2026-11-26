"""
Module `_geom.geom`
-------------------

The module `_geom.geom` is accessible as `geom`. It contains functions
to compute geometric properties of karst conduits (networkx graphs) with splays.
"""


import numpy as np
import networkx as nx

import karstnet as kn


# ----------------------------------------------------------------------------
def get_local_coordinate_system(
        G, 
        pos_attr='pos', 
        update_graph=True, 
        local_coordinate_system_attr='local_coordinate_system', 
        verbose=1):
    """Calculates a local coordinate system for at each node in the graph G.

    At each graph node, a local coordinate system is defined by three orthogonormal vectors: u, v, w.
    The direction `w` is determined based on the positions of neighboring nodes and the degree of the node 
    and based on the sequence of survey by the cavers. The direction is calculated as follows:
    - for nodes with degree 1, the direction is calculated as the vector from the node to its neighbor.
    - for nodes with degree 2, the direction is calculated as the mean direction of the direction from 
    the first neighbor to the node and the direction from the node to the second neighbor.
    - for nodes with degree greater than 2, the direction is is not defined, as the node is an intersection 
    of multiple conduits, and the dimension of the conduit will be calculated as the mean of the splays 
    instead of the projection in the uv plan.

    Parameters
    ----------
    G : networkx.Graph
        Input graph containing the cave network with nodes and edges.
        Nodes must have `pos_attr` attribute for 3D coordinates.

    pos_attr : str, optional
        The name of the node attribute that contains the 3D coordinates, by default 'pos'
    
    update_graph : bool, default True
        - if `True`: the local coordinate system is set as a node attribute (named \
        `local_coordinate_system_attr` (see below))
        - if `False`: the graph is not updated (`local_coordinate_system_attr` is not used)
    
    local_coordinate_system_attr : str, default 'local_coordinate_system'
        The name of the node attribute that will contain the local coordinate system
        (if `update_graph` is `True`), by default 'local_coordinate_system'
        
    verbose : int, default: 1
        verbose mode, higher implies more printing (info)

    Returns
    -------
    dict
        Dictionary `local_coordinate_system` with keys: node id, values: numpy.ndarray of shape (3, 3)
        representing the local coordinate system at each node, where
        loc_coor_sys[0] = u vector, loc_coor_sys[1] = v vector, loc_coor_sys[2] = w vector, 
        with w being the conduit direction, u and v being orthogonal vectors in the uv plan 
        perpendicular to w, and u having its 3rd component equal to 0 (horizontal).
    """
    local_coordinate_system = {}

    for key in G.nodes():
        # if the node exists. some nodes were survey duplicate, removed from the main graph
        neighbors = list(G.neighbors(key))
        if len(neighbors) == 0:
            if verbose > 0:
                print(f'Node {key} has no neighbors. Conduit direction cannot be defined.')
            continue

        elif len(neighbors) == 1:
            # the node is the start or end of the cave (node of degree 1)
            w = np.asarray(G.nodes(pos_attr)[key]) - np.asarray(G.nodes(pos_attr)[neighbors[0]])   
            w_norm = np.linalg.norm(w)
            if np.isclose(w_norm, 0.0):
                if verbose > 0:
                    print(f'Conduit direction for node {key} (degree 1) cannot be defined (zero vector)')
                continue
            w = w / w_norm

        elif len(neighbors) == 2:
            # node of degree 2
            # Get direction from the first neighbor to the considered node and from the considered node to the second neighbor
            pos = np.asarray(G.nodes(pos_attr)[key])
            w1 = pos - np.asarray(G.nodes(pos_attr)[neighbors[0]])
            w1_norm = np.linalg.norm(w1)
            if not np.isclose(w1_norm, 0.0):
                w1 = w1 / w1_norm
            # else:
            #     w1 = np.zeros(3)
            w2 = np.asarray(G.nodes(pos_attr)[neighbors[1]]) - pos
            w2_norm = np.linalg.norm(w2)
            if not np.isclose(w2_norm, 0.0):
                w2 = w2 / w2_norm
            # else:
            #     w2 = np.zeros(3)
            
            # Get the average direction of the two vectors
            w = w1 + w2 
            w_norm = np.linalg.norm(w)
            if np.isclose(w_norm, 0.0):
                if verbose > 0:
                    print(f'Conduit direction for node {key} (degree 2) cannot be defined (zero vector)')
                continue
            w = w / w_norm
            
        else:
            # No projection for conduit intersections
            # The node should not contain any direction informations
            if verbose > 0:
                print(f'Node {key} has degree {len(neighbors)}. Conduit direction cannot be defined.')
            continue

        
        # Calculate the u and v vector for each node, to be used for the projection of the splay in the uv plan
        # * u, v are orthogonal vectors calculated, perpendicular to the conduit direction w
        # * u has its 3rd component equal to 0 (horizontal)
        # * v = w ^ u

        # Define u : normalized vector in the horizontal plan, orthogonal to vector w
        # if np.allclose(w[:2], 0.0):
        #     u = np.array([0.0, 1.0, 0.0])    
        # else:
        #     u = np.array([-w[1], w[0], 0.0])/np.sqrt(w[0]**2+w[1]**2)
        if np.allclose(w[:2], 0.0):
            if verbose > 0:
                print(f'Node {key} has a vertical conduit direction. Conduit direction is not computed.')
            continue

        u = np.array([-w[1], w[0], 0.0])/np.sqrt(w[0]**2+w[1]**2)
        v = np.cross(w, u)
  
        local_coordinate_system[key] = np.vstack((u, v, w))

    if update_graph:
        nx.set_node_attributes(G, local_coordinate_system, local_coordinate_system_attr)

    return local_coordinate_system
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def get_splay_vectors(
        G,
        pos_attr='pos',
        splay_attr='splays',
        update_graph=True,
        splay_vector_attr='splay_vectors',
        verbose=1):
    """The splay vectors are calculated by subtracting the node position from the splay position, 
    resulting in a vector that points from the node to the splay.


    Parameters
    ----------
    G : networkx.Graph
        Input graph containing the cave network with nodes and edges.
        Nodes must have `pos_attr` attribute for 3D coordinates and `splay_attr` attribute for splay measurements.

    pos_attr : str, optional
        The name of the node attribute that contains the 3D coordinates, by default 'pos'

    splay_attr : str, optional
        The name of the node attribute that contains the splay measurements, by default 'splays'
        The splay measurements are expected to be a list of one or more sets of 3D coordinates for each node.
        ex: {node_id:[[x1,y1,z1], [x2,y2,z2], ...], ...}

    update_graph : bool, default True
        - if `True`: the splay vectors are set as a node attribute (named \
        `splay_vector_attr` (see below))
        - if `False`: the graph is not updated (`splay_vector_attr` is not used)

    splay_vector_attr : str, default 'splay_vectors'
        The name of the node attribute that will contain the splay vectors,
        (if `update_graph` is `True`), by default 'splay_vectors'

    verbose : int, default: 1
        verbose mode, higher implies more printing (info)

    Returns
    -------
    dict
        Dictionnary with keys: node id, values: list of splay vectors in (x, y, z) coordinates, 
        each key: value item of the dictionary is of the form:
        id: [array([xi, yi, zi]), ...]
    """
    splay_vectors = {} 
    splays = dict(G.nodes(splay_attr))
    keys = list(splays.keys())
    pos3d = nx.get_node_attributes(G, pos_attr)
    for key in keys:
        if splays[key] is None:
            if verbose > 0:
                print(f"Node {key} has no splay data. Skipping.")
            continue

        if len(splays[key]) == 0:
            if verbose > 0:
                print(f"Node {key} has empty splay data. Skipping.")
            continue

        if key not in pos3d:
            if verbose > 0:
                print(f"Node {key} has no position data. Skipping.")
            continue

        splay_vectors[key] = [np.array(splay) - np.array(pos3d[key]) for splay in splays[key]]
 
    if update_graph:
        nx.set_node_attributes(G, splay_vectors, splay_vector_attr)

    return splay_vectors
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def get_splay_vectors_in_local_coordinate_system(
        G, 
        pos_attr='pos', 
        splay_attr='splays',
        update_graph=True,
        local_coordinate_system_attr='local_coordinate_system',
        local_splay_vector_attr='local_splay_vectors',
        verbose=1):
    """Compute the splay vectors in the local coordinate system.

    The local coordinate system is defined by three orthonormal vectors: u, v, w, where w is the conduit direction,
    and u and v are orthogonal vectors in the uv plan perpendicular to w, with u having its 3rd component equal to 0 (horizontal)
    (see `kn.geom.get_local_coordinate_system` function).

    Parameters
    ----------
    G : networkx.Graph
        Input graph containing the cave network with nodes and edges.
        Nodes must have `pos_attr` attribute for 3D coordinates and `splay_attr` attribute for splay measurements.

    pos_attr : str, optional
        The name of the node attribute that contains the 3D coordinates, by default 'pos'

    splay_attr : str, optional
        The name of the node attribute that contains the splay measurements, by default 'splays'
        The splay measurements are expected to be a list of one or more sets of 3D coordinates for each node.
        ex: {node_id:[[x1,y1,z1], [x2,y2,z2], ...], ...}

    update_graph : bool, default True
        - if `True`: the splay vectors are set as a node attribute (named \
        `splay_vector_attr` (see below))
        - if `False`: the graph is not updated (`splay_vector_attr` is not used)

    local_coordinate_system_attr : str, default 'local_coordinate_system'
        The name of the node attribute that will contain the local coordinate system
        (if `update_graph` is `True`), by default 'local_coordinate_system'

    local_splay_vector_attr : str, default 'local_splay_vectors'
        The name of the node attribute that will contain the splay vectors,
        (if `update_graph` is `True`), by default 'splay_vectors'

    verbose : int, default: 1
        verbose mode, higher implies more printing (info)

    Returns
    -------
    local_coordinate_system : dictionary
        Dictionary with keys: node id, values: numpy.ndarray of shape (3, 3)
        representing the local coordinate system at each node, where
        loc_coor_sys[0] = u vector, loc_coor_sys[1] = v vector, loc_coor_sys[2] = w vector, 
        with w being the conduit direction, u and v being orthogonal vectors in the uv plan 
        perpendicular to w, and u having its 3rd component equal to 0 (horizontal).
    
    local_splay_vector : dictionary 
        Dictionnary with keys: node id, values: list of splay vectors in (u, v, w) coordinates, 
        each key: value item of the dictionary is of the form:
        id: [array([ui, vi, wi]), ...]
    """
    # Compute local coordinate systems and splay vectors (dictionaries)
    local_coordinate_system = kn.geom.get_local_coordinate_system(G, pos_attr=pos_attr, update_graph=False, verbose=verbose)
    splay_vectors = kn.geom.get_splay_vectors(G, pos_attr=pos_attr, splay_attr=splay_attr, update_graph=False)

    # Compute splay vectors in the local coordinate system
    local_splay_vectors = {}
    for key in G.nodes():
        if key in local_coordinate_system and key in splay_vectors:
            u, v, w = local_coordinate_system[key]
            local_splay_vectors[key] = [np.array([np.dot(splay_vector, u), np.dot(splay_vector, v), np.dot(splay_vector, w)]) for splay_vector in splay_vectors[key]]

    if update_graph:
        nx.set_node_attributes(G, local_coordinate_system, local_coordinate_system_attr)
        nx.set_node_attributes(G, local_splay_vectors, local_splay_vector_attr)

    return local_coordinate_system, local_splay_vectors
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def calc_conduit_dimensions(
        G, 
        pos_attr='pos',
        splay_attr='splays',
        max_dist_to_plane=None,
        max_angle_with_plane=None,
        set_mean_splay_length_for_intersections=True,
        output='csdim',
        update_graph=True,
        verbose=1):
    """Estimate the cross-section dimensions of the conduits in the cave network represented by the graph G.
    The cross-sectional dimensions are calculated based on the projected splay lengths in the u and v directions, 
    which are orthogonal to the conduit direction at each node (projection onto the uv plane containing the node). 
    
    If the node is an intersection (degree > 2), the cross-sectional dimensions are calculated as the mean of the splay lengths, 
    otherwise, they are calculated as the difference between the maximum and minimum projected splay lengths in the u and v directions.
    
    Different outputs can be requested, including the geometric mean of the cross-section radius, the cross-section width and height. 
    The function returns a dictionary with the requested geometric properties for each node in the graph. 
    The cross-sectional dimensions are calculated as follows:

    - The cross-sectional width is calculated as the difference between the maximum and minimum projected splay lengths in the u direction.
    - The cross-sectional height is calculated as the difference between the maximum and minimum projected splay lengths in the v direction.
    - The geometric mean of the cross-section radius is calculated as the square root of the product of half the cross-section width and half the cross-section height.  

    Some caveats to consider:
    -------------------------
    - Depending on the survey data, some nodes may have no splays or only one splay, which can result in zero cross-sectional dimensions or biased estimates. 
    - The splay measurements may not be uniformly distributed around the node, which can also affect the accuracy of the cross-sectional dimensions. 
    - Some cavers measure splays in the front and back of the station, sometime far away from the station and the projection of those splays may not \
    be representative of the actual cross-section.
    - For nodes with degree greater than 2 (intersections), the cross-sectional dimensions are calculated as \
    the mean of the projected splay lengths, rather than the projections in the u and v directions. 

    Parameters
    ----------
    G : networkx.Graph
        Input graph containing the cave network with nodes and edges.
        Nodes must have `pos_attr` attribute for 3D coordinates and `splay_attr` attribute for splay measurements.

    pos_attr : str, optional
        The name of the node attribute that contains the 3D coordinates, by default 'pos'

    splay_attr : str, optional
        The name of the node attribute that contains the splay measurements, by default 'splays'
        The splay measurements are expected to be a list of one or more sets of 3D coordinates for each node.
        ex: {node_id:[[x1,y1,z1], [x2,y2,z2], ...], ...}

    max_dist_to_plane : float, optional
        A splay is ignored if the distance between its endpoint and the uv plane is greater than this value.
        By defaut (`None`), all splays are considered, regardless of their distance to the uv plane.

    max_angle_with_plane : float, optional
        A splay is ignored if the angle between its vector and the uv plane is greater than this value (in degrees).
        By defaut (`None`), all splays are considered, regardless of their angle with the uv plane.

    set_mean_splay_length_for_intersections : bool, default True
        If `True`, for nodes with degree greater than 2 (intersections), the conduit dimensions are defined by 
        using the mean of splay vectors length in every direction.
        If `False`, for nodes with degree greater than 2 (intersections), the conduit dimensions are not defined.

    output : str {'csdim' (default), 'cs_width', 'cs_height', 'geommean_radius', 'plotting_radius'}
        The type of geometric property to compute and return. Options are:
        
        - 'csdim': cross-sectional dimensions (width and height)
        - 'cs_width': cross-sectional width
        - 'cs_height': cross-sectional height
        - 'geommean_radius': geometric mean of the cross-section radius
        - 'plotting_radius': u_min, u_max, v_min, v_max (max and min projected splay lengths in the u and v directions)

    update_graph : bool, default True
        - if `True`: the computed attribute (see `output`) is set as a node attribute (named `output`)
        - if `False`: the graph is not updated

    verbose : int, default: 1
        verbose mode, higher implies more printing (info)

    Returns
    -------
    dictionnary
        Depending on the 'output' parameter, returns different geometric properties of the conduits:
        
        'csdim': dictionnary 
            Dictionnary with keys: node id, values: [cs_width, cs_height]
            Cross-sectional dimensions, width and height, calculated as the difference 
            between the maximum and minimum projected splay lengths in the u and v directions, respectively.
            ex: csdim = kn.geom.calc_conduit_dimensions(G, output = 'csdim')

        'geommean_radius': dictionnary 
            Dictionnary with keys: node id, values: geometric mean of the cross-section radius calculated as
            np.sqrt((cs_width/2)*(cs_height/2))

        'cs_width': dictionnary
            Cross-sectional width calculated as the difference between the maximum and minimum projected splay lengths in the u direction.
            Dictionnary with keys: node id, values: cross-section width calculated as max value of projected splay in u direction - min value of projected splay in u direction

        'cs_height': dictionnary
            Cross-sectional height calculated as the difference between the maximum and minimum projected splay lengths in the v direction.
            Dictionnary with keys: node id, values: cross-section height calculated as max value of projected splay in v direction - min value of projected splay in v direction

        'plotting_radius': dictionnary
            Dictionnary with keys: node id, values: [u_min, u_max, v_min, v_max] (max and min projected splay lengths in the u and v directions)
    """
    if output not in ['csdim', 'cs_width', 'cs_height', 'geommean_radius', 'plotting_radius']:
        raise ValueError(f"Invalid output value: {output}. Must be one of 'csdim', 'geommean_radius', 'cs_width', 'cs_height'.")

    # Get splay vectors in usual 3D coordinates
    splay_vectors = kn.geom.get_splay_vectors(G, pos_attr=pos_attr, splay_attr=splay_attr, update_graph=False)

    # Compute local coordinate systems and local splay vectors (dictionaries)
    local_coordinate_system, local_splay_vectors = kn.geom.get_splay_vectors_in_local_coordinate_system(
            G, pos_attr=pos_attr, splay_attr=splay_attr, update_graph=False, verbose=verbose)

    if max_dist_to_plane is None:
        max_dist_to_plane = np.inf

    if max_angle_with_plane is None:
        max_angle_with_plane = np.inf
    else:
        max_angle_with_plane = np.deg2rad(max_angle_with_plane) # set angle in radians

    output_dict = {}

    for key in G.nodes():
        if key not in splay_vectors or len(splay_vectors[key]) == 0:
            if verbose > 0:
                print(f'Node {key} has no splays. Conduit dimensions cannot be defined.')
            continue

        if key in local_splay_vectors and key in local_coordinate_system:
            # Get the local splay vectors for the current node
            local_splays = local_splay_vectors[key]
            if len(local_splays) == 0:
                if verbose > 0:
                    print(f'Node {key} has no splays. Conduit dimensions cannot be defined.')
                continue

            # Compute geometric properties based on the local splay vectors
            local_splays_u, local_splays_v, local_splays_w = np.asarray(local_splays).T

            # Filter out splays that are too far from the plane or at too large an angle
            valid_splays = np.all([
                                np.abs(local_splays_w) <= max_dist_to_plane, 
                                np.arctan2(np.abs(local_splays_w), np.sqrt(local_splays_u**2 + local_splays_v**2)) <= max_angle_with_plane
                                ], axis=0)
            if not np.any(valid_splays):
                if verbose > 0:
                    print(f'Node {key} has no valid splays after filtering by distance to plane and angle. Conduit dimensions cannot be defined.')
                continue

            u_min = min(0.0, min(local_splays_u[valid_splays]))
            u_max = max(0.0, max(local_splays_u[valid_splays]))
            v_min = min(0.0, min(local_splays_v[valid_splays]))
            v_max = max(0.0, max(local_splays_v[valid_splays]))
        
        else:
            if G.degree[key] > 2 and set_mean_splay_length_for_intersections:
                if verbose > 0:
                    print(f'Node {key} at intersection (degree {G.degree[key]}). Conduit dimensions are defined by using the mean of splay vectors length in every direction.')
            
                # Compute the mean length of the splay vectors and set dimensions
                mean_splay_length = float(np.round(np.linalg.norm(np.asarray(splay_vectors[key]), axis=1).mean(), 2))

                u_min = - mean_splay_length
                u_max = mean_splay_length
                v_min = - mean_splay_length
                v_max = mean_splay_length

            else:
                if verbose > 0:
                    print(f'Node {key} has no local coordinate system or splay vectors in the local coordinate system. Conduit dimensions cannot be defined.')
                continue

        cs_width = float(np.round(u_max - u_min, 2))
        cs_height = float(np.round(v_max - v_min, 2))

        geommean_radius = 0.5 * np.sqrt(cs_width * cs_height)

        # Set the output dictionary based on the requested output type
        if output == 'csdim':
            output_dict[key] = [cs_width, cs_height]
        elif output == 'cs_width':
            output_dict[key] = cs_width
        elif output == 'cs_height':
            output_dict[key] = cs_height
        elif output == 'geommean_radius':
            output_dict[key] = geommean_radius
        elif output == 'plotting_radius':
            output_dict[key] = [u_min, u_max, v_min, v_max]
                
   
    if update_graph:
        nx.set_node_attributes(G, output_dict, output)

    return output_dict
# ----------------------------------------------------------------------------




#%LRUD TO SPLAYS 
#################################



def calc_splay_from_lrud(station_coordinates, normal_vector_dir,lrud, other_normal=None):
    """Transform left-right-up-down caving survey measurements to splays.
    Calculation made by Nina Egli

    This assumes that:
    1. the normal vector corresponds to the main direction of the conduit
    2. the normal vector is pointing in the direction of survey
    3. left and right measurements are always made in the horizontal plan
    4. up and down are made perpendicular to the direction of the conduit and not strait up or down

    Parameters
    ----------
    station_coordinates : list or array
        x,y,z coordinates of the station in meters or feet. example: [x,y,z]
    normal_vector_dir : list or array
        direction of the conduit. x,y,z vector calculated with .... in the direction of the survey.
    lrud : list or array
        left, right, up, down measurements in the same unit of measurements as the coordinates.

    Returns
    -------
    list of array
        left, right, up, down splays list of coordinates of the position measured on the wall

    """
    #unit vector pointing north
    unit_vector_z = np.array([0,0,1])
    #we normalize the normal vector by its norm    
    normal_vector_plan_normalized =  (normal_vector_dir/np.linalg.norm(normal_vector_dir))
    #extract the left right up down scalar values
    left_scalar = lrud[0]
    right_scalar = lrud[1]
    up_scalar = lrud[2]
    down_scalar = lrud[3]

    #calculate the unit vector to the right and down by taking the vectoriel product
    unit_vector_right = np.cross(normal_vector_plan_normalized,unit_vector_z) 
    #in case of vertical conduit where the normal is perfectly vertical, we use the normal of the previous or next point  
    if other_normal is not None and np.linalg.norm(unit_vector_right) < 0.01:
        unit_vector_right = (np.cross(other_normal,normal_vector_plan_normalized)*np.sign(normal_vector_plan_normalized[2]))

    assert np.linalg.norm(unit_vector_right) > 0.01, f'unit vector right is null: {unit_vector_right}, station coordinates:{station_coordinates}, normal_vector_plan_normalized:{normal_vector_plan_normalized}, normal_vector_dir:{normal_vector_dir}, other_normal:{other_normal}, lrud:{lrud}'

    unit_vector_right = unit_vector_right / np.linalg.norm(unit_vector_right)

    # assert np.abs(1-np.linalg.norm(unit_vector_right))<0.000001, f'unit vector right normalized is not unit vector, instead{np.linalg.norm(unit_vector_right)}'



    unit_vector_down = np.cross(normal_vector_plan_normalized,unit_vector_right) 
    unit_vector_down = unit_vector_down / np.linalg.norm(unit_vector_down)
    # assert np.abs(1-np.linalg.norm(unit_vector_down))<0.000001, f'unit vector down is not unit vector, instead{np.linalg.norm(unit_vector_down)}'

    #calculate position
    left_vector = station_coordinates - unit_vector_right * left_scalar
    right_vector = station_coordinates + unit_vector_right * right_scalar
    up_vector = station_coordinates - unit_vector_down * up_scalar
    down_vector = station_coordinates + unit_vector_down * down_scalar

    return [left_vector.tolist(), right_vector.tolist(), up_vector.tolist(), down_vector.tolist()]
    
#---------------------------------------------------------------


def calc_and_add_lrud_splays(G,df_lrud, direction, other_normal=None, pos_attr='pos', splay_attr='splays', direction_attr='direction'):
    """_summary_

    Parameters
    ----------
    G : networkx graph
        _description_
    df_lrud : pandas dataframe (with columns in any order):
        - 'fulladdress': ...
        - 'previous_fulladdress': name of the station from which the surveyors came from
        - 'following_fulladdres': name of the next station to which the surveyors is going
        - 'left'
        - 'right'
        - 'up'
        - 'down'

    direction : dictionnary
        keys: local node ids 
        values: [x,y,z] direction of the conduit at the point, normalized.
    """
    splays_list = []
    node_key_list = []
    for i,row in df_lrud.iterrows():
        # node_id = find_key_from_fulladdress(G, row.fulladdress)

    
        #CACULATE THE SPLAYS
        #######################
        node_pos = np.array(G.nodes(pos_attr)[row.node_id])  
        lrud = [row.left, row.right, row.up, row.down]
        
        splays_list = splays_list + list(calc_splay_from_lrud(node_pos, direction[row.node_id], lrud, other_normal=other_normal[row.node_id]))
        node_key_list = node_key_list + [row.node_id,row.node_id,row.node_id,row.node_id] #there is always four splays
        
    dict_splays_lrud = kn.utils.list2dict(node_key_list, splays_list)
    nx.set_node_attributes(G, dict_splays_lrud,splay_attr)
    nx.set_node_attributes(G, direction, direction_attr )




def calc_and_add_csdim(G,df_lrud):
    #CALCULATE CROSS-SECTIONAL WITH AND HEIGHT
    ##############################################
    #there is more than one group of lrud per node (since we regroupted them by coordinates)
    #this calculate the average values for left, right, up, down for each unique node.
    df_lrud['cs_width'] = df_lrud['left'] + df_lrud['right']
    df_lrud['cs_height'] = df_lrud['up'] + df_lrud['down']
    # for cases where there is more than one goupe of lrud per nodes, this groupes it by node and make the mean of the width
    #this is a personal choice. we could also take the max or the min...
    df_mean = df_lrud.groupby('node_id')[['cs_width','cs_height']].mean().round(2)

    dict_csdim = dict(zip(list(df_mean.index.values),[list(a) for a in (zip(df_mean['cs_width'],df_mean['cs_height']))]))
    nx.set_node_attributes(G, dict_csdim,'csdim')


#---------------------------------------------


def calc_csdim_from_th(G, therion_file, units='meter'):
    from collections import defaultdict

    lines = open(therion_file).readlines()

    # read and extract lrud values into a dictionnary
    list_dimension = {}
    current_path = []
    in_dimension_section = False

    for line in lines:
        line = line.strip()
        if line.startswith('survey'):
            survey_name = line.split()[1]
            current_path.append(survey_name)
            in_dimension_section = False
        elif line.startswith('data dimensions'):
            in_dimension_section = True
        elif line.startswith('endcentreline'):
            in_dimension_section = False
        elif line.startswith('endsurvey'):
            current_path.pop()
            in_dimension_section = False
        else:
            parts = line.strip().replace('\t',' ').split()
            if in_dimension_section:
                if parts and parts[1][0].isdigit():
                    index = parts[0]
                    values = list(map(float, parts[1:5]))
                    if units == 'feet':
                        values = list(np.array(values)/3.281)
                    path = '.'.join(current_path + [index])
                    list_dimension[path] = values
                else:
                    in_dimension_section = False
    

    #reverse dictionnary to have unique values for each dictionnary entry, since there can be more than one fulladdress per node id
    dict_fulladdress_reverse = { v: k for k, l in nx.get_node_attributes(G,'fulladdress').items() for v in l } 

    # Replace keys in dict1 with values from dict2
    result = defaultdict(list)
    for k, v in list_dimension.items():
        new_key = dict_fulladdress_reverse.get(k, k)
        result[new_key].append(v)

    # Convert defaultdict back to a regular dictionary if needed
    lrud_dict = dict(result)

    #sum values left+ right and up+down, and mean the results for stations with multiple lrud

    csdim = {}
    for key, value in lrud_dict.items():
        if len(value)==1:
            #check that the values are not zero
            if all(item == 0 for item in value)==False:
                cs_width = value[0][0]+value[0][1]
                cs_height = value[0][2]+value[0][3]
                csdim[key] = [cs_width,cs_height]
            else:
                cs_width_list = []
                cs_height_list = []
                for lrud in value:
                    if all(item == 0 for item in value)==False:
                        cs_width_list.append(lrud[0]+lrud[1])
                        cs_height_list.append(lrud[2]+lrud[3])
                cs_width = np.mean(cs_width_list)
                cs_height = np.mean(cs_height_list)
                csdim[key] = [cs_width,cs_height]

    #attach values to the graph    
    return csdim   