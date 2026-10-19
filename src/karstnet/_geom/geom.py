import numpy as np
import networkx as nx
import karstnet as kn


#% SPLAYS TO CONDUIT DIMENSIONS


    

def get_conduit_direction(G, pos_attr='pos'):
    """Calculates a vector direction of the conduit at each node in the graph G. 
    The direction is determined based on the positions of neighboring nodes and the degree of the node 
    and based on the sequence of survey by the cavers. The direction is calculated as follows:
    - for nodes with degree 1, the direction is calculated as the vector from the node to its neighbor.
    - for nodes with degree 2, the direction is calculated as the difference between the vectors from the node to its two neighbors.
    - for nodes with degree greater than 2, the direction is set to an empty list. In this case, the direction is not defined, 
    as the node is an intersection of multiple conduits, and the dimension of the conduit will be calculated as the mean of the splays 
    instead of the projection in the uv plan.


    Parameters
    ----------
    G : networkx.Graph
        Input graph containing the cave network with nodes and edges.
        Nodes must have 'pos' attribute for 3D coordinates.

    pos_attr : str, optional
        The name of the node attribute that contains the 3D coordinates, by default 'pos'
        
    Returns
    -------
    dictionnary
        Dictionnary with keys: node id, values: array[x,y,z]) direction
        conduit_directions (numpy.ndarray, list): list or array of [x,y,z] values 
            of the direction vector of the conduit at the point of interest.
            vector origins has to be (0,0,0).
            use function ... to determin the right direction of the vector.
    
    """
    keys = list(G.nodes()) #list(dict(G.degree()).keys())
    conduit_directions = {}
    u_vectors = {}
    v_vectors = {}
    
    # get conduit direction
    pos3d = nx.get_node_attributes(G,pos_attr)

        
    for i,key in enumerate(keys):
        #if the node exists. some nodes were survey duplicate, removed from the main graph
        degree = G.degree()[key]
        # print(i)
        # if i==695:
        #     print(i)
        # calc w at the beginning  and end of a conduit
        if degree == 1:           
            #if node is the start or end of the cave
            neighbors = list(G.neighbors(key))#[n for n in G.neighbors(key)]
            conduit_directions[key] = np.array(G.nodes(pos_attr)[key]) - np.array(G.nodes(pos_attr)[neighbors[0]])   
            if all(v == 0 for v in conduit_directions[key]): print(f'Conduit direction for node {key} is zero: {conduit_directions[key]}')

        #calc w in the middle of a conduit
        elif degree == 2:
            #if the node $P_{i}$ is degree 2, $w=P_{i+1}-P_{i-1}$ 
            #conduit_directions[key] = np.array(pos3d[keys[i+1]])-np.array(pos3d[keys[i-1]])
            neighbors = [n for n in G.neighbors(key)]
            #find length of each vector, to transform in unit vector
            pos = np.array(G.nodes(pos_attr)[key])
            pos_before = np.array(G.nodes(pos_attr)[neighbors[0]])
            pos_after = np.array(G.nodes(pos_attr)[neighbors[1]])
            vector_0 = pos - pos_before
            vector_1 = pos_after - pos
            length_0 = np.linalg.norm(vector_0)
            length_1 = np.linalg.norm(vector_1)
            #normalize (to prevent a change in direction)
            conduit_directions[key] = vector_0/length_0 - vector_1/length_1
            #conduit_directions[key] = np.array(G.nodes(pos_attr)[neighbors[0]])/length_0 - np.array(G.nodes(pos_attr)[neighbors[1]])/length_1
            if all(v == 0 for v in conduit_directions[key]): print(conduit_directions[key])
    
        elif degree > 2:
            # no projection for conduit intersections
            conduit_directions[key] = []

        #calculate the u and v vector for each node, to be used for the projection of the splay in the uv plan
        #orthogonal vectors u and v are calculated in the horizontal plan, perpendicular to the conduit direction w.
        #------------------------------------------------------------------------------------------------------------

        if list(conduit_directions[key]) != []:
            #normalize vector w
            w = np.array(conduit_directions[key])
            w /= np.linalg.norm(w)
            #find normalized vector u in the horizontal plan, orthogonal to vector w
            if ((np.abs(w[0])>np.abs(w[1])) and (np.abs(w[0])!=0)):
                #checks that x is not equal to zero and that its longer than y
                u = np.array([-w[1]/w[0],1,0])
            elif ((np.abs(w[1])>np.abs(w[0])) and (np.abs(w[1])!=0)):
                u = np.array([1,-w[0]/w[1],0])

            else:
                print(f'{w[0]} and {w[1]} are equal to zero. Normalize vector cannot be calculated')
                return np.array([0,0])    

            u /= np.linalg.norm(u)
            u_vectors[key] = u
            #find normalized vector v, orthogonal to u, in the uv plan
            v_vectors[key] = np.cross(u,w)
        

    return conduit_directions, u_vectors, v_vectors

def get_splay_vectors(G,pos_attr='pos',splay_attr='splays'):
    """The projected splay direction is calculated by subtracting the node position from the splay position, 
    resulting in a vector that points from the node to the splay.


    Parameters
    ----------
    G : networkx.Graph
        Input graph containing the cave network with nodes and edges.
        Nodes must have 'pos' attribute for 3D coordinates and 'splays' attribute for splay measurements.

    pos_attr : str, optional
        The name of the node attribute that contains the 3D coordinates, by default 'pos'

    splay_attr : str, optional
        The name of the node attribute that contains the splay measurements, by default 'splays'
        The splay measurements are expected to be a list of one or more sets of 3D coordinates for each node.
        ex: {node_id:[[x1,y1,z1], [x2,y2,z2], ...], ...}

    Returns
    -------
    dictionnary
        Dictionnary with keys: node id, values: list of projected splay [x,y,z] direction in the uv plan.
        id: [[xi,yi,zi], ...]
    """


    splay_vectors = {} 
    splays = dict(G.nodes(splay_attr))
    keys = list(splays.keys())
    pos3d = nx.get_node_attributes(G,pos_attr)
    for key in keys:
    #calculate the splay direction relative to its node
        #if the node exists. some nodes were survey duplicate, removed from the main graph
        if key in pos3d: #key in pos3d.keys():
            list_temp = []
            if splays[key]:
                for j in np.arange(len(splays[key])):
                    list_temp.append(np.array(splays[key][j])-np.array(pos3d[key]))
                splay_vectors[key] = list_temp

    return splay_vectors

def calc_projected_splay_lengths(G, pos_attr='pos', splay_attr='splays'):
    """Projects a splay (disto measurement from a station to the surrounding wall) 
    in the u,v plan perpendicular to the cave conduit direction.


    Args:
        u (numpy.ndarray): The u vector for the uv plan.
        v (numpy.ndarray): The v vector for the uv plan.
        splay_vectors (numpy.ndarray, list): list or array of [x,y,z] values 
            of the direction of splay. Splay values - node.
            vector origins has to be (0,0,0).

    Returns:
        numpy.ndarray: projected splay length in the u and v direction. 
            (splay vector projectedin the uv plan)
            array of two values [u,v]. 
            u is in horizontal x,y plane, and v is orthogonal to u in the uv plan.
    """
    # #normalize vector w
    # w = np.array(conduit_directions)
    # w /= np.linalg.norm(w)
    # #find normalized vector u in the horizontal plan, orthogonal to vector w
    # if ((np.abs(w[0])>np.abs(w[1])) and (np.abs(w[0])!=0)):
    #     #checks that x is not equal to zero and that its longer than y
    #     u = np.array([-w[1]/w[0],1,0])
    # elif ((np.abs(w[1])>np.abs(w[0])) and (np.abs(w[1])!=0)):
    #     u = np.array([1,-w[0]/w[1],0])

    # else:
    #     print(f'{w[0]} and {w[1]} are equal to zero. Normalize vector cannot be calculated')
    #     return np.array([0,0])    

    # u /= np.linalg.norm(u)
    # #find normalized vector v, orthogonal to u, in the uv plan
    # v = np.cross(u,w)
    #calculate the length of the projected vector w in the u and v direction
    conduit_directions, u_vectors, v_vectors = kn.geom.get_conduit_direction(G, pos_attr)
    splay_vectors = kn.geom.get_splay_vectors(G, pos_attr, splay_attr)


    projected_splay_lengths_u = {}
    projected_splay_lengths_v = {}

    # keys = list(splay_vectors.keys())
    for key, splay_vector_list in splay_vectors.items():
    #calculate the projected splay length component in u and v direction
        #sometimes there is splays that were duplicates. 
        # if key in conduit_directions:

        for splay_vector in splay_vector_list: 
            su_len=[]
            sv_len=[]
            if len(conduit_directions[key])>0:
                  
                # s = calc_projected_splay_lengths(conduit_directions[key],splay_vectors)
                su = float(np.dot(u_vectors[key],splay_vector))
                sv = float(np.dot(v_vectors[key],splay_vector))
                # print('s:', s)
                if su is None:
                    pass
                else:
                    su_len.append(float(su))
                    sv_len.append(float(sv))
            else:
                #here we are at an intersection, so we calculate the splay lengths instead of the projections in the uv plan
                su_len.append(np.linalg.norm(splay_vector))
                sv_len.append(np.linalg.norm(splay_vector))      

            # # store su and sv in a dictionnary
            # # this is only interesting for plots...
            projected_splay_lengths_u[key] = su_len
            projected_splay_lengths_v[key] = sv_len    


    # s_u = float(np.dot(u,splay_vectors))
    # s_v = float(np.dot(v,splay_vectors))

    # return np.array([s_u,s_v])
    return projected_splay_lengths_u, projected_splay_lengths_v



def calc_conduit_dimensions(G, output = 'cs_dimensions',pos_attr='pos',splay_attr='splays'):#conduit_directions,splays_direction):
    """Estimate the cross-section dimensions of the conduits in the cave network represented by the graph G.
    The cross-sectional dimensions are calculated based on the projected splay lengths in the u and v directions, 
    which are orthogonal to the conduit direction at each node. 
    
    If the node is an intersection (degree > 2), the cross-sectional dimensions are calculated as the mean of the projected splay lengths, 
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
    - Some cavers measure splays in the front and back of the station, sometime far away from the station and the projection of those splays may not 
    be representative of the actual cross-section.
    - For nodes with degree greater than 2 (intersections), the cross-sectional dimensions are  calculated as 
    the mean of the projected splay lengths, rather than the projections in the u and v directions. 


    Parameters
    ----------
    G : networkx.Graph
        Input graph containing the cave network with nodes and edges.
        Nodes must have 'pos' attribute for 3D coordinates and 'splays' attribute for splay measurements.

    Returns
    -------
    dictionnary
        Depending on the 'output' parameter, returns different geometric properties of the conduits:
        
        'cs_dimensions': dictionnary 
            Dictionnary with keys: node id, values: [cs_width, cs_height]
            Cross-sectional dimensions, width and height, calculated as the difference 
            between the maximum and minimum projected splay lengths in the u and v directions, respectively.
            ex: cs_dimensions = kn.geom.calc_conduit_dimensions(G, output = 'cs_dimensions')

        'geometric mean': dictionnary 
            Dictionnary with keys: node id, values: geometric mean of the cross-section radius calculated as
            np.sqrt((cs_width/2)*(cs_height/2))

        'cs_width': dictionnary
            Cross-sectional width calculated as the difference between the maximum and minimum projected splay lengths in the u direction.
            Dictionnary with keys: node id, values: cross-section width calculated as max value of projected splay in u direction - min value of projected splay in u direction

        'cs_height': dictionnary
            Cross-sectional height calculated as the difference between the maximum and minimum projected splay lengths in the v direction.
            Dictionnary with keys: node id, values: cross-section height calculated as max value of projected splay in v direction - min value of projected splay in v direction


    """
    

    conduit_directions, u_vectors, v_vectors = get_conduit_direction(G, pos_attr)
    splay_vectors = get_splay_vectors(G,pos_attr,splay_attr)
    projected_splay_lengths_u, projected_splay_lengths_v = calc_projected_splay_lengths(G, pos_attr,splay_attr)


    cs_width = {}
    cs_height = {}
    csdim = {}
    geommean_radius = {}


    for key in projected_splay_lengths_u:

        #sometimes there is splays that were duplicates. Verify than the node is in the list of splay vectors, 
        # since some nodes were survey duplicates and removed from the main graph
        if key in splay_vectors.keys():
            if key == 24:
                print(key)
            su_len=projected_splay_lengths_u[key]
            sv_len=projected_splay_lengths_v[key]
            if len(conduit_directions[key])>0:

                # add zero to the array, in case all splays have been taken in the same direction
                su_len.append(float(0))
                sv_len.append(float(0))
        
                # print(key, '---', su, sv)
                # # extract the futherst coordinate on each side
                # print(key, '---', max(su), min(su), max(sv), min(sv))
                            
                cs_width[key]  = float(np.round(max(su_len)-min(su_len),2))
                cs_height[key]  = float(np.round(max(sv_len)-min(sv_len),2))
                csdim[key] = [cs_width[key],cs_height[key]]
                geommean_radius[key] = np.sqrt((cs_width[key]/2)*(cs_height[key]/2))
                
            else:
                #for intersections, calculate the mean of the splays, instead of the projections
                cs_width[key] = float(np.round(np.mean(su_len),2))
                cs_height[key] = float(np.round(np.mean(sv_len),2))
                csdim[key] = [cs_width[key],cs_height[key]]
                geommean_radius[key] = np.sqrt((cs_width[key]/2)*(cs_height[key]/2))

    if output == 'csdim':
        return csdim
    
    elif output == 'geometric_mean':
        return geommean_radius

    elif output == 'cs_width':
        return cs_width

    elif output == 'cs_height':
        return cs_height




    
    

#---------------------------------------------------------------


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
