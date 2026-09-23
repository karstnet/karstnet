import numpy as np
import networkx as nx
import karstnet as kn

def calc_projected_splay(conduit_direction,splay_direction):
    """Projects a splay (disto measurement from a station to the surrounding wall) 
    in the u,v plan perpendicular to the cave conduit direction.


    Args:
        conduit_direction (numpy.ndarray, list): list or array of [x,y,z] values 
            of the direction vector of the conduit at the point of interest.
            vector origins has to be (0,0,0).
            use function ... to determin the right direction of the vector.
        splay_direction (numpy.ndarray, list): list or array of [x,y,z] values 
            of the direction of splay. Splay values - node.
            vector origins has to be (0,0,0).

    Returns:
        numpy.ndarray: projected splay length in the u and v direction. 
            (splay vector projectedin the uv plan)
            array of two values [u,v]. 
            u is in horizontal x,y plane, and v is orthogonal to u in the uv plan.
    """
    #normalize vector w
    w = np.array(conduit_direction)
    w /= np.linalg.norm(w)
    #find normalized vector u in the horizontal plan, orthogonal to vector w
    if ((np.abs(w[0])>np.abs(w[1])) and (np.abs(w[0])!=0)):
        #checks that x is not equal to zero and that its longer than y
        u = np.array([-w[1]/w[0],1,0])
    elif ((np.abs(w[1])>np.abs(w[0])) and (np.abs(w[1])!=0)):
        u = np.array([1,-w[0]/w[1],0])

    else:
        print('%s and %s are equal to zero. Normalize vector cannot be calculated'%(w[0],w[1]))
        return np.array([0,0])    

    u /= np.linalg.norm(u)
    #find normalized vector v, orthogonal to u, in the uv plan
    v = np.cross(u,w)
    #calculate the length of the projected vector w in the u and v direction
    s_u = np.dot(u,splay_direction)
    s_v = np.dot(v,splay_direction)

    return np.array([s_u,s_v])
    
    



def get_conduit_direction(G):
    keys = list(G.nodes()) #list(dict(G.degree()).keys())
    conduit_directions = {}
    
    # get conduit direction
    
    for i,key in enumerate(keys):
        degree = G.degree()[key]
        # print(i)
        # if i==695:
        #     print(i)
        # calc w at the beginning  and end of a conduit
        if degree == 1:           
            #if node is the start or end of the cave
            neighbors = list(G.neighbors(key))#[n for n in G.neighbors(key)]
            conduit_directions[key] = np.array(G.nodes('pos')[key]) - np.array(G.nodes('pos')[neighbors[0]])   
            if all(v == 0 for v in conduit_directions[key]): print(conduit_directions[key])        

        #calc w in the middle of a conduit
        elif degree == 2:
            #if the node $P_{i}$ is degree 2, $w=P_{i+1}-P_{i-1}$ 
            #conduit_directions[key] = np.array(pos3d[keys[i+1]])-np.array(pos3d[keys[i-1]])
            neighbors = [n for n in G.neighbors(key)]
            #find length of each vector, to transform in unit vector
            pos = np.array(G.nodes('pos')[key])
            pos_before = np.array(G.nodes('pos')[neighbors[0]])
            pos_after = np.array(G.nodes('pos')[neighbors[1]])
            vector_0 = pos - pos_before
            vector_1 = pos_after - pos
            length_0 = np.linalg.norm(vector_0)
            length_1 = np.linalg.norm(vector_1)
            #normalize (to prevent a change in direction)
            conduit_directions[key] = vector_0/length_0 - vector_1/length_1
            #conduit_directions[key] = np.array(G.nodes('pos')[neighbors[0]])/length_0 - np.array(G.nodes('pos')[neighbors[1]])/length_1
            if all(v == 0 for v in conduit_directions[key]): print(conduit_directions[key])
    
        elif degree > 2:
            # no projection for conduit intersections
            conduit_directions[key] = []
        

    return conduit_directions


def get_splay_direction(G):
    splays_direction = {} 
    splays = dict(G.nodes('splays'))
    keys = list(splays.keys())
    pos3d = nx.get_node_attributes(G,'pos')
    for key in keys:
    #calculate the splay direction relative to its node
        #if the node exists. some nodes were survey duplicate, removed from the main graph
        if key in pos3d.keys():
            list_temp = []
            if splays[key]:
                for j in np.arange(len(splays[key])):
                    list_temp.append(np.array(splays[key][j])-np.array(pos3d[key]))
                splays_direction[key] = list_temp

    return splays_direction

def calc_conduit_dimensions(G, output = 'all'):#conduit_directions,splays_direction):
    '''
    FROM SPLAYS
    #possible outputs: 
        'geometric mean': np.sqrt((cs_width/2)*(cs_height/2))
        'cs_dimensions': {[cs_width, cs_height]}
        'projections': project splay in u, projected splay in v
        'all': geometric mean, cs_width, cs_height, project splay in u, projected splay in v
    
    '''
    
    
    
    conduit_directions = get_conduit_direction(G)
    splays_direction = get_splay_direction(G)

    projected_splays_u = {}
    projected_splays_v = {}

    cs_width = {}
    cs_height = {}
    geommean_radius = {}

    keys = list(splays_direction.keys())
    for key in keys:
    #calculate the projected splay length component in u and v direction
        #sometimes there is splays that were duplicates. 
        if key in conduit_directions:
            su=[]
            sv=[]
            if len(conduit_directions[key])>0:
                for splay_direction in splays_direction[key]:       
                    s = calc_projected_splay(conduit_directions[key],splay_direction)
                    # print('s:', s)
                    if s is None:
                        pass
                    else:
                        su.append(s[0])
                        sv.append(s[1])
                # # store su and vs in a dictionnary
                # # this is only interesting for plots...
                projected_splays_u[key] = su
                projected_splays_v[key] = sv   

                # add zero to the array, in case all splays have been taken in the same direction
                su.append(0)
                sv.append(0)
                # print(key, '---', su, sv)
                # # extract the futherst coordinate on each side
                # print(key, '---', max(su), min(su), max(sv), min(sv))
            
                
                cs_width[key]  = float(np.round(max(su)-min(su),2))
                cs_height[key]  = float(np.round(max(sv)-min(sv),2))
                geommean_radius[key] = np.sqrt((cs_width[key]/2)*(cs_height[key]/2))
                
            else:
                #for intersections, calculate the mean of the splays, instead of the projections
                for splay_direction in splays_direction[key]: 
                    su.append(np.linalg.norm(splay_direction))
                    sv.append(np.linalg.norm(splay_direction))
                cs_width[key] = float(np.round(np.mean(su),2))
                cs_height[key] = float(np.round(np.mean(sv),2))
                geommean_radius[key] = np.sqrt((cs_width[key]/2)*(cs_height[key]/2))
    

    if output=='all':
        return geommean_radius, cs_width, cs_height, projected_splays_u, projected_splays_v
    
    elif output == 'geometric mean':
        return geommean_radius
    
    elif output == 'csdim':
        csdim={}
        # print(cs_width)
        for k in cs_width.keys():
            # csdim[k] = list(d[k] for d in [cs_width,cs_height])
            # print(type(cs_width[k]))
            csdim[k] = [cs_width[k], cs_height[k]]
        return csdim
    
    elif output == 'projections':
        return projected_splays_u, projected_splays_v
    

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


def calc_and_add_lrud_splays(G,df_lrud, direction, other_normal=None):
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
        node_pos = np.array(G.nodes('pos')[row.node_id])  
        lrud = [row.left, row.right, row.up, row.down]
        
        splays_list = splays_list + list(calc_splay_from_lrud(node_pos, direction[row.node_id], lrud, other_normal=other_normal[row.node_id]))
        node_key_list = node_key_list + [row.node_id,row.node_id,row.node_id,row.node_id] #there is always four splays
        
    dict_splays_lrud = kn.utils.list2dict(node_key_list, splays_list)
    nx.set_node_attributes(G, dict_splays_lrud,'splays')
    nx.set_node_attributes(G, direction, 'direction' )




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
