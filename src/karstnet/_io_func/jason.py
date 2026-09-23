# #EXPORT TO JSON
# def nx2json(G,outputpath):
#     """Export edges to a json file to Gis visualization. Requires the field 'pos' as an attribute of the graph. pos must be a list of three values x,y,z. 

#     Parameters
#     ----------
#     G : netowrkx graph created with therion import
#         _description_
#     outputpath : string
#         full path to folder where the files should be saved. the name of file will be the name of the cave stored in the graph
#     """
    
#     json_text = []
#     #EXPORT TO JSON
#     #create lines of text
#     for i,edge in enumerate(G.edges()):
#         json_text.append('{ "type": "Feature", "properties": { "name": "%s" }, "geometry": { "type": "MultiLineString", "coordinates": [ [ %s, %s ] ] } }'
#                         %(i, G.nodes('pos')[edge[0]], G.nodes('pos')[edge[1]]))

    
#     #write file
#     filepath = (outputpath,'json')
#     with open(filepath + G.graph['cavename'] + '.json', 'w') as f:
#         f.write('\n'.join(json_text))