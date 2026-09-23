#%%
import geopandas as gpd
import os

#%%

def shpCompass2th(cavename=str, 
                  shapefile_3dsta=str, 
                  shapefile_3dshots=str, 
                  export_folder=str, 
                  unit_dimension='meter',
                  crs_in=None,
                  crs_out=None
                  ):
    # create fix points
    """_summary_

    Parameters
    ----------
    cavename : str
        name of the cave that will be used to create the therion file
    shapefile_3dsta : str
        path to the ...3Dsta.shp file exported from Compass
    shapefile_3dshots : str
        path to the ...3Dshots.shp file exported from Compass
    export_folder : str
        path to the folder where the therion files will be exported
    unit_dimension : str, optional
        units of the left right up down dimensions. by default, Compass exports in feet. 
        make sure you know how it was exported, 
        by default 'feet'
    crs : str, optional
        if it is known, the coordinate system can be defined here. Example: 'EPSG:4326'
        by default, no crs will be added.
    """


    stations = gpd.read_file(shapefile_3dsta)
    shots = gpd.read_file(shapefile_3dshots)


    # fix syntax issue in compass file that won't compile in therion
    stations.STATION = stations.STATION.str.replace(' ','')
    shots.FR_STATION = shots.FR_STATION.str.replace(' ','')
    shots.TO_STATION = shots.TO_STATION.str.replace(' ','')

    stations.STATION = stations.STATION.str.replace('`','-')
    shots.FR_STATION = shots.FR_STATION.str.replace('`','-')
    shots.TO_STATION = shots.TO_STATION.str.replace('`','-')    

    stations.STATION = stations.STATION.str.replace('*',"_star_")
    shots.FR_STATION = shots.FR_STATION.str.replace('*',"_star_")
    shots.TO_STATION = shots.TO_STATION.str.replace('*',"_star_")

    stations.STATION = stations.STATION.str.replace('?',"_question_")
    shots.FR_STATION = shots.FR_STATION.str.replace('?',"_question_")
    shots.TO_STATION = shots.TO_STATION.str.replace('?',"_question_")

    stations.STATION = stations.STATION.str.replace('!',"_exclamation_")
    shots.FR_STATION = shots.FR_STATION.str.replace('!',"_exclamation_")
    shots.TO_STATION = shots.TO_STATION.str.replace('!',"_exclamation_")

    stations.STATION = stations.STATION.str.replace('$','_dollar_')
    shots.FR_STATION = shots.FR_STATION.str.replace('$','_dollar_')
    shots.TO_STATION = shots.TO_STATION.str.replace('$','_dollar_')

    stations.STATION = stations.STATION.str.replace('+','_plus_')
    shots.FR_STATION = shots.FR_STATION.str.replace('+','_plus_')
    shots.TO_STATION = shots.TO_STATION.str.replace('+','_plus_')

    stations.STATION = stations.STATION.str.replace('%','_percent_')
    shots.FR_STATION = shots.FR_STATION.str.replace('%','_percent_')
    shots.TO_STATION = shots.TO_STATION.str.replace('%','_percent_')


    #export fix points
    stations.COMPASS_X = stations.COMPASS_X.astype(str)
    stations.COMPASS_Y = stations.COMPASS_Y.astype(str)
    stations.COMPASS_Z = stations.COMPASS_Z.astype(str)

    #export left right up down
    stations.LEFT = stations.LEFT.astype(str)
    stations.RIGHT = stations.RIGHT.astype(str)
    stations.UP = stations.UP.astype(str)
    stations.DOWN = stations.DOWN.astype(str)
    stations.COMMENT = stations.COMMENT.astype(str)

    stations['fixes'] = 'fix ' + stations[["STATION","COMPASS_X", "COMPASS_Y",  "COMPASS_Z"]].apply(" ".join, axis=1) # y,x,z in this order for OxBelHa
    if crs_in:
        fix_points_list = ['centreline'] + [f'cs {crs_out}'] + stations[["fixes","COMMENT"]].apply(" # ".join, axis=1).str.replace(' # nan','').tolist() + ['endcentreline']
    else:
        fix_points_list = ['centreline'] + stations[["fixes","COMMENT"]].apply(" # ".join, axis=1).str.replace(' # nan','').tolist() + ['endcentreline']
    stations['data_dimensions'] = stations[["STATION", "LEFT", "RIGHT", "UP", "DOWN"]].apply(" ".join, axis=1).str.replace('nan', '0')
    data_dimensions_list = ['centreline',f'units dimension {unit_dimension}','data dimensions station left right up down'] + stations[["data_dimensions","COMMENT"]].apply(" # ".join, axis=1).str.replace(' # nan','').tolist() + ['endcentreline']


    centrelines = []
    #loop through the survey's
    for survey in shots.SURVEY.unique():
        # print(survey)
        subset = shots.loc[shots['SURVEY'] == survey]
        date_string = f'date {subset.YEAR.unique()[0]}.{subset.MONTH.unique()[0]}.{subset.DAY.unique()[0]}'
        comment_string = f'# {subset.SECTION.unique()[0]}'
        nosurvey_list = subset[["FR_STATION", "TO_STATION"]].apply(" ".join, axis=1).tolist()
        centrelines = centrelines + [''] + ['centreline'] + [date_string] + [comment_string] + ['data nosurvey from to'] + nosurvey_list + ['endcentreline'] + ['']
        

    all_lines = [f'survey {cavename}']  + [''] + fix_points_list + [''] +  centrelines + [''] +  data_dimensions_list + [''] + ['endsurvey']

    # Write to a .th file
    # therion_filepath = os.path.join(new_directory,filepath)
    # print(therion_filepath)
    with open(os.path.join(export_folder,f'{cavename}_convert.th'), 'w', encoding='utf-8') as file:
        for line in all_lines:
            file.write(line + '\n')

    if crs_out:
        config_list = [ f'source {cavename}_convert.th',
                        f'cs {crs_out}',
                        f'export database -o {cavename}.sql',
                        f'export model -o {cavename}.3d',
                        f'export model -o {cavename}.lox',
                        f'export model -o {cavename} -fmt esri']
    else:
        config_list = [ f'source {cavename}_convert.th',
                        f'export database -o {cavename}.sql',
                        f'export model -o {cavename}.3d',
                        f'export model -o {cavename}.lox',
                        f'export model -o {cavename} -fmt esri']


    with open(os.path.join(export_folder,f'config.thconfig'), 'w', encoding='utf-8') as file:
        for line in config_list:
            file.write(line + '\n')


#%%


# %%
