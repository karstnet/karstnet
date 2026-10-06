

===================
Karstnet functions 
===================

-------------------------
Import - Export functions
-------------------------
Import/Export functions can be used to load cave survey data into a networkx graph object or directly into a Karsnet object.

Caving survey data exist in various file formats. We only provide import/export tools for certain formats.


1. Therion
----------
Therion is an open sources caving program (<https://therion.speleo.sk>). Two common format from Therion are supported by Karstnet:

- **SQL format by Therion [.th]:** Therion archive format. --> link to the notebook
- **Aven format by Survex [.3d]:** Survex 3D visualization format used in Therion for visualization.  --> link to the notebook


2. KNdata-public
----------------
The KNData-public database archive format are text files organised in two specific structures used to store cave survey data specifically for graph analysis. The two formats are: 

- **yaml format [.yaml]:** This format is comprehensive and can store most of the information about a cave survey. It is the most complete format we provide for the KNdata-public database. Yaml file are structed into dictionaries and lists.
- **csv file format [.csv]:** This format is a simple text format that can be used to store cave survey data in a tabular form. It is less comprehensive than the yaml format but can be easily read and processed by various software tools.

text_files

3. Gocad
--------
- format readable by Gocad (.pl)
- pline

4. GIS
------
- Shapefile



Jason - not working yet



-------------------------
Cleaning functions
-------------------------


-------------------------
Plotting functions
-------------------------








