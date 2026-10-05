

=================
Karstnet functions 
=================

-------------------------
Import - Export functions
-------------------------
Those functions can be used to load cave survey data into a networkx graph object or alternatively directly into a Karsnet object.

Caving survey data exist in various file formats. We only provide import/export tools for certain formats.


1. Therion
Therion is an open sources caving program (<https://therion.speleo.sk>). Two common format have been translated:
- SQL format by Therion is their archive format. --> link to the notebook
- Aven format by Survex is the 3D visualization format. --> link to the notebook


1. KNdata-public
We use homemade formats developped to store cave survey data in the KNData-public database. Those formats are text formats and are used to store cave survey data specifically for graph analysis. The formats are: 
- yaml format. This format is comprehensive and can store most of the information about a cave survey. It is the most complete format we provide for the KNdata-public database. Yaml file are structed into dictionaries and lists.
- csv file format. This format is a simple text format that can be used to store cave survey data in a tabular form. It is less comprehensive than the yaml format but can be easily read and processed by various software tools.

text_files

3. Gocad
- format readable by Gocad (.pl)
- pline

4. GIS
- Shapefile



Jason - not working yet



-------------------------
Cleaning functions
-------------------------


-------------------------
Plotting functions
-------------------------








