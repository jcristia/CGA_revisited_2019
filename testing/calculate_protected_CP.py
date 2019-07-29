#!/usr/bin/env python
"""Calculates amount of 'Conservation Priorities' currently protected by existing conservation areas in the North Shelf Bioregion"""
# Karin Bodtker edited this version (0827) to correct the calculation of spatial overlap to consider if a HU or CP
# should be counted in an MPA
# further edits January 2019: codes used in inclusion matrix are changed; these are hard-coded into modules
#
## Built in modules ##

import os, sys, csv, re

## Third party modules ##

import arcpy
import pandas as pd
#################################
### Configure these variables ###
#################################

### print_status & detailed_status ###
#
# If True various messages on the status of the script will be printed to the standard output
# while it is processing. Nothing will be printed to standard output if this is disabled.
#
# However, if an error occurs or exception is thrown that will still go to stderr. 
#
##

print_status = True
detailed_status = True

### source_mxd ###
#
# Path as string pointing to the mxd containing all of the MPAT layers to be used in the analysis
#
##

source_mxd = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\spatial\CP_HU_MPA_Layers.mxd'

### scaling_attribute & scaling_attribute_file ###
#
# Defines which attributes to use (if necessary) for scaling a layers area based on importance.
# If all layers have one consistent attribute name then set scaling_attribute. If layers have
# varying attribute names then set scaling_attribute_file. You must not set both.
#
# scaling_attribute should be the name of a scaling attribute used by every relevant layer as a
# string. If a layer does not have that attribute present it is assumed that no scaling is
# necessary. If using scaling_attribute_file then set this value to None
#
# scaling_attribute_file should be the path to a CSV file containing two columns. The first
# should be the feature class name or the name of the shapefile without the file extension and
# the second column should be the corresponding scaling attribute. If using scaling_attribute
# (i.e. all relevent layers have one name for scaling layer) then set this value to None.
#
# The scaling attribute used by a layer should be of a numeric data type. If both
# scaling_attribute and scaling_attribute_file are None then it is assumed that no scaling
# should occur.
#
### Karin assumes this is no longer relevant?

scaling_attribute = None   # previously: 'RI'

scaling_attribute_file = None

### working_gdb_folder ###
#
# The path to a folder where temporary files are to be stored. Creates a new geodatabase with
# (hopefully) a unique name. The script deletes the database at the end but if an error is
# is thrown then one might linger. I would suggest cleaning out that folder every once in a
# while. 
#
##

working_gdb_folder = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\spatial\working_TEMP'

### sr_code ###
#
# The wkid for the spatial reference that all data sets will be projected to and calculations
# performed in.
#
##

sr_code = 3005 # NAD_1983_BC_Environment_Albers

### cp_presence_threshold, hu_presence_threshold, & layer_presence_threshold_file ###
#
# These define the threshold percent over which a layer is considerd present within an MPA
# This is defined as a decimal number from 0 - 1 (eg. 0.05 => 5%)
#
# Optionally layer_presence_threshold_file can be defined as the path to a CSV file defining
# these thresholds for each layer. The first column is the name of the layer and the second
# is the threshold percentage as a decimal number between 0 and 1. If a layer is not present
# in the CSV then the appropriate value will be used from cp_presence_threshold or
# hu_presence_threshold. If layer_presence_threshold_file is set to None then only those
# defined variables are used.
#
# This CSV should NOT have a header. Including one may result in an error or undefined
# behaviour.
#
##

cp_presence_threshold = 0.0
hu_presence_threshold = 0.01

layer_presence_threshold_file = None

### mpa_name_fields ###
#
# A list of strings that are possible field names containing MPA names
#
##

mpa_name_fields = ['UID']

### imatrix_path ###
#
# The path to the interaction matrix as a CSV. This CSV must have the the type of
# human use in the first column and conservation priorities in the third column
# such that when non-alphabetic characters (including white space) are stripped out
# and all characters are converted to lower case it matches the third element in
# HU dataset names and the fourth element in CP dataset names, respectively.
#
# For example a row reading like this:
#
# [Bottom Trawl, Demersal sharks/skates, Black skate / sandpaper skate, VERY HIGH]
#
# Would correspond to any HU feature beginning with the following:
#
#   mpatt_hu_bottomtrawl_...
#
# And any CP feature layer beginning with the following:
#
#   mpatt_eco_fish_blackskatesandpaperskate_...
#
# Finally any the values in the fourth column should be any of the following:
#
#   LOW, MEDIUM, HIGH, VERY HIGH
#
# Representing the consequence score.
#
# Finally, this CSV should have a header row. It will be discarded in reading but if it
# doesn't exist then you will miss the first row of your data.
##

imatrix_path = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\input\interactionmatrix_MgmtF2_20190124_withIH.csv'

### output1_path & output2_path ###
#
# Paths to CSV files to be output from this script
#
##

output1_path = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\output\table1.csv'

output2_path = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\output\table2.csv'

### output3_path ###
#
# Path to CSV file to be output from this script
# Used for sliver threshold table
#
##

output3_path = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\output\table3.csv'

### complexFeatureClasses ###
#
# List of strings representing the feature class names of those features that
# are too complex to process as is. These get split into single part features
# which hopefully will solve most of those issues.
#
# For features classes that split by attribute prior to this script, you can just
# include the name up until the "_{value}". You do not need to enter in every instance
# of it.
#
##

complexFeatureClasses = ['eco_coarse_bottompatches_polygons_d', 'eco_coarse_geomorphicunits_polygons_d', 'eco_coarse_coastalclasses_lines_d']

### cleanUpTempData ###
#
# If True then temporary data is deleted after it is used
# set to False if you want to inspect the temp data after the script has run
# This might not be a good idea if you are doing a full run of the script
#
##

cleanUpTempData = True

### inclusion_matrix_path ###
#
# Path to the inclusion matrix file as CSV. The first row starting with the
# second cell should be MPATT HU feature class names. The first column starting
# with the second row should contain MPA names as found in the MPA feature
# classes.
#
# The where these MPAs and HU names intersect should be one of the following:
#
#   Y: This feature should be included in the MPA
#   N: This feature should NOT be included in the MPA
#   U: We are uncertain if the feature should be included in the MPA
#
# Any other value will be treated as blank and the script will determine inclusion
# based on spatial calculations
#
# The behaviour for what the script does when encountering these values can be
# further configured by setting the override_y, override_n, and override_u
# variables below
#

inclusion_matrix_path = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\input\SpatialI_AssessI_20190326.csv'

### override_y & _n & _u ###
#
# These variables set the behaviour for what should be done when testing for
# inclusion and there is a value for the given HU-MPA combo in the inclusion matrix.
#
# If override_y is False then an HU-MPA combo in the inclusion matrix with the value
# Y will be included in an MPA. If True it will allow the script to decide inclusion
# based on the spatial calculations
#
# If override_n is False then an HU-MPA combo in the inclusion matrix with the value
# N will NOT be included in the MPA. If True it will allow the script to decide
# inclusion based on the spatial calculations
#
# /!\ override_u BEHAVES DIFFERENTLY /!\
# If override_u is set to True then an HU-MPA combo in the inclusion matrix with the
# value U will be included in the MPA. If False then it will NOT be included in the
# MPA. If None then it will allow the script to decide inclusion based on the spatial
# calculations
#

override_y = False
override_n = True

override_u = True # This one behaves differently from _y and _n please read above and
                  # be careful when setting

### Conservation Priority area overlap ###
#
# Each conservation priority is intersected with each subgregion and ecosection to get the total area
# of the CP that falls in each. A dictionary is created to hold these values.
# Set to true if you want a completely new dictionary, otherwise 'False' will check if feature exists in dict,
# and if not it will do an interesect and add it.
# For repetitive runs of the CGA where the CP inputs don't change, it will save time to not do an intersect and dissolve
# on every CP on every run.
#

cpOverlap_newDict = False
cpOverlap_DictPath = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\input\cpOverlap_rev.csv'

### Join eco UIDs to table 1 output ###
#
# This requires pandas
#

ecoUIDs_path = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\input\mpatt_eco_UID-simple_20180823.csv'
output1join_path = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\output\table1_joined.csv'

### Table 4 output ###
#
# A list of every cp and hu interaction by mpa
#

output4_path = r'C:\Users\jcristia\Documents\GIS\DFO\Python_Script\CGA_revisted_2019\testing\output\table4.csv'

### Special Cases ###
#
# A dictionary of human uses that are split by activity but use the same spatial dataset
# For instance, there is just one spatial dataset for seine fishing, but in the inclusion matrix
# we would like to distinguish between targeted species allowed in each mpa (e.g. salmon, herring, sardine).
# If an mpa allows more than 1 type, we only want to count this once.
#

hu_multiple = {'hu_co_demersalfishing_trapcom_d': {'variants':
                                                   ['hu_co_demersalfishing_traprec_d_PRAWN',
                                                    'hu_co_demersalfishing_traprec_d_CRAB']},
               'hu_co_pelagicfishing_purseseine_d': {'variants':
                                                   ['hu_co_pelagicfishing_purseseine_d_SALMON',
                                                    'hu_co_pelagicfishing_purseseine_d_HERRING',
                                                    'hu_co_pelagicfishing_purseseine_d_SARDINE']}}


######################
### Implementation ###
######################

###         ###
## Functions ##
###         ###

## fieldExists ##
#
# Checks if a field with a given name exists in the given feature class
#
# Returns True or False
#

def fieldExists(layer, field):
    for field in arcpy.ListFields(layer):
        if field.name == field:
            return True

    return False

## calculateArea ##
#
# Creates a double field in a feature class with the given name and
# populates it with the area of the associated feature
#

def calculateArea(layer, area_field):
    arcpy.AddField_management(layer, area_field, 'DOUBLE')
    arcpy.CalculateField_management(layer, area_field, '!shape.area!', 'PYTHON_9.3')

## calculateTotalArea ##
#
# Returns the total area of a feature class by summing up the values
# contained in the passed field
#

def calculateTotalArea(layer, area_field):
    summed_total = 0
    with arcpy.da.SearchCursor(layer, area_field) as cursor:
        for row in cursor:
            summed_total = summed_total + row[0]

    return summed_total
    
## createMPAdict ##
#
# creates a dictionary of MPA attributes to use as a lookup reference
# this was easier than carrying all these attributes through in the MPA dataset and every subsequent dictionary
# NOTE: this means that for all future data prep, any MPA dataset must have the same schema
#
        
def createMPAdict(source_mxd, merged_name_field):
    # Get MPA layers
    mxd = arcpy.mapping.MapDocument(source_mxd)
    layers = arcpy.mapping.ListLayers(mxd)
    mpa_layers = [lyr for lyr in layers if lyr.isFeatureLayer and lyr.datasetName.startswith('mpatt_mpa')]

    mpa_dict = {}

    for layer in mpa_layers:
        with arcpy.da.SearchCursor(layer, ["UID","NAME_E"]) as mpa_cursor:
            for row in mpa_cursor:
                mpa_dict[row[0]] = {}
                mpa_dict[row[0]] = {'name' : row[1]}

    return mpa_dict

## prepareMPAs ##
#
# Gets MPA layers from mxd by reading the dataset name and merges them all
# into one file (mpa_merged_name). Also projects them into a consistent
# spatial reference (sr_code), puts all the MPA names into one column
# (mpa_name_field), and calculates the area of each (mpa_area_field)
#
        
def prepareMPAs(source_mxd, sr_code, mpa_area_field, mpa_area_attribute_section, final_mpa_fc_name, merged_name_field, mpa_name_fields, mpa_subregion_field, subregions_ALL, ecosections_layer, mpa_marine_area):
    # Get MPA layers
    mxd = arcpy.mapping.MapDocument(source_mxd)
    layers = arcpy.mapping.ListLayers(mxd)
    mpa_layers = [lyr for lyr in layers if lyr.isFeatureLayer and lyr.datasetName.startswith('mpatt_mpa')]

    # Load layers into workspace (and project)
    working_layers = []
    for lyr in mpa_layers:
        arcpy.Project_management(lyr.dataSource, lyr.datasetName,
                                 arcpy.SpatialReference(sr_code))
        working_layers.append(lyr.datasetName)

    # Set up field mappings (need a single consistent name field)
    fm = arcpy.FieldMappings()
    for lyr in working_layers:
        fm.addTable(lyr)

    fmap = arcpy.FieldMap()
    for lyr in working_layers:
        name_field = None
        for field in arcpy.ListFields(lyr):
            if field.name in mpa_name_fields:
                name_field = field.name
                break

        if name_field is None:
            raise ValueError('MPA Layer: {0} does not have field name in mpa_name_fields'.format(lyr))
                
        fmap.addInputField(lyr, name_field)

    nf = fmap.outputField
    nf.name = merged_name_field
    fmap.outputField = nf
    fm.addFieldMap(fmap)

    for field in fm.fields:
        if field.name != merged_name_field and field.name != mpa_marine_area:
            fm.removeFieldMap(fm.findFieldMapIndex(field.name))

    # Perform merge and calculate area field
    arcpy.Merge_management(working_layers, "mpas_merged", fm)

    # Clean up the individual mpa files
    if cleanUpTempData:
        for layer in working_layers:
            arcpy.Delete_management(layer)

    # Determine which subregion each MPA is in
    arcpy.AddField_management("mpas_merged", mpa_subregion_field,"TEXT")
    arcpy.Intersect_analysis(["mpas_merged",subregions_ALL], "mpa_sub_intersect", "NO_FID")
    with arcpy.da.UpdateCursor("mpas_merged", [merged_name_field, mpa_subregion_field]) as cursor_mpa:
        for mpa in cursor_mpa:
            mpa_name = (mpa[0].replace("'", "''")).encode('utf8') # the where clause requires double apostrophes
            where = "{0} = '{1}'".format(merged_name_field, mpa_name)
            with arcpy.da.SearchCursor("mpa_sub_intersect", [merged_name_field, "subregion", "Shape_Area"], where) as cursor_mpasub:
                shpArea = 0.0
                subr = None
                for row in cursor_mpasub:
                    if row[2] > shpArea:
                        shpArea = row[2]
                        subr = row[1]
                mpa[1] = subr
                cursor_mpa.updateRow(mpa)
    arcpy.Delete_management("mpa_sub_intersect")
    
    # changed field name to _TOTAL so this needs to be done before the intersect with ecosections
    # so that the area of the total mpa gets carried forward
    calculateArea("mpas_merged", mpa_area_field)
    # now that we are using just the marine area of the protected area, we should just calculate the total as being the marine area
    arcpy.CalculateField_management("mpas_merged", mpa_area_field, '!' + mpa_marine_area + '!', 'PYTHON_9.3') 

    # intersect mpas and ecosections, then dissolve by mpa and ecosection
    arcpy.Intersect_analysis(["mpas_merged",ecosections_layer], "mpa_ecosect_intersect")
    arcpy.Dissolve_management("mpa_ecosect_intersect", final_mpa_fc_name, [merged_name_field, "ecosection"],[[mpa_subregion_field, "FIRST"],[mpa_area_field, "FIRST"]],"MULTI_PART")
    # Rename fields to get rid of the labels the dissolve appends to the beginning
    renameField(final_mpa_fc_name, 'FIRST_' + mpa_subregion_field, mpa_subregion_field, 'TEXT')
    renameField(final_mpa_fc_name, 'FIRST_' + mpa_area_field, mpa_area_field, 'DOUBLE')   
    # JC: add area field to calculate the area of each piece of an mpa in overlapping ecosections
    calculateArea(final_mpa_fc_name, mpa_area_attribute_section)
    # clean up merge and intersect datasets
    arcpy.Delete_management("mpas_merged")
    arcpy.Delete_management("mpa_ecosect_intersect")

    return final_mpa_fc_name

## buildScalingDict ##
#
# Reads a two column CSV where the first column is
# an mpatt dataset name and the second is the name
# of the field that contains the scaling factor used
# to scale a features area in that dataset before it
# is determined whether or not it exists within an MPA
#
# Returns a dict where the dataset name points to the
# scaling factor
#

def buildScalingDict(scaling_attribute_file):
    scaling_fields = {}

    with open(scaling_attribute_file, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            fc_name = row[0]
            scaling_attribute = row[1]

            scaling_fields[fc_name] = scaling_attribute

    return scaling_fields

## loadLayer ##
#
# Copies a layer into the temporary workspace, reprojects it, calculates the area
# and creates a consistently named scaling attribute
#
# If is_complex == True the layer is split to single part features before other
# calculations are performed. This will often make processes work that would fail
# otherwise
#
# Returns the output layer name which can be different
#

def loadLayer(source_mxd, layer_name, sr_code, new_bc_area_field, new_bc_total_area_field,
              scaling_dict, scaling_attribute, new_scaling_field, is_complex, density_field, value_type):
    if detailed_status:
        print 'Loading ' + layer_name
    
    # Find layer in mxd
    mxd = arcpy.mapping.MapDocument(source_mxd)
    layer = arcpy.mapping.ListLayers(mxd, layer_name)[0]

    # Load layer into workspace and reproject
    working_layer = layer.datasetName
    orig_name = working_layer

    arcpy.env.XYResolution = "0.0001 Meters" # project will fail if a resolution is set too low
    arcpy.Project_management(layer.dataSource, layer.datasetName,
                             arcpy.SpatialReference(sr_code))
    
    if is_complex:
        if detailed_status:
            print '...Exploding complex feature class'
        arcpy.MultipartToSinglepart_management(working_layer, working_layer+'_c')

        # Clean up
        if cleanUpTempData:
            arcpy.Delete_management(working_layer)
        
        working_layer = working_layer+'_c'

    # Create consistent scaling field
    arcpy.AddField_management(working_layer, new_scaling_field, 'DOUBLE')

    # If layer is in scaling_dict use that attribute for scaling_attribute
    if scaling_dict is not None and working_layer in scaling_dict:
        scaling_attribute = scaling_dict[working_layer]

    # If scaling_attribute is set and it exists in the feature class
    # copy into new_scaling_field
    if scaling_attribute is not None and len(
            arcpy.ListFields(working_layer, scaling_attribute)) >= 1:
        arcpy.CalculateField_management(working_layer, new_scaling_field,
                                        '!{0}!'.format(scaling_attribute), 'PYTHON_9.3')
    else: # If no scaling then just use 1 for scaling
        arcpy.CalculateField_management(working_layer, new_scaling_field,
                                        '1', 'PYTHON_9.3')
    
    keep_fields = [new_scaling_field, 'ecosection', density_field]

    # Delete fields that aren't important
    for field in arcpy.ListFields(working_layer):
        # Don't delete Object ID or Geometry 
        if field.type in ['OID','Geometry']:
            continue

        # Don't try to delete required fields
        if field.required:
            continue

        # Don't delete fields in keep_fields
        if field.name not in keep_fields:
            arcpy.DeleteField_management(working_layer, field.name)

    # Add new area fields
    arcpy.AddField_management(working_layer, new_bc_area_field, "DOUBLE")
    arcpy.AddField_management(working_layer, new_bc_total_area_field, "DOUBLE")

    # Calculate feature area and total area
    calculateArea(working_layer, new_bc_area_field)
    # calculate total based on if it is a density or area based feature
    # we still need new_bc_area_field to be area at this point, but the total needs to be changed
    if value_type == 'density':
        total_area = calculateTotalArea(working_layer, density_field)
    else:
        total_area = calculateTotalArea(working_layer, new_bc_area_field)
    arcpy.CalculateField_management(working_layer, new_bc_total_area_field, total_area, 'PYTHON_9.3')

    if working_layer != orig_name:
        if arcpy.Exists(orig_name):
            arcpy.Rename_management(orig_name, orig_name + '_original')
        arcpy.Rename_management(working_layer, orig_name)
        
    return orig_name

## loadRegionLayer ##
#
# Similar to loadLayer above but doesn't add a scaling attribute (not needed)
# 

def loadRegionLayer(source_mxd, layer_name, sr_code, new_bc_area_field, new_bc_total_area_field):
    if detailed_status:
        print 'Loading ' + layer_name
        
    # Find layer in mxd
    mxd = arcpy.mapping.MapDocument(source_mxd)
    layer = arcpy.mapping.ListLayers(mxd, layer_name)[0]

    # Load layer into workspace and reproject
    working_layer = layer.datasetName
    arcpy.Project_management(layer.dataSource, layer.datasetName,
                             arcpy.SpatialReference(sr_code))

    # Add new area fields
    arcpy.AddField_management(working_layer, new_bc_area_field, "DOUBLE")
    arcpy.AddField_management(working_layer, new_bc_total_area_field, "DOUBLE")

    # Calculate feature area and total area
    calculateArea(working_layer, new_bc_area_field)
    total_area = calculateTotalArea(working_layer, new_bc_area_field)
    arcpy.CalculateField_management(working_layer, new_bc_total_area_field, total_area, 'PYTHON_9.3')

    return working_layer

## buildThresholdDict ##
#
# Reads a csv where the first column is an mpatt dataset name
# and the second is the threshold used to determine whether it
# is present in an MPA or not
#
# Returns a dict where the key is an mpatt dataset name that
# points to the threshold for that layer
#

def buildThresholdDict(layer_presence_threshold_file):
    thresholds = {}

    with open(layer_presence_threshold_file, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            fc_name = row[0]
            layer_threshold = row[1]

            thresholds[fc_name] = layer_threshold

    return thresholds


## renameField ##
#
# Adds a new field with the desired name, copies values from an old
# field with an undesireable name, and deletes that old field
#

def renameField(working_layer, ifield, ofield, ftype):
    arcpy.AddField_management(working_layer, ofield, ftype)
    arcpy.CalculateField_management(working_layer, ofield, '!{0}!'.format(ifield), 'PYTHON_9.3')
    arcpy.DeleteField_management(working_layer, ifield)

## readMPAInclusionMatrix ##
#
# Reads a CSV structured like a matrix where the first row contains mpatt feature class names
# and the first column contains MPA names matching those in the MPA feature classes. The values
# are any of the following: (the values below are outdated)
#
# Y: This feature should be included in the MPA
# N: This feature should NOT be included in the MPA
# U: We are uncertain if the feature should be included in the MPA
#
# Any other value will be treated as if it were left a blank
#
# Should only be used for HU features

def readMPAInclusionMatrix(mpath):
    inclusion_matrix = {}
    with open(mpath, 'rb') as csvfile:
        reader = csv.reader(csvfile)

        # Read header first
        header = reader.next()

        for row in reader:
            # Get mpa from first column
            mpa = row[0]
            inclusion_matrix[mpa] = {}

            # Iterate through row skipping first col
            for i in range(1, len(row)):
                # Get feature class name from header
                fc_name = header[i]
                # Set inclusion value to whatever value is in the file unless its blank
                inclusion_matrix[mpa][fc_name] = row[i].strip() if row[i].strip() in ('O', 'X', 'C','na') else None
                
    return inclusion_matrix

## shouldInclude ##
#
# Tests to see if a feature belongs in an MPA guided by the inclusion matrix
#

def shouldInclude(pct_in_mpa, threshold, im, fc, mpa):
    # If mpa is not in inclusion matrix use conventional inclusion test
    if mpa not in im:
        return pct_in_mpa > threshold

    # If fc not in inclusion matrix use conventional inclusion test
    if fc not in im[mpa] and fc not in hu_multiple:
        return pct_in_mpa > threshold

    # if fc has a variation, then do a separate set of tests
    if fc in hu_multiple:  #JC 20190703: TO DO. This is just to get the ival, then do condition testing.
        #print fc + " is a hu_multiple"
        ival_list = []
        for variant in hu_multiple[fc]['variants']:
            ivaltemp = im[mpa][variant]
            ival_list.append(ivaltemp)
        #print ival_list
        if not ival_list: #if list is emtpy
            return pct_in_mpa > threshold

        # If hu is permitted or not restricted and override is disabled then include it
        if ('O' in ival_list or 'C' in ival_list) and not override_y:
            return True
        # See above but for restricted
        if ('X' in ival_list or 'na' in ival_list) and not override_n:
            return False
        # Otherwise override whatever value is in i_val with conventional test
        return pct_in_mpa > threshold

    # Get inclusion value
    i_val = im[mpa][fc]

    # If inclusion value was blank use conventional inclusion test
    if i_val is None:
        return pct_in_mpa > threshold

    # If inclusion value is permitted or not restricted and override is disabled then include it
    if i_val in ['O', 'C'] and not override_y:
        return True

    # See above but for restricted
    if i_val in ['X', 'na'] and not override_n:
        return False

    # If value is U do what override asks
 #   if i_val in ['U', '?', 'Y?'] and override_u is not None:
  #      return override_u

    # Otherwise override whatever value is in i_val with conventional test
    return pct_in_mpa > threshold

## process_geometry ##
#
# Does the geometric heavy lifting. Determines what area of a feature falls in each MPA.
# Also performs a subregional analysis if a region_layer is provided.
#
# Returns a feature class with the adjusted and clipped area of the feature class
# and the original area of the feature class in addition to information about the enclosing
# MPA and percentages comparing the clipped/adjusted size to the MPA size and the original
# size of the fc
#
    
def process_geometry(base_layer, final_mpa_fc_name, clipped_adjusted_area, scaling_attribute,
                     mpa_name_attribute, mpa_area_attribute, new_bc_total_area_field,
                     pct_of_mpa_field, pct_of_total_field, mpa_subregion_field, mpa_area_attribute_section,
                     clipped_adj_area_mpaTotal, pct_of_mpa_field_Total, density_field, value_type):
            
    working_intersect = base_layer + '_Intersect'

    # Intersect with by MPAs and explode to singlepart
    if detailed_status:
        print '...Intersecting ' + base_layer
    arcpy.Intersect_analysis([base_layer, final_mpa_fc_name], working_intersect)

    # Calculate area after clipping and adjust it by the scaling factor
    # Do this before dissolving because otherwise you can't capture scaling factors
    # or overlapping area
    arcpy.AddField_management(working_intersect, clipped_adjusted_area, 'DOUBLE')
    arcpy.CalculateField_management(working_intersect, clipped_adjusted_area,
                                    '!shape.area!*!{0}!'.format(scaling_attribute), 'PYTHON_9.3')

    # if it is a density/diversity based feature, check if cell was clipped and rescale density value
    if value_type == 'density':
        with arcpy.da.UpdateCursor(working_intersect, [density_field, clipped_adjusted_area, new_bc_area_field]) as cursor:
            for row in cursor:
                if row[1] != row[2]:
                    newValue = (row[1]/row[2]) * row[0]  # newValue = (newarea/oldarea) * value
                    row[1] = newValue
                    row[2] = row[0] # make the feature area the original density value
                else:
                    row[1] = row[0]
                    row[2] = row[0]
                cursor.updateRow(row)

    # 2080507 This is where the feature count functionality was.
    # It was taking way too long to process. To separate out by mpa and ecosection it requires reseting the cursor many times
    #  - almost as many times as there are features, so for datasets with 19,000 features, this becomes very slow.
    # The previous version maintains this code. The version before the previous version has the code where it only splits
    # out by mpa. If it ends up where I still need feature count, then I should revert back to the verison from two versions ago.

    # Dissolve by mpa_name_attribute field summing adjusted area
    working_dissolved = base_layer + '_Dissolved'
    
    if detailed_status:
        print '...Dissolving ' + base_layer
    arcpy.Dissolve_management(working_intersect, working_dissolved, [mpa_name_attribute, "ecosection"],
                              [[clipped_adjusted_area, 'SUM'], [new_bc_total_area_field, 'FIRST'],
                               [mpa_area_attribute, 'FIRST'], [mpa_area_attribute_section, 'FIRST'],
                              [mpa_subregion_field, 'FIRST']])

    # Rename fields for simplicities sake
    renameField(working_dissolved, 'SUM_' + clipped_adjusted_area, clipped_adjusted_area, 'DOUBLE')
    renameField(working_dissolved, 'FIRST_' + new_bc_total_area_field, new_bc_total_area_field, 'DOUBLE')
    renameField(working_dissolved, 'FIRST_' + mpa_area_attribute, mpa_area_attribute, 'DOUBLE')
    renameField(working_dissolved, 'FIRST_' + mpa_area_attribute_section, mpa_area_attribute_section, 'DOUBLE')
    renameField(working_dissolved, 'FIRST_' + mpa_subregion_field, mpa_subregion_field, 'TEXT')

    
    arcpy.AddField_management(working_dissolved, pct_of_mpa_field, 'DOUBLE')
    arcpy.AddField_management(working_dissolved, clipped_adj_area_mpaTotal, 'DOUBLE')
    arcpy.AddField_management(working_dissolved, pct_of_mpa_field_Total, 'DOUBLE')
    arcpy.AddField_management(working_dissolved, pct_of_total_field, 'DOUBLE')

    # Calculate percentages
    arcpy.CalculateField_management(working_dissolved, pct_of_mpa_field,
                                    '!{0}!/!{1}!'.format(clipped_adjusted_area,mpa_area_attribute),
                                    'PYTHON_9.3')
    arcpy.CalculateField_management(working_dissolved, pct_of_total_field,
                                    '!{0}!/!{1}!'.format(clipped_adjusted_area,new_bc_total_area_field),
                                    'PYTHON_9.3')
    
    # calculate new total fields
    # get unique list of mpas
    mpa_list = []
    with arcpy.da.SearchCursor(working_dissolved, [mpa_name_attribute]) as cursor:
        for row in cursor:
            if row[0] not in mpa_list:
                mpa_list.append(row[0])
    # add up feature spatial areas by mpa (Parts in ecosections?)
    for mpa in mpa_list:
        mpa_name = (mpa.replace("'", "''")).encode('utf8') # the where clause requires double apostrophes
        where = "{0} = '{1}'".format(mpa_name_attribute, mpa_name)
        # karin changed the calc here to shape area to get areal overlap stat
        # JC followup notes: this makes sense since the mpaTotal fields (clipped_adj_area_mpaTotal and pct_of_mpa_field_Total)
        # should be area-based for what they are used for. Only the pct_of_mpa_field_total gets referenced again
        # when dealing with slivers (which needs to be area based). The "adj" part of the name is a bit misleading
        # but we will leave it as is for now.
        #with arcpy.da.UpdateCursor(working_dissolved, [clipped_adjusted_area, clipped_adj_area_mpaTotal], where) as cursor:
        with arcpy.da.UpdateCursor(working_dissolved, ['SHAPE@AREA', clipped_adj_area_mpaTotal],
                                   where) as cursor:
            sum_area = 0.0
            for row in cursor:
                sum_area += row[0]
            cursor.reset()
            for row in cursor:
                row[1] = sum_area
                cursor.updateRow(row)
# karin changed the calc here to shape area to get areal overlap stat
    #arcpy.CalculateField_management(working_dissolved, pct_of_mpa_field_Total,
     #                               '!shape.area!/!{0}!'.format(mpa_area_attribute),
      #                              'PYTHON_9.3')
    arcpy.CalculateField_management(working_dissolved, pct_of_mpa_field_Total,
                                    '!{0}!/!{1}!'.format(clipped_adj_area_mpaTotal, mpa_area_attribute),
                                    'PYTHON_9.3')


    # Clean up
    if cleanUpTempData:
        for layer in arcpy.ListFeatureClasses(base_layer + '_*'):
            if layer != working_dissolved:
                arcpy.Delete_management(layer)
        
    return working_dissolved

## calculate_presence ##
#
# Returns a dict with various info about the passed layers presence within each
# MPA. The first key in the dict is the mpa name which points to another dict.
# This dict has the following keys with the following values:
#
#     'clip_area'     -> The adjusted and clipped area of the layer in the mpa
#     'orig_area'     -> The original area of the layer
#     'mpa_area'      -> The area of the MPA
#     'region_area'   -> The area of the clipped region (if applicable otherwise None)
#     'pct_in_mpa'    -> clip_area / mpa_area
#     'pct_of_region' -> clip_area / region_area (if applicable otherwise None)
#     'pct_of_total'  -> clip_area / orig_area

def calculate_presence(working_layer, final_mpa_fc_name, clipped_adjusted_area,
                       pct_of_total_field, pct_of_mpa_field, mpa_name_attribute,
                       scaling_attribute, threshold, subregion, imatrix, mpa_subregion_field, mpa_area_attribute_section,
                       clipped_adj_area_mpaTotal, pct_of_mpa_field_Total, density_field, value_type):
    mpas = {}
    sliver_freq = {} # to get sliver frequencies

    # JC: this line doesn't make sense. For one, it will never be none, since if there is no
    # subregion code then 'region' gets passed to subregion.
    # Also, why is it recalculating this? We already know the total area. Is it just easier than
    # grabbing the value of the first field? (this is probably it - it needs that value to write
    # to the dictionary).
    # Also, it should not pass 'new_bc_total_are_field'
    # Lastly, even though region_area is used in some calculations, ultimately, none of those results
    # are ever outputted anywhere.
    # JC: changing the second variable from 'new_bc_total_area_field' to 'new_bc_area_field'
    region_area = calculateTotalArea(working_layer,
                                     new_bc_area_field) if subregion is not None else None

    # Crunch the geometry for the whole region
    processed_layer = process_geometry(working_layer, final_mpa_fc_name, clipped_adjusted_area,
                                       scaling_attribute, mpa_name_attribute, mpa_area_attribute,
                                       new_bc_total_area_field, pct_of_mpa_field, pct_of_total_field,
                                       mpa_subregion_field, mpa_area_attribute_section, clipped_adj_area_mpaTotal,
                                      pct_of_mpa_field_Total, density_field, value_type)

    # Read the statistics for the whole region into a dict
    with arcpy.da.SearchCursor(
            processed_layer,
            [mpa_name_attribute,mpa_area_attribute,pct_of_mpa_field,pct_of_total_field,
             new_bc_total_area_field,clipped_adjusted_area, "ecosection", mpa_subregion_field,
             mpa_area_attribute_section, clipped_adj_area_mpaTotal, pct_of_mpa_field_Total]
    ) as cursor:
        # i.e. for each mpa which technically has hu/cp in it
        for row in cursor:
            mpa_name, mpa_area, pct_of_mpa = row[0], row[1], row[2]
            pct_of_total, hucp_og_area, hucp_clip_area = row[3], row[4], row[5]
            ecosect, subreg_mpa, mpa_area_ecosect  = row[6], row[7], row[8]
            hucp_clip_area_mpaTotal, pct_of_mpa_Total = row[9], row[10]

            # each HU/CP is a dict w/ info on its name, clipped area, and total area of
            # the original layer
            #
            # Checks if hu/cp makes up greater than 5% (or whatever) of mpa
            datasetname = working_layer if subregion is not None else '_'.join(working_layer.split('_')[:-1])
            
            if mpa_name not in sliver_freq:
                sliver_freq[mpa_name] = {'pct_overlap_cphu_mpa': pct_of_mpa_Total}
                # this should only need to be written once, even if there are multiple features for each mpa

            if shouldInclude(pct_of_mpa_Total, threshold, imatrix, datasetname, mpa_name):
                pct_of_region = (hucp_clip_area / region_area) if region_area is not None else None
                if mpa_name not in mpas:
                    mpas[mpa_name] = {}
                mpas[mpa_name][ecosect] = {'subregion': subreg_mpa,
                                  'clip_area': hucp_clip_area,
                                  'orig_area': hucp_og_area,
                                  'mpa_area': mpa_area,
                                  'region_area': region_area,
                                  'pct_in_mpa': pct_of_mpa,
                                  'pct_of_region': pct_of_region,
                                  'pct_of_total': pct_of_total}
              
    # Clean up workspace
    if cleanUpTempData:
        arcpy.Delete_management(working_layer)
        arcpy.Delete_management(processed_layer)
            
    return mpas, sliver_freq



##
# Build CP area overlap dictionary from csv file
##

def buildOverlapDict(cpOverlap_DictPath, cp_area_overlap_dict):

    with open(cpOverlap_DictPath, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        reader.next()
        for row in reader:
            fc_name = row[0]
            section = row[1]
            area = float(row[2])

            if fc_name not in cp_area_overlap_dict:
                cp_area_overlap_dict[fc_name] = {}
            if section not in cp_area_overlap_dict[fc_name]:
                cp_area_overlap_dict[fc_name][section] = {}

            cp_area_overlap_dict[fc_name][section] = {'Area': area}

    return cp_area_overlap_dict

##
# Intersect cp layers with subregions/ecosections
# This is done to get the total area OR VALUE of a cp within each subregion/ecosection
##

def calcCPlyrOverlap(cp_area_overlap_dict, working_layer, ecosections_layer, subregions_ALL, density_field, value_type):
    
    # intersect
    subr_union = working_layer + '_subUnion'
    ecos_union = working_layer + '_ecoUnion'

    # Union with subregions-ecosections
    # We don't want overlapping parts to be combined so only a union works
    if detailed_status:
        print '...Union ' + working_layer
    arcpy.Union_analysis([working_layer, subregions_ALL], subr_union)
    arcpy.Union_analysis([working_layer, ecosections_layer], ecos_union)

    # delete records that do not overlap
    FID_wlyr = "FID_" + working_layer
    if len(FID_wlyr) > 64:
        cut = len(FID_wlyr) - 64
        FID_wlyr = FID_wlyr[:-cut]
    FID_subr = "FID_" + subregions_ALL
    FID_ecos = "FID_" + ecosections_layer
    with arcpy.da.UpdateCursor(subr_union, [FID_wlyr, FID_subr]) as cursor:
        for row in cursor:
            if row[0] == -1 or row[1] == -1:
                cursor.deleteRow()
    with arcpy.da.UpdateCursor(ecos_union, [FID_wlyr, FID_ecos]) as cursor:
        for row in cursor:
            if row[0] == -1 or row[1] == -1:
                cursor.deleteRow()

    # add area field that will be summed when dissolving
    ecosub_area_field = 'ecosub_area'
    arcpy.AddField_management(subr_union, ecosub_area_field, "DOUBLE")
    arcpy.AddField_management(ecos_union, ecosub_area_field, "DOUBLE")
    arcpy.CalculateField_management(subr_union, ecosub_area_field, '!shape.area!', 'PYTHON_9.3')
    arcpy.CalculateField_management(ecos_union, ecosub_area_field, '!shape.area!', 'PYTHON_9.3')

    # if it is a density/diversity based feature, recalculate area field in case it was clipped
    if value_type == 'density':
        with arcpy.da.UpdateCursor(subr_union, [density_field, ecosub_area_field, new_bc_area_field]) as cursor:
            for row in cursor:
                if row[1] != row[2]:
                    newValue = (row[1]/row[2]) * row[0]  # newValue = (newarea/oldarea) * value
                    row[1] = newValue
                else:
                    row[1] = row[0]
                cursor.updateRow(row)

    if value_type == 'density':
        with arcpy.da.UpdateCursor(ecos_union, [density_field, ecosub_area_field, new_bc_area_field]) as cursor:
            for row in cursor:
                if row[1] != row[2]:
                    newValue = (row[1]/row[2]) * row[0]  # newValue = (newarea/oldarea) * value
                    row[1] = newValue
                else:
                    row[1] = row[0]
                cursor.updateRow(row)

    # Dissolve by ecosection/subregion field summing ecosub_area_field
    subr_dissolved = working_layer + '_subDissolved'
    ecos_dissolved = working_layer + '_ecoDissolved'
    
    if detailed_status:
        print '...Dissolving ' + working_layer
    arcpy.Dissolve_management(subr_union, subr_dissolved, ['subregion'], [[ecosub_area_field, 'SUM']])
    arcpy.Dissolve_management(ecos_union, ecos_dissolved, ['ecosection'], [[ecosub_area_field, 'SUM']])

    cp_area_overlap_dict[working_layer] = {}
    with arcpy.da.SearchCursor(subr_dissolved, ['subregion', "SUM_" + ecosub_area_field]) as cursor:
        for row in cursor:
            cp_area_overlap_dict[working_layer][row[0]] = {'Area' : row[1]}
    with arcpy.da.SearchCursor(ecos_dissolved, ['ecosection', "SUM_" + ecosub_area_field]) as cursor:
        for row in cursor:
            cp_area_overlap_dict[working_layer][row[0]] = {'Area' : row[1]}

    # delete
    if cleanUpTempData:
        arcpy.Delete_management(subr_union)
        arcpy.Delete_management(ecos_union)
        arcpy.Delete_management(subr_dissolved)
        arcpy.Delete_management(ecos_dissolved)

    return cp_area_overlap_dict




## calcEffectivenessScore ##
#
# Takes the number of interactions for a cp broken down by severity and spits out an effectiveness score
#

def calcEffectivenessScore(num_high, num_mod, num_low):
    if num_high > 0 or num_mod > 4: # High
        return 0.0
    elif num_mod == 4: # Moderate-High
        return 0.24
    elif num_mod == 3: # Moderate
        return 0.6
    elif num_mod > 0 and num_mod < 3: # Low impact
        return 0.85
    else: # Negligible
        return 1.0


## countInteractions ##
#
# Counts up the interactions for easy reading into calcEffectivenessScore
#

def countInteractions(i_list):
    num_high = 0
    num_mod = 0
    num_low = 0
    
    for interaction in i_list:
        if interaction == 'HIGH':
            num_high = num_high + 1
        elif interaction == 'MODERATE':
            num_mod = num_mod + 1
#        elif interaction == 'LOW':
#           num_low = num_low + 1
            
    return (num_high, num_mod, num_low)

## loadInteractionsMatrix ##
#
# Reads the interactions matrix from a csv file. An external library could be used
# to read it from the Excel workbook but that seemed unnecessary for this task.
#
# Returns one of the more straightforward dicts of this script. The first key is
# the CP dataset name and the second key is the HU dataset name. The value is the
# interaction severity
# 

def loadInteractionsMatrix(imatrix_path):
    imatrix = {}
    
    with open(imatrix_path, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        reader.next()

        regex = re.compile('[^a-zA-Z]')

        for row in reader:
            # Convert matrix entries to lower case and remove spaces
            # and non-alphabetic ccharacters. In theory this will make
            # them match up to the layer dataset names
            #
            # eg. Tufted puffins -> tuftedpuffins unfortunately the
            # mpatt dataset is called mpatt_eco_birds_tuftedpuffin_colonies_seasketch
            # without the s
            #

            #cp = regex.sub('', row[2]).lower()
            cp = '_'.join(row[1].split('_')[2:4])  # we are now referencing the 3rd and 4th parts of the UID
            hu = regex.sub('', row[3]).lower()
            interaction = row[5]

            # The file uses a different conventions than the docs
            # I was working off of
            if interaction == 'VERY HIGH' or interaction == 'Major Negative':
                interaction = 'HIGH'
                
            if interaction == 'MEDIUM' or interaction == 'Minor Negative':
                interaction = 'MODERATE'

#            if interaction == 'Negligible':
#                interaction = 'LOW'

            if cp not in imatrix:
                imatrix[cp] = {}

            imatrix[cp][hu] = interaction         

    return imatrix

## determineInteraction ##
#
# Get the relevent part of the dataset names and return the interaction
# from the imatrix
#

def determineInteraction(imatrix, cp, hu):
    cp = '_'.join(cp.split('_')[2:4])
    length = len(cp)            #karin added code to truncate ri1, etc, off colonies UIDs so that they match the interaction matrix
    if cp[length-3:-1] == 'ri':
        cp = cp[:-3]
        length = 0
        
    if cp in imatrix:
        if hu in hu_multiple:
            hu_orig = hu
            scores = []
            for variant in hu_multiple[hu_orig]['variants']:
                hu = variant.split('_')
                hu = (''.join(hu[3] + hu[-1])).lower()
                if hu in imatrix[cp]:
                    if detailed_status:
                        print "processing " + cp + " in determineInteraction"
                    scores.append(imatrix[cp][hu])
            if 'HIGH' in scores:
                return 'HIGH'
            elif 'MODERATE' in scores:
                return 'MODERATE'
        else:
            hu = hu.split('_')[3]
            if hu in imatrix[cp]:
                if detailed_status:
                    print "processing " + cp + " in determineInteraction"
                return imatrix[cp][hu]          
    #else:
        #print cp + " not in imatrix"
        # I don't want to make this an error since it is possible that a cp has no interactions
        # Therefore it's very important that names match between files and the imatrix
    return None

## identifyInteractions ##
#
# Finds HUs and CPs in the same MPA and checks if they have an interaction
# Builds a dict with interaction scaling factors for each CP in each MPA
#

def identifyInteractions(hu_in_mpas, cp_in_mpas, imatrix):
    cp_in_mpa_i = {}

    for mpa in hu_in_mpas:
        if mpa not in cp_in_mpas:
            continue

        if mpa not in cp_in_mpa_i:
            cp_in_mpa_i[mpa] = {}

        for ecosection in cp_in_mpas[mpa]:
            # we only need to know the interaction once. The ecosection doesn't matter here.
            for cp in cp_in_mpas[mpa][ecosection]:
                if cp not in cp_in_mpa_i[mpa]:
                        cp_in_mpa_i[mpa][cp] = {'interactions': [],
                                                'eff_score': None}
                else:
                    continue  # we only need the cp once per mpa, so if we have already encountered it then skip
                        
                for hu in hu_in_mpas[mpa]:
                    interaction = determineInteraction(imatrix, cp, hu)
                    if interaction is not None:
                        cp_in_mpa_i[mpa][cp]['interactions'].append(interaction)

    # add in any mpas that were not in hu_in_mpas but were in cp_in_mpas. These need to be carried forward, even if they do not have any interactions.
    for mpa in cp_in_mpas:
        if mpa not in cp_in_mpa_i:
            cp_in_mpa_i[mpa] = {}
            for ecosection in cp_in_mpas[mpa]:
                for cp in cp_in_mpas[mpa][ecosection]:
                    if cp not in cp_in_mpa_i[mpa]:
                        cp_in_mpa_i[mpa][cp] = {'interactions': [],
                                                'eff_score': None}
                else:
                    continue  # we only need the cp once per mpa, so if we have already encountered it then skip


    return cp_in_mpa_i

## prepareOutputTable1 ##
#
# Pulls all the data together to make a dict that can be used to create the final output
# Calculates effectiveness scores and uses it to scale the CPs
#
# The returned dict is structured like so: o_table_1[ MPA NAME ][ REGION ][ CP ]
# this points to a dict that has the following keys:
#
# 'mpa_area'     -> area of the enclosing mpa
# 'og_area'      -> total area of the original CP
# 'uscaled_area' -> area of the clipped and adjusted CP in the MPA
# 'scaled_area'  -> uscaled_area but multiplied by the CPs effectiveness score
# 'pct_of_mpa'   -> scaled_area / mpa_area
# 'pct_of_og'    -> scaled_area / og_area


def prepareOutputTable1(cp_in_mpa_i, cp_in_mpas):
    o_table_1 = {}

    ecosections = ["Johnstone Strait", "Continental Slope", "Dixon Entrance", "Strait of Georgia", "Juan de Fuca Strait", "Queen Charlotte Strait", "North Coast Fjords", "Hecate Strait", "Queen Charlotte Sound", "Vancouver Island Shelf", "Transitional Pacific", "Subarctic Pacific"]
    
    for mpa in cp_in_mpa_i:
        if mpa not in o_table_1:
            o_table_1[mpa] = {}

        for cp in cp_in_mpa_i[mpa]:
            for ecosection in ecosections:
                if ecosection in cp_in_mpas[mpa] and cp in cp_in_mpas[mpa][ecosection]:                    
                    if ecosection not in o_table_1[mpa]:
                        o_table_1[mpa][ecosection] = {}

                    # Pre populate dict with already known values
                    if cp not in o_table_1[mpa][ecosection]:
                        o_table_1[mpa][ecosection][cp] = {'mpa_area': cp_in_mpas[mpa][ecosection][cp]['mpa_area'],
                                                      'og_area': cp_in_mpas[mpa][ecosection][cp]['orig_area'],
                                                      'unscaled_area': cp_in_mpas[mpa][ecosection][cp]['clip_area'],
                                                      'scaled_area': None,
                                                      'pct_of_mpa': None,
                                                      'pct_of_og': None,
                                                      'pct_of_og_unscaled': None,
                                                      'subregion': cp_in_mpas[mpa][ecosection][cp]['subregion']}

                    # Calculate effectiveness
                    num_high, num_mod, num_low = countInteractions(cp_in_mpa_i[mpa][cp]['interactions'])

                    cp_in_mpa_i[mpa][cp]['eff_score'] = calcEffectivenessScore(num_high,
                                                                               num_mod,
                                                                               num_low)

                    # Rescale areas and calculate new percentages
                    eff_score = cp_in_mpa_i[mpa][cp]['eff_score']
                    unscaled_area = o_table_1[mpa][ecosection][cp]['unscaled_area']
                    o_table_1[mpa][ecosection][cp]['scaled_area'] = eff_score * unscaled_area

                    scaled_area = o_table_1[mpa][ecosection][cp]['scaled_area']
                    mpa_area = o_table_1[mpa][ecosection][cp]['mpa_area']
                    og_area = o_table_1[mpa][ecosection][cp]['og_area']

                    o_table_1[mpa][ecosection][cp]['pct_of_mpa'] = scaled_area / mpa_area
                    o_table_1[mpa][ecosection][cp]['pct_of_og'] = scaled_area / og_area
                    o_table_1[mpa][ecosection][cp]['pct_of_og_unscaled'] = unscaled_area / og_area
                
    return o_table_1

## writeOutputTable1 ##
#
# Writes the final output to disk
#

def writeOutputTable1(otable, opath, mpa_dict):

    with open(opath, 'wb') as f:
        w = csv.writer(f)

        w.writerow(['UID', 'name', 'subregion', 'ecosection', 'CP', 'proportion_scaled', 'proportion_unscaled', 'value_scaled','value_unscaled', 'total_value'])

        for mpa in otable:
            name = mpa_dict[mpa]['name']
            for ecosection in otable[mpa]:
                for cp in otable[mpa][ecosection]:
                    pct_of_og = otable[mpa][ecosection][cp]['pct_of_og']
                    subregion = otable[mpa][ecosection][cp]['subregion']
                    unscaled_area = otable[mpa][ecosection][cp]['unscaled_area']
                    scaled_area = otable[mpa][ecosection][cp]['scaled_area']
                    total_area = otable[mpa][ecosection][cp]['og_area']
                    pct_of_og_unscaled = otable[mpa][ecosection][cp]['pct_of_og_unscaled']
                    w.writerow([mpa.encode('utf8'), name.encode('utf8'), subregion, ecosection, cp, pct_of_og, pct_of_og_unscaled, scaled_area, unscaled_area, total_area])

def createOutputTable2(o_table_1, cp_area_overlap_dict):
    table2 = {}

    fields = ['original', 'protected', 'pct']

    for mpa in o_table_1:
        for ecosection in o_table_1[mpa]:
            for cp in o_table_1[mpa][ecosection]:
                cp_data = o_table_1[mpa][ecosection][cp]

                if cp not in table2:
                    table2[cp] = {}

                if ecosection not in table2[cp]:
                    table2[cp][ecosection] = {}

                subregion = cp_data['subregion']
                # I had to wrap any mention of subregion in this function in the below IF statement
                # The issue: a few rare cases where the mpa is assigned a subregion because some of it overlaps that subregion,
                #  but some of it is also overlaps a blank area or a different subregion, AND the cp touches the mpa, BUT the cp
                # does not touch the subregion. Therefore there ends up being a discrepancy where the subregion is not in the
                # cp_area_overlap_dict of that cp, but it is in the o_table_1 (since this one is based on the mpa
                # I'm thinking the best way to handle this: if the cp is not in the subregion, then it is not in the subregion.
                #  So perhaps it is just an if statement before to just not write it for that one.
                # However, this now results in any MPA that has a subregion of NONE not being carried forward from this point forward.
                # This does not matter for now, since I don't write None out in table2. If for some reason I need to know that
                # areas that don't fall in a subregion then I will need to consider this.
                # Also, I considered the scenario where an mpa overlaps two subregions but it is only assigned one,
                # and the cp overlaps the part of the unlisted subregion. This should be fine because as of right now this
                #  only affects the middle piece of the hecate strait sponge reef mpa, which we can just split out if needed.
                # This if statement can help identify those areas:
                #if subregion not in cp_area_overlap_dict[cp] and subregion is not None:
                #    print cp
                #    print cp_data
                #    print mpa
                #    print cp_area_overlap_dict[cp]

                if subregion in cp_area_overlap_dict[cp]:
                    if subregion not in table2[cp]:
                        table2[cp][subregion] = {}

                for field in fields:
                    if field not in table2[cp][ecosection]:
                        table2[cp][ecosection][field] = 0.0
                
                if subregion in cp_area_overlap_dict[cp]:
                    for field in fields:
                        if field not in table2[cp][subregion]:
                            table2[cp][subregion][field] = 0.0
                
                # get total area in ecosection
                orig_area_eco = cp_area_overlap_dict[cp][ecosection]['Area']

                if subregion in cp_area_overlap_dict[cp]:
                    if subregion is not None:
                        orig_area_sub = cp_area_overlap_dict[cp][subregion]['Area']
                    else:
                        # I think the script will never get here now that it is wrapped in the if statement above this one.
                        #  I will leave this in though in case I revert back.
                        orig_area_sub = 1 # this shouldn't matter since we won't write the pct of subregion-None out to table 2 anyways

                # Sum up protected area from all MPAs for CP
                table2[cp][ecosection]['original'] = orig_area_eco
                table2[cp][ecosection]['protected'] = table2[cp][ecosection]['protected'] \
                                                 + cp_data['scaled_area']

                if subregion in cp_area_overlap_dict[cp]:
                    table2[cp][subregion]['original'] = orig_area_sub
                    table2[cp][subregion]['protected'] = table2[cp][subregion]['protected'] \
                                                     + cp_data['scaled_area']
    # Calculate percentages
    for cp in table2:
        for ecosub in table2[cp]:
            if table2[cp][ecosub]['original'] != 0:
                table2[cp][ecosub]['pct'] = table2[cp][ecosub]['protected']/table2[cp][ecosub]['original']

    return table2


def writeOutputTable2(o_table_2, ofile):

    with open(ofile, 'wb') as f:
        w = csv.writer(f)

        # Write header
        w.writerow(['cp', 'ecosection_subregion', 'proportion', 'original_area', 'protected_area'])

        for cp in o_table_2:
            for eco_sub in o_table_2[cp]:
                if eco_sub is not None:
                    pct_of_ecosub = o_table_2[cp][eco_sub]['pct']
                    orig_area = o_table_2[cp][eco_sub]['original']
                    protected = o_table_2[cp][eco_sub]['protected']
                    w.writerow([cp, eco_sub, pct_of_ecosub, orig_area, protected])



def writeOutputTable3(percent_overlap, output3_path):
    cols = ['mpa','type','cp_hu','percent_area_overlap']

    with open(output3_path, 'wb') as f:
        w = csv.writer(f)

        # Write header
        w.writerow(cols)

        for mpa in percent_overlap:
            for layer_type in percent_overlap[mpa]:
                for cphu in percent_overlap[mpa][layer_type]:
                        pct_o = percent_overlap[mpa][layer_type][cphu]['pct_overlap_cphu_mpa']
                        w.writerow([mpa.encode('utf8'), layer_type, cphu, pct_o])

def joinUIDtoTable1(output1_path, ecoUIDs_path, output1join_path):
    
    a = pd.read_csv(output1_path)
    b = pd.read_csv(ecoUIDs_path)

    joined = a.merge(b, left_on = 'CP', right_on = 'Desktop_UID', how = 'left')

    joined.to_csv(output1join_path, index = False)

# output table 4: list of cp and hus that interact    
def createOutputTable4(hu_in_mpas, cp_in_mpas, imatrix, output4_path):
    cphu_int = {}

    for mpa in hu_in_mpas:
        if mpa not in cp_in_mpas:
            continue

        if mpa not in cphu_int:
            cphu_int[mpa] = {}

        for ecosection in cp_in_mpas[mpa]:
            # we only need to know the interaction once. The ecosection doesn't matter here.
            for cp in cp_in_mpas[mpa][ecosection]:

                for hu in hu_in_mpas[mpa]:
                    interaction = determineInteraction(imatrix, cp, hu)
                    if interaction is not None:

                        if cp not in cphu_int[mpa]:
                            cphu_int[mpa][cp] = {}
                        cphu_int[mpa][cp][hu] = interaction
                                          
    # write to output csv
    cols = ['mpa','cp','hu','score']

    with open(output4_path, 'wb') as f:
        w = csv.writer(f)

        # Write header
        w.writerow(cols)

        for mpa in cphu_int:
            for cp in cphu_int[mpa]:
                for hu in cphu_int[mpa][cp]:
                        score = cphu_int[mpa][cp][hu]
                        w.writerow([mpa.encode('utf8'), cp, hu, score])


  ##               ##
###  Program start  ###
  ##               ##
        
# Generate unique name for temp gdb and make it the workspace
i = 0;
while arcpy.Exists(os.path.join(working_gdb_folder, 'temp{0}.gdb'.format(str(i)))):
    i = i + 1
    
working_gdb = os.path.join(working_gdb_folder, 'temp{0}.gdb'.format(str(i)))

arcpy.CreateFileGDB_management(os.path.dirname(working_gdb), os.path.basename(working_gdb))
arcpy.env.workspace = working_gdb

#####
### Load Ecosection layer into workspace
#####

if print_status:
    print "Preparing Ecosections"

new_bc_area_field = 'etp_bc_area'
new_bc_total_area_field = 'etp_bc_total_area'
new_scaling_field = 'etp_scaling'

layer_list = arcpy.mapping.ListLayers(arcpy.mapping.MapDocument(source_mxd))
for lyr in layer_list:
    if lyr.isFeatureLayer and (lyr.datasetName == 'eco_coarse_ecosections_polygons_d'):
        ecosections = lyr
ecosections_layer = loadLayer(source_mxd, ecosections.name, sr_code,
                              new_bc_area_field, new_bc_total_area_field,
                              None, scaling_attribute, new_scaling_field,
                              None, "value", "area")

#####
### Load subregional layer into workspace
### This is the one layer that has all the subregions in it.
### It is used to determine which subregion each MPA is in
#####

layer_list = arcpy.mapping.ListLayers(arcpy.mapping.MapDocument(source_mxd))
for lyr in layer_list:
   if lyr.isFeatureLayer and (lyr.datasetName.startswith('rgn_subregions')):
       subregions_ALL = loadRegionLayer(source_mxd, lyr.name,
                                                sr_code, new_bc_area_field,
                                                new_bc_total_area_field)


#####
### Load MPA layers into workspace, create consistent name attribute, and merge together
#####

if print_status:
    print "Preparing MPAs"

mpa_area_attribute = 'etp_mpa_area_TOTAL'
merged_name_field = 'NAME_UID'   # make sure this is unique and does not exist in any of the input mpa datasets. If it is not unique it screws up the field mapping. It may not throw an error and is hard to detect.
final_mpa_fc_name = 'mpas'
mpa_subregion_field = 'subregion_mpa'
mpa_area_attribute_section = 'etp_mpa_area_SECTION'
mpa_marine_area = 'marine_m2' # this is now required: we have combined terrestrial and marine portions of a protected area,
#  but we are only concerned with the calculation of the marine area

# create the mpa dictionary lookup
mpa_dict = createMPAdict(source_mxd, merged_name_field)

final_mpa_fc_name = prepareMPAs(source_mxd, sr_code, mpa_area_attribute, mpa_area_attribute_section,
                                final_mpa_fc_name, merged_name_field, mpa_name_fields, mpa_subregion_field, subregions_ALL, ecosections_layer, mpa_marine_area)

#####
### Load subregional layers into workspace
#####

layer_list = arcpy.mapping.ListLayers(arcpy.mapping.MapDocument(source_mxd))
layer_list = [lyr for lyr in layer_list if lyr.isFeatureLayer \
              and (lyr.datasetName.startswith('rgn_subregion_'))]

rlayers = {}

for layer in layer_list:
    subregion = layer.datasetName.split('_')[3]
    rlayers[subregion] = loadRegionLayer(source_mxd, layer.name,
                                         sr_code, new_bc_area_field,
                                         new_bc_total_area_field)


#####
### Load HU/CP layers into workspace, calculate areas etc
#####

# Load CP area overlap dictionary
cp_area_overlap_dict = {}
if cpOverlap_newDict is False:
    cp_area_overlap_dict = buildOverlapDict(cpOverlap_DictPath, cp_area_overlap_dict)

# Load attribute scaling file if necessary
scaling_dict = None
if scaling_attribute_file is not None:
    scaling_dict = buildScalingDict(scaling_attribute_file)

threshold_dict = None
if layer_presence_threshold_file is not None:
    threshold_dict = buildThresholdDict(layer_presence_threshold_file)

inclusion_matrix = readMPAInclusionMatrix(inclusion_matrix_path)

# Generate layer list based on dataset names (Karin adds that layer list includes only the HU and ECO layers
layer_list = arcpy.mapping.ListLayers(arcpy.mapping.MapDocument(source_mxd))
layer_list = [lyr for lyr in layer_list if lyr.isFeatureLayer \
              and (lyr.datasetName.startswith('eco_') \
                   or lyr.datasetName.startswith('hu_'))]

arcpy.env.overwriteOutput = True

clipped_adjusted_area = 'etp_ac_area_adj'
pct_of_total_field = 'pct_of_total'
pct_of_mpa_field = 'pct_of_mpa'
# JC fields added; note that a density feature needs to have a field named "value"
clipped_adj_area_mpaTotal = 'etp_ac_area_adj_mpaTotal'
pct_of_mpa_field_Total = 'pct_of_mpa_Total'
density_field = "value"

hu_in_mpas,cp_in_mpas = {}, {}
percent_overlap = {}
for lyr in layer_list:
    if print_status:
        print "Processing " + lyr.name + " for presence in MPAs"

    # Load layer into memory, reprojecting to albers, calculate areas etc
    # If a layer is complex and causes processing to fail populate
    # complexFeatureClasses above and hopefully that fixes it
    is_complex = False
    for fc in complexFeatureClasses:
        if lyr.datasetName.startswith(fc):
            is_complex = True
            break 

    # determine if layer values are based on area or density/diversity; value_type now holds info on what attribute to look for
    value_type = 'area'
    for field in arcpy.ListFields(lyr.dataSource):
        if field.name == density_field:
            value_type = 'density'
            break

    working_layer = loadLayer(source_mxd, lyr.name, sr_code,
                              new_bc_area_field, new_bc_total_area_field,
                              scaling_dict, scaling_attribute, new_scaling_field,
                              is_complex, density_field, value_type)

    layer_type = 'cp' if working_layer.startswith('eco_') else 'hu'

    # Set to default hu/cp presence threshold and overwrite with value in threshold_dict
    # if possible
    threshold = hu_presence_threshold if layer_type == 'hu' else cp_presence_threshold
    if threshold_dict is not None and working_layer in threshold_dict:
        threshold = threshold_dict[working_layer]

    # Check last element in layer dataset name to see if has been subregionally clipped
    # and get that subregion layer - Karin asks WHY???
    subregion = working_layer.split('_')[-1]
    rlayer = rlayers[subregion] if subregion in rlayers else None
    subregion = 'region' if rlayer is None else subregion
    

    # find area OR VALUE of cp in each ecosection and subregion
    if layer_type == 'cp' and working_layer not in cp_area_overlap_dict:
        cp_area_overlap_dict = calcCPlyrOverlap(cp_area_overlap_dict, working_layer,
                         ecosections_layer, subregions_ALL, density_field, value_type)
    
    # Determine if in which MPAs and calculate statistics
    mpa_presence, sliver_freq = calculate_presence(working_layer, final_mpa_fc_name, clipped_adjusted_area,
                                      pct_of_total_field, pct_of_mpa_field, merged_name_field,
                                      new_scaling_field, threshold, subregion, inclusion_matrix,
                                      mpa_subregion_field, mpa_area_attribute_section, clipped_adj_area_mpaTotal,
                                      pct_of_mpa_field_Total, density_field, value_type)


    # If subregion fc split off that subregion tag on the fc name
    if subregion != 'region':
        working_layer = '_'.join(working_layer.split('_')[:-1])

    if layer_type == 'hu':
        for mpa in mpa_presence:            
            if mpa not in hu_in_mpas:
                hu_in_mpas[mpa] = {}
            hu_in_mpas[mpa][working_layer] = mpa_presence[mpa]
            # all we need to know is if an hu occurs in an mpa. We don't care about its area measurements at this point, so I can just keep this as is.
    else:        
        for mpa in mpa_presence:
            if mpa not in cp_in_mpas:
                cp_in_mpas[mpa] = {}
            for ecosection in mpa_presence[mpa]:
                if ecosection not in cp_in_mpas[mpa]:
                    cp_in_mpas[mpa][ecosection] = {}
                cp_in_mpas[mpa][ecosection][working_layer] = mpa_presence[mpa][ecosection]

    # Populate percent_overlap dictionary
    for mpa in sliver_freq:
        if mpa not in percent_overlap:
            percent_overlap[mpa] = {}
        if layer_type not in percent_overlap[mpa]:
            percent_overlap[mpa][layer_type] = {}
        percent_overlap[mpa][layer_type][working_layer] = sliver_freq[mpa]


# cpoverlap dict to csv to be used in future sessions
cols = ['cp','section_region','area_overlap']
with open(cpOverlap_DictPath, 'wb') as f:
    w = csv.writer(f)
    w.writerow(cols)
    for cp in cp_area_overlap_dict:
        for sec_reg in cp_area_overlap_dict[cp]:
                area_o = cp_area_overlap_dict[cp][sec_reg]['Area']
                w.writerow([cp, sec_reg, area_o])


# Tack in dummy data for HU that should be in each MPA according to the inclusion matrix
# but didn't have spatial data that sufficiently intersected
if not override_y:
    for mpa in inclusion_matrix:
        for hu in inclusion_matrix[mpa]:
            if inclusion_matrix[mpa][hu] in ['O', 'C']:
                if mpa not in hu_in_mpas:
                    hu_in_mpas[mpa] = {}

                for fc in hu_multiple:
                    if hu in hu_multiple[fc]['variants']:
                        hu = fc
                        break

                if hu not in hu_in_mpas[mpa]:
                    hu_in_mpas[mpa][hu] = {'ecosect_placeholder': # I dont think these values or the ecosection name matters here since all we need to know is if an hu occurs in an mpa.
                                           {'clip_area': 1,
                                           'orig_area': 1,
                                           'mpa_area': 1,
                                           'region_area': 1,
                                           'pct_in_mpa': 1,
                                           'pct_of_region': 1,
                                           'pct_of_total': 1}}
# Clean up
if cleanUpTempData:
    arcpy.Delete_management(working_gdb)
    
#####
### Find HU-CP interactions within MPAs and write output table 1
#####

imatrix = loadInteractionsMatrix(imatrix_path)

cp_in_mpa_i = identifyInteractions(hu_in_mpas, cp_in_mpas, imatrix)

o_table_1 = prepareOutputTable1(cp_in_mpa_i, cp_in_mpas)

writeOutputTable1(o_table_1, output1_path, mpa_dict)

#####
### Create and write output table 2
#####

o_table_2 = createOutputTable2(o_table_1, cp_area_overlap_dict)

writeOutputTable2(o_table_2, output2_path)

#####
### Write percent overlap (sliver) table
#####

writeOutputTable3(percent_overlap, output3_path)

#####
### Join eco UID table to table1
#####

joinUIDtoTable1(output1_path, ecoUIDs_path, output1join_path)

#####
### Create and write output table 4 (cp-hu list)
#####

createOutputTable4(hu_in_mpas, cp_in_mpas, imatrix, output4_path)