#!/usr/bin/env python

"""MiMB_MetabFBA_Implementation.py - Example code for Methods in Molecular Biology
GitHub Repo: https://github.com/michaelwitting/MiMB_CompMetabo

Requires (use pip to install):

cobra
ecsher
"""

__author__ = "Jake Hattwell, Janna Hastings"
__copyright__ = "None"
__credits__ = ["Jake Hattwell, Janna Hastings"]
__license__ = "CCO"
__version__ = "1"
__maintainer__ = "Jake Hattwell"
__email__ = "j.hattwell@uq.edu.au"
__status__ = "Live"


import csv 
import cobra
import cobra.flux_analysis
from cobra.core import Reaction

mapping_file = "PATH_TO_MAPPING_TABLE"
metabolite_file = "PATH_TO_DIFFERENCES_TABLE"
model_file = "PATH_TO_MODEL_FILE"

with open(mapping_file,"r") as file:
    mapping_data_input = list(csv.reader(file,delimiter="\t"))
    mapping_data = {entry[0]:entry[1] for entry in mapping_data_input}


with open(metabolite_file,'r') as file:
    metabolite_data = list(csv.reader(file,delimiter="\t"))
    metabolites = [row[0] for row in metabolite_data][1:]

###############################################################################

def addMetabolomicsReaction(model, metabolites, react_name, coefficient_str=1):
   print(react_name,":",metabolites)
   # Build a reaction
   coefficient_str=str(coefficient_str)
   reaction_str = coefficient_str + (" + "+coefficient_str).join(metabolites) + " <==> " + react_name+"_c"
   reaction = Reaction(react_name)
   reaction.name = react_name
   reaction.subsystem = 'Metabolomics integration' 
   # Add the reaction to the model
   model.add_reactions([reaction])
   reaction.build_reaction_from_string(reaction_str)
   #print(reaction,":::::",reaction_str)
   return reaction

###############################################################################

def integrate_metabolomics(model,group):
    ## generate a list of the UP, DOWN and NODIFF assignments to metabolites for the column "group"
    metabo_diff_values = [row[metabolite_data[0].index(group)] for row in metabolite_data][1:]

    # filter the increased and decreased metabolites and translate the identifiers to model identifiers using the mapping data
    up_mets_codes = [metabolites[i] for i in range(len(metabolites)) if metabo_diff_values[i] == "UP"]
    up_mets_codes = [mapping_data[i] for i in up_mets_codes if i in mapping_data and mapping_data[i] != ""]
    down_mets_codes = [metabolites[i] for i in range(len(metabolites)) if metabo_diff_values[i] == "DOWN"]
    down_mets_codes = [mapping_data[i] for i in down_mets_codes if i in mapping_data and mapping_data[i] != ""]
    # arbitrary threshold
    threshold = 2.5

    # If there are changed metabolites, add a sink reaction for them.
    if (len(up_mets_codes)>0):
        up_mets_react = addMetabolomicsReaction(model,up_mets_codes, "up_metabolites")
        model.metabolites.up_metabolites_c.compartment = 'c'
        model.add_boundary(model.metabolites.up_metabolites_c, type='sink',lb=0,ub=threshold)
    if (len(down_mets_codes)>0):
        down_mets_react = addMetabolomicsReaction(model,down_mets_codes, "down_metabolites")
        model.metabolites.down_metabolites_c.compartment = 'c'
        model.add_boundary(model.metabolites.down_metabolites_c, type='sink',lb=-1*threshold,ub=0)

    # Add the demand reactions to the model objective
    if (len(up_mets_codes)>0 and len(down_mets_codes)>0):
        model.objective = model.reactions.BIO0100.flux_expression + model.reactions.up_metabolites.flux_expression - model.reactions.down_metabolites.flux_expression
    elif (len(up_mets_codes)>0):
        model.objective = model.reactions.BIO0100.flux_expression + model.reactions.up_metabolites.flux_expression
    elif (len(down_mets_codes)>0):
        model.objective = model.reactions.BIO0100.flux_expression - model.reactions.down_metabolites.flux_expression
    return model

###############################################################################

model = cobra.io.read_sbml_model(model_file)
model = integrate_metabolomics(model,"INSERT_COLUMN_NAME")
solution = cobra.flux_analysis.pfba(model)
print(solution)
