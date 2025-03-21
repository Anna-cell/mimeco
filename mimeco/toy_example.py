
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 14:38:53 2025

@author: e158401a
"""

import os
#os.chdir("/home/e158401a/Documents/mimeco")
import cobra
from mimeco import analysis
import pandas as pd
from importlib.resources import files


#test bacteria - bacteria interaction inference

#BiGGs namespace

Western_diet = pd.read_csv(files("mimeco.resources").joinpath("Western_diet_BiGG.csv"), index_col = 0)

model1 = cobra.io.read_sbml_model("/home/e158401a/Documents/models/embl_gems/models/b/bifidobacterium/Bifidobacterium_adolescentis_ATCC_15703.xml")
model2 = cobra.io.read_sbml_model("/home/e158401a/Documents/models/embl_gems/models/l/lactobacillus/Lactobacillus_rhamnosus_GG_GG_ATCC_53103.xml")
model1.solver = "cplex"
model2.solver = "cplex"

#score and type
import time
start_time = time.time()
int_score, int_type = analysis.interaction_score_and_type(model1, model2, Western_diet, 
                                                        undescribed_metabolites_constraint="as_is", plot=True)
print("--- %s seconds ---" % (time.time() - start_time))

int_score, int_type = analysis.interaction_score_and_type(model1, model2, Western_diet, 
                                                        undescribed_metabolites_constraint="partially_constrained", plot=True)

# crossfed metabolites
model1_biomass_id = "Growth"
model2_biomass_id = "Growth"

potential_exchange, data = analysis.crossfed_metabolites_plotdata(model1 = model1, model2 = model2, medium = Western_diet, undescribed_metabolites_constraint = "partially_constrained",
                               solver = "cplex", model1_biomass_id = model1_biomass_id, model2_biomass_id = model2_biomass_id, plot = True, retrieve_data = "selection")


potential_exchange2, data2 = analysis.crossfed_metabolites_plotdata(model1 = model2, model2 = model1, medium = Western_diet, undescribed_metabolites_constraint = "partially_constrained",
                               solver = "cplex", model1_biomass_id = model2_biomass_id, model2_biomass_id = model1_biomass_id, plot=True, retrieve_data = "all")
#with enterocyte
start_time = time.time()
int_score, int_type = analysis.enterocyte_interaction_score_and_type(model1, Western_diet, undescribed_metabolites_constraint="as_is", 
                                                                   namespace="BIGG", plot=True)
print("--- %s seconds ---" % (time.time() - start_time))

potential_crossfeeding = analysis.enterocyte_crossfed_metabolites(model = model1, medium = Western_diet, undescribed_metabolites_constraint = "as_is", 
                                                                solver = "cplex", model_biomass_id = model1_biomass_id, namespace = "BIGG", 
                                                                plot = True, sample_size = 1000)

"""
# AGORA namespace

bact_agora2 = cobra.io.load_matlab_model("/home/e158401a/Documents/CH2_FB/AGORA2_bacteria/Lactobacillus_rhamnosus_GG_ATCC_53103.mat")
bact_agora2.solver = "cplex"

#host = cobra.io.read_sbml_model("/home/e158401a/Documents/mimeco/mimeco/resources/enterocyte_VMH_v3.xml")
#host2 = cobra.io.read_sbml_model("/home/e158401a/Documents/CH2_FB/enterocyte_VMH.xml")

bact_agora2_biomass_id = "Growth" 
potential_crossfeeding = mimeco.enterocyte_crossfed_metabolites(model = bact_agora2, medium = Western_diet, undescribed_metabolites_constraint = "as_is", 
                                                                solver = "cplex", model_biomass_id = bact_agora2, namespace = "AGORA", 
                                                                plot = True, sample_size = 1000)
"""