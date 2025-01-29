
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 14:38:53 2025

@author: e158401a
"""

import cobra
import pickle
import mimeco

resources_path = "/home/e158401a/Documents/pareto_analysis/resources/"
with open(resources_path+"mediums_lumen_v2.pickle", "rb") as fp:   # Unpickling
    mediums_lumen = pickle.load(fp)

Western_diet = mediums_lumen["Western_diet"]
model1 = cobra.io.read_sbml_model("/home/e158401a/Documents/models/embl_gems/models/l/lactobacillus/Lactobacillus_plantarum_WCFS1.xml")
model2 = cobra.io.read_sbml_model("/home/e158401a/Documents/models/embl_gems/models/a/akkermansia/Akkermansia_muciniphila_ATCC_BAA_835.xml")
model1.solver = "cplex"
model2.solver = "cplex"


model1_biomass_id = "Growth"
model2_biomass_id = "Growth"

int_score, int_type = mimeco.interaction_score_and_type(model1, model2, Western_diet, 
                                                        undescribed_metabolites_constraint="partially_constrained")

potential_exchange = mimeco.exchanged_metabolites(model1 = model1, model2 = model2, medium = Western_diet, undescribed_metabolites_constraint = "partially_constrained",
                               solver = "cplex", model1_biomass_id = model1_biomass_id, model2_biomass_id = model2_biomass_id)

potential_exchange2to1 = mimeco.exchanged_metabolites(model1 = model2, model2 = model1, medium = Western_diet, undescribed_metabolites_constraint = "partially_constrained",
                               solver = "cplex", model1_biomass_id = model1_biomass_id, model2_biomass_id = model2_biomass_id)


with open("/home/e158401a/Documents/MIMECO/tests/resources/Lactobacillus_plantarum_WCFS1_Akkermansia_muciniphila_ATCC_BAA_835_WD_PC_potential_exchange.pickle", "wb") as fp:
    pickle.dump(potential_exchange, fp)
    
Western_diet.to_csv("resources/Western_diet.csv", index = True)

import pandas as pd

WD = pd.read_csv("resources/Western_diet.csv", index_col = 0)
