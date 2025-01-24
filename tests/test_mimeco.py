#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:32:14 2025

@author: e158401a
"""
import pickle
import cobra
import utils
import mimeco



model1 = cobra.io.read_sbml_model("/home/e158401a/Documents/MIMECO/tests/resources/Lactobacillus_plantarum_WCFS1.xml")
model2 = cobra.io.read_sbml_model("/home/e158401a/Documents/MIMECO/tests/resources/Akkermansia_muciniphila_ATCC_BAA_835.xml")
model1.solver = "cplex"
model2.solver = "cplex"

with open("/home/e158401a/Documents/MIMECO/tests/resources/Western_diet.pickle", "rb") as fp:   # Unpickling
    Western_diet = pickle.load(fp)

with open("/home/e158401a/Documents/MIMECO/tests/resources/Lactobacillus_plantarum_WCFS1_Akkermansia_muciniphila_ATCC_BAA_835_WD_PC_interaction_type.pickle", "rb") as fp:   # Unpickling
    interaction_type_LPAM = pickle.load(fp)

with open("/home/e158401a/Documents/MIMECO/tests/resources/Lactobacillus_plantarum_WCFS1_Akkermansia_muciniphila_ATCC_BAA_835_WD_PC_interaction_score.pickle", "rb") as fp:   # Unpickling
    interaction_score_LPAM = pickle.load(fp)

with open("/home/e158401a/Documents/MIMECO/tests/resources/Lactobacillus_plantarum_WCFS1_Akkermansia_muciniphila_ATCC_BAA_835_WD_PC_medium_dict.pickle", "rb") as fp:   # Unpickling
    medium_dict_LPAM = pickle.load(fp)

with open("/home/e158401a/Documents/MIMECO/tests/resources/Lactobacillus_plantarum_WCFS1_Akkermansia_muciniphila_ATCC_BAA_835_WD_PC_potential_exchange.pickle", "rb") as fp:   # Unpickling
    exchange_potential_LPAM = pickle.load(fp)

model1.solver = "cplex"
model2.solver = "cplex"

model1_biomass_id = "Growth"
model2_biomass_id = "Growth"

def decomposed_interaction_score_and_type(model1, model2, medium, undescribed_metabolites_constraint):
    metabolic_dict = utils.create_ecosystem_metabolic_dict(model1, model2)
    with model1:
        model1, constrained_medium_dict1 = utils.restrain_medium(model1, medium, undescribed_metabolites_constraint)
        solo_growth_model1 = model1.optimize().objective_value
    with model2:
        model2, constrained_medium_dict2 = utils.restrain_medium(model2, medium, undescribed_metabolites_constraint)
        solo_growth_model2 = model2.optimize().objective_value
    print(solo_growth_model1, solo_growth_model2)
    medium_dict = {**constrained_medium_dict1, **constrained_medium_dict2}
    model1 = utils.unrestrain_medium(model1)
    model2 = utils.unrestrain_medium(model2)
    sol_mofba = utils.mo_fba(model1, model2, metabolic_dict, medium_dict)
    xy, maxi_model1, maxi_model2 = utils.pareto_parsing(sol_mofba, solo_growth_model1, solo_growth_model2)
    interaction_score = utils.infer_interaction_score(xy)
    interaction_type = utils.infer_interaction_type(interaction_score, maxi_model1, maxi_model2,solo_growth_model1, solo_growth_model2)
    test_results = {"constrained_medium_dict":medium_dict, "interaction_score":interaction_score, "interaction_type":interaction_type}
    return test_results

def test_constrained_medium_dict():
    """Test the content of constrained_medium_dict created by restrain_medium"""
    expected = medium_dict_LPAM
    actual = decomposed_interaction_score_and_type(model1, model2, Western_diet, undescribed_metabolites_constraint="partially_constrained")["constrained_medium_dict"]
    assert actual == expected, "medium constraint not built correctly"

def test_interaction_score():
    expected = interaction_score_LPAM
    actual = decomposed_interaction_score_and_type(model1, model2, Western_diet, undescribed_metabolites_constraint="partially_constrained")["interaction_score"]
    assert actual == expected, "interaction_score incorrect"

def test_interaction_type():
    expected = interaction_type_LPAM
    actual = decomposed_interaction_score_and_type(model1, model2, Western_diet, undescribed_metabolites_constraint="partially_constrained")["interaction_type"]
    assert actual == expected, "interaction_type incorrect"

def test_potential_exchange():
    expected = exchange_potential_LPAM
    actual = mimeco.exchanged_metabolites(model1 = model1, model2 = model2, medium = Western_diet, undescribed_metabolites_constraint = "partially_constrained",
                                   solver = "cplex", model1_biomass_id = model1_biomass_id, model2_biomass_id = model2_biomass_id)
    assert actual == expected, "potential_exchange incorrect"