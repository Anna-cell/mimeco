#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 09:46:31 2025

@author: e158401a
"""

import utils
import pickle
import warnings

def interaction_score_and_type(model1, model2, medium, undescribed_metabolites_constraint):
    """
    A function that, given 2 models in the same namespace and a defined medium, analyses the interaction between the two models,
    infering qualitative (interaction_type) and quantitative (interaction_score) information on their metabolic interaction. 
    
    Parameters
    ----------
    model1 : cobra.Model
    model2 : cobra.Model
    medium : pandas series
        Index : metabolites names
        values  : Availability of corresponding metabolite in the medium as a positive flux value. 
    undescribed_metabolites_constraint : string ("blocked" or "partially_constrained"). 
        How strictly constrained are the medium metabolites for which the flux is not described in the medium dataframe.
        "blocked" : They are not available in the medium at all (can result in model unable to grow)
        "partially_constrained" : They are made available with an influx in the medium of 1 mmol.gDW^-1.h^-1

    Returns
    -------
    interaction_score : float
        Predicts the nature of the interaction between model1 and model 2. 
        Score < 0 predicts a competitive interaction,
        Score = 0 predicts a neutral interaction
        Score > 0 predicts a positive interaction
    interaction_type : string
        Qualitative description of the interaction.
    """

    metabolic_dict = utils.create_ecosystem_metabolic_dict(model1, model2)
    #Infers maximal objective value of both models seperately, in the given mdeium.
    with model1:
        model1, constrained_medium_dict1 = utils.restrain_medium(model1, medium, undescribed_metabolites_constraint)
        solo_growth_model1 = model1.optimize().objective_value
    with model2:
        model2, constrained_medium_dict2 = utils.restrain_medium(model2, medium, undescribed_metabolites_constraint)
        solo_growth_model2 = model2.optimize().objective_value
    if solo_growth_model1 == solo_growth_model2 == 0:
        raise RuntimeError("Both models had a null objective value when modeled alone in the given medium."+
                           " To enable this analysis, you need to adjust the medium or models. You can also"+
                           " try to lighten the medium constraint by using the \"partially_constrained\""+
                           " option for the undescribed_metabolites_constraint argument.") 
    elif solo_growth_model1 == 0 or solo_growth_model2 == 0:
        warnings.warn("One model had a null objective value when modeled alone in the given medium."+
                      " If this is not an expected result, you might want to use the \"partially_constrained\""+
                      " option for the undescribed_metabolites_constraint argument, or redefine your medium or model.")
    medium_dict = {**constrained_medium_dict1, **constrained_medium_dict2} #Translate medium constraint for mocba
    # mocba will create new exchange reaction exterior to both models. the original exchange reactions, if restrained, will prevent 
    #exchanges between organisms. Here we unconstrain them.
    model1 = utils.unrestrain_medium(model1)
    model2 = utils.unrestrain_medium(model2)
    sol_mofba = utils.mo_fba(model1, model2, metabolic_dict, medium_dict)[0] #get multi-objective solution (pareto front)
    xy, maxi_model1, maxi_model2 = utils.pareto_parsing(sol_mofba, solo_growth_model1, solo_growth_model2) #parse and normalize pareto front
    interaction_score = utils.infer_interaction_score(xy) #measure AUC of Pareto front, translates it in quantitative interaction prediction
    interaction_type = utils.infer_interaction_type(interaction_score, maxi_model1, maxi_model2,solo_growth_model1, solo_growth_model2) #infers interaction type from pareto front's shape.
    return interaction_score, interaction_type



def crossfed_metabolites(model1, model2, medium, undescribed_metabolites_constraint, solver, model1_biomass_id, model2_biomass_id):
    """
    A function that, given 2 models in the same namespace and a defined medium, predicts metabolic exchanges that
    are correlated with the increase of model2 objective value.

    Parameters
    ----------
    model1 : cobra.Model
    model2 : cobra.Model
    medium : pandas series
        Index : metabolites names
        values  : Availability of corresponding metabolite in the medium as a positive flux value. 
    undescribed_metabolites_constraint : string ("blocked" or "partially_constrained"). 
        How strictly constrained are the medium metabolites for which the flux is not described in the medium dataframe.
        "blocked" : They are not available in the medium at all (can result in model unable to grow)
        "partially_constrained" : They are made available with an influx in the medium of 1 mmol.gDW^-1.h^-1
    solver : string
        solver supported by the cobra toolbox. "cplex" or "gurobi" are recommended but require prior installation.
    model1_biomass_id : string
        id of the reaction used as objective in model1 (if the objective coefficient is not null for several reactions, 
                                                        then a new reaction must be built to constrain the model to a given 
                                                        objective value through its flux)
    model2_biomass_id : string
        id of the reaction used as objective in model2 (if the objective coefficient is not null for several reactions, 
                                                        then a new reaction must be built to constrain the model to a given 
                                                        objective value through its flux)
    Returns
    -------
    potential_crossfeeding : dictionnary
        keys : metabolites id
        values : [proportion of samples featuring inverse secretion/uptake for a same metabolite, 
        proportion of samples with metabolite exchange from model1 to model2, 
        proportion of samples with metabolite exchange from model2 to model1]
    """
    
    metabolic_dict = utils.create_ecosystem_metabolic_dict(model1, model2)
    #Infers maximal objective value of both models seperately, in the given mdeium.
    with model1:
        model1, constrained_medium_dict1 = utils.restrain_medium(model1, medium, undescribed_metabolites_constraint)
        solo_growth_model1 = model1.optimize().objective_value
    with model2:
        model2, constrained_medium_dict2 = utils.restrain_medium(model2, medium, undescribed_metabolites_constraint)
        solo_growth_model2 = model2.optimize().objective_value
    if solo_growth_model1 == solo_growth_model2 == 0:
        raise RuntimeError("Both models had a null objective value when modeled alone in the given medium."+
                           " To enable this analysis, you need to adjust the medium or models. You can also"+
                           " try to lighten the medium constraint by using the \"partially_constrained\""+
                           " option for the undescribed_metabolites_constraint argument.") 
    elif solo_growth_model1 == 0 or solo_growth_model2 == 0:
        warnings.warn("One model had a null objective value when modeled alone in the given medium."+
                      " If this is not an expected result, you might want to use the \"partially_constrained\""+
                      " option for the undescribed_metabolites_constraint argument, or redefine your medium or model.")
    medium_dict = {**constrained_medium_dict1, **constrained_medium_dict2} #Translate medium constraint for mocba
    # mocba will create new exchange reaction exterior to both models. the original exchange reactions, if restrained, will prevent 
    #exchanges between organisms. Here we unconstrain them.
    model1 = utils.unrestrain_medium(model1)
    model2 = utils.unrestrain_medium(model2)
    sol_mofba, ecosys = utils.mo_fba(model1, model2, metabolic_dict, medium_dict) #get multi-objective solution (pareto front), and ecosystem model from mocba
    xy, maxi_model1, maxi_model2 = utils.pareto_parsing(sol_mofba, solo_growth_model1, solo_growth_model2) #parse and normalize pareto front
    cobra_ecosys = utils.mocba_to_cobra(ecosys) #Translate mocba ecosystem into cobra.Model
    cobra_ecosys.solver = solver 
    model1_id = model1.id
    model2_id = model2.id
    sampling = utils.pareto_sampling(cobra_ecosys, xy, solo_growth_model1, solo_growth_model2, model1_id, model2_id, model1_biomass_id, model2_biomass_id, sample_size = 1000)
    correlation_reactions = utils.correlation(sampling)
    potential_crossfeeding = utils.crossfed_mets(model1 = model1, model1_id = model1_id, sampling = sampling, correlation_reactions = correlation_reactions, 
                                        model2_id = model2_id, model2_biomass_id=model2_biomass_id)
    return potential_crossfeeding