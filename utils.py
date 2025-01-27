#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 10:37:38 2025

@author: e158401a
"""

import pandas as pd
from sklearn import metrics
import mocbapy
from mocbapy.EcosystemModel import create_model, bensolve_default_options
import mocbapy.analysis
import math


def create_ecosystem_metabolic_dict(model1, model2):
    """
    Builds a dictionnary for how to name the metabolites of each inputted model in the ecosystem model. 
    Here, we keep the same ids. 

    Parameters
    ----------
    model1 : cobra.Model
    model2 : cobra.Model

    Returns
    -------
    metabolic_dict : dictionnary {(met_id, model): met_id}
        Dictionnary of all metabolites of the inputted models, and how to name them in the ecosystem model. 
    """

    model1_metabolites = []
    for met in model1.metabolites:
        model1_metabolites.append((met.id, model1))
    model2_metabolites = []
    for met in model2.metabolites:
        model2_metabolites.append((met.id, model2))
    all_metabolites = list(set(model1_metabolites + model2_metabolites))
    metabolic_dict = {x:x[0] for x in all_metabolites}
    return metabolic_dict

def restrain_medium(model, medium, undescribed_metabolites_constraint):
    """
    Builds the dictionnary used for constraining the medium of the ecosystem model based on inputted medium data.

    Parameters
    ----------
    model : cobra.Model
    medium : pandas series
        Index : metabolites names
        values  : Availability of corresponding metabolite in the medium as a positive flux value. 
    undescribed_metabolites_constraint : string ("blocked" or "partially_constrained"). 
        How strictly constrained are the medium metabolites for which the flux is not described in the medium dataframe.
        "blocked" : They are not available in the medium at all (can result in model unable to grow)
        "partially_constrained" : They are made available with an influx in the medium of 1 mmol.gDW^-1.h^-1

    Returns
    -------
    model : cobra.Model
    constrained_medium_dict : dictionnary {met_id:(lower_bound,upper_bound)}
        Guide constraint of the ecosystem medium based on inputted medium data. 
        Controls fluxes of metabolites entering the emodels external environment.
    """

    constrained_medium_dict = {}
    for reac in model.exchanges:
        met_ex, suffixe = no_compartment_id(list(reac.metabolites.keys())[0].id)
        if met_ex in list(medium.index):
            constrained_medium_dict[met_ex+suffixe] = (-medium.loc[met_ex], reac._upper_bound)
            reac.lower_bound = -medium.loc[met_ex]
        elif undescribed_metabolites_constraint == "blocked":
            constrained_medium_dict[met_ex+suffixe] = (0, reac._upper_bound)
            reac.lower_bound = 0
        elif undescribed_metabolites_constraint == "partially_constrained":
            constrained_medium_dict[met_ex+suffixe] = (-1, reac._upper_bound)
            reac.lower_bound = -1
    if constrained_medium_dict == {}:
        warnings.warn("The inputted medium constraint does not match the model's namespace. The medium could not be applied to the ecosystem model")
    return model, constrained_medium_dict

def unrestrain_medium(model):
    """
    Opens the exchange reaction of individual models to enable exchange with paired model. 
    Environment constraints will be reaplied at the scale of the ecosystem in the construction
    of the ecosystem model.

    Parameters
    ----------
    model : cobra.Model

    Returns
    -------
    model : cobra.Model
    """

    for reac in model.exchanges:
        reac.bounds = (-1000, 1000)
    return model

def mo_fba(model1, model2, metabolic_dict, constrained_medium_dict):
    """
    Compute multi-objective FBA between the two given models
    
    Parameters
    ----------
    model1 : cobra.Model
    model2 : cobra.Model
    metabolic_dict : dictionnary {(met_id, model): met_id}
        Dictionnary of all metabolites of the inputted models, and how to name them in the ecosystem model. 
    constrained_medium_dict : dictionnary {met_id:(lower_bound,upper_bound)}
        Guide constraint of the ecosystem medium based on inputted medium data. 
        Controls fluxes of metabolites entering the emodels external environment.

    Returns
    -------
    sol_mofba : 
        Multi-objective solution (Pareto front) of the ecosystem model
    """

    model1 = unrestrain_medium(model1)
    model2 = unrestrain_medium(model2)
    ecosys = create_model(model_array=[model1, model2], metabolic_dict=metabolic_dict, diet = constrained_medium_dict)
    bensolve_opts = bensolve_default_options()
    bensolve_opts['message_level'] = 0
    sol_mofba = mocbapy.analysis.mo_fba(ecosys, options=bensolve_opts)
    return sol_mofba, ecosys

def pareto_parsing(sol_mofba, solo_growth_model1, solo_growth_model2):
    """
    Parse the Pareto front and returns its points normalized by the optimal 
    growth of each model optimized in isolation

    Parameters
    ----------
    sol_mofba : 
        Multi-objective solution (Pareto front) of the ecosystem model
    solo_growth_model1 : float
        Objective value of model1 when optimized alone in the described medium
    solo_growth_model2 : float
        Objective value of model2 when optimized alone in the described medium

    Returns
    -------
    xy : pandas dataframe
    Normalized pareto points
    maxi_model1 : float
        maximal objective value of model1 in the ecosystem model
    maxi_model2 : float
        maximal objective value of model1 in the ecosystem model
    """

    x = []
    y = []
    #Initialize analysis variables
    maxi_model1 = (solo_growth_model1, 0)
    maxi_model2 = (0, solo_growth_model2)
    for i in range(len(sol_mofba.Primal.vertex_value)): #Parse Pareto front points
        if sol_mofba.Primal.vertex_type[i] == 1: #Select points (vs vectors) of the Pareto front
            #Normalize points based on model solo growth
            x.append(sol_mofba.Primal.vertex_value[i][0]/solo_growth_model1) 
            y.append(sol_mofba.Primal.vertex_value[i][1]/solo_growth_model2) 
            
            #Evaluate if model grows better in ecosystem compared to solo model
            if sol_mofba.Primal.vertex_value[i][0] > maxi_model1[0]: 
                maxi_model1 = sol_mofba.Primal.vertex_value[i]
            if sol_mofba.Primal.vertex_value[i][1] > maxi_model2[1]:
                maxi_model2 = sol_mofba.Primal.vertex_value[i]

    #Add infered ecosystem Pareto points
    xy = pd.DataFrame({'x': x, 'y': y})
    xy.sort_values('x', inplace=True)
               
    #Add initial points, corresponding to solo models optimal objective values
    #Values added are slightly dfferent than 1 and 0 to make sure the serie of coordinate continuous
    if not ((xy['x'] == 1) & (xy['y'] == 0)).any():
        xy = xy.append({'x' : 1.00001, 'y' : -0.00001}, ignore_index=True) 
    if not ((xy['x'] == 0) & (xy['y'] == 1)).any():
        xy.loc[-1] = [-0.00001, 1.00001]
        xy.index = xy.index + 1  # shifting index
        xy = xy.sort_index()  # sorting by index
    
    return xy, maxi_model1, maxi_model2

def infer_interaction_score(xy):
    """
    Calculates pareto front's area under the score to determine the interaction score of the ecosystem.
    
    Parameters
    ----------
    xy : pandas dataframe
    Normalized pareto points

    Returns
    -------
    interaction_score : float
        Predicts the nature of the interaction between model1 and model 2. 
        Score < 0 predicts a competitive interaction,
        Score = 0 predicts a neutral interaction
        Score > 0 predicts a positive interaction
    """

    try:
        AUC = metrics.auc(x = xy['x'], y = xy['y'])
    except: 
    #If the growth alone is lower than with paired model, x is not monotonous since it goes up and then down. AUC determination takes one more step.
    #We calculate the AUC where x is monotonous, the AUC of the inverted part of the Pareto, and substract the last part from the first.
        AUC1 = metrics.auc(x = xy['x'][0:-1], y = xy['y'][0:-1])
        AUC2 = metrics.auc(x = xy['x'][-2:], y = xy['y'][-2:])
        AUC = AUC1 - AUC2
    interaction_score = AUC-1
    return interaction_score
    

def infer_interaction_type(interaction_score, maxi_model1, maxi_model2, solo_growth_model1, solo_growth_model2):
    """
    Infers inetraction type from the models solo maximal objective values, ecosystem maximal objective value and interaction score.

    Parameters
    ----------
    interaction_score : float
        Predicts the nature of the interaction between model1 and model 2. 
        Score < 0 predicts a competitive interaction,
        Score = 0 predicts a neutral interaction
        Score > 0 predicts a positive interaction
    maxi_model1 : float
        maximal objective value of model1 in the ecosystem model
    maxi_model2 : float
        maximal objective value of model2 in the ecosystem model
    solo_growth_model1 : float
        Objective value of model1 when optimized alone in the described medium
    solo_growth_model2 : float
        Objective value of model2 when optimized alone in the described medium

    Returns
    -------
    interaction_type : string
        Qualitative description (competition, neutrality, favors model1, favors model2,
        limited mutualism or mutualism) of the metabolic interaction between models
    """

    interaction_type_code = ["0","0","0"]
    #Evaluate parameters defining interaction category:
        #interaction_type[0] : model1 growth better in ecosystem than alone
        #interaction_type[1] : model2 growth better in ecosystem than alone
        #interaction_type[2] : Both models share an optimal solution
    if maxi_model1[0] > solo_growth_model1+(0.01*solo_growth_model1): #solution is an approximation, so slightly varies between instances, 
                                                                      #we make sure the difference in growth is not an artefact
        interaction_type_code[0] = "1"        
    if maxi_model2[1] > solo_growth_model2+(0.01*solo_growth_model2):
        interaction_type_code[1] = "1"
    if tuple(maxi_model1) == tuple(maxi_model2) :
        interaction_type_code[2] = "1"
    interaction_type_code=''.join(interaction_type_code)
    if interaction_score < -0.0001 and interaction_type_code == "000": #Make sure score is negative despite approximation : competition
        interaction_type_code = "-000"
    elif interaction_score > -0.0001 and interaction_score < 0.0001: #Make sure score is equivalent to 0 : neutrality
        interaction_type_code = "=000" 
    if interaction_type_code not in ["-000", "=000", "100","010","110", "111"]:
        raise RuntimeError("There was a problem while infering interaction_type. It is probably in the definition of the model or medium.")
    interaction_type_translation = {"-000":"competition", "=000": "neutrality", 
                                    "100":"favors model1", "010":"favors model2",
                                    "110":"limited mutualism", "111":"mutualism"}
    interaction_type = interaction_type_translation[interaction_type_code]
    return interaction_type

def mocba_to_cobra(ecosys):
    """
    Converts the ecosystem model built through mocbapy to a cobra.Model

    Parameters
    ----------
    ecosys : mocbapy.EcosystemModel

    Returns
    -------
    cobra_model : cobra.Model
    """

    cobra_model = cobra.Model('ecosys')
    for m in range(len(ecosys.sysmetabolites)):
        cobra_model.add_metabolites(cobra.Metabolite(ecosys.sysmetabolites[m]))
    for r in range(len(ecosys.sysreactions)):
        reaction = cobra.Reaction(ecosys.sysreactions[r])
        reaction.lower_bound = ecosys.lb[r]
        reaction.upper_bound = ecosys.ub[r]
        dict_metabolites = {}
        cobra_model.add_reactions([reaction])
        for m in range(len(ecosys.sysmetabolites)):
            if ecosys.Ssigma[m,r] != 0:
                dict_metabolites[ecosys.sysmetabolites[m]] = ecosys.Ssigma[m,r]
        reaction.add_metabolites(dict_metabolites)
    return cobra_model

def pareto_sampling(cobra_ecosys, xy, solo_growth_model1, solo_growth_model2, model1_id, model2_id, model1_biomass_id, model2_biomass_id, sample_size = 1000):
    """
    Samples the Pareto front, infering a solution for <sample_size> points on the pareto front.

    Parameters
    ----------
    cobra_ecosys : mocbapy.EcosystemModel
    xy : pandas dataframe
    Normalized pareto points
    solo_growth_model1 : float
        Objective value of model1 when optimized alone in the described medium
    solo_growth_model2 : float
        Objective value of model2 when optimized alone in the described medium
    model1_id : string
        Model denomination in the cobra.Model of model1
    model2_id : string
        Model denomination in the cobra.Model of model2
    model1_biomass_id : string
        id of the reaction used as objective in model1 (if the objective coefficient is not null for several reactions, 
                                                        then a new reaction must be built to constrain the model to a given 
                                                        objective value through its flux)
    model2_biomass_id : string
        id of the reaction used as objective in model2 (if the objective coefficient is not null for several reactions, 
                                                        then a new reaction must be built to constrain the model to a given 
                                                        objective value through its flux)
    sample_size : int, optional
        Number of samples sampled from the Pareto front to infer correlation between exchange reactions and biomass. The default is 1000.

    Returns
    -------
    sampling : pandas.dataframe
        columns : reactions_id
        rows : string(objective-value-model1_objective-value-model2) for a given sample
    """

    fba = cobra_ecosys.optimize() #just to get the reactions names for index of sampling
    sampling_dict = {}
    sampling_dict["reactions"] = fba.fluxes.index
    sumdist = 0    
    # Dividing pareto front into segments between two extreme points. Denormalize pareto point so that it can become objective_value constraints.
    # Measure total distance of the pareto front, to determine constant distance between samples points to represent homogeneous representation of pareto front metabolic phenotypes.
    for p in xy.index[1:]:
        x1 = xy.loc[p-1, 'x'] * solo_growth_model1
        y1 = xy.loc[p-1, 'y'] * solo_growth_model2
        x2 = xy.loc[p, 'x'] * solo_growth_model1
        y2 = xy.loc[p, 'y'] * solo_growth_model2
        dist = math.hypot(x2-x1, y2-y1)
        sumdist = sumdist + dist
    dist_1pt = sample_size/sumdist
    #Then browse the pareto again to do the sampling.
    for p in xy.index[1:]:
        x1 = xy.loc[p-1, 'x'] * solo_growth_model1
        y1 = xy.loc[p-1, 'y'] * solo_growth_model2
        x2 = xy.loc[p, 'x'] * solo_growth_model1
        y2 = xy.loc[p, 'y'] * solo_growth_model2
        dist = math.hypot(x2-x1, y2-y1)
        nbpnts_xy = round(dist*dist_1pt)
        if nbpnts_xy>0:
            dx = (x2-x1)/nbpnts_xy
            dy = (y2-y1)/nbpnts_xy
            # Sample a segment of the pareto front.
            # Constraint biomass reaction flux value based on pareto front points coordinates (pix, piy).
            for i in range (0, nbpnts_xy):
                pix = x1+dx*i
                piy = y1+dy*i
                cobra_ecosys.reactions.get_by_id(model1_biomass_id+":"+model1_id).bounds = (pix, pix)
                cobra_ecosys.reactions.get_by_id(model2_biomass_id+":"+model2_id).bounds = (piy, piy)
                fba = cobra_ecosys.optimize() #One sample
                sampling_dict[str(pix)+"_"+str(piy)] = list(fba.fluxes) #keep objective values information in sampling dataframe      
    sampling = pd.DataFrame(data=sampling_dict)
    sampling.index = sampling["reactions"]
    sampling = sampling.T
    sampling = sampling.drop("reactions")
    return sampling

def correlation(sampling):
    """
    Measure correlation between all reactions of the ecosystem model. 
    Correlation type : spearman, because we are interested in common evolution, not proportionnality.

    Parameters
    ----------
    sampling : pandas.dataframe
        columns : reactions_id
        rows : string(objective-value-model1_objective-value-model2) for a given sample

    Returns
    -------
    correlation_reactions : pandas dataframe
        a correlation matrix featuring all reactions of the ecosystem model
    """

    sampling = sampling.astype(float)
    correlation_reactions = sampling.corr(method ='spearman')
    correlation_reactions = correlation_reactions.fillna(0)
    return correlation_reactions

def oppositeSigns(x, y): #returns boolean. 1 if x y opposite sign, 0 otherwhise
    #return (y >= 0) if (x < 0) else (y < 0) 
    if x < 0:
        return (y > 0) # Because I want a net import 
    elif x > 0: # Because I want a net import 
        return (y < 0)
    else:
        return False
 
def reac_to_met_id(reac, model_id):
    """
    Returns exchanged metabolic id from an exchange reaction. 

    Parameters
    ----------
    reac : cobra.Reaction.id
    model_id : string
        Model denomination in the cobra.Model

    Returns
    -------
    met : string
        metabolite id
    """

    met = reac.replace("_e:"+model_id, "")
    met = met.replace("EX_", "")
    return met
    
def exchanged_mets(model1, model1_id, sampling, correlation_reactions, model2_id, model2_biomass_id):
    """
    Infers metabolites that are exchanged between organisms in the ecosystem models, correlated with an increasing model2 objective value.

    Parameters
    ----------
    model1 : cobra.Model
    model1_id : cobra.Model
    sampling : pandas.dataframe
        columns : reactions_id
        rows : string(objective-value-model1_objective-value-model2) for a given sample
    correlation_reactions : pandas dataframe
        a correlation matrix featuring all reactions of the ecosystem model
    model2_id : string
        Model denomination in the cobra.Model of model2
    model2_biomass_id : string
        id of the reaction used as objective in model2 (if the objective coefficient is not null for several reactions, 
                                                        then a new reaction must be built to constrain the model to a given 
                                                        objective value through its flux)

    Returns
    -------
    potential_exchange : dictionnary
        keys : metabolites id
        values : [nb of samples featuring inverse secretion/uptake for a same metabolite, direction of the exchange]
    """

    potential_exchange = {}
    for ex_reac in model1.exchanges:      
        ecosys_reac_id_model1 = ex_reac.id+":"+model1_id
        ecosys_reac_id_model2 = ex_reac.id+":"+model2_id
        met_id = reac_to_met_id(ecosys_reac_id_model1, model1_id)      
        #if a metabolite has an exchange reaction in both models
        if ecosys_reac_id_model1 in correlation_reactions.index and ecosys_reac_id_model2 in correlation_reactions.index:
            #if both exchange reactions have at least one non-null flux value among all samples.
            #and if both reactions are inversely correlated (fluxes variation are going opposite ways, one toward secretion, the other toward uptake)
            if sum(sampling[ecosys_reac_id_model1])!=0 and sum(sampling[ecosys_reac_id_model2])!=0 and correlation_reactions.loc[ecosys_reac_id_model1, ecosys_reac_id_model2] <= -0.5:
                # If the uptake / secretion of given metabolite in model1, associated with its secretion / uptake in model2, is correlated with increased model2 objective value
                if correlation_reactions.loc[ecosys_reac_id_model1, model2_biomass_id+":"+model2_id] > 0.8:
                    exchange = 0
                    way = ["0","0"]
                    for s in sampling.index: #parse all solutions for metabolite of interest
                        # if metabolite is secreted in one model, and uptaken in the other
                        if oppositeSigns(sampling.loc[s, ecosys_reac_id_model1], sampling.loc[s, ecosys_reac_id_model2]) :
                            exchange = exchange+1
                            if sampling.loc[s, ecosys_reac_id_model1] <0 :
                                way[0] = "1"
                            elif sampling.loc[s, ecosys_reac_id_model2]<0:
                                way[1] = "1"
                    if exchange >10 and met_id not in potential_exchange.keys():
                        potential_exchange[met_id] = [exchange, way]
    return potential_exchange


    
def no_compartment_id(metabolite):
    """
    Separates metabolic id from the compartment information. Aknowledges BiGG and AGORA namespaces.

    Parameters
    ----------
    metabolite : cobra.Metabolite

    Returns
    -------
    metabolite : string
        metabolite_id, without its compartment information
    suffixe : string
        compartment information of the inputted metabolite
    """
    
    a = metabolite
    done = False
    if metabolite[-2:] == "_e":
        metabolite = metabolite.replace('_e','')
        if metabolite != a:
            suffixe='_e'
            done = True
    metabolite = metabolite.replace(('(e)'),'')
    if metabolite != a and not done:
        suffixe='(e)'
        done = True
    metabolite = metabolite.replace(('_c'),'')
    if metabolite != a and not done:
        suffixe='_c'
        done = True
    metabolite = metabolite.replace(('[c]'),'')
    if metabolite != a and not done:
        suffixe='[c]'
        done = True
    metabolite = metabolite.replace(('[m]'),'')
    if metabolite != a and not done:
        suffixe='[m]'
        done = True
    metabolite = metabolite.replace(('[x]'),'')
    if metabolite != a and not done:
        suffixe='[x]'
        done = True
    metabolite = metabolite.replace(('[r]'),'')
    if metabolite != a and not done:
        suffixe='[r]'
        done = True
    metabolite = metabolite.replace(('[e]'),'')
    if metabolite != a and not done:
        suffixe='[e]'
        done = True
    metabolite = metabolite.replace(('[n]'),'')
    if metabolite != a and not done:
        suffixe='[n]'
        done = True
    metabolite = metabolite.replace(('[u]'),'')
    if metabolite != a and not done:
        suffixe='[u]'
        done = True
    metabolite = metabolite.replace(('[luC]'),'')
    if metabolite != a and not done:
        suffixe='[luC]'
        done = True
    metabolite = metabolite.replace(('[bpC]'),'')
    if metabolite != a and not done:
        suffixe='[bpC]'
        done = True
    metabolite = metabolite.replace(('[g]'),'')
    if metabolite != a and not done:
        suffixe='[g]'
        done = True
    metabolite = metabolite.replace(('[l]'),'')
    if metabolite != a and not done:
        suffixe='[l]'
        done = True
    if done == False: 
        print("no suffixe found for", metabolite)
        suffixe = ""
    return metabolite, suffixe