#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 16:11:15 2025

@author: Anna Lambert
"""
import pandas as pd
import mimeco.utils as utils

blood_reactions = pd.read_csv("mimeco/resources/blood_reactions_id_enterocyte.csv")
blood_reactions = blood_reactions["blood_exchanges"].tolist()

# AAD Must be lb ub because blood is different... will need different format. (two columns)
#!! blood exchange reaction in the other way as lumen, and general convention. 


def restrain_blood_exchange_enterocyte(model, namespace, medium_blood = "AAD"):
    """
    Restrains exchanges with the blood compartment for the enterocyte. 
    Default constraint : Average American diet (AAD) from https://doi.org/10.1093/hmg/ddt119
    Function and medium built for cobra.Model of small intestinal epithelial cell adapted from https://doi.org/10.1093/hmg/ddt119

    Parameters
    ----------
    model : cobra.Model
        should be small intestinal epithelial cell adapted from https://doi.org/10.1093/hmg/ddt119
    namespace : string
        "BIGG" : enterocyte and medium in the BiGG namespace. Compatible with CarveMe.
        "AGORA" : enterocyte and medium in the Agora namespace: Compatible with Agora and VMH models. (Built with Model SEED / Kbase)
    medium_blood : pandas.DataFrame, optional
        A pandas.DataFrame defining blood exchange constraints for the enterocyte
        Index : Exchanged metabolites with the blood (excpet default AAD where it is exchange reactions)
        column 1 : header = "lb", lower_bound to constrain the reaction with
        column 2 : header = "ub", upper_bound to constrain the reaction with
        NOTE : blood exchange reaction are written in the opposite direction as usual exchange reaction. 
        As a result, a negative flux is an export flux (from the cell to the blood) and a positive flux is an import flux (from the blood to the cell). 
        default : None; In this case, applies the Average American diet (AAD) from https://doi.org/10.1093/hmg/ddt119
    Returns
    -------
    model : cobra.Model
        sIEC with blood exchanges constrained
    """
    if namespace == "BIGG":
        AAD_medium_blood = pd.read_csv("mimeco/resources/AAD_BiGG.tsv", sep="\t", index_col = 0)
        #blood_suffixe = "(e)"
    elif namespace == "AGORA":
        AAD_medium_blood = pd.read_csv("mimeco/resources/AAD_VMH.tsv", sep="\t", index_col = 0)
        #blood_suffixe = "_b"
    #Constrain exchanges with blood compartment
    if isinstance(medium_blood, str):
        if medium_blood == "AAD": #default medium AAD
            medium_blood = AAD_medium_blood
            for reac in medium_blood.index:
                model.reactions.get_by_id(reac).bounds = (AAD_medium_blood.loc[reac,'lb'],AAD_medium_blood.loc[reac,'ub'])
        else:
            raise TypeError("The inputted blood medium is a string when it should be a pandas.DataFrame with exchanged metabolites as index"+ 
                             "and lb and ub as columns. You can also keep the default value \"AAD\" to use an average american diet as constraints.")
    else: #custom medium
        for blood_reac in blood_reactions:
            ex_met = utils.no_compartment_id(model.reactions.get_by_id(blood_reac).reactants[0].id)[0]
            if ex_met in medium_blood.index:
                ex_reac = model.reactions.get_by_id(blood_reac)
                ex_reac.bounds = (medium_blood.loc[ex_met,'lb'], medium_blood.loc[ex_met,'ub'])
            else:
                ex_reac = model.reactions.get_by_id(blood_reac)
                ex_reac.bounds = (0,1000) #undefined metabolites can be secreted in the blood but not absorbed by the enterocyte.
    return model