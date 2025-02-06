#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 16:11:15 2025

@author: Anna Lambert
"""
import pandas as pd

blood_reactions = pd.read_csv("resources/blood_reactions_id_enterocyte.csv")
AAD_medium_blood = pd.read_csv("resources/AAD_notebook.tsv", sep="\t", index_col = 0)

# AAD Must be lb ub because blood is different... will need different format. (two columns)
#!! blood exchange reaction in the other way as lumen, and general convention. 


def restrain_blood_exchange_enterocyte(model, medium_blood = None):
    """
    Restrains exchanges with the blood compartment for the enterocyte. 
    Default constraint : Average American diet (AAD) from https://doi.org/10.1093/hmg/ddt119
    Function and medium built for cobra.Model of small intestinal epithelial cell adapted from https://doi.org/10.1093/hmg/ddt119

    Parameters
    ----------
    model : cobra.Model
        should be small intestinal epithelial cell adapted from https://doi.org/10.1093/hmg/ddt119
    medium_blood : pandas.DataFrame, optional
        A pandas.DataFrame defining blood exchange constraints for the enterocyte
        Index : Exchange reactions with the blood
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
    
    #Constrain exchanges with blood compartment
    if medium_blood == None:
        medium_blood = AAD_medium_blood
        for reac in medium_blood.index:
            model.reactions.get_by_id(reac).bounds = (AAD_medium_blood.loc[reac,'lb'],AAD_medium_blood.loc[reac,'ub'])
    else:
        for met in medium_blood.index:
            ex_reac_id = "EX_"+met+"(e)"
            if ex_reac_id in blood_reactions:
                ex_reac = model.reactions.get_by_id(ex_reac_id)
                ex_reac.bounds = (medium_blood.loc[ex_reac,'lb'], medium_blood.loc[ex_reac,'ub'])
    return model

#restrain medium works the same with enterocyte and bacteria. 