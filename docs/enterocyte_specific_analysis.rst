Enterocyte specific analysis
============================

If you are interested in the interaction of an organism with the human small intestine, the global usage has been adapted in two specific functions. 
In this case, one inputted model is joined to a small intestinal epithelial cell (sIEC) in a pairwise ecosystem. 

Enterocyte model
----------------
The sIEC is adapted from the one built in "Predicting the impact of diet and enzymopathies on human small intestinal epithelial cells", Sahoo et al, 2013,
It has been translated from matlab to sbml, and its namespace has been adapted to the **BiGG's namespace**, which is **compatible with CarveMe reconstruction**. Therefore, the model inputted should be built in the BIGG's namespace to be compatible. An option enabling the analysis in the AGORA / Virtual Metabolic Human namespace is in development, but requires further testing to be reliable.

Infering interaction score and type
------------------------------------
The arguments model, medium and undescribed_metabolites_constraint are the same as in global usage. 
The sIEC has two external compartments : the lumen and the blood. The constraints inputted in the argument "medium" are applied to the lumen compartment, as it is the medium shared between the enterocyte and the other organism. The blood medium is by default constrained based on an "Average American Diet (AAD) extracted from from https://doi.org/10.1093/hmg/ddt119. 

.. code-block:: python

  int_score, int_type = mimeco.analysis.enterocyte_interaction_score_and_type(model1, Western_diet,
                                                                    undescribed_metabolites_constraint="as_is", solver = "gurobi")

* The optional argument ``medium_blood`` allows you to input your custom blood medium. It should be a pandas.DataFrame of the following format:

   - Index : Exchanged metabolites with the blood (except default AAD where it is exchange reactions)
   - column 1 : header = "lb", lower_bound to constrain the reaction with
   - column 2 : header = "ub", upper_bound to constrain the reaction with

* The optional argument ``plot`` is set to "False" by default. When set to "True", the function will show a matplotlib plot of the Pareto front to ease the analysis. 

.. code-block:: python

  int_score, int_type = mimeco.enterocyte_interaction_score_and_type(model1, Western_diet, 
                                                                    undescribed_metabolites_constraint="as_is",
                                                                    solver = "gurobi" 
                                                                    namespace="BIGG", plot = False)

Predicting cross-feeding between models
---------------------------------------

This function predicts exchanged metabolites that favor enterocyte maintenance. Its options are the same as in the corresponding functions in global usage. 

.. code-block:: python
  
  potential_crossfeeding = mimeco.analysis.enterocyte_crossfed_metabolites(model = model1, medium = Western_diet, undescribed_metabolites_constraint = "as_is", 
                                                                solver = "gurobi", model_biomass_id = model1_biomass_id, namespace = "BIGG", 
                                                                plot = False, sample_size = 1000)

