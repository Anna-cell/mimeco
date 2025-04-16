Global usage
============

Models and medium formats
-------------------------

**MIMEco** infers the interaction between two organisms in a given medium. These **two organisms** must be inputted as **cobrapy models**.
Cobrapy allows easy import of models of sbml, json, yaml, matlab and pickle formats, as described in `their documentation  <https://cobrapy.readthedocs.io/en/latest/io.html>`_.

The **medium** must be inputted in the form of a **pandas Series** where the index is the medium metabolite identifier, 
and the corresponding value is the availability of corresponding metabolite in the medium as a positive flux value (float).

.. list-table:: Example of medium
   :widths: 25 25
   :header-rows: 1
   
   * - index
     - Influx (mmol.gDW-1.h-1)
   * - ala__L
     - 0.3
   * - arachd
     - 0.00386729
   * - chsterol
     - 0.0575975
   * - glc__D
     - 0.294868
   * - ...
     - ...


Infering interaction score and type
------------------------------------

First, import your models and medium in the format described previously. It is highly advised to define their solver as "cplex" or "gurobi" 
rather than the default glpk. 

.. code-block:: python

    model1 = cobra.io.read_sbml_model("resources/model1.xml")
    model1.solver = "gurobi"
    model2 = cobra.io.read_sbml_model("resources/model2.xml")
    model2.solver = "gurobi"

The medium is an important parameter, it will define the metabolic environment of the organisms. The interaction will change depending on the available nutrients.

.. code-block:: python

    medium = pd.read_csv("resources/Western_diet.csv", index_col = 0)

It is also important to precise the **level of constraint for influx of metabolites that are not described in the medium**. 

* ``undescribed_metabolites_constraint = "blocked"``, any exchange reaction for a metabolite that is not described in the medium will have its **lower bound set to 0**, preventing the metabolite to be made available in the medium if not from organism's secretion. When playing with various models and / or medium, it often prevents the models from having a non-null objective value.

* ``undescribed_metabolites_constraint = "partially_constrained"`` sets the **lower_bound of undescribed metabolites to -1**, allowing a limited influx in the medium. Adequately restraining important metabolites in a medium creates limiting metabolites that will restrain organisms' growth even with an imperfectly constrained medium. However, this can impact predicted metabolic pathways, including interactions between models. Ideally, the medium should be complete enough to enable the modeled organisms to grow in a "blocked" context.

* ``undescribed_metabolites_constraint = "as_is"`` keeps the originally set bounds for undefined metabolites exchange reactions.

:py:func:`interaction_score_and_type()` function takes two models and a medium as variables, and an ``undescribed_metabolites_constraint`` level. You must also state the solver you use: 

.. code-block:: python

    interaction_score, interaction_type = 
    mimeco.analysis.interaction_score_and_type(model1, model2, medium, 
    undescribed_metabolites_constraint="partially_unconstrained", solver="gurobi")

Options
~~~~~~~
The optional argument ``plot`` is set to "False" by default. When set to "True", the function will show a matplotlib plot of the Pareto front to ease the analysis. 

.. code-block:: python

    interaction_score, interaction_type = 
    mimeco.analysis.interaction_score_and_type(model1, model2, medium, 
    undescribed_metabolites_constraint="partially_unconstrained", plot = True)


Predicting cross-feeding between models
----------------------------------------

The function :py:func:`crossfed_metabolites()` predicts the metabolites exchanged between the two models, which are correlated with the increase of model2 objective value.
In other words, **exchanged metabolites favoring model2's objective** (usually, growth). This can help identify cross-feeding.

In addition to the precedently described inputs, this function necessitates the following elements:

* ``solver``: solver that you use (advised : "cplex" or "gurobi")

* ``model1_biomass_id``: id (str) of the reaction used as objective in model1 (if the objective coefficient is not null for several
  reactions then a new reaction must be built to constrain the model to a given objective value through its flux)

* ``model2_biomass_id``: id (str) of the reaction used as objective in model2 (if the objective coefficient is not null for several
  reactions then a new reaction must be built to constrain the model to a given objective value through its flux)

.. code-block:: python

    potential_crossfeeding = mimeco.analysis.crossfed_metabolites(model1, model2, 
    medium, undescribed_metabolites_constraint, solver, 
    model1_biomass_id, model2_biomass_id)

The output is a dictionnary formatted as :

.. code-block:: python

    {metabolic id : 
    [proportion of samples featuring inverse secretion/ uptake for given metabolite,
    proportion of samples with metabolite exchange from model1 to model2,
    proportion of samples with metabolite exchange from model2 to model1]}

As the selected metabolites are the ones favoring model2, it is interesting to run the function twice while inversing models positions.

Options
~~~~~~~

* The optional argument ``plot`` is set to **"False"** by default. When set to **"True"**, the function will show matplotlib plots of the exchanges of crossfed metabolites along the Pareto front. See <Practical example> for illustration. 

* The optional argument ``sample_size`` is set to **1000** by default. It is the amount of solutions sampled along the Pareto front, on which the crossfeeding analysis depends. 

* The optional argument ``retrieve_data`` is set to **"no"** by default. 
  
  When set to **"selection"**, the function 
  returns two variables: the potential_crossfeeding dictionnary and relevant data in the form of a pandas.DataFrame. 
  This dataFrame contains the flux of exchange reactions of interest in each sampled solution on the Pareto front. 
  Reactions of interest are exchange reaction for a metabolite predicted as crossfed in both organisms.
  
  When set to **"all"**, the function returns two variables: the potential_crossfeeding dictionnary and a Dataframe containing the sampling results for every reactions of the ecosystem model.

* The optional argument ``exchange_correlation`` is set to 0.5 by default. Defines the threshold for a correlation between secretion and uptake of a same metabolite by paired models for this metabolite to be considered exchanged between models.

* The optional argument ``biomass_correlation`` is set to 0.8 by default. Defines the correlation threshhold between the exchange of the metabolite and the biomass production of model2 for its selection as crossfed.

.. code-block:: python

   potential_crossfeeding = crossfed_metabolites(model1, model2, 
   medium, undescribed_metabolites_constraint, solver, 
   model1_biomass_id, model2_biomass_id, 
   plot = True, sample_size = "10000", retrieve_data = True)


See <Practical example> for an application of both function and interpretation of results.
