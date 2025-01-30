Usage
=======

Models and medium formats
-------------------------

**MIMEco** infers the interaction between two organisms in a given medium. These **two organisms** must be inputted as **cobrapy models**.
Cobrapy allows easy import of models of sbml, json, yaml, matlab and pickle formats. 

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

First, import your models and medium in the format described precedently. It is highly advised to define their solver as "cplex" or "gurobi" 
rather than the default glpk.

.. code-block:: python

    model1 = cobra.io.read_sbml_model("resources/model1.xml")
    model1.solcer = "cplex"
    model2 = cobra.io.read_sbml_model("resources/model2.xml")
    model2.solver = "cplex"

The medium is an important parameter, it will define the metabolic environment of the organisms. Interaction will change depending on the available nutrients.

.. code-block:: python

    medium = pd.read_csv("resources/Western_diet.csv", index_col = 0)

It is also important to precise the level of constraint for influx of metabolites that are not described in medium. 
When ``undescribed_metabolites_constraint = "blocked"``, any exchange reaction for a metabolite that is not in medium 
will have its **lower bound set to 0**, preventing the metabolite to be made available in the medium if not from organism's secretion.
When playing with various models and / or medium, it often prevents the models to have a non-null objective value.

``undescribed_metabolites_constraint = "partially_constrained"`` sets the **lower_bound of undescribed metabolites to -1**, allowing a limited influx in the medium.
Adequatly restraining important metabolites in medium creates limiting metabolites that will restrain organisms growth even with an imperfectly constrained medium.
However, this can impact predicted metabolic pathways, including interactions between models. Ideally, the medium should be 
complete enough to enable the modeled organisms to grow in a "blocked" context.

:py:func:`interaction_score_and_type()` function takes two models and a medium as variables, and an ``undescribed_metabolites_constraint`` level: 

.. code-block:: python

    interaction_score, interaction_type = 
    mimeco.interaction_score_and_type(model1, model2, medium, 
    undescribed_metabolites_constraint="partially_unconstrained")

Predicting cross-feeding between models
----------------------------------------

The function :py:func:`crossfed_metabolites()` predicts the metabolites exchanged between the two models, that are correlated with the increase of model2 objective value.
In other words, **exchanged metabolites favoring model2's objective** (usually, growth). This can help identify cross-feeding.

In addition to the precedently described inputs, this function necessitate the following elements :
* ``solver`` : solver that you use (advised : "cplex" or "gurobi")
* ``model1_biomass_id`` : id (str) of the reaction used as objective in model1 (if the objective coefficient 
is not null for several reactions then a new reaction must be built to constrain the model to a given 
objective value through its flux)
* ``model2_biomass_id`` : id (str) of the reaction used as objective in model2 (if the objective coefficient 
is not null for several reactions then a new reaction must be built to constrain the model to a given 
objective value through its flux)

.. code-block:: python

    potential_crossfeeding = crossfed_metabolites(model1, model2, 
    medium, undescribed_metabolites_constraint, solver, 
    model1_biomass_id, model2_biomass_id)

The output is a dictionnary formatted as :
``{metabolic id : [proportion of samples featuring inverse secretion/ uptake for given metabolite, 
                  proportion of samples with metabolite exchange from model1 to model2,
                  proportion of samples with metabolite exchange from model2 to model1]}``

As the selected metabolites are the one favoring model2, it is interesting to run the function twice while inversing models position.

See <Practical example> for an application of both function and interprtation of results.
