MIMEco : Metabolic Interaction Modeling in Ecosystems
=======================================================

MIMEco is python package that explores communities metabolic interactions using multi-objective linear programming on GEnome-scale Metabolic models (GEMs). 
This tool automates the inference of interaction type, interaction score and exchanged metablites between two models in a given condition.

The concept of the methodology is described in Lambert, A., Budinich, M., Mahé, M., Chaffron, S., & Eveillard, D. (2024). Community metabolic modeling of host-microbiota interactions through multi-objective optimization. Iscience, 27(6). (http://doi.org/10.1016/j.isci.2024.110092)

Note : A technical note dedicated to this package will be published

.. image:: mimeco/resources/MIMEco_logo.png
  :width: 500
  :alt: MIMEco logo

Documentation
~~~~~~~~~~~~~

For detailed documentation, please go to `readthedocs : mimeco <https://mimeco.readthedocs.io/en/latest/>`_

Dependancies
~~~~~~~~~~~~~

**GLPK**:
MIMEco depends on benpy, which needs glpk to function. Its installation is clearly described in `benpy's pyPI page <https://pypi.org/project/benpy/#annex-installing-glpk>`_

**Efficient solver**
To use MIMEco, you will need to download a solver, preferably CPLEX or gurobi. Both are free for academics, but require to get a license to be used at full capacity.
A tutorial on how to make gurobi work with mimeco is written in the `documentation <https://mimeco.readthedocs.io/en/latest/installation.html#installing-gurobi>`_

Installation
~~~~~~~~~~~~~

MIMEco is available on pyPI. You can install it with :code:`pip install mimeco`