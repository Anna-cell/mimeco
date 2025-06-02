Installation
============

MIMEco is available in `pyPI <https://pypi.org/project/mimeco/>`_::
    
    pip install mimeco

Dependencies 
-------------

**GLPK**:
MIMEco depends on benpy, which needs glpk to function. Its installation is clearly described in `benpy's pyPI page <https://pypi.org/project/benpy/#annex-installing-glpk>`_

**Efficient solver**:
To infer reliable results, it is important to have an efficient solver, such as gurobi or CPLEX.
While an unlicensed version of gurobi exists on pyPI and conda, it is limited to small models and won't work with typical metabolic models.
Gurobi at its full capacity is free for academics, but necesitates to create an account to retrieve a license. 

Installing gurobi
~~~~~~~~~~~~~~~~~

1. Create an account on `gurobi.com <https://www.gurobi.com>`_
    To get a free license, you must be identified on the website
2. Download a licence `here <https://portal.gurobi.com/iam/licenses/request>`_          
    Dowload the named-user or the WLS licence, depending on if you need to use mimeco on a single machine or more.
3. Download the Gurobi Optimizer `here <https://www.gurobi.com/downloads/gurobi-software/>`_
4. Choose a directory for installation. **Move your license** (gurobi.lic) in it and **unzip the Gurobi Optimizer** compressed file there with :code:`tar xvfz`. 

   For example::

    tar xvfz /home/user/Downloads/gurobi12.0.1_linux64.tar.gz /home/installation/gurobi

5. **Set your environment variables** to find gurobi. This process depends on your operating system. A complete explanation can be found `here <https://support.gurobi.com/hc/en-us/articles/13443862111761-How-do-I-set-system-environment-variables-for-Gurobi>`_.
   To follow our precedent example on ubuntu, you should edit your .bashrc to add::

    export GRB_LICENSE_FILE=/home/installation/gurobi/gurobi.lic
    export GUROBI_HOME="/home/installation/gurobi/gurobi1201/linux64"
    export PATH="${PATH}:${GUROBI_HOME}/bin"
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

6. **Test your License.** For that, type :code:`gurobi_cl --license` in a terminal, and verify that your license is indeed found. 
7. **Test your installation** by importing gurobi in your python script with :code:`import gurobipy`

If steps 6 and 7 work, you should be able to use "gurobi" as solver option when using mimeco, as long as you import gurobipy in your script.