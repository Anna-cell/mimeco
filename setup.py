from setuptools import setup, find_packages

setup(
name='mimeco',
version='0.1',
author='Anna Lambert',
author_email='anna.lambert@univ-nantes.fr',
description='Multi-objective GEMs metabolic interaction inference',
packages=find_packages(),
classifiers=[
'Programming Language :: Python :: 3',
#'License :: OSI Approved :: MIT License',
'Operating System :: OS Independent',
],
python_requires='>=3.7',
install_requires=
pandas >= 0.22.1 #should be updated from cobra anyway
cobra >= 1.3.5
matplotlib >= 3.5.2
scikit-learn >= 1.0.2 
matplotlib >= 3.5.2)
#mocpaby, solver
#warnings, math : builtin