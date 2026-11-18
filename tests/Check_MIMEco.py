#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_mimeco.py
===============

Diagnostic of mimeco.analysis before a push: every function of analysis.py is run on the usual
models, and the results are compared with those of a reference version of mimeco, e.g. 0.1.20 (PyPI).
The reference file also stores where it comes from (mimeco version, solver, sample size, date), which
is printed before the comparison.

Usage
-----
1. Once, in an environment where the reference mimeco is installed (e.g. pip install mimeco==0.1.20):
       python check_mimeco.py reference      -> writes reference_results.pickle next to this script
2. After each change in mimeco_source:
       python check_mimeco.py local          -> runs mimeco_source and prints every error, then every
                                                difference with the reference
   ('python check_mimeco.py' alone also runs the local mode)

Cases
-----
  carveme     L. plantarum x A. oryzae, CarveMe models, YPD (BiGG namespace)
  gapseq      A. oryzae x L. plantarum, gapseq models, YPD (ModelSEED namespace)
  enterocyte  CarveMe L. plantarum with the enterocyte, Western diet (BiGG namespace)
  agora2      L. plantarum x L. casei, gap-filled AGORA2 JSON models, MRS+ (local mode only: mimeco
              0.1.20 fails on these models, which is the bug being fixed)
"""
import os
import sys
import pickle
import datetime
from importlib.metadata import version
import traceback
import warnings

import matplotlib
matplotlib.use("Agg")                 # plots are drawn (to test the plot functions) but not shown
import matplotlib.pyplot as plt
import pandas as pd
import cobra
import logging
# cobra warns about every exchange reaction it adds or skips when reading a CarveMe model
# ('Adding exchange reaction ...', 'Ignoring reaction ... since it already exists'): harmless, hidden here.
logging.getLogger("cobra").setLevel(logging.ERROR)


# =============================================================================
# SETTINGS
# =============================================================================
# Everything is found from the location of this script (tests/), so the test runs on any machine:
# models and media in tests/resources/, the mimeco code to test in the folder above tests/.
TESTS = os.path.dirname(os.path.abspath(__file__)) + "/"           # .../mimeco_source/mimeco/tests/
RESOURCES = TESTS + "resources/"
MIMECO_SOURCE = os.path.dirname(os.path.abspath(TESTS)) + "/"      # .../mimeco_source/mimeco/
REFERENCE_FILE = TESTS + "reference_results.pickle"

SOLVER = "gurobi"
# Full academic license; the gurobi.lic files inside site-packages/gurobipy are size-limited trial
# licenses, too small for ecosystem models. A GRB_LICENSE_FILE already set in the terminal is kept.
os.environ.setdefault("GRB_LICENSE_FILE", os.path.expanduser("~/Installations/gurobi1303/gurobi.lic"))
SAMPLE_SIZE = 200                     # Pareto samples used by the test in crossfed_metabolites. MIMEco's
                                      # default is 1000; 200 is faster, and both versions use the same value
TOLERANCE = 1e-6                      # relative tolerance when comparing numbers

# "reference" (once, with the reference mimeco) or "local" (mimeco_source), given on the command line, e.g.
# 'python check_mimeco.py reference'. Without argument (e.g. run from Spyder): "local".
MODE = sys.argv[1] if len(sys.argv) > 1 else "local"


# =============================================================================
# IMPORT THE RIGHT MIMECO
# =============================================================================
if MODE == "local":
    sys.path.insert(0, MIMECO_SOURCE)             # mimeco_source is found before the conda mimeco
try:
    import mimeco
    from mimeco import analysis
except Exception:                                 # e.g. a SyntaxError in analysis.py
    print("!! mimeco cannot be imported:")
    print(traceback.format_exc(limit=0))
    sys.exit(1)

print(f"Mode: {MODE}, mimeco imported from {os.path.dirname(mimeco.__file__)}")
from_source = os.path.abspath(mimeco.__file__).startswith(os.path.abspath(MIMECO_SOURCE))
if MODE == "local" and not from_source:
    sys.exit("!! local mode, but mimeco does not come from mimeco_source")
if MODE == "reference" and from_source:
    print("   warning: the reference is made with mimeco_source, not with an installed version")


# =============================================================================
# FUNCTIONS
# =============================================================================
def load(path, model_id):
    """Load an SBML or JSON model with a short id, e.g. 'L_pla', used in the ecosystem ids."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if path.endswith(".json"):
            model = cobra.io.load_json_model(path)
        else:
            model = cobra.io.read_sbml_model(path)
    model.id = model_id
    model.solver = SOLVER
    return model


def copy_model(model):
    """Copy of a model for one call.
    With Gurobi, cobra copies the solver through an LP file, which Gurobi cannot read back when a
    reaction id is an LP keyword, e.g. reaction 'St' in A_ory_gurobi.xml ('st' = 'subject to').
    The copy is therefore made with glpk, then the solver is set back to SOLVER."""
    model.solver = "glpk"
    new = model.copy()
    model.solver = SOLVER
    new.solver = SOLVER
    return new


def run(name, function, **kwargs):
    """Run one mimeco function, print OK or the error (last lines of the traceback), never stop.
    The models are copied so that one call cannot change the models of the next one."""
    for k in kwargs:
        if isinstance(kwargs[k], cobra.Model):
            kwargs[k] = copy_model(kwargs[k])
    print(f"\n--- {name}")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = function(**kwargs)
        plt.close("all")
        print(f"OK {name}")
        return result
    except Exception:
        plt.close("all")
        print(f"!! ERROR {name}")
        print("".join(traceback.format_exc().splitlines(keepends=True)[-6:]))
        return None


def to_compare(result):
    """Keep only what can be pickled and compared, e.g. an ecosystem model -> its size and growth."""
    if isinstance(result, cobra.Model):
        return {"n_reactions": len(result.reactions), "n_metabolites": len(result.metabolites),
                "objective_value": result.slim_optimize(error_value=0.0)}
    if isinstance(result, tuple):
        return tuple(to_compare(r) for r in result)
    return result


def compare(name, ref, new):
    """Print each difference between the reference and the local result. Returns the number found."""
    if isinstance(ref, tuple):
        return sum(compare(f"{name}[{i}]", r, n) for i, (r, n) in enumerate(zip(ref, new)))
    if isinstance(ref, dict):
        return sum(compare(f"{name}['{k}']", ref[k], new.get(k)) for k in ref)
    if isinstance(ref, pd.DataFrame):
        n_diff = 0
        added = [c for c in new.columns if c not in ref.columns]
        if added:
            print(f"   info {name}: new columns {added}")
        for c in ref.columns:
            if c not in new.columns:
                print(f"!! {name}: column '{c}' missing")
                n_diff += 1
            elif list(new[c]) != list(ref[c]):
                same = (len(new[c]) == len(ref[c]) and pd.api.types.is_numeric_dtype(ref[c])
                        and (abs(new[c].values - ref[c].values) <= TOLERANCE * (1 + abs(ref[c].values))).all())
                if not same:
                    print(f"!! {name}: column '{c}' differs\n   reference: {list(ref[c])[:10]}"
                          f"\n   local:     {list(new[c])[:10]}")
                    n_diff += 1
        return n_diff
    if isinstance(ref, float):
        if new is None or abs(new - ref) > TOLERANCE * (1 + abs(ref)):
            print(f"!! {name}: reference {ref}, local {new}")
            return 1
        return 0
    if ref != new:
        print(f"!! {name}: reference {ref}, local {new}")
        return 1
    return 0


# =============================================================================
# MODELS AND MEDIA
# =============================================================================
print("\nLoading models ...")
# Copies of the usual validation models (mimeco_validation/models, media) and of the Lpla_Lcas models
lpla_carveme = load(RESOURCES + "L_pla_gurobi.xml", "L_pla")
aory_carveme = load(RESOURCES + "A_ory_gurobi.xml", "A_ory")
lpla_gapseq = load(RESOURCES + "L_pla_gapseq.xml", "L_pla")
aory_gapseq = load(RESOURCES + "A_ory_gapseq.xml", "A_ory")
ypd_bigg = pd.read_csv(RESOURCES + "YPD_media_BiGG.csv", index_col=0)
ypd_modelseed = pd.read_csv(RESOURCES + "YPD_media_ModelSEED.csv", index_col=0)
western_diet = pd.read_csv(RESOURCES + "Western_diet_BiGG.csv", index_col=0)

# Same arguments as in the validation scripts (Condition1.py, Condition_LGG_e.py).
# Biomass ids: 'Growth' in CarveMe, the biomass exchange 'EX_cpd11416_c0' in gapseq.
pairs = {
    "carveme": dict(model1=lpla_carveme, model2=aory_carveme, medium=ypd_bigg,
                    biomass1="Growth", biomass2="Growth"),
    "gapseq": dict(model1=aory_gapseq, model2=lpla_gapseq, medium=ypd_modelseed,
                   biomass1="EX_cpd11416_c0", biomass2="EX_cpd11416_c0"),
}
if MODE == "local":
    pairs["agora2"] = dict(
        model1=load(RESOURCES + "Lactobacillus_plantarum_WCFS1_AGORA2_Lambert_gapfilled.json", "Lpla"),
        model2=load(RESOURCES + "Lactobacillus_casei_casei_BL23_AGORA2_Lambert_gapfilled.json", "Lcas"),
        medium=pd.read_csv(RESOURCES + "MRS_mimeco_medium.csv", index_col=0),
        biomass1="biomass524", biomass2="biomass525")
constraint = dict(undescribed_metabolites_constraint="partially_constrained", undescribed_met_lb=-0.0001)


# =============================================================================
# RUN EVERY FUNCTION OF analysis.py
# =============================================================================
results = {}
plot = MODE == "local"                # the plot functions are only tested on the local code

for case, p in pairs.items():
    results[case + " interaction_score_and_type"] = run(
        f"{case}: interaction_score_and_type", analysis.interaction_score_and_type,
        model1=p["model1"], model2=p["model2"], medium=p["medium"], plot=plot, **constraint)

    crossfed = run(
        f"{case}: crossfed_metabolites", analysis.crossfed_metabolites,
        model1=p["model1"], model2=p["model2"], medium=p["medium"], solver=SOLVER,
        model1_biomass_id=p["biomass1"], model2_biomass_id=p["biomass2"],
        sample_size=SAMPLE_SIZE, plot=plot, retrieve_data="selection", **constraint)
    results[case + " crossfed_metabolites"] = crossfed
    # each crossfed exchange must point to a real exchange reaction of model1
    if MODE == "local" and crossfed is not None:
        exchanges = [r.id for r in p["model1"].exchanges]
        for rid in crossfed[0].get("reaction_id", []):
            if rid not in exchanges:
                print(f"!! {case}: reaction_id '{rid}' is not an exchange of model1")

    if hasattr(analysis, "extract_pairwize_model"):          # not in every version of mimeco
        results[case + " extract_pairwize_model"] = run(
            f"{case}: extract_pairwize_model", analysis.extract_pairwize_model,
            model1=p["model1"], model2=p["model2"], medium=p["medium"], solver=SOLVER,
            model1_biomass_id=p["biomass1"], model2_biomass_id=p["biomass2"], **constraint)

enterocyte = dict(model=lpla_carveme, medium=western_diet, solver=SOLVER,
                  undescribed_metabolites_constraint="as_is")
results["enterocyte enterocyte_interaction_score_and_type"] = run(
    "enterocyte: enterocyte_interaction_score_and_type",
    analysis.enterocyte_interaction_score_and_type, plot=plot, **enterocyte)
results["enterocyte enterocyte_crossfed_metabolites"] = run(
    "enterocyte: enterocyte_crossfed_metabolites", analysis.enterocyte_crossfed_metabolites,
    model_biomass_id="Growth", sample_size=SAMPLE_SIZE, plot=plot, retrieve_data="selection",
    **enterocyte)

results = {k: to_compare(v) for k, v in results.items()}


# =============================================================================
# SAVE (reference) OR COMPARE (local)
# =============================================================================
n_errors = sum(v is None for v in results.values())
print(f"\n{'=' * 70}\n{n_errors} function call(s) failed out of {len(results)}")

if MODE == "reference":
    # where the reference comes from, printed before each comparison
    try:
        mimeco_version = version("mimeco")             # installed package, e.g. '0.1.20'
    except Exception:                                  # mimeco_source is not an installed package
        mimeco_version = "not installed (source folder)"
    metadata = {"mimeco_version": mimeco_version, "mimeco_path": os.path.dirname(mimeco.__file__),
                "solver": SOLVER, "sample_size": SAMPLE_SIZE,
                "date": datetime.date.today().isoformat(), "cobra_version": cobra.__version__}
    with open(REFERENCE_FILE, "wb") as fh:
        pickle.dump({"metadata": metadata, "results": results}, fh)
    print(f"Reference written in {REFERENCE_FILE}: {metadata}")
else:
    with open(REFERENCE_FILE, "rb") as fh:
        saved = pickle.load(fh)
    metadata, reference = saved["metadata"], saved["results"]
    print("Reference: " + ", ".join(f"{k} = {v}" for k, v in metadata.items()))
    # a different solver or sample size can change the sampling, hence give false differences
    if metadata["solver"] != SOLVER or metadata["sample_size"] != SAMPLE_SIZE:
        print(f"   warning: reference made with solver {metadata['solver']} and sample_size "
              f"{metadata['sample_size']}, this run uses {SOLVER} and {SAMPLE_SIZE}")
    n_diff = 0
    for name, ref in reference.items():
        if ref is None:
            print(f"   info {name}: failed with the reference mimeco, nothing to compare")
        elif results.get(name) is None:
            print(f"!! {name}: works with the reference mimeco, fails with mimeco_source")
            n_diff += 1
        else:
            n_diff += compare(name, ref, results[name])
    print(f"{n_diff} difference(s) with the reference")
    if n_errors == 0 and n_diff == 0:
        print("-> everything works and gives the same results: OK to push")
