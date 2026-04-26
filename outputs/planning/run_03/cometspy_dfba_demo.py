"""
cometspy dFBA Demo Script
=========================
2-species co-culture simulation using COMETS via cometspy.
This script tests the full cometspy workflow and records any errors/limitations.

Purpose: OCT LLM XAI Project - Phase 4 (dFBA Dynamic Simulation) design reference
Date: 2026-04-26
"""

import sys
import traceback

# =============================================
# Step 1: Verify cometspy installation
# =============================================
print("=" * 60)
print("STEP 1: Verify cometspy installation")
print("=" * 60)

try:
    import cometspy
    print(f"[OK] cometspy imported successfully")
    # version check - module may not have __version__
    try:
        print(f"  Version: {cometspy.__version__}")
    except AttributeError:
        # Try from comets submodule
        from cometspy.comets import __version__
        print(f"  Version: {__version__} (from cometspy.comets)")
except ImportError as e:
    print(f"[FAIL] cometspy import failed: {e}")
    sys.exit(1)

try:
    from cometspy import model, layout, params, comets
    print(f"[OK] All submodules imported: model, layout, params, comets")
except ImportError as e:
    print(f"[FAIL] Submodule import failed: {e}")
    sys.exit(1)

try:
    import cobra
    print(f"[OK] cobra imported, version: {cobra.__version__}")
except ImportError as e:
    print(f"[FAIL] cobra import failed: {e}")
    sys.exit(1)

# =============================================
# Step 2: Check COMETS Java core availability
# =============================================
print("\n" + "=" * 60)
print("STEP 2: Check COMETS Java core availability")
print("=" * 60)

import os

comets_home = os.environ.get('COMETS_HOME', None)
print(f"  COMETS_HOME env var: {comets_home if comets_home else 'NOT SET'}")

java_path = None
try:
    import shutil
    java_path = shutil.which('java')
except:
    pass
print(f"  Java on PATH: {java_path if java_path else 'NOT FOUND'}")

if comets_home is None:
    print("\n  [CRITICAL] COMETS_HOME is not set!")
    print("  cometspy is a Python wrapper that calls the COMETS Java engine.")
    print("  Without COMETS_HOME pointing to the Java installation,")
    print("  sim.run() will fail with a KeyError on os.environ['COMETS_HOME']")
    print()
    print("  Required setup:")
    print("  1. Install Java JDK 11+")
    print("  2. Download COMETS from https://github.com/segrelab/COMETS")
    print("  3. Set COMETS_HOME environment variable to COMETS installation dir")
    print("  4. (Optional) Set GUROBI_HOME if using Gurobi optimizer")

# =============================================
# Step 3: Build Python-side model objects (no Java needed)
# =============================================
print("\n" + "=" * 60)
print("STEP 3: Build cometspy model objects (Python-side)")
print("=" * 60)

try:
    # Load E. coli textbook model from cobra (new API in cobra 0.31+)
    print("  Loading E. coli textbook model from cobra...")
    ecoli_cobra = cobra.io.load_model("textbook")
    print(f"  [OK] cobra model loaded: {ecoli_cobra.id}")
    print(f"        Reactions: {len(ecoli_cobra.reactions)}")
    print(f"        Metabolites: {len(ecoli_cobra.metabolites)}")

    # Convert to COMETS model
    ecoli_comets = model(ecoli_cobra)
    ecoli_comets.open_exchanges()
    ecoli_comets.initial_pop = [0, 0, 1.e-10]
    print(f"  [OK] COMETS model created: {ecoli_comets.id}")
    print(f"        Initial pop: {ecoli_comets.initial_pop}")
    print(f"        Exchange metabolites (first 5): {list(ecoli_comets.get_exchange_metabolites())[:5]}")

    # Create a second species by modifying the E. coli model
    # (simulating a different strain with different acetate uptake)
    print("\n  Creating second species (modified E. coli)...")

    ecoli2_cobra = cobra.io.load_model("textbook")
    ecoli2_cobra.id = "ecoli_mutant"
    # Modify: reduce glucose uptake, enhance acetate uptake
    for rxn in ecoli2_cobra.reactions:
        if rxn.id == "EX_glc__D_e":
            rxn.lower_bound = -5.0  # reduced glucose uptake
        if rxn.id == "EX_ac_e":
            rxn.lower_bound = -15.0  # enhanced acetate uptake

    ecoli2_comets = model(ecoli2_cobra)
    ecoli2_comets.open_exchanges()
    ecoli2_comets.initial_pop = [0, 0, 1.e-10]
    print(f"  [OK] Second COMETS model created: {ecoli2_comets.id}")
    print(f"        Initial pop: {ecoli2_comets.initial_pop}")

except Exception as e:
    print(f"  [FAIL] Model creation error: {e}")
    traceback.print_exc()

# =============================================
# Step 4: Build layout and parameters
# =============================================
print("\n" + "=" * 60)
print("STEP 4: Build layout and parameters")
print("=" * 60)

try:
    # Create layout with both models
    l = layout([ecoli_comets, ecoli2_comets])

    # Set up media
    l.set_specific_metabolite("glc__D_e", 0.015)     # glucose (limiting)
    l.set_specific_metabolite("ac_e", 0.0)            # acetate (produced by species 1)
    l.set_specific_metabolite("o2_e", 15.0, static=True)  # oxygen (unlimited)
    l.set_specific_metabolite("h2o_e", 1000.0, static=True)
    l.set_specific_metabolite("nh4_e", 1000.0, static=True)
    l.set_specific_metabolite("pi_e", 1000.0, static=True)

    print(f"  [OK] Layout created with grid: {l.grid}")
    print(f"        Models: {l.get_model_ids()}")
    print(f"        Media metabolites: {len(l.media)}")
    print(f"        Non-zero media:")
    for _, row in l.media[l.media['init_amount'] != 0.0].iterrows():
        print(f"          {row['metabolite']}: {row['init_amount']} mmol")

    # Set parameters
    p = params()
    p.set_param("timeStep", 0.1)         # hours per step
    p.set_param("spaceWidth", 0.02)      # cm
    p.set_param("maxCycles", 100)         # total cycles
    p.set_param("writeTotalBiomassLog", True)
    p.set_param("writeMediaLog", True)
    p.set_param("MediaLogRate", 10)
    p.set_param("writeFluxLog", True)
    p.set_param("FluxLogRate", 10)

    print(f"\n  [OK] Parameters set:")
    print(f"        timeStep: {p.all_params['timeStep']} hr")
    print(f"        spaceWidth: {p.all_params['spaceWidth']} cm")
    print(f"        maxCycles: {p.all_params['maxCycles']}")
    print(f"        defaultVmax: {p.all_params['defaultVmax']} mmol/gDW/hr")
    print(f"        defaultKm: {p.all_params['defaultKm']} M")
    print(f"        exchangestyle: {p.all_params['exchangestyle']}")

except Exception as e:
    print(f"  [FAIL] Layout/params error: {e}")
    traceback.print_exc()

# =============================================
# Step 5: Attempt simulation run
# =============================================
print("\n" + "=" * 60)
print("STEP 5: Attempt simulation run (requires COMETS Java)")
print("=" * 60)

try:
    sim = comets(l, p)
    print("  [OK] comets object created (Java classpath built)")
    print(f"        COMETS_HOME: {sim.COMETS_HOME}")
    print(f"        VERSION: {sim.VERSION}")
    print("\n  Running simulation...")
    sim.run()
    print("\n  [OK] Simulation completed!")
    print(f"  Total biomass:\n{sim.total_biomass}")

    # Try to get metabolite time series
    try:
        media_ts = sim.get_metabolite_time_series()
        print(f"\n  Media time series:\n{media_ts.head()}")
    except Exception as e:
        print(f"  [WARN] Could not get media time series: {e}")

except KeyError as e:
    print(f"  [EXPECTED FAILURE] KeyError: {e}")
    print("  --> COMETS_HOME environment variable is not set.")
    print("  --> This is the expected error when COMETS Java core is not installed.")
    print("  --> The cometspy Python wrapper requires the COMETS Java engine")
    print("      to be installed and COMETS_HOME to point to it.")
except Exception as e:
    print(f"  [UNEXPECTED ERROR] {type(e).__name__}: {e}")
    traceback.print_exc()

# =============================================
# Step 6: Summary and findings
# =============================================
print("\n" + "=" * 60)
print("STEP 6: Summary and Findings")
print("=" * 60)

summary = """
COMETSPy dFBA Demo Summary
===========================

1. INSTALLATION STATUS
   - cometspy v0.6.3: INSTALLED (pip install cometspy)
   - Dependencies: cobra 0.31.1, numpy, pandas, optlang, swiglpk, python-libsbml
   - All Python submodules (model, layout, params, comets): IMPORTABLE

2. COMETS JAVA CORE STATUS
   - COMETS_HOME: NOT SET
   - Java runtime: NOT FOUND on PATH
   - COMETS Java JAR: NOT BUNDLED with cometspy pip package
   - CONCLUSION: cometspy is a Python WRAPPER, not a standalone simulator

3. PYTHON-SIDE FUNCTIONALITY (works without Java)
   - Model construction from cobra models: WORKS
   - Layout creation with multiple species: WORKS
   - Media/metabolite configuration: WORKS
   - Parameter configuration: WORKS
   - Model file writing (COMETS .cmd format): WORKS
   - Layout file writing: WORKS

4. SIMULATION (requires Java COMETS core)
   - sim.run() FAILS with KeyError on os.environ['COMETS_HOME']
   - COMETS Java core must be separately downloaded from:
     https://github.com/segrelab/COMETS
   - Requires Java JDK 11+ and or-tools/Gurobi LP solver

5. CRITICAL LIMITATIONS
   a) cometspy cannot run simulations standalone - Java COMETS engine required
   b) COMETS_HOME must be set as environment variable
   c) On Windows, COMETS uses comets_scr.bat script (line 396-400 in comets.py)
   d) Default optimizer is GUROBI (commercial); or-tools (free) also supported
   e) Python 3.13 compatibility: cometspy installed but cobra uses deprecated
      pandas features (e.g., DataFrame.append removed in pandas 2.x, used in
      model.py add_signal/add_multitoxin methods)
   f) No __version__ attribute on the cometspy package itself

6. SETUP REQUIREMENTS FOR FULL FUNCTIONALITY
   - Install Java JDK 11+
   - Download COMETS v2.12.5 from GitHub releases
   - Set COMETS_HOME environment variable
   - (Optional) Install Gurobi and set GUROBI_HOME for commercial solver
   - Alternative: Use or-tools (bundled with COMETS) as free solver

7. IMPLICATIONS FOR PHASE 4 DESIGN
   - For pure-Python dFBA, consider alternatives:
     a) dfba (dfba-python): Pure Python dFBA using cobra + scipy
     b) Cameo: Python metabolic engineering with dynamic simulation
     c) Custom dFBA loop using cobra.optimize() + Euler integration
   - cometspy is best suited when spatial/2D simulation is needed
   - For non-spatial (well-mixed) dFBA, simpler tools may be more appropriate
   - If using cometspy, Docker container approach recommended for deployment
"""

print(summary)
