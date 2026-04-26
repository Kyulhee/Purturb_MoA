"""
Step 8: dFBA Demo - Fixed version mismatch and encoding issues
Uses progress=False to avoid ascii decode error
Removes version-mismatch parameters
"""
import os
import sys
import subprocess
import traceback

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo_run_log.txt")

class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            try:
                f.write(obj)
                f.flush()
            except:
                pass
    def flush(self):
        for f in self.files:
            try:
                f.flush()
            except:
                pass

log_fh = open(LOG_FILE, "w", encoding="utf-8")
sys.stdout = Tee(sys.__stdout__, log_fh)
sys.stderr = Tee(sys.__stderr__, log_fh)

COMETS_HOME = r"C:\Users\hgh97\comets_2.12.4"
JDK_HOME = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8"

os.environ["JAVA_HOME"] = JDK_HOME
os.environ["COMETS_HOME"] = COMETS_HOME
os.environ["PATH"] = os.path.join(JDK_HOME, "bin") + ";" + os.environ.get("PATH", "")

print("=" * 60)
print("dFBA Demo - Two-species E. coli co-culture (Fixed)")
print("=" * 60)

import cometspy as c
import cobra

print(f"cometspy: 0.6.3")
print(f"cobra: {cobra.__version__}")

# Load E. coli core model
ecoli = cobra.io.load_model("textbook")
print(f"Model: {ecoli.id}, {len(ecoli.reactions)} reactions, {len(ecoli.metabolites)} metabolites")

# Create two species for co-culture
model1 = c.model(ecoli)
model1.id = "ecoli_wt"
model1.initial_pop = [0, 0, 1e-6]
model1.open_exchanges()

model2 = c.model(ecoli)
model2.id = "ecoli_mutant"
model2.initial_pop = [0, 0, 1e-6]
model2.open_exchanges()

# Layout
lay = c.layout([model1, model2])
lay.set_specific_metabolite("glc__D_e", 0.01)
lay.set_specific_metabolite("nh4_e", 1000.0)
lay.set_specific_metabolite("pi_e", 1000.0)

# Parameters
p = c.params()
p.set_param("maxCycles", 30)       # Shorter run
p.set_param("timeStep", 0.1)
p.set_param("writeTotalBiomassLog", True)
p.set_param("writeMediaLog", True)
p.set_param("BiomassLogRate", 1)
p.set_param("MediaLogRate", 1)
p.set_param("deathRate", 0.0)

# Remove problematic parameters that don't exist in COMETS 2.10.0
# The version mismatch causes "Unknown parameter" errors
print("\nRemoving version-mismatch parameters...")
skip_keys = ['velocitymulticonvlogformat', 'velocityMultiConvLogName',
             'velocityMultiConvLogRate', 'writeVelocityMultiConvLog',
             'numCyclesPerDeletion']
for key in list(p.all_params.keys()):
    if any(sk.lower() in key.lower() for sk in skip_keys):
        print(f"  Removing: {key} = {p.all_params[key]}")
        del p.all_params[key]

# Create simulation object
print("\nCreating COMETS simulation...")
sim = c.comets(lay, p)
print(f"COMETS VERSION: {sim.VERSION}")

# Also remove from sim parameters after object creation
for key in list(sim.parameters.all_params.keys()):
    if any(sk.lower() in key.lower() for sk in skip_keys):
        del sim.parameters.all_params[key]

# Run without progress bar (avoids ascii decode error)
print("\nRunning simulation (maxCycles=30, timeStep=0.1h)...")
print("Using progress=False to avoid encoding issues")

try:
    sim.run(delete_files=False, progress=False)
    print("\n=== SIMULATION COMPLETED ===")

    # Results
    if hasattr(sim, 'total_biomass') and sim.total_biomass is not None:
        print("\n--- Total Biomass ---")
        print(sim.total_biomass.to_string())
        csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "total_biomass.csv")
        sim.total_biomass.to_csv(csv_path, index=False)
        print(f"\nSaved to {csv_path}")
    else:
        print("\nNo total_biomass data available")

    if hasattr(sim, 'media') and sim.media is not None:
        print("\n--- Media (first 30 rows) ---")
        print(sim.media.head(30).to_string())
        csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "media.csv")
        sim.media.to_csv(csv_path, index=False)
        print(f"\nSaved to {csv_path}")
    else:
        print("\nNo media data available")

    # Print run output
    if sim.run_output:
        print("\n--- Simulation Output (last 3000 chars) ---")
        # Decode with error handling
        try:
            print(sim.run_output[-3000:])
        except:
            print(sim.run_output[-3000:].encode('utf-8', errors='replace').decode('utf-8'))

    if sim.run_errors and sim.run_errors != "STDERR empty.":
        print("\n--- Simulation Errors ---")
        print(sim.run_errors[:3000])

except UnicodeDecodeError as e:
    print(f"\nUnicodeDecodeError during run: {e}")
    print("The simulation may have actually run. Checking for output files...")

    # Try to read the output files directly
    working_dir = sim.working_dir
    print(f"Working dir: {working_dir}")

    # Look for output log files
    import glob
    for pattern in ["totalbiomass*", "biomasslog*", "medialog*", "fluxlog*"]:
        matches = glob.glob(os.path.join(working_dir, pattern))
        for m in matches:
            print(f"\nFound: {m}")
            try:
                with open(m, 'r', encoding='utf-8', errors='replace') as f:
                    content = f.read()
                    print(content[:1000])
            except:
                print("  Could not read file")

except Exception as e:
    print(f"\nERROR: {type(e).__name__}: {e}")
    traceback.print_exc()

    # Try to read partial output
    try:
        if hasattr(sim, 'run_output') and sim.run_output:
            print(f"\nPartial output (first 5000 chars):\n{sim.run_output[:5000]}")
    except:
        pass

# Also try running COMETS directly via command line for comparison
print("\n" + "=" * 60)
print("Alternative: Direct Java execution test")
print("=" * 60)

try:
    # Build classpath
    import glob as globmod
    classpath_jars = []
    for root, dirs, files in os.walk(os.path.join(COMETS_HOME, "lib")):
        for f in files:
            if f.endswith(".jar"):
                classpath_jars.append(os.path.join(root, f))
    for f in os.listdir(os.path.join(COMETS_HOME, "bin")):
        if f.endswith(".jar"):
            classpath_jars.append(os.path.join(COMETS_HOME, "bin", f))

    cp = ";".join(classpath_jars)
    java_exe = os.path.join(JDK_HOME, "bin", "java.exe")

    result = subprocess.run(
        [java_exe, "-classpath", cp,
         "edu.bu.segrelab.comets.Comets",
         "-loader", "edu.bu.segrelab.comets.fba.FBACometsLoader"],
        capture_output=True, text=True, timeout=15,
        cwd=COMETS_HOME
    )
    print(f"Direct Java execution return code: {result.returncode}")
    print(f"stdout: {result.stdout[:2000]}")
    if result.stderr:
        print(f"stderr: {result.stderr[:1000]}")
except Exception as e:
    print(f"Direct execution failed: {e}")

log_fh.close()
