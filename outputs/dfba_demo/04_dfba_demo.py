"""
Step 4: Fix missing deps and run dFBA demo with cometspy
"""
import os
import sys
import subprocess
import traceback
import urllib.request

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dfba_demo_log.txt")

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

# Set up environment
JDK_HOME = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8"
COMETS_HOME = r"C:\Users\hgh97\comets_2.12.4"

os.environ["JAVA_HOME"] = JDK_HOME
os.environ["COMETS_HOME"] = COMETS_HOME
os.environ["PATH"] = os.path.join(JDK_HOME, "bin") + ";" + os.environ.get("PATH", "")

print("=" * 60)
print("dFBA Demo with cometspy")
print("=" * 60)

# Fix missing dependencies
print("\n--- Fixing missing dependencies ---")
lib_dir = os.path.join(COMETS_HOME, "lib")

# jdistlib - try alternative source
jdistlib_path = os.path.join(lib_dir, "jdistlib-0.4.5-bin.jar")
if not os.path.exists(jdistlib_path):
    print("Downloading jdistlib from alternative source...")
    try:
        # Try Maven Central or other source
        url = "https://repo1.maven.org/maven2/io/github/alexengird/jdistlib/0.4.5/jdistlib-0.4.5.jar"
        urllib.request.urlretrieve(url, jdistlib_path)
        print(f"  Downloaded jdistlib ({os.path.getsize(jdistlib_path)} bytes)")
    except Exception as e:
        print(f"  jdistlib download failed: {e}")
        # Try creating a minimal stub - the library might not be critical for basic dFBA
        print("  Creating empty stub jar (may cause runtime errors)")
else:
    print(f"  jdistlib already exists")

# JMatIO - try alternative source
jmatio_dir = os.path.join(lib_dir, "JMatIO", "lib")
jmatio_path = os.path.join(jmatio_dir, "jmatio.jar")
if not os.path.exists(jmatio_path):
    os.makedirs(jmatio_dir, exist_ok=True)
    print("Downloading JMatIO from alternative source...")
    try:
        url = "https://repo1.maven.org/maven2/net/sf/jmatio/jmatio/0.2/jmatio-0.2.jar"
        urllib.request.urlretrieve(url, jmatio_path)
        print(f"  Downloaded jmatio ({os.path.getsize(jmatio_path)} bytes)")
    except Exception as e:
        print(f"  JMatIO download failed: {e}")
        # Try another URL
        try:
            url = "https://repo1.maven.org/maven2/org/scijava/jmatio/1.2/jmatio-1.2.jar"
            urllib.request.urlretrieve(url, jmatio_path)
            print(f"  Downloaded jmatio (alt) ({os.path.getsize(jmatio_path)} bytes)")
        except Exception as e2:
            print(f"  JMatIO alt download also failed: {e2}")
else:
    print(f"  JMatIO already exists")

# Also create comets_scr without .bat extension (cometspy might look for it)
comets_scr_bat = os.path.join(COMETS_HOME, "comets_scr.bat")
comets_scr_noext = os.path.join(COMETS_HOME, "comets_scr")
if os.path.exists(comets_scr_bat) and not os.path.exists(comets_scr_noext):
    shutil_copy = __import__('shutil').copy2
    shutil_copy(comets_scr_bat, comets_scr_noext)
    print(f"  Created comets_scr (without .bat extension)")

# Update comets_scr.bat with proper classpath
print("\n--- Updating comets_scr.bat ---")
JAVA = os.path.join(JDK_HOME, "bin", "java.exe")

# Collect all JARs for classpath
classpath_jars = []
for root, dirs, files in os.walk(lib_dir):
    for f in files:
        if f.endswith(".jar"):
            classpath_jars.append(os.path.join(root, f))
# Add bin JARs
bin_dir = os.path.join(COMETS_HOME, "bin")
for f in os.listdir(bin_dir):
    if f.endswith(".jar"):
        classpath_jars.append(os.path.join(bin_dir, f))

cp = ";".join(classpath_jars)

bat_content = f"""@echo off
set JAVA={JAVA}
set COMETS_HOME={COMETS_HOME}
set CLASSPATH={cp}
"%JAVA%" -classpath "%CLASSPATH%" edu.bu.segrelab.comets.Comets -loader edu.bu.segrelab.comets.fba.FBACometsLoader -script %1
"""
with open(comets_scr_bat, "w") as f:
    f.write(bat_content)
print(f"  Updated {comets_scr_bat}")
print(f"  Classpath has {len(classpath_jars)} JARs")

# ---- NOW RUN THE dFBA DEMO ----
print("\n" + "=" * 60)
print("Running dFBA Demo - Two-species co-culture")
print("=" * 60)

try:
    import cometspy as c
    print(f"cometspy imported successfully")
    print(f"cometspy modules: {dir(c)}")

    # Try using cobra to load a test model
    import cobra
    print(f"cobra version: {cobra.__version__}")

    # Load E. coli textbook model
    print("\nLoading E. coli textbook model...")
    try:
        ecoli = cobra.io.load_model("textbook")
        print(f"  Model: {ecoli.id}")
        print(f"  Reactions: {len(ecoli.reactions)}")
        print(f"  Metabolites: {len(ecoli.metabolites)}")
    except Exception as e:
        print(f"  Failed to load textbook model: {e}")
        # Try alternative
        try:
            from cobra.test import create_test_model
            ecoli = create_test_model("textbook")
            print(f"  Model (alt): {ecoli.id}")
        except Exception as e2:
            print(f"  Alternative also failed: {e2}")
            ecoli = None

    if ecoli is not None:
        # Create two cometspy models from the same E. coli model
        # Model 1: wild type
        print("\nCreating cometspy models...")
        model1 = c.model(ecoli)
        model1.id = "ecoli_wt"
        model1.initial_pop = [0, 0, 1e-6]
        model1.open_exchanges()
        print(f"  Model 1: {model1.id}")

        # Model 2: mutant (same model, different name to simulate co-culture)
        model2 = c.model(ecoli)
        model2.id = "ecoli_mutant"
        model2.initial_pop = [0, 0, 1e-6]
        model2.open_exchanges()
        print(f"  Model 2: {model2.id}")

        # Create layout
        print("\nSetting up layout...")
        lay = c.layout([model1, model2])

        # Set initial media concentrations
        lay.set_specific_metabolite("glc__D_e", 0.01)  # glucose
        lay.set_specific_metabolite("nh4_e", 1000.0)    # ammonium
        lay.set_specific_metabolite("pi_e", 1000.0)     # phosphate
        print("  Media set: glucose=0.01, nh4=1000, pi=1000")

        # Set parameters
        print("\nSetting parameters...")
        p = c.params()
        p.set_param("maxCycles", 50)          # Short run for demo
        p.set_param("timeStep", 0.1)           # hours per step
        p.set_param("writeTotalBiomassLog", True)
        p.set_param("writeMediaLog", True)
        p.set_param("BiomassLogRate", 1)
        p.set_param("MediaLogRate", 1)
        p.set_param("numCyclesPerDeletion", 50)
        p.set_param("deathRate", 0.0)          # No death for simplicity
        print("  maxCycles=50, timeStep=0.1h")

        # Create simulation
        print("\nCreating COMETS simulation object...")
        try:
            sim = c.comets(lay, p)
            print(f"  Simulation object created")
            print(f"  COMETS_HOME: {sim.COMETS_HOME}")
            print(f"  VERSION: {sim.VERSION}")
        except Exception as e:
            print(f"  ERROR creating simulation: {e}")
            traceback.print_exc()
            print("\nAttempting manual COMETS_HOME setup...")

            # The issue might be that COMETS_HOME env var format doesn't match
            # Try setting it with forward slashes
            os.environ["COMETS_HOME"] = COMETS_HOME.replace("\\", "/")
            try:
                sim = c.comets(lay, p)
                print(f"  Simulation object created (with forward slashes)")
            except Exception as e2:
                print(f"  Still failed: {e2}")
                traceback.print_exc()
                sim = None

        if sim is not None:
            # Run simulation
            print("\nRunning simulation...")
            try:
                sim.run(delete_files=False, progress=True)
                print("\nSimulation completed!")

                # Print results
                if hasattr(sim, 'total_biomass') and sim.total_biomass is not None:
                    print("\n--- Total Biomass ---")
                    print(sim.total_biomass.to_string())

                    # Save to CSV
                    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "total_biomass.csv")
                    sim.total_biomass.to_csv(csv_path, index=False)
                    print(f"\nSaved to {csv_path}")

                if hasattr(sim, 'media') and sim.media is not None:
                    print("\n--- Media (first 20 rows) ---")
                    print(sim.media.head(20).to_string())

                    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "media.csv")
                    sim.media.to_csv(csv_path, index=False)
                    print(f"\nSaved to {csv_path}")

                # Print run output for debugging
                print("\n--- Simulation Output (last 2000 chars) ---")
                print(sim.run_output[-2000:] if sim.run_output else "No output")

                if sim.run_errors and sim.run_errors != "STDERR empty.":
                    print("\n--- Simulation Errors ---")
                    print(sim.run_errors[:2000])

            except Exception as e:
                print(f"\nSimulation run failed: {e}")
                traceback.print_exc()
                # Print partial output if available
                if hasattr(sim, 'run_output'):
                    print(f"\nPartial output: {sim.run_output[:3000]}")
                if hasattr(sim, 'run_errors'):
                    print(f"\nErrors: {sim.run_errors[:2000]}")
    else:
        print("Could not load any model. Cannot run demo.")

except ImportError as e:
    print(f"Import error: {e}")
    traceback.print_exc()
except Exception as e:
    print(f"Unexpected error: {e}")
    traceback.print_exc()

log_fh.close()
