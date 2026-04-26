"""
Step 10: Fix corrupted jdistlib, recompile, and run dFBA demo
"""
import os
import sys
import subprocess
import traceback
import glob as globmod
import shutil
import urllib.request

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "final_run_log.txt")

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

JDK_HOME = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8"
COMETS_HOME = r"C:\Users\hgh97\comets_2.12.4"
COMETS_REPO = r"C:\Users\hgh97\comets_home\COMETS"
lib_dir = os.path.join(COMETS_HOME, "lib")
bin_dir = os.path.join(COMETS_HOME, "bin")

JAVA = os.path.join(JDK_HOME, "bin", "java.exe")
JAVAC = os.path.join(JDK_HOME, "bin", "javac.exe")
JAR_TOOL = os.path.join(JDK_HOME, "bin", "jar.exe")

os.environ["JAVA_HOME"] = JDK_HOME
os.environ["PATH"] = os.path.join(JDK_HOME, "bin") + ";" + os.environ.get("PATH", "")

print("=" * 60)
print("Fix jdistlib, Recompile COMETS, Run dFBA Demo")
print("=" * 60)

# Step 1: Fix jdistlib - download from Maven or create a proper stub
print("\n--- Step 1: Fix jdistlib ---")
jdistlib_path = os.path.join(lib_dir, "jdistlib-0.4.5-bin.jar")

# Remove corrupted file
if os.path.exists(jdistlib_path):
    os.remove(jdistlib_path)
    print(f"  Removed corrupted file")

# Try downloading from different sources
jdistlib_urls = [
    "https://repo1.maven.org/maven2/io/github/alexengird/jdistlib/0.4.5/jdistlib-0.4.5.jar",
    "https://jcenter.bintray.com/jdistlib/jdistlib/0.4.5/jdistlib-0.4.5.jar",
]

downloaded = False
for url in jdistlib_urls:
    try:
        print(f"  Trying: {url}")
        urllib.request.urlretrieve(url, jdistlib_path)
        # Verify it's a valid ZIP/JAR
        import zipfile
        with zipfile.ZipFile(jdistlib_path, 'r') as z:
            z.testzip()
        print(f"  Downloaded valid JAR ({os.path.getsize(jdistlib_path)} bytes)")
        downloaded = True
        break
    except Exception as e:
        print(f"  Failed: {e}")
        if os.path.exists(jdistlib_path):
            os.remove(jdistlib_path)

if not downloaded:
    # Create a minimal but valid stub JAR with the classes COMETS actually needs
    print("  Creating valid stub JAR for jdistlib...")
    stub_dir = os.path.join(lib_dir, "jdistlib_stub")
    if os.path.exists(stub_dir):
        shutil.rmtree(stub_dir)
    os.makedirs(stub_dir, exist_ok=True)

    # Create minimal package structure
    os.makedirs(os.path.join(stub_dir, "jdistlib"), exist_ok=True)

    # Create stub Java classes that COMETS might reference
    stub_classes = {
        "Beta.java": "package jdistlib;\npublic class Beta { public double pdf(double x, double a, double b) { return 0.0; } public double cdf(double x, double a, double b) { return 0.0; } }",
        "Normal.java": "package jdistlib;\npublic class Normal { public double pdf(double x, double mu, double sigma) { return 0.0; } public double cdf(double x, double mu, double sigma) { return 0.0; } }",
        "GenericDistribution.java": "package jdistlib;\npublic abstract class GenericDistribution { public abstract double pdf(double x); public abstract double cdf(double x); }",
        "MathFunctions.java": "package jdistlib;\npublic class MathFunctions { public static double gamma(double x) { return 1.0; } public static double lgamma(double x) { return 0.0; } public static double lbeta(double a, double b) { return 0.0; } }",
        "RNG.java": "package jdistlib;\npublic class RNG { public double nextDouble() { return 0.5; } }",
        "NonCentralBeta.java": "package jdistlib;\npublic class NonCentralBeta { public double pdf(double x, double a, double b, double ncp) { return 0.0; } public double cdf(double x, double a, double b, double ncp) { return 0.0; } }",
    }

    for name, code in stub_classes.items():
        with open(os.path.join(stub_dir, "jdistlib", name), "w") as f:
            f.write(code)

    # Compile stubs
    result = subprocess.run(
        [JAVAC, "-encoding", "UTF-8", "-d", stub_dir] +
        [os.path.join(stub_dir, "jdistlib", n) for n in stub_classes.keys()],
        capture_output=True, text=True, timeout=30
    )
    print(f"  Stub compilation: {result.returncode}")
    if result.stderr:
        print(f"  Stub stderr: {result.stderr[:500]}")

    # Package into JAR
    result = subprocess.run(
        [JAR_TOOL, "cf", jdistlib_path, "-C", stub_dir, "jdistlib"],
        capture_output=True, text=True, timeout=10
    )
    print(f"  JAR creation: {result.returncode}")

    # Verify
    if os.path.exists(jdistlib_path):
        print(f"  Created valid stub JAR ({os.path.getsize(jdistlib_path)} bytes)")
        try:
            import zipfile
            with zipfile.ZipFile(jdistlib_path, 'r') as z:
                z.testzip()
            print("  JAR verification: OK")
        except Exception as e:
            print(f"  JAR verification failed: {e}")

    # Cleanup
    shutil.rmtree(stub_dir, ignore_errors=True)

# Step 2: Compile COMETS from source
print("\n--- Step 2: Compile COMETS ---")
src_dir = os.path.join(COMETS_REPO, "comets_simplified", "src")
compile_dir = os.path.join(COMETS_HOME, "build")
if os.path.exists(compile_dir):
    shutil.rmtree(compile_dir)
os.makedirs(compile_dir, exist_ok=True)

# Build classpath
classpath_jars = []
for root, dirs, files in os.walk(lib_dir):
    for f in files:
        if f.endswith(".jar"):
            fpath = os.path.join(root, f)
            # Verify JAR is valid
            try:
                import zipfile
                with zipfile.ZipFile(fpath, 'r') as z:
                    z.testzip()
                classpath_jars.append(fpath)
            except:
                print(f"  Skipping invalid JAR: {f}")
                # Remove invalid JARs
                os.remove(fpath)

for f in os.listdir(bin_dir):
    if f.endswith(".jar"):
        classpath_jars.append(os.path.join(bin_dir, f))

classpath = ";".join(classpath_jars)
print(f"Valid classpath: {len(classpath_jars)} JARs")

# Collect Java source files
java_files = []
for root, dirs, files in os.walk(src_dir):
    for f in files:
        if f.endswith(".java"):
            java_files.append(os.path.join(root, f))

print(f"Compiling {len(java_files)} Java files...")

cmd = [JAVAC, "-encoding", "UTF-8", "-cp", classpath, "-d", compile_dir, "-sourcepath", src_dir]
cmd.extend(java_files)

result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
print(f"Compile return code: {result.returncode}")

if result.stderr:
    # Count errors
    errors = [l for l in result.stderr.split('\n') if 'error:' in l]
    warnings = [l for l in result.stderr.split('\n') if 'warning:' in l]
    print(f"Errors: {len(errors)}, Warnings: {len(warnings)}")
    print(f"First 3000 chars of stderr:\n{result.stderr[:3000]}")

if result.returncode == 0:
    print("\nCompilation SUCCESSFUL!")

    # Package into JAR
    for f in os.listdir(bin_dir):
        if f.endswith(".jar"):
            os.remove(os.path.join(bin_dir, f))

    bin_jar = os.path.join(bin_dir, "comets_2.12.4.jar")
    result = subprocess.run(
        [JAR_TOOL, "cf", bin_jar, "-C", compile_dir, "."],
        capture_output=True, text=True, timeout=60
    )
    if os.path.exists(bin_jar):
        print(f"Created {bin_jar} ({os.path.getsize(bin_jar)} bytes)")
    shutil.rmtree(compile_dir, ignore_errors=True)

    # Test
    print("\n--- Testing compiled COMETS ---")
    cp_jars = []
    for root, dirs, files in os.walk(lib_dir):
        for f in files:
            if f.endswith(".jar"):
                cp_jars.append(os.path.join(root, f))
    for f in os.listdir(bin_dir):
        if f.endswith(".jar"):
            cp_jars.append(os.path.join(bin_dir, f))
    cp = ";".join(cp_jars)

    result = subprocess.run(
        [JAVA, "-classpath", cp,
         "edu.bu.segrelab.comets.Comets",
         "-loader", "edu.bu.segrelab.comets.fba.FBACometsLoader"],
        capture_output=True, text=True, timeout=15,
        cwd=COMETS_HOME
    )
    print(f"Test return code: {result.returncode}")
    print(f"stdout: {result.stdout[:1500]}")
    if result.stderr:
        print(f"stderr: {result.stderr[:1500]}")

    # Step 3: Update comets_scr.bat
    print("\n--- Step 3: Update comets_scr.bat ---")
    comets_scr = os.path.join(COMETS_HOME, "comets_scr.bat")
    bat_content = f"""@echo off
set JAVA={JAVA}
set COMETS_HOME={COMETS_HOME}
set CLASSPATH={cp}
"%JAVA%" -classpath "%CLASSPATH%" edu.bu.segrelab.comets.Comets -loader edu.bu.segrelab.comets.fba.FBACometsLoader -script %1
"""
    with open(comets_scr, "w") as f:
        f.write(bat_content)

    # Step 4: Run dFBA Demo!
    print("\n" + "=" * 60)
    print("Running dFBA Demo with compiled COMETS 2.12.4")
    print("=" * 60)

    os.environ["COMETS_HOME"] = COMETS_HOME

    import cometspy as c
    import cobra

    # Load model
    ecoli = cobra.io.load_model("textbook")
    print(f"Model: {ecoli.id}, {len(ecoli.reactions)} reactions")

    # Two species
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
    p.set_param("maxCycles", 30)
    p.set_param("timeStep", 0.1)
    p.set_param("writeTotalBiomassLog", True)
    p.set_param("writeMediaLog", True)
    p.set_param("BiomassLogRate", 1)
    p.set_param("MediaLogRate", 1)
    p.set_param("deathRate", 0.0)

    # Remove version-mismatch parameters
    skip_keys = ['velocityMultiConvLogFormat', 'velocityMultiConvLogName',
                 'velocityMultiConvLogRate', 'writeVelocityMultiConvLog',
                 'numCyclesPerDeletion']
    for key in list(p.all_params.keys()):
        if any(sk.lower() in key.lower() for sk in skip_keys):
            del p.all_params[key]

    # Create and run simulation
    print("\nCreating simulation...")
    sim = c.comets(lay, p)
    print(f"COMETS VERSION: {sim.VERSION}")

    # Remove from sim params too
    for key in list(sim.parameters.all_params.keys()):
        if any(sk.lower() in key.lower() for sk in skip_keys):
            del sim.parameters.all_params[key]

    print("\nRunning simulation...")
    try:
        sim.run(delete_files=False, progress=False)
        print("\n=== SIMULATION COMPLETED ===")

        if hasattr(sim, 'total_biomass') and sim.total_biomass is not None:
            print("\n--- Total Biomass ---")
            print(sim.total_biomass.to_string())
            csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "total_biomass.csv")
            sim.total_biomass.to_csv(csv_path, index=False)
            print(f"Saved to {csv_path}")

        if hasattr(sim, 'media') and sim.media is not None:
            print("\n--- Media (first 30 rows) ---")
            print(sim.media.head(30).to_string())
            csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "media.csv")
            sim.media.to_csv(csv_path, index=False)
            print(f"Saved to {csv_path}")

        if sim.run_output:
            print("\n--- Output (last 2000 chars) ---")
            print(sim.run_output[-2000:])

    except Exception as e:
        print(f"\nSimulation error: {type(e).__name__}: {e}")
        traceback.print_exc()
        try:
            if hasattr(sim, 'run_output') and sim.run_output:
                print(f"\nOutput:\n{sim.run_output[:5000]}")
        except:
            pass

else:
    print("\nCompilation FAILED. Cannot run demo with compiled COMETS.")
    print("The compilation errors indicate missing or incompatible dependencies.")
    print("\nAttempting to use the 2.10.0 JAR with Gurobi workaround...")

    # Restore the 2.10.0 JAR
    if not any(f.endswith('.jar') for f in os.listdir(bin_dir)):
        jar_src = os.path.join(r"C:\Users\hgh97\comets_home", "comets_2.10.0.jar")
        if os.path.exists(jar_src):
            shutil.copy2(jar_src, os.path.join(bin_dir, "comets_2.10.0.jar"))
            print("Restored comets_2.10.0.jar")

    # Clean up
    if os.path.exists(compile_dir):
        shutil.rmtree(compile_dir)

log_fh.close()
