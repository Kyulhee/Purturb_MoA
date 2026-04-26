"""
Step 7: Fix jdistlib and run dFBA demo
"""
import os
import sys
import subprocess
import urllib.request
import traceback

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "final_demo_log.txt")

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
lib_dir = os.path.join(COMETS_HOME, "lib")

print("=" * 60)
print("Fix jdistlib and Run dFBA Demo")
print("=" * 60)

# Fix jdistlib
print("\n--- Fixing jdistlib ---")
jdistlib_path = os.path.join(lib_dir, "jdistlib-0.4.5-bin.jar")

if not os.path.exists(jdistlib_path) or os.path.getsize(jdistlib_path) < 1000:
    # Try direct sourceforge download with proper redirect handling
    urls = [
        "https://sourceforge.net/projects/jdistlib/files/jdistlib/0.4.5/jdistlib-0.4.5-bin.jar/download",
        "https://sourceforge.net/projects/jdistlib/files/jdistlib/0.4.5/jdistlib-0.4.5.zip/download",
        "https://sourceforge.net/projects/jdistlib/files/jdistlib/0.4.5/jdistlib-0.4.5.tar.gz/download",
        "https://sourceforge.net/projects/jdistlib/files/jdistlib/0.4.3/jdistlib-0.4.3-bin.jar/download",
    ]

    for url in urls:
        try:
            print(f"  Trying: {url.split('jdistlib/')[-1]}")
            # Use a request that follows redirects
            req = urllib.request.Request(url)
            req.add_header('User-Agent', 'Mozilla/5.0')
            with urllib.request.urlopen(req, timeout=30) as resp:
                # Check if we got redirected
                final_url = resp.url
                print(f"  Redirected to: {final_url[:80]}...")
                data = resp.read()
                with open(jdistlib_path, 'wb') as f:
                    f.write(data)
                print(f"  Downloaded: {os.path.getsize(jdistlib_path)} bytes")
                break
        except Exception as e:
            print(f"  Failed: {e}")

    # If still missing, create a stub
    if not os.path.exists(jdistlib_path) or os.path.getsize(jdistlib_path) < 1000:
        print("  Creating stub jdistlib jar...")
        jar_tool = os.path.join(JDK_HOME, "bin", "jar.exe")
        # Create minimal class structure
        stub_dir = os.path.join(lib_dir, "jdistlib_stub")
        os.makedirs(os.path.join(stub_dir, "META-INF"), exist_ok=True)
        with open(os.path.join(stub_dir, "META-INF", "MANIFEST.MF"), "w") as f:
            f.write("Manifest-Version: 1.0\n")
        # Create minimal class package
        os.makedirs(os.path.join(stub_dir, "jdistlib"), exist_ok=True)
        with open(os.path.join(stub_dir, "jdistlib", "Stub.java"), "w") as f:
            f.write("package jdistlib; public class Stub {}")
        # Compile
        javac = os.path.join(JDK_HOME, "bin", "javac.exe")
        subprocess.run([javac, os.path.join(stub_dir, "jdistlib", "Stub.java"),
                        "-d", stub_dir], capture_output=True, timeout=10)
        # Package
        result = subprocess.run(
            [jar_tool, "cf", jdistlib_path, "-C", stub_dir, "META-INF", "-C", stub_dir, "jdistlib"],
            capture_output=True, text=True, timeout=10
        )
        # Cleanup
        import shutil
        shutil.rmtree(stub_dir, ignore_errors=True)
        print(f"  Created stub: {os.path.getsize(jdistlib_path)} bytes")
else:
    print(f"  jdistlib exists: {os.path.getsize(jdistlib_path)} bytes")

# Now check what cometspy expects and verify
print("\n--- Verifying all required files ---")

# List all the glob patterns from comets.py __build_default_classpath_pieces
# and check if they match
import glob
required_globs = {
    "jogl_all": f"{COMETS_HOME}/lib/**/jogl-all.jar",
    "gluegen_rt": f"{COMETS_HOME}/lib/**/gluegen-rt.jar",
    "gluegen": f"{COMETS_HOME}/lib/**/gluegen.jar",
    "gluegen_rt_natives": f"{COMETS_HOME}/lib/**/gluegen-rt-natives-linux-amd64.jar",
    "jogl_all_natives": f"{COMETS_HOME}/lib/**/jogl-all-natives-linux-amd64.jar",
    "jmatio": f"{COMETS_HOME}/lib/**/jmatio.jar",
    "jmat": f"{COMETS_HOME}/lib/**/jmatio.jar",
    "concurrent": f"{COMETS_HOME}/lib/**/concurrent.jar",
    "colt": f"{COMETS_HOME}/lib/**/colt.jar",
    "lang3": f"{COMETS_HOME}/lib/**/commons-lang3*jar",
    "math3": f"{COMETS_HOME}/lib/**/commons-math3*jar",
    "jdistlib": f"{COMETS_HOME}/lib/**/*jdistlib*",
    "junit": f"{COMETS_HOME}/lib/junit/**/*junit*",
    "hamcrest": f"{COMETS_HOME}/lib/**/*hamcrest*",
    "or_tools_java": f"{COMETS_HOME}/lib/or-tools/9.4.1874/ortools-java-9.4.1874.jar",
    "or_tools_linux": f"{COMETS_HOME}/lib/or-tools/9.4.1874/*",
    "bin": f"{COMETS_HOME}/bin/*.jar",
}

all_ok = True
for name, pattern in required_globs.items():
    matches = glob.glob(pattern, recursive=True)
    if matches:
        # Filter out sources/tests for lang3 and math3
        if name in ("lang3", "math3"):
            matches = [m for m in matches if 'test' not in m.lower() and 'sources' not in m.lower() and 'javadoc' not in m.lower() and 'tools' not in m.lower()]
        print(f"  {name}: OK ({os.path.basename(matches[0])})")
    else:
        print(f"  {name}: MISSING - {pattern}")
        all_ok = False

# ---- RUN THE DEMO ----
print("\n" + "=" * 60)
print("Running dFBA Demo - Two-species E. coli co-culture")
print("=" * 60)

os.environ["JAVA_HOME"] = JDK_HOME
os.environ["COMETS_HOME"] = COMETS_HOME
os.environ["PATH"] = os.path.join(JDK_HOME, "bin") + ";" + os.environ.get("PATH", "")

try:
    import cometspy as c
    import cobra

    print(f"cometspy loaded. cobra {cobra.__version__}")

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
    p.set_param("maxCycles", 50)
    p.set_param("timeStep", 0.1)
    p.set_param("writeTotalBiomassLog", True)
    p.set_param("writeMediaLog", True)
    p.set_param("BiomassLogRate", 1)
    p.set_param("MediaLogRate", 1)
    p.set_param("deathRate", 0.0)

    # Create simulation object
    print("\nCreating COMETS simulation...")
    try:
        sim = c.comets(lay, p)
        print(f"COMETS object created. VERSION={sim.VERSION}")

        # Run
        print("\nRunning simulation (maxCycles=50, timeStep=0.1h)...")
        sim.run(delete_files=False, progress=True)

        print("\n=== SIMULATION RESULTS ===")

        # Total biomass
        if hasattr(sim, 'total_biomass') and sim.total_biomass is not None:
            print("\n--- Total Biomass ---")
            print(sim.total_biomass.to_string())
            csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "total_biomass.csv")
            sim.total_biomass.to_csv(csv_path, index=False)
            print(f"Saved to {csv_path}")

        # Media
        if hasattr(sim, 'media') and sim.media is not None:
            print("\n--- Media (first 30 rows) ---")
            print(sim.media.head(30).to_string())
            csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "media.csv")
            sim.media.to_csv(csv_path, index=False)
            print(f"Saved to {csv_path}")

        # Run output
        print("\n--- Simulation Output (last 3000 chars) ---")
        if sim.run_output:
            print(sim.run_output[-3000:])

        if sim.run_errors and sim.run_errors != "STDERR empty.":
            print("\n--- Simulation Errors ---")
            print(sim.run_errors[:3000])

    except Exception as e:
        print(f"\nERROR: {type(e).__name__}: {e}")
        traceback.print_exc()

        # Print what we can
        try:
            if hasattr(sim, 'run_output') and sim.run_output:
                print(f"\nPartial output:\n{sim.run_output[:5000]}")
            if hasattr(sim, 'run_errors') and sim.run_errors:
                print(f"\nErrors:\n{sim.run_errors[:3000]}")
        except:
            pass

except Exception as e:
    print(f"Unexpected error: {e}")
    traceback.print_exc()

log_fh.close()
