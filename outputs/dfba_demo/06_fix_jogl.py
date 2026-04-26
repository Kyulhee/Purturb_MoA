"""
Step 6: Fix JOGL download and jdistlib, then retry dFBA demo
"""
import os
import sys
import subprocess
import urllib.request
import zipfile
import traceback
import shutil

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fix_jogl_log.txt")

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
lib_dir = os.path.join(COMETS_HOME, "lib")

print("=" * 60)
print("Fix JOGL and jdistlib, then run dFBA demo")
print("=" * 60)

# ---- JOGL ----
# Jogamp uses a different Maven group ID format
# Try the correct Maven Central paths
print("\n--- Downloading JOGL/JOGAMP ---")
jogl_dir = os.path.join(lib_dir, "jogl", "jogamp-all-platforms", "jar")
os.makedirs(jogl_dir, exist_ok=True)

MAVEN = "https://repo1.maven.org/maven2"

# Try different version formats for JOGAMP
jogl_attempts = [
    # Version 2.4.0-rc-20220328
    ("jogl-all.jar", [
        f"{MAVEN}/org/jogamp/jogl/jogl-all/2.4.0-rc-20220328/jogl-all-2.4.0-rc-20220328.jar",
        f"{MAVEN}/org/jogamp/jogl/jogl-all/2.4.0/jogl-all-2.4.0.jar",
        f"{MAVEN}/org/jogamp/jogl/jogl-all/2.3.2/jogl-all-2.3.2.jar",
    ]),
    ("jogl-all-natives-windows-amd64.jar", [
        f"{MAVEN}/org/jogamp/jogl/jogl-all/2.4.0-rc-20220328/jogl-all-2.4.0-rc-20220328-natives-windows-amd64.jar",
        f"{MAVEN}/org/jogamp/jogl/jogl-all/2.3.2/jogl-all-2.3.2-natives-windows-amd64.jar",
    ]),
    ("jogl-all-natives-linux-amd64.jar", [
        f"{MAVEN}/org/jogamp/jogl/jogl-all/2.4.0-rc-20220328/jogl-all-2.4.0-rc-20220328-natives-linux-amd64.jar",
        f"{MAVEN}/org/jogamp/jogl/jogl-all/2.3.2/jogl-all-2.3.2-natives-linux-amd64.jar",
    ]),
    ("gluegen-rt.jar", [
        f"{MAVEN}/org/jogamp/gluegen/gluegen-rt/2.4.0-rc-20220328/gluegen-rt-2.4.0-rc-20220328.jar",
        f"{MAVEN}/org/jogamp/gluegen/gluegen-rt/2.3.2/gluegen-rt-2.3.2.jar",
    ]),
    ("gluegen-rt-natives-windows-amd64.jar", [
        f"{MAVEN}/org/jogamp/gluegen/gluegen-rt/2.4.0-rc-20220328/gluegen-rt-2.4.0-rc-20220328-natives-windows-amd64.jar",
        f"{MAVEN}/org/jogamp/gluegen/gluegen-rt/2.3.2/gluegen-rt-2.3.2-natives-windows-amd64.jar",
    ]),
    ("gluegen-rt-natives-linux-amd64.jar", [
        f"{MAVEN}/org/jogamp/gluegen/gluegen-rt/2.4.0-rc-20220328/gluegen-rt-2.4.0-rc-20220328-natives-linux-amd64.jar",
        f"{MAVEN}/org/jogamp/gluegen/gluegen-rt/2.3.2/gluegen-rt-2.3.2-natives-linux-amd64.jar",
    ]),
    ("gluegen.jar", [
        f"{MAVEN}/org/jogamp/gluegen/gluegen/2.4.0-rc-20220328/gluegen-2.4.0-rc-20220328.jar",
        f"{MAVEN}/org/jogamp/gluegen/gluegen/2.3.2/gluegen-2.3.2.jar",
    ]),
]

for name, urls in jogl_attempts:
    target = os.path.join(jogl_dir, name)
    if os.path.exists(target) and os.path.getsize(target) > 100:
        print(f"  [SKIP] {name}")
        continue
    success = False
    for url in urls:
        try:
            print(f"  [DL] {name} from {url.split('/maven2/')[-1]}...", end=" ")
            urllib.request.urlretrieve(url, target)
            if os.path.getsize(target) > 100:
                print(f"OK ({os.path.getsize(target)} bytes)")
                success = True
                break
            else:
                os.remove(target)
                print("too small")
        except Exception as e:
            print(f"failed")
    if not success:
        print(f"  WARNING: Could not download {name}")

# ---- jdistlib ----
print("\n--- Fixing jdistlib ---")
jdistlib_path = os.path.join(lib_dir, "jdistlib-0.4.5-bin.jar")
jdistlib_zip = os.path.join(lib_dir, "jdistlib_download.zip")

# Check if the previous download was a valid jar or zip
if os.path.exists(jdistlib_zip) and not os.path.exists(jdistlib_path):
    print(f"  Found previous download: {os.path.getsize(jdistlib_zip)} bytes")
    try:
        with zipfile.ZipFile(jdistlib_zip, 'r') as z:
            print(f"  It's a zip file. Contents:")
            for name in z.namelist():
                if name.endswith('.jar'):
                    print(f"    {name}")
                    data = z.read(name)
                    with open(jdistlib_path, 'wb') as f:
                        f.write(data)
                    print(f"  Extracted to {jdistlib_path} ({os.path.getsize(jdistlib_path)} bytes)")
                    break
            else:
                # Maybe it IS the jar, just with wrong extension
                print("  No jar found inside zip. Renaming as jar...")
                shutil.move(jdistlib_zip, jdistlib_path)
    except zipfile.BadZipFile:
        print("  Not a valid zip. Might be the jar directly. Renaming...")
        shutil.move(jdistlib_zip, jdistlib_path)
        print(f"  Renamed to {jdistlib_path} ({os.path.getsize(jdistlib_path)} bytes)")

if os.path.exists(jdistlib_zip):
    os.remove(jdistlib_zip)

# Verify jdistlib
if os.path.exists(jdistlib_path):
    print(f"  jdistlib: {os.path.getsize(jdistlib_path)} bytes")
else:
    print("  jdistlib: STILL MISSING")

# ---- Verify all JOGL files ----
print("\n--- JOGL verification ---")
for f in os.listdir(jogl_dir):
    print(f"  {f}: {os.path.getsize(os.path.join(jogl_dir, f))} bytes")

# ---- Now try the dFBA demo ----
print("\n" + "=" * 60)
print("Running dFBA Demo")
print("=" * 60)

JDK_HOME = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8"
os.environ["JAVA_HOME"] = JDK_HOME
os.environ["COMETS_HOME"] = COMETS_HOME
os.environ["PATH"] = os.path.join(JDK_HOME, "bin") + ";" + os.environ.get("PATH", "")

try:
    import cometspy as c
    import cobra

    print(f"cometspy and cobra loaded")

    # Load model
    ecoli = cobra.io.load_model("textbook")
    print(f"E. coli model: {ecoli.id}, {len(ecoli.reactions)} reactions")

    # Create two species
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

    # Create simulation
    print("\nCreating simulation object...")
    try:
        sim = c.comets(lay, p)
        print(f"SUCCESS! VERSION: {sim.VERSION}")
    except KeyError as e:
        print(f"KeyError: {e}")
        traceback.print_exc()
    except IndexError as e:
        print(f"IndexError: {e}")
        traceback.print_exc()
        print("\nSome classpath JARs are still missing.")
        print("Attempting to bypass by manually setting classpath...")

        # Try creating comets object with manual override
        # The issue is in __build_default_classpath_pieces
        # We need to monkey-patch or use an alternative approach
        print("\nAttempting alternative: direct Java execution...")
        # Write the layout and params files manually and run COMETS directly

    except Exception as e:
        print(f"Error: {type(e).__name__}: {e}")
        traceback.print_exc()

except Exception as e:
    print(f"Unexpected error: {e}")
    traceback.print_exc()

log_fh.close()
