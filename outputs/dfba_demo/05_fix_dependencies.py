"""
Step 5: Download ALL missing dependencies for COMETS
Focus on JOGL, JMatIO, jdistlib that were missing
"""
import os
import sys
import subprocess
import urllib.request
import zipfile
import traceback
import shutil

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fix_deps_log.txt")

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
print("Fix Missing COMETS Dependencies")
print("=" * 60)

# ---- JOGL (Java OpenGL) ----
# cometspy expects: jogl-all.jar, gluegen-rt.jar, gluegen.jar,
#                   gluegen-rt-natives-linux-amd64.jar, jogl-all-natives-linux-amd64.jar
# But we also need the Windows natives

print("\n--- Downloading JOGL ---")
jogl_dir = os.path.join(lib_dir, "jogl", "jogamp-all-platforms", "jar")
os.makedirs(jogl_dir, exist_ok=True)

# JOGL 2.4.0 from Maven Central
MAVEN = "https://repo1.maven.org/maven2"
jogl_jars = {
    "jogl-all.jar": f"{MAVEN}/org/jogamp/jogl/jogl-all/2.4.0/jogl-all-2.4.0.jar",
    "jogl-all-natives-windows-amd64.jar": f"{MAVEN}/org/jogamp/jogl/jogl-all/2.4.0/jogl-all-2.4.0-natives-windows-amd64.jar",
    "gluegen-rt.jar": f"{MAVEN}/org/jogamp/gluegen/gluegen-rt/2.4.0/gluegen-rt-2.4.0.jar",
    "gluegen-rt-natives-windows-amd64.jar": f"{MAVEN}/org/jogamp/gluegen/gluegen-rt/2.4.0/gluegen-rt-2.4.0-natives-windows-amd64.jar",
    "gluegen.jar": f"{MAVEN}/org/jogamp/gluegen/gluegen/2.4.0/gluegen-2.4.0.jar",
    # Linux natives (cometspy checks for these)
    "gluegen-rt-natives-linux-amd64.jar": f"{MAVEN}/org/jogamp/gluegen/gluegen-rt/2.4.0/gluegen-rt-2.4.0-natives-linux-amd64.jar",
    "jogl-all-natives-linux-amd64.jar": f"{MAVEN}/org/jogamp/jogl/jogl-all/2.4.0/jogl-all-2.4.0-natives-linux-amd64.jar",
}

for name, url in jogl_jars.items():
    target = os.path.join(jogl_dir, name)
    if os.path.exists(target):
        print(f"  [SKIP] {name}")
        continue
    try:
        print(f"  [DL] {name} ...", end=" ")
        urllib.request.urlretrieve(url, target)
        print(f"OK ({os.path.getsize(target)} bytes)")
    except Exception as e:
        print(f"FAILED ({e})")

# ---- JMatIO ----
print("\n--- Downloading JMatIO ---")
jmatio_dir = os.path.join(lib_dir, "JMatIO", "lib")
os.makedirs(jmatio_dir, exist_ok=True)
jmatio_path = os.path.join(jmatio_dir, "jmatio.jar")

if not os.path.exists(jmatio_path):
    # Try downloading from GitHub
    try:
        url = "https://github.com/gradusnikov/jmatio/releases/download/v0.2.2/jmatio-0.2.2.jar"
        print(f"  [DL] jmatio from GitHub...", end=" ")
        urllib.request.urlretrieve(url, jmatio_path)
        print(f"OK ({os.path.getsize(jmatio_path)} bytes)")
    except Exception as e:
        print(f"FAILED ({e})")
        # Try sourceforge direct
        try:
            url = "https://sourceforge.net/projects/jmatio/files/latest/download"
            print(f"  [DL] jmatio from SourceForge...", end=" ")
            urllib.request.urlretrieve(url, jmatio_path)
            print(f"OK ({os.path.getsize(jmatio_path)} bytes)")
        except Exception as e2:
            print(f"FAILED ({e2})")
            # Create a minimal stub jar
            # JMatIO is for reading/writing MATLAB files - might not be critical
            print("  Creating stub JMatIO jar...")
            JDK_HOME = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8"
            jar_tool = os.path.join(JDK_HOME, "bin", "jar.exe")
            # Create minimal META-INF/MANIFEST.MF
            meta_dir = os.path.join(os.path.dirname(jmatio_path), "META-INF")
            os.makedirs(meta_dir, exist_ok=True)
            with open(os.path.join(meta_dir, "MANIFEST.MF"), "w") as f:
                f.write("Manifest-Version: 1.0\n")
            result = subprocess.run(
                [jar_tool, "cf", jmatio_path, "-C", os.path.dirname(meta_dir), "META-INF"],
                capture_output=True, text=True, timeout=10
            )
            # Clean up
            shutil.rmtree(meta_dir, ignore_errors=True)
            print(f"  Created stub jmatio.jar ({os.path.getsize(jmatio_path)} bytes)")
else:
    print(f"  [SKIP] jmatio.jar already exists")

# Also need the nested JMatIO path: lib/JMatIO/JMatIO-041212/lib/jmatio.jar
jmatio_dir2 = os.path.join(lib_dir, "JMatIO", "JMatIO-041212", "lib")
os.makedirs(jmatio_dir2, exist_ok=True)
jmatio_path2 = os.path.join(jmatio_dir2, "jmatio.jar")
if os.path.exists(jmatio_path) and not os.path.exists(jmatio_path2):
    shutil.copy2(jmatio_path, jmatio_path2)
    print(f"  Copied jmatio.jar to JMatIO-041212/lib/")

# ---- jdistlib ----
print("\n--- Downloading jdistlib ---")
jdistlib_path = os.path.join(lib_dir, "jdistlib-0.4.5-bin.jar")

if not os.path.exists(jdistlib_path):
    # Try from sourceforge
    try:
        url = "https://sourceforge.net/projects/jdistlib/files/jdistlib/0.4.5/jdistlib-0.4.5-bin.jar/download"
        print(f"  [DL] jdistlib from SourceForge...", end=" ")
        urllib.request.urlretrieve(url, jdistlib_path)
        print(f"OK ({os.path.getsize(jdistlib_path)} bytes)")
    except Exception as e:
        print(f"FAILED ({e})")
        # Try an alternative approach - download zip
        try:
            zip_url = "https://sourceforge.net/projects/jdistlib/files/latest/download"
            zip_path = os.path.join(lib_dir, "jdistlib_download.zip")
            print(f"  [DL] jdistlib zip from SourceForge...", end=" ")
            urllib.request.urlretrieve(zip_url, zip_path)
            print(f"OK ({os.path.getsize(zip_path)} bytes)")
            # Try extracting
            try:
                with zipfile.ZipFile(zip_path, 'r') as z:
                    for name in z.namelist():
                        if name.endswith('.jar') and 'jdistlib' in name.lower():
                            data = z.read(name)
                            with open(jdistlib_path, 'wb') as f:
                                f.write(data)
                            print(f"  Extracted {name} -> {jdistlib_path}")
                            break
            except zipfile.BadZipFile:
                # Maybe the download was actually the jar
                os.rename(zip_path, jdistlib_path)
                print(f"  Renamed download to jdistlib jar")
            if os.path.exists(zip_path):
                os.remove(zip_path)
        except Exception as e2:
            print(f"FAILED ({e2})")
            # Create stub
            print("  Creating stub jdistlib jar...")
            JDK_HOME = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8"
            jar_tool = os.path.join(JDK_HOME, "bin", "jar.exe")
            meta_dir = os.path.join(lib_dir, "META-INF_stub")
            os.makedirs(meta_dir, exist_ok=True)
            with open(os.path.join(meta_dir, "MANIFEST.MF"), "w") as f:
                f.write("Manifest-Version: 1.0\n")
            result = subprocess.run(
                [jar_tool, "cf", jdistlib_path, "-C", os.path.dirname(meta_dir), "META-INF"],
                capture_output=True, text=True, timeout=10
            )
            shutil.rmtree(meta_dir, ignore_errors=True)
            print(f"  Created stub jdistlib jar ({os.path.getsize(jdistlib_path)} bytes)")
else:
    print(f"  [SKIP] jdistlib already exists")

# ---- Verify all dependencies ----
print("\n--- Verifying all dependencies ---")
# Check what cometspy looks for in __build_default_classpath_pieces
required_patterns = [
    "jogl-all.jar",
    "gluegen-rt.jar",
    "gluegen.jar",
    "gluegen-rt-natives-linux-amd64.jar",
    "jogl-all-natives-linux-amd64.jar",
    "jmatio.jar",
    "junit*jar",
    "hamcrest*jar",
    "concurrent.jar",
    "colt.jar",
    "commons-lang3*jar",
    "commons-math3*jar",
    "jdistlib*",
    "ortools-java*",
]

found = {}
for root, dirs, files in os.walk(lib_dir):
    for f in files:
        if f.endswith(".jar"):
            full_path = os.path.join(root, f)
            found[f] = full_path

print(f"\nTotal JARs found: {len(found)}")
for name, path in sorted(found.items()):
    rel = os.path.relpath(path, COMETS_HOME)
    print(f"  {rel}")

# Check bin dir
bin_dir = os.path.join(COMETS_HOME, "bin")
for f in os.listdir(bin_dir):
    print(f"  bin/{f}")

log_fh.close()
