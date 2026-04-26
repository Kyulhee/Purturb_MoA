"""
Step 9: Compile COMETS from source with UTF-8 encoding
The 2.10.0 JAR requires Gurobi - we need to compile 2.12.4 source with GLOP support
"""
import os
import sys
import subprocess
import traceback
import glob as globmod
import shutil

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compile_log.txt")

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
print("Compile COMETS 2.12.4 from Source with UTF-8 Encoding")
print("=" * 60)

# Step 1: Fix encoding issues in Java source
print("\n--- Step 1: Fix encoding issues ---")
src_dir = os.path.join(COMETS_REPO, "comets_simplified", "src")

# Find the problematic file
for root, dirs, files in os.walk(src_dir):
    for f in files:
        if f.endswith(".java"):
            fpath = os.path.join(root, f)
            try:
                with open(fpath, 'r', encoding='utf-8') as fh:
                    fh.read()
            except UnicodeDecodeError:
                print(f"  Encoding issue: {fpath}")
                # Read with latin-1 and fix
                with open(fpath, 'r', encoding='latin-1') as fh:
                    content = fh.read()
                # Replace non-ASCII characters in comments
                # The problematic line: z[i] >= v[i] and z[i] >= -v[i]
                # The special chars are likely >= (greater-than-or-equal) symbols
                import re
                # Replace unicode math symbols with ASCII equivalents
                replacements = {
                    '\u2265': '>=',  # >=
                    '\u2264': '<=',  # <=
                    '\u2260': '!=',  # !=
                    '\u00b1': '+/-',  # +/-
                }
                for old, new in replacements.items():
                    content = content.replace(old, new)
                with open(fpath, 'w', encoding='utf-8') as fh:
                    fh.write(content)
                print(f"  Fixed: {fpath}")

# Step 2: Collect all JARs for classpath
print("\n--- Step 2: Build classpath ---")
classpath_jars = []
for root, dirs, files in os.walk(lib_dir):
    for f in files:
        if f.endswith(".jar"):
            classpath_jars.append(os.path.join(root, f))

# Remove the old 2.10.0 JAR from bin if it exists
for f in os.listdir(bin_dir):
    if f.endswith(".jar"):
        fpath = os.path.join(bin_dir, f)
        classpath_jars.append(fpath)

classpath = ";".join(classpath_jars)
print(f"Classpath: {len(classpath_jars)} JARs")

# Step 3: Compile
print("\n--- Step 3: Compile Java source ---")
compile_dir = os.path.join(COMETS_HOME, "build")
if os.path.exists(compile_dir):
    shutil.rmtree(compile_dir)
os.makedirs(compile_dir, exist_ok=True)

# Collect Java source files
java_files = []
for root, dirs, files in os.walk(src_dir):
    for f in files:
        if f.endswith(".java"):
            java_files.append(os.path.join(root, f))

print(f"Found {len(java_files)} Java source files")

# Compile with UTF-8 encoding
cmd = [JAVAC, "-encoding", "UTF-8", "-cp", classpath, "-d", compile_dir, "-sourcepath", src_dir]
cmd.extend(java_files)

print(f"Compiling {len(java_files)} files with -encoding UTF-8...")
result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
print(f"Compile return code: {result.returncode}")

if result.stdout:
    print(f"stdout: {result.stdout[:2000]}")
if result.stderr:
    # Print only first 3000 chars of errors
    stderr_text = result.stderr
    print(f"stderr ({len(stderr_text)} chars):")
    print(stderr_text[:3000])

if result.returncode == 0:
    print("\nCompilation successful!")

    # Step 4: Package into JAR
    print("\n--- Step 4: Package into JAR ---")
    # Remove old JAR
    for f in os.listdir(bin_dir):
        if f.endswith(".jar"):
            os.remove(os.path.join(bin_dir, f))

    bin_jar = os.path.join(bin_dir, "comets_2.12.4.jar")
    result = subprocess.run(
        [JAR_TOOL, "cf", bin_jar, "-C", compile_dir, "."],
        capture_output=True, text=True, timeout=60
    )
    print(f"JAR creation return code: {result.returncode}")
    if os.path.exists(bin_jar):
        print(f"Created {bin_jar} ({os.path.getsize(bin_jar)} bytes)")
    else:
        print("JAR creation FAILED")

    # Clean up build dir
    shutil.rmtree(compile_dir, ignore_errors=True)

    # Step 5: Test the compiled JAR
    print("\n--- Step 5: Test compiled COMETS ---")

    # Rebuild classpath with new JAR
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
    print(f"stdout: {result.stdout[:1000]}")
    if result.stderr:
        print(f"stderr: {result.stderr[:1000]}")

else:
    print("\nCompilation FAILED. Checking specific errors...")

    # Try to identify which files have issues
    error_files = set()
    for line in result.stderr.split('\n'):
        if '.java:' in line:
            parts = line.split(':')
            if len(parts) >= 1:
                error_files.add(parts[0])

    print(f"\nFiles with errors: {len(error_files)}")
    for f in sorted(error_files):
        print(f"  {f}")

    # Try compiling without the problematic files
    print("\nAttempting compilation without problematic files...")
    good_files = [f for f in java_files if f not in error_files]
    if good_files:
        cmd = [JAVAC, "-encoding", "UTF-8", "-cp", classpath, "-d", compile_dir, "-sourcepath", src_dir]
        cmd.extend(good_files)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        print(f"Partial compile return code: {result.returncode}")
        if result.stderr:
            print(f"stderr: {result.stderr[:2000]}")

# Step 6: Update comets_scr.bat
print("\n--- Step 6: Update comets_scr.bat ---")
comets_scr = os.path.join(COMETS_HOME, "comets_scr.bat")

cp_jars = []
for root, dirs, files in os.walk(lib_dir):
    for f in files:
        if f.endswith(".jar"):
            cp_jars.append(os.path.join(root, f))
for f in os.listdir(bin_dir):
    if f.endswith(".jar"):
        cp_jars.append(os.path.join(bin_dir, f))
cp = ";".join(cp_jars)

bat_content = f"""@echo off
set JAVA={JAVA}
set COMETS_HOME={COMETS_HOME}
set CLASSPATH={cp}
"%JAVA%" -classpath "%CLASSPATH%" edu.bu.segrelab.comets.Comets -loader edu.bu.segrelab.comets.fba.FBACometsLoader -script %1
"""
with open(comets_scr, "w") as f:
    f.write(bat_content)
print(f"Updated {comets_scr}")

log_fh.close()
