"""
Step 2: Explore COMETS source structure and build COMETS_HOME
"""
import os
import sys
import shutil
import subprocess
import traceback

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "explore_log.txt")

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

repo = r"C:\Users\hgh97\comets_home\COMETS"

print("=== COMETS Source Structure ===")
for root, dirs, files in os.walk(repo):
    depth = root.replace(repo, "").count(os.sep)
    if depth > 3:
        dirs.clear()
        continue
    indent = "  " * depth
    print(f"{indent}{os.path.basename(root)}/")
    subindent = "  " * (depth + 1)
    for f in files[:20]:
        print(f"{subindent}{f}")
    if len(files) > 20:
        print(f"{subindent}... and {len(files)-20} more")

# Check .classpath
cp_file = os.path.join(repo, ".classpath")
if os.path.exists(cp_file):
    print("\n=== .classpath ===")
    with open(cp_file, "r") as f:
        print(f.read())

# Check comets_simplified
simple = os.path.join(repo, "comets_simplified")
if os.path.exists(simple):
    print("\n=== comets_simplified structure ===")
    for root, dirs, files in os.walk(simple):
        depth = root.replace(simple, "").count(os.sep)
        if depth > 3:
            dirs.clear()
            continue
        indent = "  " * depth
        print(f"{indent}{os.path.basename(root)}/")
        subindent = "  " * (depth + 1)
        for f in files[:20]:
            print(f"{subindent}{f}")

    # Read README.txt in comets_simplified
    readme = os.path.join(simple, "README.txt")
    if os.path.exists(readme):
        print("\n=== comets_simplified/README.txt ===")
        with open(readme, "r") as f:
            print(f.read())

# Check .travis.yml for build instructions
travis = os.path.join(repo, ".travis.yml")
if os.path.exists(travis):
    print("\n=== .travis.yml ===")
    with open(travis, "r") as f:
        print(f.read())

# Check the full dependencies.md
dep_file = os.path.join(repo, "dependencies.md")
if os.path.exists(dep_file):
    print("\n=== dependencies.md (full) ===")
    with open(dep_file, "r") as f:
        print(f.read())

# Now try to build COMETS_HOME structure from the downloaded jar
# The comets_2.10.0.jar can serve as the main executable
COMETS_HOME = r"C:\Users\hgh97\comets_home"

# Create bin directory with the JAR
bin_dir = os.path.join(COMETS_HOME, "bin")
os.makedirs(bin_dir, exist_ok=True)

# Copy the JAR as the main bin
jar_src = os.path.join(COMETS_HOME, "comets_2.10.0.jar")
jar_dst = os.path.join(bin_dir, "comets_2.10.0.jar")
if os.path.exists(jar_src) and not os.path.exists(jar_dst):
    shutil.copy2(jar_src, jar_dst)
    print(f"\nCopied JAR to bin/: {jar_dst}")

# For Windows, cometspy expects a comets_scr script in COMETS_HOME
# Let's check what the script should look like
# From the comets.py source, on Windows it runs:
#   COMETS_HOME\comets_scr  "script_file"
# This is likely a batch file or shell script

# Create comets_scr.bat for Windows
scr_path = os.path.join(COMETS_HOME, "comets_scr.bat")
java_exe = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8\bin\java.exe"

# First, let's check what's in the JAR
print("\n=== JAR contents (comets_2.10.0.jar) ===")
try:
    result = subprocess.run(
        [java_exe, "-jar", jar_dst, "--help"],
        capture_output=True, text=True, timeout=15,
        cwd=COMETS_HOME
    )
    print(f"Return code: {result.returncode}")
    print(f"stdout: {result.stdout[:2000]}")
    print(f"stderr: {result.stderr[:2000]}")
except Exception as e:
    print(f"JAR execution test: {e}")

# List JAR contents
try:
    result = subprocess.run(
        [java_exe, "-jar", jar_dst],
        capture_output=True, text=True, timeout=15,
        cwd=COMETS_HOME
    )
    print(f"\nJAR run (no args) return code: {result.returncode}")
    print(f"stdout: {result.stdout[:2000]}")
    print(f"stderr: {result.stderr[:2000]}")
except Exception as e:
    print(f"JAR execution: {e}")

# Try listing jar contents with jar tool
jar_tool = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8\bin\jar.exe"
if os.path.exists(jar_tool):
    try:
        result = subprocess.run(
            [jar_tool, "tf", jar_dst],
            capture_output=True, text=True, timeout=15
        )
        lines = result.stdout.strip().split("\n")
        print(f"\nJAR contains {len(lines)} entries. Key entries:")
        for line in lines:
            if any(kw in line for kw in ["Comets", "comets", "FBA", "fba", "Loader", "MANIFEST", "Main"]):
                print(f"  {line}")
    except Exception as e:
        print(f"jar tf failed: {e}")

log_fh.close()
