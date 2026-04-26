"""
Step 1: Download and install COMETS Java core
Tries multiple approaches:
  1. Download from runcomets.org
  2. Download from GitHub releases (v2.10.0)
  3. Clone and build from source
"""
import urllib.request
import os
import json
import sys
import traceback

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "install_log.txt")

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

print("=" * 60)
print("COMETS Java Core Installation Script")
print("=" * 60)

# Check Java
print("\n--- Java Check ---")
jdk_path = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8"
if os.path.exists(jdk_path):
    java_exe = os.path.join(jdk_path, "bin", "java.exe")
    print(f"JDK found at: {jdk_path}")
    print(f"java.exe exists: {os.path.exists(java_exe)}")
    os.environ["JAVA_HOME"] = jdk_path
    os.environ["PATH"] = os.path.join(jdk_path, "bin") + ";" + os.environ.get("PATH", "")
    print(f"JAVA_HOME set to: {jdk_path}")
else:
    print("JDK not found at expected location")

# Test java
import subprocess
try:
    result = subprocess.run(["java", "-version"], capture_output=True, text=True, timeout=10)
    print(f"Java version check: returncode={result.returncode}")
    print(f"  stdout: {result.stdout[:200]}")
    print(f"  stderr: {result.stderr[:200]}")
except Exception as e:
    print(f"Java check failed: {e}")

# Try downloading COMETS
COMETS_HOME = r"C:\Users\hgh97\comets_home"
os.makedirs(COMETS_HOME, exist_ok=True)

success = False

# Approach 1: Try downloading from runcomets.org
print("\n--- Approach 1: runcomets.org ---")
try:
    # Try to find download URL
    base_url = "https://www.runcomets.org"
    req = urllib.request.Request(base_url)
    req.add_header('User-Agent', 'Mozilla/5.0')
    with urllib.request.urlopen(req, timeout=15) as resp:
        content = resp.read().decode('utf-8', errors='ignore')
        print(f"Page length: {len(content)}")
        # Look for download links
        import re
        links = re.findall(r'href=["\']([^"\']*(?:download|zip|jar|comets)[^"\']*)["\']', content, re.I)
        print(f"Potential download links: {links[:10]}")
except Exception as e:
    print(f"runcomets.org failed: {e}")

# Approach 2: GitHub API - check all releases and their assets
print("\n--- Approach 2: GitHub Releases ---")
try:
    url = 'https://api.github.com/repos/segrelab/COMETS/releases'
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Python')
    with urllib.request.urlopen(req, timeout=15) as resp:
        releases = json.loads(resp.read().decode())
        for r in releases[:5]:
            print(f'Release: {r["tag_name"]} - {r["name"]}')
            assets = r.get('assets', [])
            if assets:
                for a in assets:
                    print(f'  Asset: {a["name"]} - {a["browser_download_url"]}')
            else:
                print(f'  No assets attached')
            # Check for zipball
            zipball = r.get('zipball_url', '')
            print(f'  Zipball: {zipball}')
except Exception as e:
    print(f"GitHub API failed: {e}")

# Approach 3: Clone the repo and look for pre-built distribution
print("\n--- Approach 3: Clone COMETS repo ---")
comets_repo = os.path.join(COMETS_HOME, "COMETS")
try:
    if not os.path.exists(comets_repo):
        print("Cloning COMETS repository...")
        result = subprocess.run(
            ["git", "clone", "--depth", "1", "https://github.com/segrelab/COMETS.git", comets_repo],
            capture_output=True, text=True, timeout=120
        )
        print(f"Clone returncode: {result.returncode}")
        if result.returncode != 0:
            print(f"  stderr: {result.stderr[:500]}")
        else:
            print("Clone successful!")
            # List structure
            for item in os.listdir(comets_repo):
                print(f"  {item}")
    else:
        print(f"COMETS repo already cloned at {comets_repo}")
        for item in os.listdir(comets_repo):
            print(f"  {item}")
    success = True
except Exception as e:
    print(f"Clone failed: {e}")
    traceback.print_exc()

# Approach 4: Download v2.10.0 JAR (last version with assets)
print("\n--- Approach 4: Download comets_2.10.0.jar ---")
jar_path = os.path.join(COMETS_HOME, "comets_2.10.0.jar")
try:
    jar_url = "https://github.com/segrelab/comets/releases/download/v2.10.0/comets_2.10.0.jar"
    if not os.path.exists(jar_path):
        print(f"Downloading {jar_url}...")
        urllib.request.urlretrieve(jar_url, jar_path)
        print(f"Downloaded to {jar_path} ({os.path.getsize(jar_path)} bytes)")
    else:
        print(f"JAR already exists at {jar_path} ({os.path.getsize(jar_path)} bytes)")
except Exception as e:
    print(f"JAR download failed: {e}")

# Now check if the cloned repo has a build system
print("\n--- Checking for build system ---")
if os.path.exists(comets_repo):
    # Look for build files
    for root, dirs, files in os.walk(comets_repo):
        depth = root.replace(comets_repo, "").count(os.sep)
        if depth > 2:
            continue
        for f in files:
            if f in ['build.xml', 'build.gradle', 'pom.xml', 'Makefile', '.classpath']:
                print(f"  Found: {os.path.join(root, f)}")

    # Check dependencies.md
    dep_file = os.path.join(comets_repo, "dependencies.md")
    if os.path.exists(dep_file):
        print("\n--- dependencies.md content ---")
        with open(dep_file, 'r') as f:
            content = f.read()
            print(content[:2000])

    # Check for lib/ directory
    lib_dir = os.path.join(comets_repo, "lib")
    if os.path.exists(lib_dir):
        print(f"\nlib/ directory contents:")
        for f in os.listdir(lib_dir):
            print(f"  {f}")

    # Check for bin/ directory
    bin_dir = os.path.join(comets_repo, "bin")
    if os.path.exists(bin_dir):
        print(f"\nbin/ directory contents:")
        for f in os.listdir(bin_dir):
            print(f"  {f}")

    # Check for comets_simplified
    simple_dir = os.path.join(comets_repo, "comets_simplified")
    if os.path.exists(simple_dir):
        print(f"\ncomets_simplified/ contents:")
        for item in os.listdir(simple_dir):
            print(f"  {item}")

# Summary
print("\n" + "=" * 60)
print("INSTALLATION SUMMARY")
print("=" * 60)
print(f"COMETS_HOME target: {COMETS_HOME}")
print(f"COMETS repo cloned: {os.path.exists(comets_repo)}")
print(f"JAR downloaded: {os.path.exists(jar_path)}")
if os.path.exists(jar_path):
    print(f"JAR size: {os.path.getsize(jar_path)} bytes")

# Check if we have enough for a working COMETS_HOME
has_bin = os.path.exists(os.path.join(COMETS_HOME, "bin"))
has_lib = os.path.exists(os.path.join(COMETS_HOME, "lib"))
print(f"Has bin/: {has_bin}")
print(f"Has lib/: {has_lib}")
print(f"Ready for cometspy: {has_bin and has_lib}")

log_fh.close()
