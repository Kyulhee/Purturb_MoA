"""
Step 3: Set up COMETS_HOME with all required Java dependencies
Downloads all dependency JARs and creates the proper directory structure.
Then compiles the COMETS Java source code.
"""
import os
import sys
import subprocess
import urllib.request
import traceback
import shutil

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "setup_comets_log.txt")

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

# Set up Java
JDK_HOME = r"C:\Users\hgh97\.jdk\jdk-17.0.18+8"
os.environ["JAVA_HOME"] = JDK_HOME
os.environ["PATH"] = os.path.join(JDK_HOME, "bin") + ";" + os.environ.get("PATH", "")

JAVA = os.path.join(JDK_HOME, "bin", "java.exe")
JAVAC = os.path.join(JDK_HOME, "bin", "javac.exe")
JAR_TOOL = os.path.join(JDK_HOME, "bin", "jar.exe")

COMETS_HOME = r"C:\Users\hgh97\comets_2.12.4"
COMETS_REPO = r"C:\Users\hgh97\comets_home\COMETS"

print("=" * 60)
print("COMETS_HOME Setup Script")
print("=" * 60)

# Verify Java
print(f"\nJava: {os.path.exists(JAVA)} - {JAVA}")
print(f"Javac: {os.path.exists(JAVAC)} - {JAVAC}")

# Check if COMETS_HOME already exists and is complete
bin_dir = os.path.join(COMETS_HOME, "bin")
lib_dir = os.path.join(COMETS_HOME, "lib")
comets_scr = os.path.join(COMETS_HOME, "comets_scr.bat")

if os.path.exists(bin_dir) and os.path.exists(lib_dir) and os.path.exists(comets_scr):
    print(f"\nCOMETS_HOME already set up at {COMETS_HOME}")
    print("Skipping setup.")
else:
    # Create directory structure
    os.makedirs(bin_dir, exist_ok=True)
    os.makedirs(lib_dir, exist_ok=True)

    # ---- Download dependencies ----
    # Maven Central base URL
    MAVEN = "https://repo1.maven.org/maven2"

    dependencies = {
        # colt and concurrent
        "colt/lib/colt.jar": f"{MAVEN}/colt/colt/1.2.0/colt-1.2.0.jar",
        "colt/lib/concurrent.jar": f"{MAVEN}/concurrent/concurrent/1.3.4/concurrent-1.3.4.jar",

        # jdistlib
        "jdistlib-0.4.5-bin.jar": "https://sourceforge.net/projects/jdistlib/files/jdistlib/0.4.5/jdistlib-0.4.5-bin.jar/download",

        # Apache commons-lang3
        "commons-lang3-3.9/commons-lang3-3.9.jar": f"{MAVEN}/org/apache/commons/commons-lang3/3.9/commons-lang3-3.9.jar",
        "commons-lang3-3.9/commons-lang3-3.9-sources.jar": f"{MAVEN}/org/apache/commons/commons-lang3/3.9/commons-lang3-3.9-sources.jar",

        # Apache commons-math3
        "commons-math3-3.6.1/commons-math3-3.6.1.jar": f"{MAVEN}/org/apache/commons/commons-math3/3.6.1/commons-math3-3.6.1.jar",
        "commons-math3-3.6.1/commons-math3-3.6.1-tools.jar": f"{MAVEN}/org/apache/commons/commons-math3/3.6.1/commons-math3-3.6.1-tools.jar",

        # Apache commons-rng
        "commons-rng-1.0/commons-rng-client-api-1.0.jar": f"{MAVEN}/org/apache/commons/commons-rng-client-api/1.0/commons-rng-client-api-1.0.jar",
        "commons-rng-1.0/commons-rng-core-1.0.jar": f"{MAVEN}/org/apache/commons/commons-rng-core/1.0/commons-rng-core-1.0.jar",
        "commons-rng-1.0/commons-rng-simple-1.0.jar": f"{MAVEN}/org/apache/commons/commons-rng-simple/1.0/commons-rng-simple-1.0.jar",
        "commons-rng-1.0/commons-rng-sampling-1.0.jar": f"{MAVEN}/org/apache/commons/commons-rng-sampling/1.0/commons-rng-sampling-1.0.jar",
        "commons-rng-1.0/commons-rng-jmh-1.0.jar": f"{MAVEN}/org/apache/commons/commons-rng-jmh/1.0/commons-rng-jmh-1.0.jar",

        # JUnit
        "junit/junit-4.12.jar": f"{MAVEN}/junit/junit/4.12/junit-4.12.jar",
        "junit/hamcrest-core-1.3.jar": f"{MAVEN}/org/hamcrest/hamcrest-core/1.3/hamcrest-core-1.3.jar",

        # JMatIO
        "JMatIO/lib/jmatio.jar": "https://sourceforge.net/projects/jmatio/files/jmatio/0.2/jmatio-0.2.zip/download",

        # or-tools (for GLOP optimizer)
        "or-tools/9.4.1874/ortools-java-9.4.1874.jar": f"{MAVEN}/com/google/ortools/ortools-java/9.4.1874/ortools-java-9.4.1874.jar",
        "or-tools/9.4.1874/ortools-win32-x86-64-9.4.1874.jar": f"{MAVEN}/com/google/ortools/ortools-win32-x86-64/9.4.1874/ortools-win32-x86-64-9.4.1874.jar",
    }

    print(f"\nDownloading {len(dependencies)} dependencies...")
    success_count = 0
    fail_count = 0
    for rel_path, url in dependencies.items():
        target = os.path.join(lib_dir, rel_path)
        if os.path.exists(target):
            print(f"  [SKIP] {rel_path} (already exists)")
            success_count += 1
            continue
        os.makedirs(os.path.dirname(target), exist_ok=True)
        try:
            print(f"  [DL] {rel_path} ...", end=" ")
            urllib.request.urlretrieve(url, target)
            size = os.path.getsize(target)
            print(f"OK ({size} bytes)")
            success_count += 1
        except Exception as e:
            print(f"FAILED ({e})")
            fail_count += 1
            # Remove partial downloads
            if os.path.exists(target):
                os.remove(target)

    print(f"\nDependencies: {success_count} OK, {fail_count} FAILED")

    # ---- Compile COMETS from source ----
    print("\n--- Compiling COMETS Java source ---")
    src_dir = os.path.join(COMETS_REPO, "comets_simplified", "src")

    if os.path.exists(src_dir):
        # Collect all Java source files
        java_files = []
        for root, dirs, files in os.walk(src_dir):
            for f in files:
                if f.endswith(".java"):
                    java_files.append(os.path.join(root, f))

        print(f"Found {len(java_files)} Java source files")

        # Build classpath from all downloaded JARs
        classpath_jars = []
        for root, dirs, files in os.walk(lib_dir):
            for f in files:
                if f.endswith(".jar"):
                    classpath_jars.append(os.path.join(root, f))

        classpath = ";".join(classpath_jars)
        print(f"Classpath has {len(classpath_jars)} JARs")

        # Compile
        compile_dir = os.path.join(COMETS_HOME, "build")
        os.makedirs(compile_dir, exist_ok=True)

        if java_files:
            cmd = [JAVAC, "-cp", classpath, "-d", compile_dir, "-sourcepath", src_dir]
            cmd.extend(java_files)
            print(f"Compiling {len(java_files)} files...")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            print(f"Compile return code: {result.returncode}")
            if result.stdout:
                print(f"stdout: {result.stdout[:1000]}")
            if result.stderr:
                print(f"stderr: {result.stderr[:2000]}")

            if result.returncode == 0:
                # Package into JAR
                bin_jar = os.path.join(bin_dir, "comets_2.12.4.jar")
                cmd = [JAR_TOOL, "cf", bin_jar, "-C", compile_dir, "."]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
                print(f"JAR creation return code: {result.returncode}")
                if os.path.exists(bin_jar):
                    print(f"Created {bin_jar} ({os.path.getsize(bin_jar)} bytes)")
            else:
                print("Compilation failed. Trying with downloaded JAR instead...")
                # Copy the downloaded JAR
                jar_src = os.path.join(r"C:\Users\hgh97\comets_home", "comets_2.10.0.jar")
                jar_dst = os.path.join(bin_dir, "comets_2.10.0.jar")
                if os.path.exists(jar_src):
                    shutil.copy2(jar_src, jar_dst)
                    print(f"Copied {jar_src} to {jar_dst}")
    else:
        print("Source directory not found, using downloaded JAR")
        jar_src = os.path.join(r"C:\Users\hgh97\comets_home", "comets_2.10.0.jar")
        jar_dst = os.path.join(bin_dir, "comets_2.10.0.jar")
        if os.path.exists(jar_src):
            shutil.copy2(jar_src, jar_dst)

    # ---- Create comets_scr.bat for Windows ----
    print("\n--- Creating comets_scr.bat ---")
    # Build classpath for the bat file
    classpath_jars = []
    for root, dirs, files in os.walk(lib_dir):
        for f in files:
            if f.endswith(".jar"):
                classpath_jars.append(os.path.join(root, f))
    # Add bin JAR
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
    with open(comets_scr, "w") as f:
        f.write(bat_content)
    print(f"Created {comets_scr}")
    print(f"Classpath has {len(classpath_jars)} JARs")

# ---- Verify COMETS_HOME ----
print("\n" + "=" * 60)
print("COMETS_HOME VERIFICATION")
print("=" * 60)
print(f"COMETS_HOME: {COMETS_HOME}")
print(f"bin/ exists: {os.path.exists(bin_dir)}")
if os.path.exists(bin_dir):
    for f in os.listdir(bin_dir):
        print(f"  {f}")
print(f"lib/ exists: {os.path.exists(lib_dir)}")
if os.path.exists(lib_dir):
    for item in os.listdir(lib_dir):
        print(f"  {item}/")
print(f"comets_scr.bat exists: {os.path.exists(comets_scr)}")

# Try running COMETS
print("\n--- Test COMETS execution ---")
try:
    result = subprocess.run(
        [comets_scr],
        capture_output=True, text=True, timeout=15,
        cwd=COMETS_HOME
    )
    print(f"Return code: {result.returncode}")
    print(f"stdout: {result.stdout[:1000]}")
    print(f"stderr: {result.stderr[:1000]}")
except Exception as e:
    print(f"Test failed: {e}")
    traceback.print_exc()

# Set COMETS_HOME env var for current session
os.environ["COMETS_HOME"] = COMETS_HOME
print(f"\nCOMETS_HOME environment variable set to: {COMETS_HOME}")
print("NOTE: You may need to set COMETS_HOME as a system environment variable.")

log_fh.close()
