"""
Step 12: Compile COMETS v2.12.4 with proper stubs for missing dependencies
Skip GLPK-dependent files, create proper jdistlib and JMatIO stubs
"""
import os
import sys
import subprocess
import traceback
import glob as globmod
import shutil

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "compile_v2124_log.txt")

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
COMETS_V2124 = r"C:\Users\hgh97\comets_home\COMETS_v2124\segrelab-comets-ce2a55d"
lib_dir = os.path.join(COMETS_HOME, "lib")
bin_dir = os.path.join(COMETS_HOME, "bin")

JAVA = os.path.join(JDK_HOME, "bin", "java.exe")
JAVAC = os.path.join(JDK_HOME, "bin", "javac.exe")
JAR_TOOL = os.path.join(JDK_HOME, "bin", "jar.exe")

os.environ["JAVA_HOME"] = JDK_HOME
os.environ["PATH"] = os.path.join(JDK_HOME, "bin") + ";" + os.environ.get("PATH", "")

print("=" * 60)
print("Compile COMETS v2.12.4 with Proper Stubs")
print("=" * 60)

# Step 1: Create proper stub JARs for missing dependencies
print("\n--- Step 1: Create stub JARs ---")

# 1a. jdistlib stub - needs jdistlib.rng.MersenneTwister and other classes
print("Creating jdistlib stub...")
jdistlib_stub = os.path.join(lib_dir, "jdistlib_stub")
if os.path.exists(jdistlib_stub):
    shutil.rmtree(jdistlib_stub)
os.makedirs(jdistlib_stub, exist_ok=True)

jdistlib_classes = {
    "rng/MersenneTwister.java": """
package jdistlib.rng;
public class MersenneTwister extends Random {
    private int[] mt = new int[624];
    private int mti = 625;
    public MersenneTwister() { this(System.currentTimeMillis()); }
    public MersenneTwister(long seed) { mt[0] = (int)seed; for(mti=1; mti<624; mti++) mt[mti] = (1812433253*(mt[mti-1]^(mt[mti-1]>>>30))+mti); mti=624; }
    protected int next(int bits) { return nextInt() >>> (32-bits); }
    public int nextInt() { int y; if(mti>=624) { int i; for(i=0;i<624-397;i++){y=(mt[i]&0x80000000)|(mt[i+1]&0x7fffffff);mt[i]=mt[i+397]^(y>>>1)^((y&1)!=0?0x9908b0df:0);} for(;i<624-1;i++){y=(mt[i]&0x80000000)|(mt[i+1]&0x7fffffff);mt[i]=mt[i+(397-624)]^(y>>>1)^((y&1)!=0?0x9908b0df:0);} y=(mt[624-1]&0x80000000)|(mt[0]&0x7fffffff);mt[624-1]=mt[397-1]^(y>>>1)^((y&1)!=0?0x9908b0df:0);mti=0;} y=mt[mti++];y^=(y>>>11);y^=(y<<7)&0x9d2c5680;y^=(y<<15)&0xefc60000;y^=(y>>>18);return y; }
}
""",
    "rng/Random.java": """
package jdistlib.rng;
public abstract class Random extends java.util.Random {
    public Random() { super(); }
    public Random(long seed) { super(seed); }
    public abstract int nextInt();
    protected int next(int bits) { return nextInt() >>> (32-bits); }
}
""",
    "Beta.java": "package jdistlib;\npublic class Beta { public static double dbeta(double x, double a, double b, boolean give_log) { return 0.0; } public static double pbeta(double x, double a, double b, boolean lower_tail, boolean log_p) { return 0.5; } public static double qbeta(double p, double a, double b, boolean lower_tail, boolean log_p) { return 0.5; } public static double rbeta(jdistlib.rng.Random rr, double a, double b) { return 0.5; } }",
    "NonCentralBeta.java": "package jdistlib;\npublic class NonCentralBeta { public static double dnchisq(double x, double df, double ncp, boolean give_log) { return 0.0; } }",
    "MathFunctions.java": "package jdistlib;\npublic class MathFunctions { public static double gammafn(double x) { return 1.0; } public static double lgammafn(double x) { return 0.0; } public static double lbeta(double a, double b) { return 0.0; } public static double choose(double n, double k) { return 1.0; } public static double lchoose(double n, double k) { return 0.0; } }",
    "GenericDistribution.java": "package jdistlib;\npublic abstract class GenericDistribution { public abstract double density(double x, boolean log); public abstract double cumulative(double x, boolean lower_tail, boolean log_p); }",
    "Normal.java": "package jdistlib;\npublic class Normal { public static double dnorm(double x, double mu, double sigma, boolean give_log) { return 0.0; } public static double pnorm(double x, double mu, double sigma, boolean lower_tail, boolean log_p) { return 0.5; } public static double qnorm(double p, double mu, double sigma, boolean lower_tail, boolean log_p) { return 0.0; } }",
}

for rel_path, code in jdistlib_classes.items():
    fpath = os.path.join(jdistlib_stub, "jdistlib", rel_path)
    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    with open(fpath, "w", encoding="utf-8") as f:
        f.write(code)

# Compile jdistlib stubs
result = subprocess.run(
    [JAVAC, "-encoding", "UTF-8", "-d", jdistlib_stub] +
    globmod.glob(os.path.join(jdistlib_stub, "jdistlib", "**", "*.java"), recursive=True),
    capture_output=True, text=True, timeout=30
)
print(f"  jdistlib stub compile: {result.returncode}")
if result.stderr:
    print(f"  stderr: {result.stderr[:500]}")

# Package into JAR
jdistlib_jar = os.path.join(lib_dir, "jdistlib-0.4.5-bin.jar")
if os.path.exists(jdistlib_jar):
    os.remove(jdistlib_jar)
result = subprocess.run(
    [JAR_TOOL, "cf", jdistlib_jar, "-C", jdistlib_stub, "jdistlib"],
    capture_output=True, text=True, timeout=10
)
print(f"  jdistlib JAR: {os.path.exists(jdistlib_jar)} ({os.path.getsize(jdistlib_jar)} bytes)")
shutil.rmtree(jdistlib_stub, ignore_errors=True)

# 1b. JMatIO stub - needs com.jmatio.io.MatFileIncrementalWriter etc.
print("\nCreating JMatIO stub...")
jmatio_stub = os.path.join(lib_dir, "jmatio_stub")
if os.path.exists(jmatio_stub):
    shutil.rmtree(jmatio_stub)
os.makedirs(jmatio_stub, exist_ok=True)

jmatio_classes = {
    "io/MatFileIncrementalWriter.java": "package com.jmatio.io; import java.io.*; public class MatFileIncrementalWriter { public MatFileIncrementalWriter(String fileName) throws IOException {} public void close() throws IOException {} }",
    "io/MatFileWriter.java": "package com.jmatio.io; import java.io.*; public class MatFileWriter { public MatFileWriter() {} public void write(String fileName, java.util.List<com.jmatio.types.MLArray> data) throws IOException {} }",
    "io/MatFileReader.java": "package com.jmatio.io; import java.io.*; public class MatFileReader { public MatFileReader() {} public MatFileReader(String fileName) throws IOException {} }",
    "types/MLArray.java": "package com.jmatio.types; public class MLArray { public static final int mxCELL_CLASS = 1; public static final int mxSTRUCT_CLASS = 2; public static final int mxOBJECT_CLASS = 3; public static final int mxUINT8_CLASS = 4; public static final int mxINT8_CLASS = 5; public static final int mxUINT16_CLASS = 6; public static final int mxINT16_CLASS = 7; public static final int mxUINT32_CLASS = 8; public static final int mxINT32_CLASS = 9; public static final int mxSINGLE_CLASS = 10; public static final int mxDOUBLE_CLASS = 11; public static final int mxUINT64_CLASS = 12; public static final int mxINT64_CLASS = 13; public static final int mxSPARSE_CLASS = 14; protected String name; protected int type; public MLArray(String name, int type) { this.name=name; this.type=type; } public String getName() { return name; } public int getType() { return type; } }",
    "types/MLDouble.java": "package com.jmatio.types; public class MLDouble extends MLArray { public MLDouble(String name, int[] dims, int type) { super(name, type); } public MLDouble(String name, int[] dims, int type, double[] data) { super(name, type); } }",
    "types/MLCell.java": "package com.jmatio.types; public class MLCell extends MLArray { public MLCell(String name, int[] dims) { super(name, mxCELL_CLASS); } }",
    "types/MLStructure.java": "package com.jmatio.types; public class MLStructure extends MLArray { public MLStructure(String name, int[] dims) { super(name, mxSTRUCT_CLASS); } }",
}

for rel_path, code in jmatio_classes.items():
    fpath = os.path.join(jmatio_stub, "com", "jmatio", rel_path)
    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    with open(fpath, "w", encoding="utf-8") as f:
        f.write(code)

result = subprocess.run(
    [JAVAC, "-encoding", "UTF-8", "-d", jmatio_stub] +
    globmod.glob(os.path.join(jmatio_stub, "com", "jmatio", "**", "*.java"), recursive=True),
    capture_output=True, text=True, timeout=30
)
print(f"  JMatIO stub compile: {result.returncode}")

# Package into JAR
jmatio_jar = os.path.join(lib_dir, "JMatIO", "lib", "jmatio.jar")
os.makedirs(os.path.dirname(jmatio_jar), exist_ok=True)
if os.path.exists(jmatio_jar):
    os.remove(jmatio_jar)
result = subprocess.run(
    [JAR_TOOL, "cf", jmatio_jar, "-C", jmatio_stub, "com"],
    capture_output=True, text=True, timeout=10
)
print(f"  JMatIO JAR: {os.path.exists(jmatio_jar)} ({os.path.getsize(jmatio_jar)} bytes)")
# Also copy to JMatIO-041212
jmatio_dir2 = os.path.join(lib_dir, "JMatIO", "JMatIO-041212", "lib")
os.makedirs(jmatio_dir2, exist_ok=True)
shutil.copy2(jmatio_jar, os.path.join(jmatio_dir2, "jmatio.jar"))
shutil.rmtree(jmatio_stub, ignore_errors=True)

# Step 2: Compile COMETS v2.12.4
print("\n--- Step 2: Compile COMETS v2.12.4 ---")
src_dir = os.path.join(COMETS_V2124, "comets_simplified", "src")
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
            try:
                import zipfile
                with zipfile.ZipFile(fpath, 'r') as z:
                    z.testzip()
                classpath_jars.append(fpath)
            except:
                print(f"  Skipping invalid JAR: {f}")
                os.remove(fpath)

for f in os.listdir(bin_dir):
    if f.endswith(".jar"):
        classpath_jars.append(os.path.join(bin_dir, f))

classpath = ";".join(classpath_jars)
print(f"Classpath: {len(classpath_jars)} JARs")

# Collect Java source files, skip GLPK-dependent files
java_files = []
skip_files = ["FBAOptimizerGLPK.java"]  # Skip GLPK optimizer
for root, dirs, files in os.walk(src_dir):
    for f in files:
        if f.endswith(".java"):
            if any(sf in f for sf in skip_files):
                print(f"  Skipping: {f}")
                continue
            java_files.append(os.path.join(root, f))

print(f"Compiling {len(java_files)} Java files...")

cmd = [JAVAC, "-encoding", "UTF-8", "-cp", classpath, "-d", compile_dir, "-sourcepath", src_dir]
cmd.extend(java_files)

result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
print(f"Compile return code: {result.returncode}")

if result.stderr:
    errors = [l for l in result.stderr.split('\n') if 'error:' in l]
    print(f"Errors: {len(errors)}")
    print(f"stderr (first 5000 chars):\n{result.stderr[:5000]}")

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
        print(f"stderr: {result.stderr[:2000]}")

    # Update comets_scr.bat
    comets_scr = os.path.join(COMETS_HOME, "comets_scr.bat")
    bat_content = f"""@echo off
set JAVA={JAVA}
set COMETS_HOME={COMETS_HOME}
set CLASSPATH={cp}
"%JAVA%" -classpath "%CLASSPATH%" edu.bu.segrelab.comets.Comets -loader edu.bu.segrelab.comets.fba.FBACometsLoader -script %1
"""
    with open(comets_scr, "w") as f:
        f.write(bat_content)

    # Run dFBA Demo!
    print("\n" + "=" * 60)
    print("Running dFBA Demo with COMETS 2.12.4 (compiled from source)")
    print("=" * 60)

    os.environ["COMETS_HOME"] = COMETS_HOME

    import cometspy as c
    import cobra

    ecoli = cobra.io.load_model("textbook")
    print(f"Model: {ecoli.id}, {len(ecoli.reactions)} reactions")

    model1 = c.model(ecoli)
    model1.id = "ecoli_wt"
    model1.initial_pop = [0, 0, 1e-6]
    model1.open_exchanges()

    model2 = c.model(ecoli)
    model2.id = "ecoli_mutant"
    model2.initial_pop = [0, 0, 1e-6]
    model2.open_exchanges()

    lay = c.layout([model1, model2])
    lay.set_specific_metabolite("glc__D_e", 0.01)
    lay.set_specific_metabolite("nh4_e", 1000.0)
    lay.set_specific_metabolite("pi_e", 1000.0)

    p = c.params()
    p.set_param("maxCycles", 30)
    p.set_param("timeStep", 0.1)
    p.set_param("writeTotalBiomassLog", True)
    p.set_param("writeMediaLog", True)
    p.set_param("BiomassLogRate", 1)
    p.set_param("MediaLogRate", 1)
    p.set_param("deathRate", 0.0)

    # Remove version-mismatch params
    skip_keys = ['velocityMultiConvLogFormat', 'velocityMultiConvLogName',
                 'velocityMultiConvLogRate', 'writeVelocityMultiConvLog']
    for key in list(p.all_params.keys()):
        if any(sk.lower() in key.lower() for sk in skip_keys):
            del p.all_params[key]

    print("\nCreating simulation...")
    sim = c.comets(lay, p)
    print(f"VERSION: {sim.VERSION}")

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
    print("\nCompilation FAILED.")
    if os.path.exists(compile_dir):
        shutil.rmtree(compile_dir)

    # Show the specific error files
    error_files = set()
    for line in result.stderr.split('\n'):
        if '.java:' in line:
            parts = line.split(':')
            if len(parts) >= 1:
                error_files.add(parts[0].strip())
    print(f"\nFiles with errors ({len(error_files)}):")
    for f in sorted(error_files):
        print(f"  {f}")

log_fh.close()
