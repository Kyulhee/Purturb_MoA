"""
Step 11: Download the pre-built COMETS Windows distribution
From the .classpath we saw references to comets_windows/comets_2.12.4
Try to find and download this distribution
"""
import os
import sys
import subprocess
import urllib.request
import json
import traceback
import shutil

LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "download_dist_log.txt")

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

os.environ["JAVA_HOME"] = JDK_HOME
os.environ["PATH"] = os.path.join(JDK_HOME, "bin") + ";" + os.environ.get("PATH", "")

print("=" * 60)
print("Download Pre-built COMETS Windows Distribution")
print("=" * 60)

# Check GitHub releases more thoroughly
print("\n--- Checking GitHub releases for distribution files ---")
try:
    url = 'https://api.github.com/repos/segrelab/COMETS/releases'
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Python')
    with urllib.request.urlopen(req, timeout=15) as resp:
        releases = json.loads(resp.read().decode())
        for r in releases[:5]:
            tag = r["tag_name"]
            assets = r.get('assets', [])
            print(f'\nRelease {tag}:')
            print(f'  tarball_url: {r.get("tarball_url", "N/A")}')
            print(f'  zipball_url: {r.get("zipball_url", "N/A")}')
            for a in assets:
                print(f'  Asset: {a["name"]} ({a["size"]} bytes) - {a["browser_download_url"]}')

            # If no assets, download zipball
            if not assets:
                zipball_url = r.get("zipball_url", "")
                if zipball_url:
                    print(f'  Downloading zipball from {zipball_url}...')
                    target = os.path.join(os.path.dirname(LOG_FILE), f"comets_{tag}.zip")
                    try:
                        urllib.request.urlretrieve(zipball_url, target)
                        print(f'  Downloaded: {os.path.getsize(target)} bytes')
                    except Exception as e:
                        print(f'  Failed: {e}')

except Exception as e:
    print(f"GitHub API failed: {e}")

# Try the v2.12.4 zipball - it should contain the source
print("\n--- Downloading COMETS v2.12.4 source zipball ---")
v2124_zip = os.path.join(os.path.dirname(LOG_FILE), "comets_v2.12.4.zip")
v2124_dir = r"C:\Users\hgh97\comets_home\COMETS_v2124"

if not os.path.exists(v2124_dir):
    try:
        zipball_url = "https://api.github.com/repos/segrelab/COMETS/zipball/v2.12.4"
        print(f"Downloading from {zipball_url}")
        urllib.request.urlretrieve(zipball_url, v2124_zip)
        print(f"Downloaded: {os.path.getsize(v2124_zip)} bytes")

        import zipfile
        with zipfile.ZipFile(v2124_zip, 'r') as z:
            # The zipball contains a top-level directory like segrelab-COMETS-xxxxx
            names = z.namelist()
            top_dir = names[0].split('/')[0]
            print(f"Top directory in zip: {top_dir}")

            # Check if this version has a pre-built distribution
            for name in names:
                if 'comets_' in name.lower() and name.endswith('.jar'):
                    print(f"  Found JAR: {name}")
                if 'lib/' in name.lower() and name.endswith('.jar'):
                    print(f"  Found lib JAR: {name}")

            # Extract
            print("Extracting...")
            z.extractall(v2124_dir)
            print("Done!")

        # Clean up
        os.remove(v2124_zip)

    except Exception as e:
        print(f"Failed: {e}")
        traceback.print_exc()
else:
    print(f"Already extracted at {v2124_dir}")

# Check what's in the extracted directory
if os.path.exists(v2124_dir):
    print("\n--- v2.12.4 source structure ---")
    for item in os.listdir(v2124_dir):
        print(f"  {item}")

    # Check the inner directory
    inner_dirs = [d for d in os.listdir(v2124_dir) if os.path.isdir(os.path.join(v2124_dir, d))]
    if inner_dirs:
        inner = os.path.join(v2124_dir, inner_dirs[0])
        print(f"\n  Inner dir: {inner_dirs[0]}")
        for item in os.listdir(inner):
            print(f"    {item}")

        # Check if comets_simplified has more files
        simple_src = os.path.join(inner, "comets_simplified", "src")
        if os.path.exists(simple_src):
            # Count files
            java_count = 0
            for root, dirs, files in os.walk(simple_src):
                for f in files:
                    if f.endswith('.java'):
                        java_count += 1
            print(f"    Java source files: {java_count}")

# Now try a different approach: use the runcomets.org download
print("\n--- Trying runcomets.org download ---")
try:
    url = "https://www.runcomets.org/get-started"
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Mozilla/5.0')
    with urllib.request.urlopen(req, timeout=15) as resp:
        content = resp.read().decode('utf-8', errors='ignore')
        # Find download links
        import re
        links = re.findall(r'href=["\']([^"\']+)["\']', content)
        download_links = [l for l in links if any(kw in l.lower() for kw in ['download', 'zip', 'comets', 'windows', 'jar'])]
        print(f"Potential download links ({len(download_links)}):")
        for l in download_links[:20]:
            print(f"  {l}")

        # Also search for direct file references
        file_refs = re.findall(r'(https?://[^\s"\'<>]+(?:zip|jar|tar\.gz))', content, re.I)
        print(f"\nDirect file references ({len(file_refs)}):")
        for r in file_refs[:20]:
            print(f"  {r}")
except Exception as e:
    print(f"Failed: {e}")

# Try the direct download from runcomets.org
print("\n--- Trying direct download links ---")
download_attempts = [
    ("https://www.runcomets.org/s/COMETS_Windows_2.12.4.zip", "COMETS_Windows_2.12.4.zip"),
    ("https://www.runcomets.org/s/comets_2.12.4.zip", "comets_2.12.4.zip"),
    ("https://www.runcomets.org/s/comets_windows.zip", "comets_windows.zip"),
]

for url, filename in download_attempts:
    target = os.path.join(r"C:\Users\hgh97", filename)
    try:
        print(f"  Trying: {url}")
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Mozilla/5.0')
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = resp.read()
            with open(target, 'wb') as f:
                f.write(data)
            print(f"  Downloaded: {os.path.getsize(target)} bytes")

            # Check if it's a valid zip
            try:
                import zipfile
                with zipfile.ZipFile(target, 'r') as z:
                    names = z.namelist()
                    print(f"  Valid zip with {len(names)} entries")
                    # Show structure
                    for n in names[:30]:
                        print(f"    {n}")
                    if len(names) > 30:
                        print(f"    ... and {len(names)-30} more")

                    # Extract to COMETS_HOME
                    print("  Extracting to COMETS_HOME...")
                    z.extractall(COMETS_HOME)
                    print("  Done!")
                break
            except zipfile.BadZipFile:
                print("  Not a valid zip file")
                os.remove(target)

    except Exception as e:
        print(f"  Failed: {e}")
        if os.path.exists(target):
            os.remove(target)

# Final check: does COMETS_HOME have the proper structure now?
print("\n--- Final COMETS_HOME check ---")
bin_dir = os.path.join(COMETS_HOME, "bin")
lib_dir = os.path.join(COMETS_HOME, "lib")
comets_scr = os.path.join(COMETS_HOME, "comets_scr.bat")

print(f"bin/ exists: {os.path.exists(bin_dir)}")
if os.path.exists(bin_dir):
    for f in os.listdir(bin_dir):
        print(f"  {f}")
print(f"lib/ exists: {os.path.exists(lib_dir)}")
if os.path.exists(lib_dir):
    for item in os.listdir(lib_dir):
        print(f"  {item}/")
print(f"comets_scr.bat exists: {os.path.exists(comets_scr)}")

log_fh.close()
