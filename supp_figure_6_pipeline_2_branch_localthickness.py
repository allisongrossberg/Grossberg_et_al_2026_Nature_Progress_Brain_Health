# =============================================================================
# Supplementary Figure 6 pipeline, step 2 - FIJI/Jython macro (HA TLR
#   Inhibitor experiment): Local Thickness map for the branch (process-
#   preserving) masks saved by step 1
# =============================================================================
# PATHS: `input_dir` below must match step 1's `output` directory. Edit it
#   to point at your own output directory before running.
#
# Runs the real Local Thickness plugin on each already-saved
# "_branch_premask.tif" mask from step 1. Does not touch raw images or
# re-run Bio-Formats.
#
# Output: <image>-<channel>_branch_LocalThickness.tif, same directory.
# =============================================================================

import os
from ij import IJ, WindowManager, macro

macro.Interpreter.batchMode = True

# EDIT: must match step 1's `output` directory
input_dir = "/Volumes/One Touch/COAST_Paper_Review_Edits_Image_Analysis/HA_TLR_Inhibitor/"
um_to_pix = "0.311"

all_files = os.listdir(input_dir)
branch_files = sorted([
    f for f in all_files
    if f.endswith("_branch_premask.tif") and not f.startswith(".")
])
print("Found {0} branch_premask files to process".format(len(branch_files)))

failed = []
for i, fname in enumerate(branch_files):
    out_name = fname.replace("_branch_premask.tif", "_branch_LocalThickness.tif")
    if os.path.exists(input_dir + out_name):
        continue
    try:
        IJ.open(input_dir + fname)
        imp = IJ.getImage()
        IJ.run(imp, "Set Scale...", "distance=1 known={0} unit=micron global".format(um_to_pix))
        titlesBefore = set(WindowManager.getImageTitles())
        IJ.run("Local Thickness (masked, calibrated, silent)", "")
        titlesAfter = set(WindowManager.getImageTitles())
        newTitles = list(titlesAfter - titlesBefore)
        if len(newTitles) != 1:
            raise RuntimeError("Expected exactly one new window, got: " + str(newTitles))
        IJ.selectWindow(newTitles[0])
        IJ.saveAs("Tiff", input_dir + out_name)
    except Exception as e:
        print("FAILED: " + fname + " | " + str(e))
        failed.append(fname)
    finally:
        try:
            IJ.run("Close All")
        except Exception:
            pass
    if i % 200 == 0:
        print("Progress: {0}/{1}".format(i, len(branch_files)))

if failed:
    with open(input_dir + "branch_localthickness_failed.txt", "w") as f:
        for x in failed:
            f.write(x + "\n")

print("Done. Processed {0} file(s), {1} failed.".format(len(branch_files), len(failed)))
