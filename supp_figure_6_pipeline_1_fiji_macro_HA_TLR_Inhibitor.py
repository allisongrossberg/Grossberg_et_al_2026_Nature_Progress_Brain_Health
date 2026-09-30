# =============================================================================
# Supplementary Figure 6 pipeline, step 1 - FIJI/Jython macro (HA TLR
#   Inhibitor experiment): raw .oir microscopy images -> DAPI count, GFAP/
#   Vimentin MFI, and pre-watershed segmentation masks
# =============================================================================
# PATHS: `input`/`output` below are hardcoded to this project's local
#   filesystem layout (raw images and per-image outputs are far too large to
#   bundle with this code deposit). Edit both before running.
#
# Segmentation uses DAPI-seeded (nucleus-anchored) watershed: each DAPI
# nucleus is a seed, and GFAP/Vimentin masks are split at those same seed
# positions downstream (skimage.segmentation.watershed(markers=...) in step
# 3), so cell identity is shared across both channels. This macro handles
# preprocessing/thresholding only and deliberately skips FIJI's own blind
# Watershed step on the area/branch masks -- the actual per-cell split
# happens in step 3, not here.
#
# Output feeds step 2 (branch Local Thickness map) and step 3 (final per-cell
# morphology CSV).
#
# Usage:
#   /Applications/Fiji/Fiji.app/Contents/MacOS/fiji-macos --run <this file>
#   (NOT --headless -- Bio-Formats does not work under ImageJ's --headless
#   flag even in macro mode; batchMode=True below suppresses visible windows
#   instead)
# =============================================================================

import os

from ij import IJ, WindowManager, macro, CompositeImage
from ij.plugin.frame import RoiManager

macro.Interpreter.batchMode = True  # headless-safe: suppress image display

# EDIT: point these at your own raw-image and output directories
input = "/Volumes/One Touch/Image Analysis Folders/Image_Analysis_Input_7_18_24_COAST_HA_40x_TLR_Exp_1_3_Comb_ADE_ReDo/"
output = "/Volumes/One Touch/COAST_Paper_Review_Edits_Image_Analysis/HA_TLR_Inhibitor/"
um_to_pix = "0.311"

all_files = os.listdir(input)
files = [f for f in all_files if f.endswith(".oir") and not f.startswith(".")]  # exclude macOS AppleDouble ._ metadata files
files = sorted(files)

# windowless=true avoids a disambiguation dialog that would otherwise crash
# headless runs (no display to draw it on).
DEFAULT_OPTIONS = "color_mode=Default view=Hyperstack windowless=true"
COMPOSITE_OPTIONS = "color_mode=Composite view=Hyperstack windowless=true"


def process_one_file(file):
    print("Now Processing File: {0}".format(file))
    IJ.run("Bio-Formats Importer", "open=[" + input + file + "]" + DEFAULT_OPTIONS)

    originalFileName = os.path.basename(file)
    originalFileNameWithoutExtension = os.path.splitext(originalFileName)[0]
    originalFilePath = input + originalFileName

    # ---- DAPI count ----
    IJ.run("Set Scale...", "distance=1 known={0} unit=micron global".format(um_to_pix))

    IJ.selectWindow(originalFileName)
    imp = IJ.getImage()
    n_channels = imp.getNChannels()
    if n_channels == 4:
        IJ.run("Arrange Channels...", "new=1243")
    IJ.run("Split Channels")

    blueChannelName = "C1-" + originalFileName
    greenChannelName = "C2-" + originalFileName
    redChannelName = "C3-" + originalFileName

    IJ.selectWindow(blueChannelName)
    duplicateTitle = "duplicate_" + blueChannelName
    IJ.run("Duplicate...", "title=" + duplicateTitle)

    IJ.selectWindow(blueChannelName)
    imp = IJ.getImage()
    IJ.run("Gaussian Blur...", "sigma=2")
    IJ.run("Gamma...", "value=0.5")
    IJ.run("Max...", "value=1500")
    IJ.setAutoThreshold(imp, "Default dark no-reset")
    IJ.runMacro("""setOption("BlackBackground", true)""")
    IJ.run("Convert to Mask")
    IJ.run("Watershed")

    IJ.run("Set Measurements...", "area mean standard min max perimeter shape limit redirect=[None] decimal=2")
    IJ.run("Analyze Particles...", "size=35.00-Infinity circularity=0.50-1.00 show=Outlines display summarize overlay composite add")
    IJ.run("Measure")

    rm = RoiManager.getInstance()
    nROIs = rm.getCount()
    indexes = range(0, nROIs + 1)
    IJ.selectWindow("Results")
    for i in indexes:
        IJ.runMacro('setResult("Image File", {0}, "{1}")'.format(i, blueChannelName))

    IJ.run("Bio-Formats Importer", "open=[" + input + file + "]" + COMPOSITE_OPTIONS)
    IJ.selectWindow(originalFileName)
    imp = IJ.getImage()
    rm = RoiManager.getInstance()
    IJ.run("Colors...", "channels=1 slices")
    rm.runCommand(imp, "Show All")
    IJ.run("Make Composite")
    if n_channels == 4:
        IJ.runMacro("""Stack.toggleChannel(3)""")
    IJ.saveAs(imp, "Tiff", output + originalFileNameWithoutExtension + "-DAPI_ROI_Composite.tif")

    IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-DAPI_Count_Results.csv")
    IJ.selectWindow("Summary")
    IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-DAPI_Count_Summary.csv")
    IJ.run("Clear Results")

    if nROIs > 0:
        rm.setSelectedIndexes(indexes)
        IJ.runMacro("""roiManager("Add")""")
        rm.runCommand("Save", output + originalFileNameWithoutExtension + "-DAPI_RoiSet.zip")

    IJ.selectWindow(blueChannelName)
    IJ.saveAs(imp, "Tiff", output + originalFileName)

    IJ.run("Set Measurements...", "area mean standard min max perimeter shape limit redirect=[duplicate_C1-" + originalFileName + "] decimal=2")
    IJ.run("Analyze Particles...", "size=35.00-Infinity circularity=0.50-1.00 show=Outlines display summarize overlay composite add")
    IJ.run("Measure")
    IJ.selectWindow("Results")
    for i in indexes:
        IJ.runMacro('setResult("Image File", {0}, "duplicate_C1" + "{1}" + "")'.format(i, originalFileName))
    IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-DAPI_Duplicate_Image_Results.csv")
    IJ.run("Clear Results")

    if n_channels == 4:
        nfkbChannelName = "C4-" + originalFileName
        IJ.selectWindow(blueChannelName)
        IJ.run("Set Measurements...", "mean standard min max limit redirect=[C4-" + originalFileName + "] decimal=2")
        IJ.run("Analyze Particles...", "size=35.00-Infinity circularity=0.50-1.00 show=Outlines display summarize overlay composite add")
        IJ.run("Measure")
        IJ.selectWindow("Results")
        for i in indexes:
            IJ.runMacro('setResult("Image File", {0}, "C4-" + "{1}" + "")'.format(i, originalFileName))
        IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-NFkB_Results.csv")
        IJ.run("Clear Results")
        IJ.run("Bio-Formats Importer", "open=[" + input + file + "]" + COMPOSITE_OPTIONS)
        IJ.selectWindow(originalFileName)
        imp = IJ.getImage()
        rm = RoiManager.getInstance()
        IJ.run("Colors...", "channels=1 slices")
        rm.runCommand(imp, "Show All")
        IJ.run("Make Composite")
        IJ.saveAs("Tiff", output + originalFileNameWithoutExtension + "-NFkB_Composite_Image_With_ROI.tif")

    rm.reset()

    # ---- GFAP channel: 6 duplicates ----
    IJ.selectWindow(greenChannelName)
    duplicateTitle = "duplicate_" + greenChannelName
    IJ.run("Duplicate...", "title=" + duplicateTitle)
    IJ.selectWindow(greenChannelName)
    duplicateTitle2 = "duplicate_2_" + greenChannelName
    IJ.run("Duplicate...", "title=" + duplicateTitle2)
    IJ.selectWindow(greenChannelName)
    duplicateTitle3 = "duplicate_3_" + greenChannelName
    IJ.run("Duplicate...", "title=" + duplicateTitle3)
    duplicateTitle4 = "duplicate_4_" + greenChannelName
    IJ.run("Duplicate...", "title=" + duplicateTitle4)
    duplicateTitle5 = "duplicate_5_" + greenChannelName
    IJ.run("Duplicate...", "title=" + duplicateTitle5)
    duplicateTitle6 = "duplicate_6_" + greenChannelName
    IJ.run("Duplicate...", "title=" + duplicateTitle6)

    # ---- Vimentin channel: same 6-duplicate pattern ----
    IJ.selectWindow(redChannelName)
    rDupTitle3 = "duplicate_3_" + redChannelName
    IJ.run("Duplicate...", "title=" + rDupTitle3)
    IJ.selectWindow(redChannelName)
    rDupTitle4 = "duplicate_4_" + redChannelName
    IJ.run("Duplicate...", "title=" + rDupTitle4)
    IJ.selectWindow(redChannelName)
    rDupTitle5 = "duplicate_5_" + redChannelName
    IJ.run("Duplicate...", "title=" + rDupTitle5)
    IJ.selectWindow(redChannelName)
    rDupTitle6 = "duplicate_6_" + redChannelName
    IJ.run("Duplicate...", "title=" + rDupTitle6)

    # ---- GFAP MFI ----
    IJ.selectWindow(duplicateTitle)
    IJ.run("Enhance Contrast...", "saturated=1 normalize")
    IJ.run("Gaussian Blur...", "sigma=2")
    IJ.run("Subtract Background...", "rolling=50")
    IJ.setAutoThreshold(imp, "Default dark no-reset")
    IJ.runMacro("""setOption("BlackBackground", true)""")
    IJ.run("Convert to Mask")
    IJ.run("Set Measurements...", "area mean standard min max limit redirect=[duplicate_2_C2-" + originalFileName + "] decimal=2")
    IJ.run("Measure")
    IJ.selectWindow("Results")
    IJ.runMacro('setResult("Image File", {0}, "duplicate_2_C2-" + "{1}" + "")'.format(0, originalFileName))
    IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-GFAP_MFI_Results.csv")
    IJ.run("Clear Results")

    # ---- Vimentin MFI ----
    IJ.selectWindow(redChannelName)
    duplicateTitleRed = "duplicate_" + redChannelName
    IJ.run("Duplicate...", "title=" + duplicateTitleRed)
    IJ.selectWindow(redChannelName)
    IJ.run("Enhance Contrast...", "saturated=1 normalize")
    IJ.run("Gaussian Blur...", "sigma=2")
    IJ.run("Subtract Background...", "rolling=50")
    IJ.setAutoThreshold(imp, "Default dark no-reset")
    IJ.runMacro("""setOption("BlackBackground", true)""")
    IJ.run("Convert to Mask")
    IJ.run("Set Measurements...", "area mean standard min max limit redirect=[duplicate_C3-" + originalFileName + "] decimal=2")
    IJ.run("Measure")
    IJ.selectWindow("Results")
    IJ.runMacro('setResult("Image File", {0}, "duplicate_C3-" + "{1}" + "")'.format(0, originalFileName))
    IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-Vimentin_MFI_Results.csv")
    IJ.run("Clear Results")

    # Cell area branch: identical algorithm for GFAP and Vimentin.
    def run_area_branch(dupTitle3, dupTitle4, channelName, channelLabel, maskTag):
        IJ.selectWindow(dupTitle3)
        cimp = IJ.getImage()
        IJ.run("Enhance Contrast...", "saturated=0 normalize")
        IJ.run("Gaussian Blur...", "sigma=8")
        IJ.run("Subtract Background...", "rolling=52")
        IJ.run("Bandpass Filter...", "filter_large=85 filter_small=20 suppress=None tolerance=5 autoscale saturate")
        IJ.run("Despeckle")
        IJ.run("Maximum...", "radius=3")
        IJ.run("Remove Outliers...", "radius=2 threshold=50 which=Bright")
        IJ.run("Gaussian Blur...", "sigma=1")
        IJ.setAutoThreshold(cimp, "MaxEntropy dark no-reset")
        IJ.runMacro("""setOption("BlackBackground", true)""")
        IJ.run("Convert to Mask")
        IJ.run("Despeckle")
        IJ.run("Close-")
        IJ.run("Dilate")
        # No blind Watershed here -- touching cells stay merged; step 3
        # performs the DAPI-seeded split on this exact pre-split mask.
        IJ.run("Fill Holes")
        IJ.run("Dilate")
        IJ.run("Dilate")
        IJ.run("Close-")
        IJ.run("Despeckle")

        maskSaveTitle = maskTag + "_premask"
        IJ.run("Duplicate...", "title=" + maskSaveTitle)
        IJ.saveAs("Tiff", output + originalFileNameWithoutExtension + "-" + maskTag + "_premask.tif")

        # Local Thickness plugin, run once on the pre-split mask (thickness
        # is a property of the mask's own geometry, independent of the
        # later per-cell split). Detect the new window by diffing titles
        # before/after rather than assuming its exact name.
        titlesBefore = set(WindowManager.getImageTitles())
        IJ.run("Local Thickness (masked, calibrated, silent)", "")
        titlesAfter = set(WindowManager.getImageTitles())
        newTitles = list(titlesAfter - titlesBefore)
        if len(newTitles) != 1:
            raise RuntimeError("Expected exactly one new window from Local Thickness, got: " + str(newTitles))
        IJ.selectWindow(newTitles[0])
        IJ.saveAs("Tiff", output + originalFileNameWithoutExtension + "-" + maskTag + "_LocalThickness.tif")
        IJ.run("Close")

        IJ.selectWindow(dupTitle3)

        IJ.run("Set Measurements...", "area mean standard min perimeter shape integrated skewness limit redirect=[None] decimal=2")
        IJ.run("Analyze Particles...", "size=85.00-Infinity circularity=0.00-1.00 show=Outlines display summarize overlay add composite")
        IJ.run("Measure")

        rmLocal = RoiManager.getInstance()
        nR = rmLocal.getCount()
        idx = range(0, nR + 1)
        IJ.selectWindow("Results")
        for i in idx:
            IJ.runMacro('setResult("Image File", {0}, "{1}")'.format(i, channelName))
        IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-" + channelLabel + "_Cell_Count_Results.csv")
        IJ.selectWindow("Summary")
        IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-" + channelLabel + "_Count_Summary.csv")
        IJ.run("Clear Results")

        if nR > 0:
            rmLocal.setSelectedIndexes(idx)
            IJ.runMacro("""roiManager("Add")""")
            rmLocal.runCommand("Save", output + originalFileNameWithoutExtension + "-" + channelLabel + "_RoiSet.zip")

        IJ.selectWindow(dupTitle3)
        IJ.saveAs("Tiff", output + originalFileNameWithoutExtension + "-" + channelLabel + "Channel_with_ROI.tif")

        IJ.run("Set Measurements...", "area mean standard min max perimeter shape limit redirect=[" + dupTitle4 + "] decimal=2")
        IJ.run("Analyze Particles...", "size=85.00-Infinity circularity=0.00-1.00 show=Outlines display summarize overlay add composite")
        IJ.run("Measure")
        if nR == 1:
            IJ.runMacro('setResult("Area", 0, 0)')
            IJ.runMacro('setResult("Perim.", 0, 0)')
            IJ.runMacro('setResult("Circ.", 0, 0)')
            IJ.runMacro('setResult("IntDen", 0, 0)')
            IJ.runMacro('setResult("Skew", 0, 0)')
            IJ.runMacro('setResult("RawIntDen", 0, 0)')
            IJ.runMacro('setResult("AR", 0, 0)')
            IJ.runMacro('setResult("Round", 0, 0)')
            IJ.runMacro('setResult("Solidity", 0, 0)')
            IJ.runMacro('setResult("MinThr", 0, 0)')
            IJ.runMacro('setResult("MaxThr", 0, 0)')
        IJ.selectWindow("Results")
        for i in idx:
            IJ.runMacro('setResult("Image File", {0}, "{1}" + "{2}"+ "")'.format(i, dupTitle4, originalFileName))
        IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-" + channelLabel + "_Duplicate_Image_Results.csv")
        IJ.run("Clear Results")
        rmLocal.reset()
        return nR

    # Cell branching branch: identical algorithm for GFAP and Vimentin.
    def run_branching_branch(dupTitle5, channelLabel):
        IJ.selectWindow(dupTitle5)
        bimp = IJ.getImage()
        IJ.run("Bandpass Filter...", "filter_large=85 filter_small=8 suppress=None tolerance=5 autoscale saturate")
        IJ.run("Despeckle")
        IJ.run("Sharpen")
        IJ.run("Grays")
        IJ.run("Unsharp Mask...", "radius=35 mask=0.75")
        IJ.run("Despeckle")
        IJ.run("Remove Outliers...", "radius=2 threshold=50 which=Bright")
        IJ.setAutoThreshold(bimp, "MaxEntropy dark no-reset")
        IJ.runMacro("""setOption("BlackBackground", true)""")
        IJ.run("Convert to Mask")
        IJ.run("Despeckle")
        IJ.run("Close-")
        IJ.run("Remove Outliers...", "radius=3 threshold=50 which=Bright")

        # Lighter-touch preprocessing than the area branch (no heavy blur,
        # no post-threshold dilation) to preserve fine process structure;
        # step 3 reads this mask for shape/skeleton/thickness.
        branchMaskSaveTitle = channelLabel + "_branch_premask"
        IJ.run("Duplicate...", "title=" + branchMaskSaveTitle)
        IJ.saveAs("Tiff", output + originalFileNameWithoutExtension + "-" + channelLabel + "_branch_premask.tif")
        IJ.selectWindow(dupTitle5)

        IJ.run("Skeletonize")
        IJ.run("Analyze Skeleton (2D/3D)", "prune=none show display")

        IJ.selectWindow("Results")
        IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-" + channelLabel + "_Skeleton_Results.csv")
        IJ.selectWindow("Branch information")
        IJ.saveAs("Results", output + originalFileNameWithoutExtension + "-" + channelLabel + "_Skeleton_Branch_Info.csv")
        IJ.selectWindow("Tagged skeleton")
        IJ.saveAs("Tiff", output + originalFileNameWithoutExtension + "-" + channelLabel + "_Tagged_Skeletons.tif")
        IJ.run("Clear Results")

    # GFAP
    nROIs_gfap = run_area_branch(duplicateTitle3, duplicateTitle4, greenChannelName, "GFAP", "GFAP")

    IJ.run("Bio-Formats Importer", "open=[" + input + file + "]" + COMPOSITE_OPTIONS)
    IJ.selectWindow(originalFileName)
    imp = IJ.getImage()
    rm = RoiManager.getInstance()
    IJ.run("Colors...", "channels=1 slices")
    rm.runCommand(imp, "Show All")
    IJ.run("Make Composite")
    if n_channels == 4:
        IJ.runMacro("""Stack.toggleChannel(3)""")
    IJ.saveAs("Tiff", output + originalFileNameWithoutExtension + "-GFAP_Composite_Image_With_ROI.tif")
    rm.reset()

    run_branching_branch(duplicateTitle5, "GFAP")

    # Vimentin
    nROIs_vim = run_area_branch(rDupTitle3, rDupTitle4, redChannelName, "Vimentin", "Vimentin")

    IJ.run("Bio-Formats Importer", "open=[" + input + file + "]" + COMPOSITE_OPTIONS)
    IJ.selectWindow(originalFileName)
    imp = IJ.getImage()
    rm = RoiManager.getInstance()
    IJ.run("Colors...", "channels=1 slices")
    rm.runCommand(imp, "Show All")
    IJ.run("Make Composite")
    if n_channels == 4:
        IJ.runMacro("""Stack.toggleChannel(3)""")
    IJ.saveAs("Tiff", output + originalFileNameWithoutExtension + "-Vimentin_Composite_Image_With_ROI.tif")
    rm.reset()

    run_branching_branch(rDupTitle5, "Vimentin")

    IJ.run("Clear Results")
    IJ.run("Close All")
    window_list = WindowManager.getNonImageWindows()
    for window in window_list:
        IJ.selectWindow(window.title)
        IJ.run("Close")


failed_files = []
for file in files:
    try:
        process_one_file(file)
    except Exception as e:
        print("FAILED on file: {0} | error: {1}".format(file, str(e)))
        failed_files.append(file + " | " + str(e))
        # best-effort cleanup so a failure mid-file doesn't cascade
        try:
            IJ.run("Clear Results")
        except Exception:
            pass
        try:
            IJ.run("Close All")
        except Exception:
            pass

if failed_files:
    with open(output + "failed_files.txt", "w") as f:
        for item in failed_files:
            f.write("%s\n" % item)

print("Done. Processed {0} file(s), {1} failed.".format(len(files), len(failed_files)))
