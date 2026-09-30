"""
Supplementary Figure 6 pipeline, step 3 - measure_morphology_dapiseeded_v2.py
===============================================================================
Standalone Python script (not run inside FIJI): reads step 1's masks and
step 2's branch Local Thickness maps, performs its own DAPI-seeded watershed
split, and writes the final per-cell morphology CSV
(morphology_v2_dapiseeded_results.csv) that supp_figure_6b.R reads directly.

Cell-body area (hypertrophy) is measured from the area-blob mask
("*-<ch>_premask.tif"); shape (circularity, solidity, aspect ratio),
skeleton length, branch count, and thickness are measured from the
process-preserving branch mask ("*-<ch>_branch_premask.tif" +
"*-<ch>_branch_LocalThickness.tif"). Both masks are split with the same
DAPI-seeded watershed, so a cell's area and shape/branching share one
identity (dapi_seed_label) and land in the same output row.

Usage:
    python measure_morphology_dapiseeded_v2.py <macro_out_dir> <qc_dir> \
        [pixel_size_um] [min_object_px] [name_filter_substring] [max_images]

    name_filter_substring / max_images are OPTIONAL and only for quick
    validation runs on a subset; omit both for the full production run.
"""
import sys
import os
import glob
import numpy as np
import pandas as pd
import tifffile
import roifile
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.morphology import skeletonize, remove_small_objects
from skan import Skeleton
from skimage.measure import regionprops
from skimage.draw import polygon
import warnings
warnings.filterwarnings("ignore")

# Defaults match HA (0.311 um/px, 85 px^2 cutoff); RA passes 0.207 / 30.
PIXEL_SIZE_UM = 0.310740284701119
PROCESS_DIAMETER_CUTOFF_UM = 3.0
MIN_OBJECT_PX = 85
SPECKLE_MIN_PX = 8  # only removes 1-7 px isolated specks; real process tips survive


def erosion_radius_px(pixel_size_um):
    return max(1, int(round((PROCESS_DIAMETER_CUTOFF_UM / 2.0) / pixel_size_um)))


def true_skeleton_length_um(skel, pixel_size_um):
    if skel.sum() < 2:
        return 0.0
    try:
        return float(Skeleton(skel, spacing=pixel_size_um).path_lengths().sum())
    except Exception:
        return float(skel.sum()) * pixel_size_um


def skeleton_branch_count(skel):
    """Per-cell branch count = number of skeleton branches (skan paths). This
    is the 'branching' axis: diffuse/resting astrocytes have many fine
    branches, hypertrophic reactive cells retract them (fewer, thicker)."""
    if skel.sum() < 2:
        return 0
    try:
        return int(Skeleton(skel).n_paths)
    except Exception:
        return 0


def dapi_rois_to_markers(roi_path, shape):
    if not os.path.exists(roi_path):
        return np.zeros(shape, dtype=np.int32)
    rois = roifile.roiread(roi_path)
    if not isinstance(rois, list):
        rois = [rois]
    markers = np.zeros(shape, dtype=np.int32)
    next_label = 1
    for roi in rois:
        coords = roi.coordinates()
        if coords is None or len(coords) == 0:
            continue
        rr, cc = polygon(coords[:, 1], coords[:, 0], shape=shape)
        if len(rr) == 0:
            continue
        markers[rr, cc] = next_label
        next_label += 1
    return markers


def _split(mask, dapi_markers):
    """DAPI-seeded watershed split of a boolean mask."""
    distance = ndi.distance_transform_edt(mask)
    return watershed(-distance, markers=dapi_markers, mask=mask)


def measure_one_channel(image_base, channel_label, macro_out_dir, dapi_markers,
                        pixel_size_um, min_object_px):
    area_mask_path = os.path.join(macro_out_dir, "{0}-{1}_premask.tif".format(image_base, channel_label))
    branch_mask_path = os.path.join(macro_out_dir, "{0}-{1}_branch_premask.tif".format(image_base, channel_label))
    branch_lt_path = os.path.join(macro_out_dir, "{0}-{1}_branch_LocalThickness.tif".format(image_base, channel_label))
    if not os.path.exists(branch_mask_path):
        return [], None

    # Area / hypertrophy: from the area-branch blob mask (per cell).
    area_mask = tifffile.imread(area_mask_path) > 0
    area_labels = _split(area_mask, dapi_markers)
    area_by_label = {r.label: r.area for r in regionprops(area_labels)}

    # Shape / skeleton / thickness: from the process-preserving branch mask.
    branch_mask = tifffile.imread(branch_mask_path) > 0
    branch_mask = remove_small_objects(branch_mask, min_size=SPECKLE_MIN_PX)
    branch_labels = _split(branch_mask, dapi_markers)

    # Read with PIL, not tifffile -- the LocalThickness map is a byte-swapped
    # 32-bit float TIFF and this env's tifffile crashes on it under NumPy 2.0.
    from PIL import Image as _PILImage
    thickness_map_full = np.array(_PILImage.open(branch_lt_path), dtype=np.float32)

    er = erosion_radius_px(pixel_size_um)
    rows = []
    skel_full = np.zeros_like(branch_mask, dtype=bool)
    thickness_full = np.zeros(branch_mask.shape, dtype=np.float32)
    soma_full = np.zeros_like(branch_mask, dtype=bool)

    for region in regionprops(branch_labels):
        branch_area_px = region.area
        if branch_area_px < min_object_px:
            continue

        # Crop per-pixel ops to the cell's own bounding box (all local
        # operations; cell never extends past its own bbox).
        minr, minc, maxr, maxc = region.bbox
        sl = (slice(minr, maxr), slice(minc, maxc))
        cell_mask = branch_labels[sl] == region.label

        skel = skeletonize(cell_mask)
        n_skeleton_px = int(skel.sum())
        if n_skeleton_px == 0:
            continue
        skel_full[sl] |= skel

        thickness_map_cell = thickness_map_full[sl]
        thickness_at_skel = thickness_map_cell[skel]
        valid = ~np.isnan(thickness_at_skel)
        thickness_full[sl][skel] = np.nan_to_num(thickness_at_skel)
        thickness_vals = thickness_at_skel[valid]
        if len(thickness_vals) == 0:
            continue

        eroded = ndi.binary_erosion(cell_mask, structure=np.ones((3, 3)), iterations=er)
        soma_full[sl] |= eroded
        is_soma_skel = eroded[skel][valid]
        soma_vals = thickness_vals[is_soma_skel]
        process_vals = thickness_vals[~is_soma_skel]

        perim = region.perimeter if region.perimeter > 0 else np.nan
        circularity = min(1.0, (4 * np.pi * branch_area_px) / (perim ** 2)) if perim and perim > 0 else np.nan
        solidity = float(region.solidity) if region.convex_area > 0 else np.nan
        aspect_ratio = float(region.major_axis_length / region.minor_axis_length) if region.minor_axis_length > 0 else np.nan

        # hypertrophy area from the area-blob mask for this same nucleus;
        # NaN if that nucleus had no area-mask object (signal only in one mask)
        hyper_area_px = area_by_label.get(region.label, np.nan)

        rows.append(dict(
            dapi_seed_label=int(region.label),
            # hypertrophy: cell-body footprint from the area-blob mask
            area_um2=(hyper_area_px * pixel_size_um ** 2) if hyper_area_px == hyper_area_px else np.nan,
            # shape/skeleton footprint from the process-preserving mask
            branch_mask_area_um2=branch_area_px * pixel_size_um ** 2,
            circularity=circularity,
            solidity=solidity,
            aspect_ratio=aspect_ratio,               # elongation
            n_skeleton_px=n_skeleton_px,
            skeleton_length_um=true_skeleton_length_um(skel, pixel_size_um),
            n_branches=skeleton_branch_count(skel),  # branching
            whole_cell_median_thickness_um=float(np.median(thickness_vals)),
            whole_cell_mean_thickness_um=float(np.mean(thickness_vals)),
            process_median_thickness_um=float(np.median(process_vals)) if len(process_vals) else np.nan,
            soma_median_thickness_um=float(np.median(soma_vals)) if len(soma_vals) else np.nan,
            soma_area_um2=float(eroded.sum() * pixel_size_um ** 2),
            soma_area_fraction=float(eroded.sum() / branch_area_px),
            has_resolved_process=bool(len(process_vals) > 0),
            channel=channel_label,
            image=image_base,
        ))

    return rows, (branch_mask, skel_full, thickness_full, soma_full)


def save_qc_png(image_base, channel_label, premask, skel_full, thickness_full, soma_full, qc_dir):
    from PIL import Image
    h, w = premask.shape
    panel_mask = np.zeros((h, w, 3), dtype=np.uint8)
    panel_mask[premask] = [80, 80, 80]
    panel_mask[soma_full] = [30, 110, 200]
    panel_mask[skel_full] = [255, 220, 0]
    tmax = thickness_full.max() if thickness_full.max() > 0 else 1.0
    heat = (np.clip(thickness_full / tmax, 0, 1) * 255).astype(np.uint8)
    panel_heat = np.stack([heat, np.zeros_like(heat), 255 - heat], axis=-1)
    panel_heat[~premask] = 0
    combined = np.concatenate([panel_mask, panel_heat], axis=1)
    Image.fromarray(combined).save(
        os.path.join(qc_dir, "{0}-{1}_v2_QC.png".format(image_base, channel_label)))


def _process_image(args):
    """Worker: one image, both channels. Pickle-safe (module-level,
    plain-tuple args) so it runs in a multiprocessing Pool."""
    (image_base, macro_out_dir, pixel_size_um, min_object_px, save_qc) = args
    dapi_roi_path = os.path.join(macro_out_dir, "{0}-DAPI_RoiSet.zip".format(image_base))
    sample = os.path.join(macro_out_dir, "{0}-GFAP_branch_premask.tif".format(image_base))
    shape = tifffile.imread(sample).shape
    dapi_markers = dapi_rois_to_markers(dapi_roi_path, shape)
    rows_out = []
    for channel_label in ["GFAP", "Vimentin"]:
        rows, viz = measure_one_channel(image_base, channel_label, macro_out_dir,
                                        dapi_markers, pixel_size_um, min_object_px)
        rows_out.extend(rows)
        if save_qc and viz is not None:
            save_qc_png(image_base, channel_label, *viz, qc_dir=save_qc)
    return rows_out


def main(macro_out_dir, qc_dir, pixel_size_um, min_object_px, name_filter=None,
         max_images=None, save_qc=True, n_workers=None, qc_every=1, out_csv=None):
    os.makedirs(qc_dir, exist_ok=True)
    branch_files = [
        f for f in glob.glob(os.path.join(macro_out_dir, "*_branch_premask.tif"))
        if not os.path.basename(f).startswith(".")
    ]
    image_bases = sorted(set(
        os.path.basename(f).split("-GFAP_branch_premask.tif")[0].split("-Vimentin_branch_premask.tif")[0]
        for f in branch_files
    ))
    if name_filter:
        image_bases = [b for b in image_bases if name_filter in b]
    if max_images:
        image_bases = image_bases[:int(max_images)]

    tasks = [
        (b, macro_out_dir, pixel_size_um, min_object_px,
         (qc_dir if (save_qc and i % qc_every == 0) else False))
        for i, b in enumerate(image_bases)
    ]

    if n_workers is None:
        n_workers = max(1, (os.cpu_count() or 2) - 1)

    all_rows = []
    if n_workers == 1:
        for i, t in enumerate(tasks):
            all_rows.extend(_process_image(t))
            if i % 200 == 0:
                print("Progress: {0}/{1}".format(i, len(tasks)))
    else:
        import multiprocessing as mp
        with mp.Pool(n_workers) as pool:
            for i, rows in enumerate(pool.imap_unordered(_process_image, tasks, chunksize=4)):
                all_rows.extend(rows)
                if i % 200 == 0:
                    print("Progress: {0}/{1}".format(i, len(tasks)))

    df = pd.DataFrame(all_rows)
    if out_csv is None:
        out_csv = os.path.join(macro_out_dir, "morphology_v2_dapiseeded_results.csv")
    df.to_csv(out_csv, index=False)
    print("Saved {0} rows to {1} ({2} workers)".format(len(df), out_csv, n_workers))
    return df


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("macro_out_dir")
    p.add_argument("qc_dir", nargs="?", default=None)
    p.add_argument("pixel_size_um", nargs="?", type=float, default=PIXEL_SIZE_UM)
    p.add_argument("min_object_px", nargs="?", type=int, default=MIN_OBJECT_PX)
    p.add_argument("--filter", default=None, help="only images whose name contains this substring")
    p.add_argument("--max-images", type=int, default=None)
    p.add_argument("--workers", type=int, default=None, help="parallel processes (default: cores-1)")
    p.add_argument("--qc-every", type=int, default=25, help="write a QC PNG for every Nth image (1=all)")
    p.add_argument("--no-qc", action="store_true")
    p.add_argument("--out", default=None,
                   help="output CSV path (default: <macro_out_dir>/morphology_v2_dapiseeded_results.csv). "
                        "ALWAYS set this for subset/test runs so they don't clobber the full results.")
    a = p.parse_args()
    qc_dir = a.qc_dir or os.path.join(a.macro_out_dir, "qc_images_v2")
    main(a.macro_out_dir, qc_dir, a.pixel_size_um, a.min_object_px,
         name_filter=a.filter, max_images=a.max_images, save_qc=not a.no_qc,
         n_workers=a.workers, qc_every=a.qc_every, out_csv=a.out)
