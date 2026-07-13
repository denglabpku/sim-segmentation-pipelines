# Overview
3D-SIM imaging analysis using the data recorded on a Multi-SIM imaging system (Multimodality Structured Illumination Microscopy; NanoInsights-Tech Co., Ltd.) equipped with a 100× NA 1.49 oil-immersion objective (Nikon CFI SR HP Apo 100× oil), an sCMOS camera (Photometrics Kinetix; 65 nm pixel size).

Written by Zuhui Wang

# Dependent Softwares/Packages
Fiji ImageJ with TrackMate installed.

Matlab R2025b (older Matlab should be fine)

# Instructions
This pipeline is designed to analyze cluster numbers and sizes captured using 3D-SIM. The workflow consists of the following steps. Skip step 0 and if only single channel 3D-SIM analysis is required.

## 0. Bead-Based Image Alignment

### 0.1 Bead Image Preprocessing

Run the section **`%% Prepare multi-color mrc beads image`** in **`script_prepareInputStacks.m`** to convert the multi-channel bead MRC images into TIF format.

If necessary (e.g, uneven background fluctation, high background), subtract the beads image background using the ImageJ macro before localizations using Fiji.ImageJ macro:

**`general_background_substract_rollingball.ijm`**

### 0.2 Generate TrackMate Localization XML

1. Create one folder for each imaging channel (e.g., `C1_localizations`, `C2_localizations`, `C3_localizations`) to store the localization results.
2. Open a representative bead image in Fiji/ImageJ.
3. Set the voxel size via:

   **Image → Properties...**

   and click **Global** so that the settings apply to all subsequently opened images.
4. Run TrackMate manually on the representative bead image.
5. Export the TrackMate settings as an XML file and save a copy into each localization folder created above.

### 0.3 Batch Run TrackMate Localization

Run **`ij_batch_TrackMate_localization.py`** in Fiji.ImageJ and select the exported TrackMate XML file to batch-process all selected bead images.

### 0.4 Image Alignment (`uncropped_aligned_imgs`)

Run the section

**`%% REGISTRATION IMAGE USING SEPARATE CHANNEL IMAGES (MRC SIM data)`**

in **`script_beads_alignment_2or3Ch_codex.m`** to align the multi-channel images.

> **Important:** Perform image alignment **before** cropping individual cells. Image registration requires both the bead images and the experimental images to share the same image coordinates. Cropping before alignment will alter the coordinate system and lead to incorrect registration.

---

## 1. Image Preprocessing

Run the section **`%% Prepare dual-color tif SIM image (usually after alignment)`** in **`script_prepareInputStacks_MRC.m`** to crop each cell.

---

## 2. Segmentation

Run **`script_3D_segmentation_2Ch_MRC.m`** to segment clusters from the aligned multi-channel image stacks of each cell.

### 2.1 Choosing `globalOtsu_raw` vs. `globalOtsu_norm`

When photobleaching is highly visible across Z-slices (e.g., in live-cell SIM), `globalOtsu` should use `Vnorm` instead of `V` to ensure effective thresholding across all slices.

When photobleaching is minimal (e.g., in fixed-cell SIM), using `V` is sufficient. Using `Vnorm` in this case may inadvertently amplify SIM reconstruction artifacts.

### 2.2 Selecting an Optimal Z-Range

> **Important:** This step is particularly critical when using `globalOtsu_norm`.

1. Use the MATLAB-generated normalized mean/median intensity plots and normalized FFT magnitude plots to determine which channel should be inspected.
2. Select the channel with both low normalized intensity and low FFT magnitude.
3. Open the selected channel in Fiji/ImageJ.

If using `globalOtsu_norm`:

- In **Fiji**, navigate to:

  **Image → Color → Threshold**

  and enable **Stack histogram** to overlay the thresholded image on the raw (or normalized) image. This allows you to determine the slice range where Global Otsu thresholding performs accurately, and only use Z-range that is free of SIM reconstruction artifacts.

> **Testing note:** For per-Z normalization, using the **mean** generally stabilizes downstream Otsu thresholding better than the **median** when processing live-cell SIM datasets with severe photobleaching.

### 2.3 Adjusting `params.LoG` and `segmentFusion`

- Choose a LoG size that closely matches your microscope's resolution limit.
- In most cases, keep `h_value = 0`.
- `DistanceThreshold = 1` corresponds to one unit of XY resolution.
- `DistanceThreshold = 3` generally allows small clusters to be separated while preventing larger clusters from being over-segmented.
- Always visually inspect the segmentation results to verify parameter selection.

---

## 3. Visualization

Import the following into Imaris:

- Raw images
- Segmentation outputs
- Weighted centroids

> **Note:** Manually configure the image properties in Imaris so that the voxel size (pixel dimensions) exactly matches the values specified in your Step 3 configuration.

# Reference
[1] Kraus F, Miron E, Demmerle J, et al. Quantitative 3D structured illumination microscopy of nuclear structures[J]. Nature Protocols, 2017, 12(5): 1011-1028.

[2] Sage D, Neumann F R, Hediger F, et al. Automatic tracking of individual fluorescence particles: application to the study of chromosome dynamics[J]. IEEE Transactions on Image Processing, 2005, 14(9): 1372-1383.

