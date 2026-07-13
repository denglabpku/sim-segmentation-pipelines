% After-segmentation analysis for script_3D_segmentation_2Ch_MRC.m outputs.
%
% This script computes:
% 1. Nearest-neighbor distances between two channels' weighted centroids,
%    with randomized centroid controls using the same detected cluster counts.
% 2. Segmented-voxel overlap between the two channels, with randomized voxel
%    permutation controls preserving foreground voxel counts.
% 3. Per-cluster and per-cell cluster number, volume, and integrated intensity.
% 4. Summary plots for the measurements above.

clear; clc; close all;

analysis_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(analysis_script_dir);
addpath(genpath(repo_root), '-end');

%% User settings
% This root_dir matches the current MRC segmentation script. If the folder is
% not found, the script will ask you to select a seg_result folder.
root_dir = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260612_DLD1_RPB1-JF549_s5p-mintbody488';
seg_parent_dir = fullfile(root_dir, 'cropped_aligned_imgs', 'seg_result');

channel_names = ["C488", "C561"];
n_permutations = 1000;
rng_seed = 1;
save_plot_formats = ["png", "pdf"];

% When true, the script reconstructs the same nuclear ROI crop used during
% segmentation if the corresponding *_NuclearMask.mat file is available.
% Random centroid and voxel controls are then sampled only inside that ROI.
use_roi_mask_when_available = true;

seg_parent_override = getenv('AFTER_SEG_ANALYSIS_SEG_PARENT_DIR');
if ~isempty(seg_parent_override)
    seg_parent_dir = seg_parent_override;
end

n_permutations_override = getenv('AFTER_SEG_ANALYSIS_N_PERMUTATIONS');
if ~isempty(n_permutations_override)
    parsed_n_permutations = str2double(n_permutations_override);
    if isfinite(parsed_n_permutations) && parsed_n_permutations >= 0
        n_permutations = round(parsed_n_permutations);
    else
        warning('Ignoring invalid AFTER_SEG_ANALYSIS_N_PERMUTATIONS value: %s', ...
            n_permutations_override);
    end
end

if ~isfolder(seg_parent_dir)
    selected_dir = uigetdir(repo_root, ...
        'Select seg_result folder, or one per-cell segmentation result folder');
    if isequal(selected_dir, 0)
        error('No segmentation result folder was selected.');
    end
    seg_parent_dir = selected_dir;
end

analysis_out_dir = fullfile(root_dir, 'cropped_aligned_imgs', 'AfterSegAnalysis');
if ~exist(analysis_out_dir, 'dir')
    mkdir(analysis_out_dir);
end

rng(rng_seed);

runs = findSegmentationRuns(seg_parent_dir, channel_names);
if isempty(runs)
    error('No completed two-channel segmentation result folders were found under: %s', seg_parent_dir);
end
fprintf('Found %d completed segmentation result folders.\n', numel(runs));

config = struct();
config.root_dir = root_dir;
config.seg_parent_dir = seg_parent_dir;
config.analysis_out_dir = analysis_out_dir;
config.channel_names = channel_names;
config.n_permutations = n_permutations;
config.rng_seed = rng_seed;
config.use_roi_mask_when_available = use_roi_mask_when_available;

all_cluster_metrics = table();
all_cluster_summary = table();
all_nearest_neighbors = table();
all_nearest_summary = table();
all_nearest_random = table();
all_overlap_summary = table();
all_overlap_random = table();

for i_run = 1:numel(runs)
    run_name = runs(i_run).name;
    run_folder = runs(i_run).folder;
    fprintf('Processing %s (%d/%d)\n', run_name, i_run, numel(runs));

    [stats_ch1, params_ch1] = readStatsTable(runs(i_run).statsFiles(1));
    [stats_ch2, params_ch2] = readStatsTable(runs(i_run).statsFiles(2));
    voxel_size_um = getVoxelSizeUm(params_ch1, params_ch2);

    labels_ch1 = readLabelVolume(runs(i_run).labelFiles(1));
    labels_ch2 = readLabelVolume(runs(i_run).labelFiles(2));
    [labels_ch1, labels_ch2] = cropToCommonVolume(labels_ch1, labels_ch2, run_name);
    bw_ch1 = labels_ch1 > 0;
    bw_ch2 = labels_ch2 > 0;

    [analysis_mask, mask_source] = loadAnalysisMask( ...
        run_folder, run_name, size3(labels_ch1), use_roi_mask_when_available);
    if nnz((bw_ch1 | bw_ch2) & ~analysis_mask) > 0
        warning(['ROI mask for %s did not cover all segmented voxels; ', ...
            'using the full cropped volume for randomized controls.'], run_name);
        analysis_mask = true(size3(labels_ch1));
        mask_source = "WholeVolume_ROIMismatch";
    end

    centroids_ch1_um = getWeightedCentroidsUm(stats_ch1, voxel_size_um);
    centroids_ch2_um = getWeightedCentroidsUm(stats_ch2, voxel_size_um);

    metrics_ch1 = makeClusterMetricTable(stats_ch1, centroids_ch1_um, ...
        voxel_size_um, run_name, channel_names(1));
    metrics_ch2 = makeClusterMetricTable(stats_ch2, centroids_ch2_um, ...
        voxel_size_um, run_name, channel_names(2));
    all_cluster_metrics = appendTable(all_cluster_metrics, metrics_ch1);
    all_cluster_metrics = appendTable(all_cluster_metrics, metrics_ch2);
    all_cluster_summary = appendTable(all_cluster_summary, ...
        summarizeClusterMetrics(metrics_ch1, run_name, channel_names(1)));
    all_cluster_summary = appendTable(all_cluster_summary, ...
        summarizeClusterMetrics(metrics_ch2, run_name, channel_names(2)));

    nn_ch1_to_ch2 = buildNearestNeighborTable(centroids_ch1_um, centroids_ch2_um, ...
        stats_ch1, stats_ch2, run_name, channel_names(1), channel_names(2));
    nn_ch2_to_ch1 = buildNearestNeighborTable(centroids_ch2_um, centroids_ch1_um, ...
        stats_ch2, stats_ch1, run_name, channel_names(2), channel_names(1));
    all_nearest_neighbors = appendTable(all_nearest_neighbors, nn_ch1_to_ch2);
    all_nearest_neighbors = appendTable(all_nearest_neighbors, nn_ch2_to_ch1);

    nn_random_ch1_to_ch2 = computeCentroidRandomControls( ...
        centroids_ch1_um, height(stats_ch2), analysis_mask, voxel_size_um, ...
        n_permutations, run_name, channel_names(1), channel_names(2), mask_source);
    nn_random_ch2_to_ch1 = computeCentroidRandomControls( ...
        centroids_ch2_um, height(stats_ch1), analysis_mask, voxel_size_um, ...
        n_permutations, run_name, channel_names(2), channel_names(1), mask_source);
    all_nearest_random = appendTable(all_nearest_random, nn_random_ch1_to_ch2);
    all_nearest_random = appendTable(all_nearest_random, nn_random_ch2_to_ch1);

    all_nearest_summary = appendTable(all_nearest_summary, summarizeNearestNeighbors( ...
        nn_ch1_to_ch2, nn_random_ch1_to_ch2, run_name, channel_names(1), channel_names(2), ...
        height(stats_ch1), height(stats_ch2), mask_source));
    all_nearest_summary = appendTable(all_nearest_summary, summarizeNearestNeighbors( ...
        nn_ch2_to_ch1, nn_random_ch2_to_ch1, run_name, channel_names(2), channel_names(1), ...
        height(stats_ch2), height(stats_ch1), mask_source));

    overlap_actual = computeVoxelOverlap(bw_ch1, bw_ch2, analysis_mask, ...
        voxel_size_um, run_name, channel_names(1), channel_names(2), mask_source);
    overlap_random = computeVoxelRandomControls(bw_ch1, bw_ch2, analysis_mask, ...
        voxel_size_um, n_permutations, run_name, channel_names(1), channel_names(2), mask_source);
    all_overlap_random = appendTable(all_overlap_random, overlap_random);
    all_overlap_summary = appendTable(all_overlap_summary, ...
        summarizeVoxelOverlap(overlap_actual, overlap_random));
end

%% Save tables and workspace
writetable(all_cluster_metrics, fullfile(analysis_out_dir, 'cluster_metrics_per_cluster.csv'));
writetable(all_cluster_summary, fullfile(analysis_out_dir, 'cluster_metrics_per_cell_summary.csv'));
writetable(all_nearest_neighbors, fullfile(analysis_out_dir, 'nearest_neighbor_distances_per_cluster.csv'));
writetable(all_nearest_summary, fullfile(analysis_out_dir, 'nearest_neighbor_summary_per_cell.csv'));
writetable(all_nearest_random, fullfile(analysis_out_dir, 'nearest_neighbor_random_controls.csv'));
writetable(all_overlap_summary, fullfile(analysis_out_dir, 'voxel_overlap_summary_per_cell.csv'));
writetable(all_overlap_random, fullfile(analysis_out_dir, 'voxel_overlap_random_controls.csv'));

save(fullfile(analysis_out_dir, 'afterSegAnalysis_workspace.mat'), ...
    'config', 'runs', 'all_cluster_metrics', 'all_cluster_summary', ...
    'all_nearest_neighbors', 'all_nearest_summary', 'all_nearest_random', ...
    'all_overlap_summary', 'all_overlap_random', '-v7.3');

makeSummaryPlots(all_nearest_summary, all_nearest_neighbors, ...
    all_overlap_summary, all_cluster_metrics, all_cluster_summary, ...
    analysis_out_dir, save_plot_formats);

fprintf('Analysis complete. Results saved to:\n%s\n', analysis_out_dir);

%% Local functions
function runs = findSegmentationRuns(seg_parent_dir, channel_names)
    candidate_dirs = string(seg_parent_dir);
    listing = dir(seg_parent_dir);
    for i = 1:numel(listing)
        if listing(i).isdir && ~ismember(listing(i).name, {'.', '..', 'AfterSegAnalysis','Archive'})
            candidate_dirs(end+1, 1) = string(fullfile(seg_parent_dir, listing(i).name)); %#ok<AGROW>
        end
    end

    runs = struct('name', {}, 'folder', {}, 'statsFiles', {}, 'labelFiles', {});
    for i_dir = 1:numel(candidate_dirs)
        run_folder = candidate_dirs(i_dir);
        stats_files = strings(1, numel(channel_names));
        label_files = strings(1, numel(channel_names));
        has_all_files = true;

        for i_ch = 1:numel(channel_names)
            ch = channel_names(i_ch);
            stats_files(i_ch) = firstExistingFile(run_folder, ...
                [ch + "_segResults.mat", ch + "_stats3D_table.csv"]);
            label_files(i_ch) = firstExistingFile(run_folder, ch + "_good_seglabels.tif");
            if stats_files(i_ch) == "" || label_files(i_ch) == ""
                has_all_files = false;
                break;
            end
        end

        if has_all_files
            runs(end+1).name = folderBaseName(run_folder); %#ok<AGROW>
            runs(end).folder = run_folder;
            runs(end).statsFiles = stats_files;
            runs(end).labelFiles = label_files;
        end
    end
end

function path_out = firstExistingFile(folder, names)
    path_out = "";
    for i = 1:numel(names)
        test_path = fullfile(char(folder), char(names(i)));
        if isfile(test_path)
            path_out = string(test_path);
            return;
        end
    end
end

function name = folderBaseName(folder)
    [~, name_char] = fileparts(char(folder));
    name = string(name_char);
end

function [stats_table, params] = readStatsTable(stats_file)
    params = struct();
    [~, ~, ext] = fileparts(char(stats_file));
    if strcmpi(ext, '.mat')
        loaded = load(char(stats_file), 'stats3D_table', 'params');
        if ~isfield(loaded, 'stats3D_table')
            error('Missing stats3D_table in %s', stats_file);
        end
        stats_table = loaded.stats3D_table;
        if isfield(loaded, 'params')
            params = loaded.params;
        end
    else
        stats_table = readtable(char(stats_file));
    end
end

function labels = readLabelVolume(label_file)
    try
        labels = tiffreadVolume(char(label_file));
    catch
        labels = readTiffStackWithTiffClass(label_file);
    end
    labels = uint32(labels);
    if ismatrix(labels)
        labels = reshape(labels, size(labels, 1), size(labels, 2), 1);
    end
end

function stack = readTiffStackWithTiffClass(label_file)
    info = imfinfo(char(label_file));
    n_slices = numel(info);
    t = Tiff(char(label_file), 'r');
    cleanup_obj = onCleanup(@() close(t));

    for k = 1:n_slices
        setDirectory(t, k);
        slice = read(t);
        if k == 1
            stack = zeros(size(slice, 1), size(slice, 2), n_slices, 'like', slice);
        end
        stack(:, :, k) = slice;
    end
end

function [a, b] = cropToCommonVolume(a, b, run_name)
    size_a = size3(a);
    size_b = size3(b);
    if isequal(size_a, size_b)
        return;
    end
    common_size = min(size_a, size_b);
    warning('Label volume sizes differ for %s; cropping both channels to [%d %d %d].', ...
        run_name, common_size(1), common_size(2), common_size(3));
    a = a(1:common_size(1), 1:common_size(2), 1:common_size(3));
    b = b(1:common_size(1), 1:common_size(2), 1:common_size(3));
end

function sz = size3(volume)
    sz = size(volume);
    if numel(sz) < 3
        sz(3) = 1;
    end
    sz = sz(1:3);
end

function voxel_size_um = getVoxelSizeUm(params_a, params_b)
    voxel_size_um = struct();
    voxel_size_um.pixelXY = readNestedField(params_a, 'imaging', 'pixel_size', NaN);
    voxel_size_um.zStep = readNestedField(params_a, 'imaging', 'z_step', NaN);
    voxel_size_um.volume = readNestedField(params_a, 'measure', 'voxelSize', NaN);

    if isnan(voxel_size_um.pixelXY)
        voxel_size_um.pixelXY = readNestedField(params_b, 'imaging', 'pixel_size', NaN);
    end
    if isnan(voxel_size_um.zStep)
        voxel_size_um.zStep = readNestedField(params_b, 'imaging', 'z_step', NaN);
    end
    if isnan(voxel_size_um.volume)
        voxel_size_um.volume = readNestedField(params_b, 'measure', 'voxelSize', NaN);
    end
    if isnan(voxel_size_um.pixelXY)
        warning('Could not read pixel_size from params; using 1 um for XY.');
        voxel_size_um.pixelXY = 1;
    end
    if isnan(voxel_size_um.zStep)
        warning('Could not read z_step from params; using 1 um for Z.');
        voxel_size_um.zStep = 1;
    end
    if isnan(voxel_size_um.volume)
        voxel_size_um.volume = voxel_size_um.pixelXY ^ 2 * voxel_size_um.zStep;
    end
end

function value = readNestedField(s, parent_name, child_name, default_value)
    value = default_value;
    if isstruct(s) && isfield(s, parent_name)
        parent = s.(parent_name);
        if isstruct(parent) && isfield(parent, child_name)
            value = parent.(child_name);
        end
    end
end

function [analysis_mask, mask_source] = loadAnalysisMask(run_folder, run_name, volume_size, use_roi_mask)
    analysis_mask = true(volume_size);
    mask_source = "WholeVolume";
    if ~use_roi_mask
        return;
    end

    cell_tokens = regexp(char(run_name), '_Cell(\d+)', 'tokens', 'once');
    if isempty(cell_tokens)
        return;
    end
    cell_idx = str2double(cell_tokens{1});
    base_name = regexprep(char(run_name), '_Cell\d+.*$', '');

    seg_result_dir = fileparts(char(run_folder));
    image_dir = fileparts(seg_result_dir);
    mask_file = fullfile(image_dir, sprintf('%s_NuclearMask.mat', base_name));
    if ~isfile(mask_file)
        return;
    end

    loaded = load(mask_file, 'roi_info_nuc', 'ImHeight', 'ImWidth');
    if ~isfield(loaded, 'roi_info_nuc') || cell_idx > numel(loaded.roi_info_nuc)
        return;
    end

    poly_verts = loaded.roi_info_nuc{cell_idx};
    x_min = floor(min(poly_verts(:, 1)));
    x_max = ceil(max(poly_verts(:, 1)));
    y_min = floor(min(poly_verts(:, 2)));
    y_max = ceil(max(poly_verts(:, 2)));

    nuclear_mask_full = poly2mask(poly_verts(:, 1), poly_verts(:, 2), ...
        loaded.ImHeight, loaded.ImWidth);
    roi_mask_2d = nuclear_mask_full(y_min:y_max, x_min:x_max);
    if ~isequal(size(roi_mask_2d), volume_size(1:2))
        warning('Reconstructed ROI mask size does not match %s label volume.', run_name);
        return;
    end

    analysis_mask = repmat(roi_mask_2d, 1, 1, volume_size(3));
    mask_source = "NuclearROI";
end

function centroids_um = getWeightedCentroidsUm(stats_table, voxel_size_um)
    n = height(stats_table);
    centroids_um = nan(n, 3);
    vars = string(stats_table.Properties.VariableNames);

    if all(ismember(["WeightedCentroid_x_um", "WeightedCentroid_y_um", "WeightedCentroid_z_um"], vars))
        centroids_um = [double(stats_table.WeightedCentroid_x_um), ...
            double(stats_table.WeightedCentroid_y_um), ...
            double(stats_table.WeightedCentroid_z_um)];
        return;
    end

    if ismember("WeightedCentroid", vars)
        centroid_pixels = stats_table.WeightedCentroid;
    elseif all(ismember(["WeightedCentroid_1", "WeightedCentroid_2", "WeightedCentroid_3"], vars))
        centroid_pixels = [stats_table.WeightedCentroid_1, ...
            stats_table.WeightedCentroid_2, stats_table.WeightedCentroid_3];
    elseif ismember("Centroid", vars)
        centroid_pixels = stats_table.Centroid;
    elseif all(ismember(["Centroid_1", "Centroid_2", "Centroid_3"], vars))
        centroid_pixels = [stats_table.Centroid_1, stats_table.Centroid_2, stats_table.Centroid_3];
    else
        warning('No weighted centroid columns were found. Centroid distance outputs will be NaN.');
        return;
    end

    if isempty(centroid_pixels)
        return;
    end
    centroid_pixels = double(centroid_pixels);
    centroids_um = [(centroid_pixels(:, 1) - 0.5) .* voxel_size_um.pixelXY, ...
        (centroid_pixels(:, 2) - 0.5) .* voxel_size_um.pixelXY, ...
        (centroid_pixels(:, 3) - 0.5) .* voxel_size_um.zStep];
end

function cluster_table = makeClusterMetricTable(stats_table, centroids_um, voxel_size_um, cell_id, channel_name)
    n = height(stats_table);
    cluster_id = getClusterIDs(stats_table);
    volume_voxels = getNumericColumn(stats_table, 'Volume', nan(n, 1));
    volume_um3 = getNumericColumn(stats_table, 'VolumePhysical', volume_voxels .* voxel_size_um.volume);
    mean_intensity = getNumericColumn(stats_table, 'MeanIntensity', nan(n, 1));
    integrated_intensity = volume_voxels .* mean_intensity;
    equi_radius_um = getNumericColumn(stats_table, 'equiRadius', ((3 .* volume_um3) ./ (4 .* pi)) .^ (1/3));

    cluster_table = table( ...
        repmat(string(cell_id), n, 1), ...
        repmat(string(channel_name), n, 1), ...
        cluster_id, volume_voxels, volume_um3, mean_intensity, integrated_intensity, ...
        equi_radius_um, centroids_um(:, 1), centroids_um(:, 2), centroids_um(:, 3), ...
        'VariableNames', {'CellID', 'Channel', 'ClusterID', 'Volume_voxels', ...
        'Volume_um3', 'MeanIntensity', 'IntegratedIntensity_VolumeTimesMean', ...
        'EquiRadius_um', 'WeightedCentroidX_um', 'WeightedCentroidY_um', ...
        'WeightedCentroidZ_um'});
end

function summary_table = summarizeClusterMetrics(cluster_table, cell_id, channel_name)
    summary_table = table( ...
        string(cell_id), string(channel_name), height(cluster_table), ...
        nanMean(cluster_table.Volume_voxels), nanMedian(cluster_table.Volume_voxels), ...
        nanMean(cluster_table.Volume_um3), nanMedian(cluster_table.Volume_um3), ...
        nanMean(cluster_table.MeanIntensity), nanMedian(cluster_table.MeanIntensity), ...
        nanMean(cluster_table.IntegratedIntensity_VolumeTimesMean), ...
        nanMedian(cluster_table.IntegratedIntensity_VolumeTimesMean), ...
        'VariableNames', {'CellID', 'Channel', 'ClusterCount', ...
        'MeanVolume_voxels', 'MedianVolume_voxels', 'MeanVolume_um3', ...
        'MedianVolume_um3', 'MeanIntensity', 'MedianIntensity', ...
        'MeanIntegratedIntensity', 'MedianIntegratedIntensity'});
end

function ids = getClusterIDs(stats_table)
    if ismember('ClusterID', stats_table.Properties.VariableNames)
        ids = double(stats_table.ClusterID);
    else
        ids = (1:height(stats_table))';
    end
    ids = ids(:);
end

function values = getNumericColumn(stats_table, var_name, default_values)
    if ismember(var_name, stats_table.Properties.VariableNames)
        values = stats_table.(var_name);
        if iscell(values)
            values = str2double(values);
        end
        values = double(values);
    else
        values = default_values;
    end
    values = values(:);
end

function nn_table = buildNearestNeighborTable(source_points_um, target_points_um, ...
        source_stats, target_stats, cell_id, source_channel, target_channel)
    n_source = size(source_points_um, 1);
    target_cluster_ids = getClusterIDs(target_stats);
    source_cluster_ids = getClusterIDs(source_stats);

    [nearest_distance_um, nearest_target_idx] = nearestNeighborDistances(source_points_um, target_points_um);
    nearest_target_cluster_id = nan(n_source, 1);
    valid_idx = isfinite(nearest_target_idx);
    if any(valid_idx)
        nearest_target_cluster_id(valid_idx) = target_cluster_ids(nearest_target_idx(valid_idx));
    end

    nn_table = table( ...
        repmat(string(cell_id), n_source, 1), ...
        repmat(string(source_channel), n_source, 1), ...
        repmat(string(target_channel), n_source, 1), ...
        repmat(string(source_channel) + "->" + string(target_channel), n_source, 1), ...
        source_cluster_ids, nearest_target_cluster_id, nearest_distance_um, ...
        'VariableNames', {'CellID', 'SourceChannel', 'TargetChannel', 'Direction', ...
        'SourceClusterID', 'NearestTargetClusterID', 'NearestDistance_um'});
end

function [nearest_distance, nearest_idx] = nearestNeighborDistances(source_points, target_points)
    n_source = size(source_points, 1);
    nearest_distance = nan(n_source, 1);
    nearest_idx = nan(n_source, 1);

    valid_source_idx = find(all(isfinite(source_points), 2));
    valid_target_idx = find(all(isfinite(target_points), 2));
    if isempty(valid_source_idx) || isempty(valid_target_idx)
        return;
    end

    target_valid = target_points(valid_target_idx, :);
    chunk_size = 2000;
    for start_idx = 1:chunk_size:numel(valid_source_idx)
        stop_idx = min(start_idx + chunk_size - 1, numel(valid_source_idx));
        source_rows = valid_source_idx(start_idx:stop_idx);
        source_chunk = source_points(source_rows, :);

        dx = source_chunk(:, 1) - target_valid(:, 1)';
        dy = source_chunk(:, 2) - target_valid(:, 2)';
        dz = source_chunk(:, 3) - target_valid(:, 3)';
        d2 = dx .^ 2 + dy .^ 2 + dz .^ 2;
        [min_d2, local_idx] = min(d2, [], 2);

        nearest_distance(source_rows) = sqrt(min_d2);
        nearest_idx(source_rows) = valid_target_idx(local_idx);
    end
end

function random_table = computeCentroidRandomControls(source_points_um, n_random_targets, ...
        analysis_mask, voxel_size_um, n_permutations, cell_id, source_channel, target_channel, mask_source)
    random_mean = nan(n_permutations, 1);
    random_median = nan(n_permutations, 1);
    random_p25 = nan(n_permutations, 1);
    random_p75 = nan(n_permutations, 1);

    for i_perm = 1:n_permutations
        random_target_points = sampleRandomPointsInMask(n_random_targets, analysis_mask, voxel_size_um);
        [d_um, ~] = nearestNeighborDistances(source_points_um, random_target_points);
        random_mean(i_perm) = nanMean(d_um);
        random_median(i_perm) = nanMedian(d_um);
        random_p25(i_perm) = nanPrctile(d_um, 25);
        random_p75(i_perm) = nanPrctile(d_um, 75);
    end

    random_table = table( ...
        repmat(string(cell_id), n_permutations, 1), ...
        repmat(string(source_channel), n_permutations, 1), ...
        repmat(string(target_channel), n_permutations, 1), ...
        repmat(string(source_channel) + "->" + string(target_channel), n_permutations, 1), ...
        (1:n_permutations)', repmat(string(mask_source), n_permutations, 1), ...
        random_mean, random_median, random_p25, random_p75, ...
        'VariableNames', {'CellID', 'SourceChannel', 'TargetChannel', 'Direction', ...
        'Permutation', 'RandomMaskSource', 'RandomMeanDistance_um', ...
        'RandomMedianDistance_um', 'RandomP25Distance_um', 'RandomP75Distance_um'});
end

function points_um = sampleRandomPointsInMask(n_points, analysis_mask, voxel_size_um)
    points_um = zeros(n_points, 3);
    if n_points == 0
        return;
    end
    valid_idx = find(analysis_mask);
    if isempty(valid_idx)
        error('The randomization mask is empty.');
    end
    sampled_idx = valid_idx(randi(numel(valid_idx), n_points, 1));
    [row, col, z] = ind2sub(size3(analysis_mask), sampled_idx);
    points_um = [(double(col) - 1 + rand(n_points, 1)) .* voxel_size_um.pixelXY, ...
        (double(row) - 1 + rand(n_points, 1)) .* voxel_size_um.pixelXY, ...
        (double(z) - 1 + rand(n_points, 1)) .* voxel_size_um.zStep];
end

function summary_table = summarizeNearestNeighbors(nn_table, random_table, cell_id, ...
        source_channel, target_channel, n_source, n_target, mask_source)
    actual_distances = nn_table.NearestDistance_um;
    random_medians = random_table.RandomMedianDistance_um;
    actual_median = nanMedian(actual_distances);
    valid_random_medians = random_medians(isfinite(random_medians));

    if isempty(valid_random_medians) || ~isfinite(actual_median)
        p_random_median_le_actual = NaN;
    else
        p_random_median_le_actual = (1 + sum(valid_random_medians <= actual_median)) ./ ...
            (numel(valid_random_medians) + 1);
    end

    summary_table = table( ...
        string(cell_id), string(source_channel), string(target_channel), ...
        string(source_channel) + "->" + string(target_channel), ...
        n_source, n_target, string(mask_source), ...
        nanMean(actual_distances), actual_median, nanPrctile(actual_distances, 25), ...
        nanPrctile(actual_distances, 75), ...
        nanMean(random_table.RandomMeanDistance_um), nanMean(random_medians), ...
        nanStd(random_medians), actual_median - nanMean(random_medians), ...
        p_random_median_le_actual, ...
        'VariableNames', {'CellID', 'SourceChannel', 'TargetChannel', 'Direction', ...
        'SourceClusterCount', 'TargetClusterCount', 'RandomMaskSource', ...
        'ActualMeanDistance_um', 'ActualMedianDistance_um', 'ActualP25Distance_um', ...
        'ActualP75Distance_um', 'RandomMeanDistanceMean_um', ...
        'RandomMedianDistanceMean_um', 'RandomMedianDistanceStd_um', ...
        'ActualMinusRandomMedianDistance_um', 'P_RandomMedianDistance_LE_Actual'});
end

function overlap_table = computeVoxelOverlap(bw_ch1, bw_ch2, analysis_mask, ...
        voxel_size_um, cell_id, channel1, channel2, mask_source)
    bw_ch1 = bw_ch1 & analysis_mask;
    bw_ch2 = bw_ch2 & analysis_mask;
    n_ch1 = nnz(bw_ch1);
    n_ch2 = nnz(bw_ch2);
    n_overlap = nnz(bw_ch1 & bw_ch2);
    n_union = nnz(bw_ch1 | bw_ch2);

    overlap_table = table( ...
        string(cell_id), string(channel1), string(channel2), string(mask_source), ...
        n_ch1, n_ch2, n_overlap, n_union, ...
        n_ch1 .* voxel_size_um.volume, n_ch2 .* voxel_size_um.volume, ...
        n_overlap .* voxel_size_um.volume, ...
        safeDivide(n_overlap, n_ch1), safeDivide(n_overlap, n_ch2), ...
        safeDivide(n_overlap, n_union), safeDivide(2 .* n_overlap, n_ch1 + n_ch2), ...
        'VariableNames', {'CellID', 'Channel1', 'Channel2', 'RandomMaskSource', ...
        'Channel1Voxels', 'Channel2Voxels', 'OverlapVoxels', 'UnionVoxels', ...
        'Channel1Volume_um3', 'Channel2Volume_um3', 'OverlapVolume_um3', ...
        'OverlapFractionOfCh1', 'OverlapFractionOfCh2', 'JaccardIndex', ...
        'DiceCoefficient'});
end

function random_table = computeVoxelRandomControls(bw_fixed, bw_randomized_source, ...
        analysis_mask, voxel_size_um, n_permutations, cell_id, fixed_channel, randomized_channel, mask_source)
    bw_fixed = bw_fixed & analysis_mask;
    bw_randomized_source = bw_randomized_source & analysis_mask;
    n_fixed = nnz(bw_fixed);
    n_randomized = nnz(bw_randomized_source);
    valid_idx = find(analysis_mask);
    if n_randomized > numel(valid_idx)
        error('Randomized foreground voxel count exceeds the number of voxels in the analysis mask.');
    end

    fixed_inside_mask = bw_fixed(valid_idx);
    random_overlap = nan(n_permutations, 1);
    random_union = nan(n_permutations, 1);
    random_frac_fixed = nan(n_permutations, 1);
    random_frac_randomized = nan(n_permutations, 1);
    random_jaccard = nan(n_permutations, 1);
    random_dice = nan(n_permutations, 1);

    for i_perm = 1:n_permutations
        if n_randomized == 0
            overlap_voxels = 0;
        else
            selected_local_idx = randperm(numel(valid_idx), n_randomized);
            overlap_voxels = sum(fixed_inside_mask(selected_local_idx));
        end
        union_voxels = n_fixed + n_randomized - overlap_voxels;
        random_overlap(i_perm) = overlap_voxels;
        random_union(i_perm) = union_voxels;
        random_frac_fixed(i_perm) = safeDivide(overlap_voxels, n_fixed);
        random_frac_randomized(i_perm) = safeDivide(overlap_voxels, n_randomized);
        random_jaccard(i_perm) = safeDivide(overlap_voxels, union_voxels);
        random_dice(i_perm) = safeDivide(2 .* overlap_voxels, n_fixed + n_randomized);
    end

    random_table = table( ...
        repmat(string(cell_id), n_permutations, 1), ...
        repmat(string(fixed_channel), n_permutations, 1), ...
        repmat(string(randomized_channel), n_permutations, 1), ...
        repmat(string(mask_source), n_permutations, 1), ...
        (1:n_permutations)', random_overlap, random_union, ...
        random_overlap .* voxel_size_um.volume, random_frac_fixed, ...
        random_frac_randomized, random_jaccard, random_dice, ...
        'VariableNames', {'CellID', 'FixedChannel', 'RandomizedChannel', ...
        'RandomMaskSource', 'Permutation', 'RandomOverlapVoxels', ...
        'RandomUnionVoxels', 'RandomOverlapVolume_um3', ...
        'RandomOverlapFractionOfFixedChannel', ...
        'RandomOverlapFractionOfRandomizedChannel', 'RandomJaccardIndex', ...
        'RandomDiceCoefficient'});
end

function summary_table = summarizeVoxelOverlap(actual_table, random_table)
    random_jaccard = random_table.RandomJaccardIndex;
    random_frac_ch1 = random_table.RandomOverlapFractionOfFixedChannel;
    random_frac_ch2 = random_table.RandomOverlapFractionOfRandomizedChannel;
    random_dice = random_table.RandomDiceCoefficient;

    p_random_jaccard_ge_actual = rightTailPermutationP(random_jaccard, actual_table.JaccardIndex);
    p_random_overlap_ch1_ge_actual = rightTailPermutationP(random_frac_ch1, actual_table.OverlapFractionOfCh1);
    p_random_overlap_ch2_ge_actual = rightTailPermutationP(random_frac_ch2, actual_table.OverlapFractionOfCh2);

    summary_table = actual_table;
    summary_table.RandomOverlapFractionOfCh1Mean = nanMean(random_frac_ch1);
    summary_table.RandomOverlapFractionOfCh1Std = nanStd(random_frac_ch1);
    summary_table.RandomOverlapFractionOfCh2Mean = nanMean(random_frac_ch2);
    summary_table.RandomOverlapFractionOfCh2Std = nanStd(random_frac_ch2);
    summary_table.RandomJaccardIndexMean = nanMean(random_jaccard);
    summary_table.RandomJaccardIndexStd = nanStd(random_jaccard);
    summary_table.RandomDiceCoefficientMean = nanMean(random_dice);
    summary_table.RandomDiceCoefficientStd = nanStd(random_dice);
    summary_table.ActualMinusRandomJaccardIndex = actual_table.JaccardIndex - nanMean(random_jaccard);
    summary_table.P_RandomJaccard_GE_Actual = p_random_jaccard_ge_actual;
    summary_table.P_RandomOverlapFractionOfCh1_GE_Actual = p_random_overlap_ch1_ge_actual;
    summary_table.P_RandomOverlapFractionOfCh2_GE_Actual = p_random_overlap_ch2_ge_actual;
end

function makeSummaryPlots(nn_summary, nn_rows, overlap_summary, ...
        cluster_metrics, cluster_summary, analysis_out_dir, save_plot_formats)
    if height(nn_summary) > 0
        plotNearestNeighborBars(nn_summary, analysis_out_dir, save_plot_formats);
    end
    if height(nn_rows) > 0
        plotNearestNeighborHistograms(nn_rows, analysis_out_dir, save_plot_formats);
    end
    if height(overlap_summary) > 0
        plotOverlapSummary(overlap_summary, analysis_out_dir, save_plot_formats);
    end
    if height(cluster_metrics) > 0
        plotClusterMetrics(cluster_metrics, cluster_summary, analysis_out_dir, save_plot_formats);
    end
end

function plotNearestNeighborBars(nn_summary, analysis_out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 900 430]);
    directions = unique(nn_summary.Direction, 'stable');
    group_labels = strings(1, numel(directions));
    mean_actual = nan(numel(directions), 1);
    mean_random = nan(numel(directions), 1);

    hold on;
    for i_dir = 1:numel(directions)
        idx = nn_summary.Direction == directions(i_dir);
        actual_vals = nn_summary.ActualMedianDistance_um(idx);
        random_vals = nn_summary.RandomMedianDistanceMean_um(idx);
        x_actual = (i_dir - 1) * 3 + 1;
        x_random = (i_dir - 1) * 3 + 2;
        mean_actual(i_dir) = nanMean(actual_vals);
        mean_random(i_dir) = nanMean(random_vals);

        bar(x_actual, mean_actual(i_dir), 0.8, 'FaceColor', [0.20 0.40 0.80], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.75);
        bar(x_random, mean_random(i_dir), 0.8, 'FaceColor', [0.75 0.75 0.75], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.9);
        scatterWithJitter(x_actual, actual_vals, [0.02 0.12 0.35]);
        scatterWithJitter(x_random, random_vals, [0.20 0.20 0.20]);
        group_labels(i_dir) = directions(i_dir);
    end
    ylabel('Median nearest-neighbor distance (um)');
    set(gca, 'XTick', (0:numel(directions)-1) * 3 + 1.5, 'XTickLabel', group_labels, ...
        'TickLabelInterpreter', 'none', 'Box', 'off', 'LineWidth', 1);
    legend({'Actual', 'Randomized'}, 'Location', 'best', 'Box', 'off');
    title('Weighted-centroid nearest-neighbor distances');
    saveFigure(fig, analysis_out_dir, 'nearest_neighbor_actual_vs_random', save_plot_formats);
    close(fig);
end

function plotNearestNeighborHistograms(nn_rows, analysis_out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 900 430]);
    directions = unique(nn_rows.Direction, 'stable');
    tiledlayout(1, numel(directions), 'TileSpacing', 'compact', 'Padding', 'compact');
    for i_dir = 1:numel(directions)
        nexttile;
        idx = nn_rows.Direction == directions(i_dir);
        vals = nn_rows.NearestDistance_um(idx);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            histogram(vals, 'Normalization', 'probability', 'FaceColor', [0.20 0.40 0.80], ...
                'EdgeColor', 'none');
        end
        xlabel('Nearest distance (um)');
        ylabel('Probability');
        title(directions(i_dir), 'Interpreter', 'none');
        set(gca, 'Box', 'off', 'LineWidth', 1);
    end
    saveFigure(fig, analysis_out_dir, 'nearest_neighbor_distance_histograms', save_plot_formats);
    close(fig);
end

function plotOverlapSummary(overlap_summary, analysis_out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 1050 380]);
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plotActualRandomPair(overlap_summary.JaccardIndex, overlap_summary.RandomJaccardIndexMean, ...
        'Jaccard index');
    title('Voxel overlap');

    nexttile;
    plotActualRandomPair(overlap_summary.OverlapFractionOfCh1, ...
        overlap_summary.RandomOverlapFractionOfCh1Mean, ...
        sprintf('Fraction of %s voxels', overlap_summary.Channel1(1)));
    title('Overlap over channel 1');

    nexttile;
    plotActualRandomPair(overlap_summary.OverlapFractionOfCh2, ...
        overlap_summary.RandomOverlapFractionOfCh2Mean, ...
        sprintf('Fraction of %s voxels', overlap_summary.Channel2(1)));
    title('Overlap over channel 2');

    saveFigure(fig, analysis_out_dir, 'voxel_overlap_actual_vs_random', save_plot_formats);
    close(fig);
end

function plotClusterMetrics(cluster_metrics, cluster_summary, analysis_out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 1050 760]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    channels = unique(cluster_metrics.Channel, 'stable');

    nexttile;
    hold on;
    for i_ch = 1:numel(channels)
        idx = cluster_summary.Channel == channels(i_ch);
        vals = cluster_summary.ClusterCount(idx);
        bar(i_ch, nanMean(vals), 0.7, 'FaceColor', channelColor(i_ch), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.75);
        scatterWithJitter(i_ch, vals, [0.05 0.05 0.05]);
    end
    set(gca, 'XTick', 1:numel(channels), 'XTickLabel', channels, ...
        'TickLabelInterpreter', 'none', 'Box', 'off', 'LineWidth', 1);
    ylabel('Cluster count per cell');
    title('Cluster number');

    nexttile;
    hold on;
    for i_ch = 1:numel(channels)
        vals = cluster_metrics.Volume_um3(cluster_metrics.Channel == channels(i_ch));
        vals = vals(isfinite(vals) & vals > 0);
        if ~isempty(vals)
            histogram(vals, 'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
                'LineWidth', 1.5, 'EdgeColor', channelColor(i_ch));
        end
    end
    xlabel('Cluster volume (um^3)');
    ylabel('Probability');
    title('Cluster size');
    legend(channels, 'Interpreter', 'none', 'Box', 'off');
    set(gca, 'Box', 'off', 'LineWidth', 1);

    nexttile;
    hold on;
    for i_ch = 1:numel(channels)
        vals = cluster_metrics.IntegratedIntensity_VolumeTimesMean(cluster_metrics.Channel == channels(i_ch));
        vals = vals(isfinite(vals) & vals > 0);
        if ~isempty(vals)
            histogram(log10(vals), 'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
                'LineWidth', 1.5, 'EdgeColor', channelColor(i_ch));
        end
    end
    xlabel('log10(Volume * MeanIntensity)');
    ylabel('Probability');
    title('Cluster integrated intensity');
    legend(channels, 'Interpreter', 'none', 'Box', 'off');
    set(gca, 'Box', 'off', 'LineWidth', 1);

    nexttile;
    hold on;
    for i_ch = 1:numel(channels)
        vals = cluster_metrics.EquiRadius_um(cluster_metrics.Channel == channels(i_ch));
        vals = vals(isfinite(vals) & vals > 0);
        if ~isempty(vals)
            histogram(vals, 'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
                'LineWidth', 1.5, 'EdgeColor', channelColor(i_ch));
        end
    end
    xlabel('Equivalent radius (um)');
    ylabel('Probability');
    title('Equivalent radius');
    legend(channels, 'Interpreter', 'none', 'Box', 'off');
    set(gca, 'Box', 'off', 'LineWidth', 1);

    saveFigure(fig, analysis_out_dir, 'cluster_metrics_summary', save_plot_formats);
    close(fig);
end

function plotActualRandomPair(actual_vals, random_vals, y_label_text)
    hold on;
    actual_vals = actual_vals(:);
    random_vals = random_vals(:);
    bar(1, nanMean(actual_vals), 0.7, 'FaceColor', [0.20 0.40 0.80], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.75);
    bar(2, nanMean(random_vals), 0.7, 'FaceColor', [0.75 0.75 0.75], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9);
    for i = 1:numel(actual_vals)
        if isfinite(actual_vals(i)) && isfinite(random_vals(i))
            plot([1 2], [actual_vals(i) random_vals(i)], '-', 'Color', [0.65 0.65 0.65]);
        end
    end
    scatterWithJitter(1, actual_vals, [0.02 0.12 0.35]);
    scatterWithJitter(2, random_vals, [0.20 0.20 0.20]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Actual', 'Randomized'}, ...
        'Box', 'off', 'LineWidth', 1);
    ylabel(y_label_text);
end

function scatterWithJitter(x_center, vals, marker_color)
    vals = vals(:);
    vals = vals(isfinite(vals));
    if isempty(vals)
        return;
    end
    if isscalar(vals)
        jitter = 0;
    else
        jitter = linspace(-0.16, 0.16, numel(vals))';
    end
    scatter(x_center + jitter, vals, 24, 'filled', ...
        'MarkerFaceColor', marker_color, 'MarkerFaceAlpha', 0.65, ...
        'MarkerEdgeColor', 'none');
end

function color = channelColor(i_ch)
    palette = [0.85 0.30 0.20; 0.20 0.40 0.80; 0.20 0.60 0.35; 0.55 0.35 0.75];
    color = palette(mod(i_ch - 1, size(palette, 1)) + 1, :);
end

function saveFigure(fig, analysis_out_dir, base_name, save_plot_formats)
    for i_fmt = 1:numel(save_plot_formats)
        fmt = lower(save_plot_formats(i_fmt));
        out_path = fullfile(analysis_out_dir, base_name + "." + fmt);
        if fmt == "fig"
            savefig(fig, char(out_path));
        elseif fmt == "pdf"
            exportgraphics(fig, char(out_path), 'ContentType', 'vector');
        else
            exportgraphics(fig, char(out_path), 'Resolution', 300);
        end
    end
end

function out = appendTable(base_table, add_table)
    if width(base_table) == 0
        out = add_table;
    else
        out = [base_table; add_table];
    end
end

function value = safeDivide(num, den)
    if den == 0
        value = NaN;
    else
        value = double(num) ./ double(den);
    end
end

function p_value = rightTailPermutationP(random_values, actual_value)
    random_values = random_values(isfinite(random_values));
    if isempty(random_values) || ~isfinite(actual_value)
        p_value = NaN;
    else
        p_value = (1 + sum(random_values >= actual_value)) ./ (numel(random_values) + 1);
    end
end

function value = nanMean(x)
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = mean(x);
    end
end

function value = nanMedian(x)
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = median(x);
    end
end

function value = nanStd(x)
    x = x(isfinite(x));
    if numel(x) < 2
        value = NaN;
    else
        value = std(x);
    end
end

function value = nanPrctile(x, p)
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = prctile(x, p);
    end
end
