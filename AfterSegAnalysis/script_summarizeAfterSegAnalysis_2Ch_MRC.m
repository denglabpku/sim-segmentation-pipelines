% Merge and summarize outputs from multiple script_afterSegAnalysis_2Ch_MRC.m runs.
%
% Use this after running script_afterSegAnalysis_2Ch_MRC.m separately on
% multiple biological replicate datasets. The script merges the exported CSV
% tables, adds replicate metadata, writes replicate-level summary CSV files,
% and generates combined plots similar to the per-dataset analysis script.

clear; clc; close all;

analysis_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(analysis_script_dir);
addpath(genpath(repo_root), '-end');

%% User settings
% Add one root_dir per biological replicate. Each root is expected to contain:
% cropped_aligned_imgs\AfterSegAnalysis\*.csv
root_dirs = [
    % "E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned\20260415_8N_JF549-Flag647",...
    % "E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned\20260417_8N_JF549-Flag647",...
    % "E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned\20260429_8N_JF549-Flag647",...
    "E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned\20260417_11_JF549-Flag647"
    ];

% Optional readable names. Leave empty to use the final folder name of each root_dir.
replicate_names = strings(0, 1);

% If your AfterSegAnalysis folders are not under root_dir\cropped_aligned_imgs,
% set analysis_dirs directly and leave root_dirs empty.
analysis_dirs = strings(0, 1);

analysis_relative_dir = fullfile('cropped_aligned_imgs', 'AfterSegAnalysis');
summary_out_dir = fullfile('E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned', 'MergedSummary_2Ch_MRC_11_JF549-Flag647');
save_plot_formats = ["png", "pdf"];
merge_random_control_tables = true;

% Optional environment overrides for batch testing or reuse:
%   AFTER_SEG_SUMMARY_ROOT_DIRS     semicolon-separated root_dir list
%   AFTER_SEG_SUMMARY_ANALYSIS_DIRS semicolon-separated AfterSegAnalysis dirs
%   AFTER_SEG_SUMMARY_OUT_DIR       output folder
root_dirs = applyStringListOverride(root_dirs, 'AFTER_SEG_SUMMARY_ROOT_DIRS');
analysis_dirs = applyStringListOverride(analysis_dirs, 'AFTER_SEG_SUMMARY_ANALYSIS_DIRS');
summary_out_override = getenv('AFTER_SEG_SUMMARY_OUT_DIR');
if ~isempty(summary_out_override)
    summary_out_dir = summary_out_override;
end

if ~exist(summary_out_dir, 'dir')
    mkdir(summary_out_dir);
end

datasets = buildDatasetList(root_dirs, replicate_names, analysis_dirs, analysis_relative_dir);
if isempty(datasets)
    error(['No completed AfterSegAnalysis folders were found. ', ...
        'Edit root_dirs or analysis_dirs in this script.']);
end
fprintf('Found %d AfterSegAnalysis folders to merge.\n', numel(datasets));

%% Merge exported tables
merged_cluster_metrics = table();
merged_cluster_summary = table();
merged_nearest_neighbors = table();
merged_nearest_summary = table();
merged_nearest_random = table();
merged_overlap_summary = table();
merged_overlap_random = table();

for i_dataset = 1:numel(datasets)
    dataset = datasets(i_dataset);
    fprintf('Merging %s (%d/%d)\n', dataset.DatasetName, i_dataset, numel(datasets));

    merged_cluster_metrics = appendCsvWithDataset(merged_cluster_metrics, dataset, ...
        'cluster_metrics_per_cluster.csv', true);
    merged_cluster_summary = appendCsvWithDataset(merged_cluster_summary, dataset, ...
        'cluster_metrics_per_cell_summary.csv', true);
    merged_nearest_neighbors = appendCsvWithDataset(merged_nearest_neighbors, dataset, ...
        'nearest_neighbor_distances_per_cluster.csv', true);
    merged_nearest_summary = appendCsvWithDataset(merged_nearest_summary, dataset, ...
        'nearest_neighbor_summary_per_cell.csv', true);
    merged_overlap_summary = appendCsvWithDataset(merged_overlap_summary, dataset, ...
        'voxel_overlap_summary_per_cell.csv', true);

    if merge_random_control_tables
        merged_nearest_random = appendCsvWithDataset(merged_nearest_random, dataset, ...
            'nearest_neighbor_random_controls.csv', false);
        merged_overlap_random = appendCsvWithDataset(merged_overlap_random, dataset, ...
            'voxel_overlap_random_controls.csv', false);
    end
end

if height(merged_cluster_metrics) == 0 && height(merged_nearest_summary) == 0 && ...
        height(merged_overlap_summary) == 0
    error('No readable CSV result tables were found in the selected AfterSegAnalysis folders.');
end

%% Replicate-level summaries
nearest_summary_by_replicate = summarizeNumericByGroup(merged_nearest_summary, ...
    ["DatasetID", "DatasetName", "Direction", "SourceChannel", "TargetChannel"], ...
    ["ActualMedianDistance_um", "ActualMeanDistance_um", ...
    "RandomMedianDistanceMean_um", "ActualMinusRandomMedianDistance_um"]);

overlap_summary_by_replicate = summarizeNumericByGroup(merged_overlap_summary, ...
    ["DatasetID", "DatasetName", "Channel1", "Channel2"], ...
    ["OverlapFractionOfCh1", "OverlapFractionOfCh2", "JaccardIndex", ...
    "DiceCoefficient", "RandomOverlapFractionOfCh1Mean", ...
    "RandomOverlapFractionOfCh2Mean", "RandomJaccardIndexMean", ...
    "RandomDiceCoefficientMean", "ActualMinusRandomJaccardIndex"]);

cluster_summary_by_replicate = summarizeNumericByGroup(merged_cluster_summary, ...
    ["DatasetID", "DatasetName", "Channel"], ...
    ["ClusterCount", "MeanVolume_um3", "MedianVolume_um3", ...
    "MeanIntensity", "MeanIntegratedIntensity", "MedianIntegratedIntensity"]);

%% Save merged tables
writeIfNotEmpty(merged_cluster_metrics, summary_out_dir, 'merged_cluster_metrics_per_cluster.csv');
writeIfNotEmpty(merged_cluster_summary, summary_out_dir, 'merged_cluster_metrics_per_cell_summary.csv');
writeIfNotEmpty(merged_nearest_neighbors, summary_out_dir, 'merged_nearest_neighbor_distances_per_cluster.csv');
writeIfNotEmpty(merged_nearest_summary, summary_out_dir, 'merged_nearest_neighbor_summary_per_cell.csv');
writeIfNotEmpty(merged_overlap_summary, summary_out_dir, 'merged_voxel_overlap_summary_per_cell.csv');
writeIfNotEmpty(nearest_summary_by_replicate, summary_out_dir, 'summary_nearest_neighbor_by_replicate.csv');
writeIfNotEmpty(overlap_summary_by_replicate, summary_out_dir, 'summary_voxel_overlap_by_replicate.csv');
writeIfNotEmpty(cluster_summary_by_replicate, summary_out_dir, 'summary_cluster_metrics_by_replicate.csv');

if merge_random_control_tables
    writeIfNotEmpty(merged_nearest_random, summary_out_dir, 'merged_nearest_neighbor_random_controls.csv');
    writeIfNotEmpty(merged_overlap_random, summary_out_dir, 'merged_voxel_overlap_random_controls.csv');
end

save(fullfile(summary_out_dir, 'merged_afterSegAnalysis_workspace.mat'), ...
    'datasets', 'merged_cluster_metrics', 'merged_cluster_summary', ...
    'merged_nearest_neighbors', 'merged_nearest_summary', 'merged_nearest_random', ...
    'merged_overlap_summary', 'merged_overlap_random', ...
    'nearest_summary_by_replicate', 'overlap_summary_by_replicate', ...
    'cluster_summary_by_replicate', '-v7.3');

makeMergedSummaryPlots(merged_nearest_summary, merged_nearest_neighbors, ...
    merged_overlap_summary, merged_cluster_metrics, merged_cluster_summary, ...
    nearest_summary_by_replicate, overlap_summary_by_replicate, ...
    cluster_summary_by_replicate, summary_out_dir, save_plot_formats);

fprintf('Merged summary complete. Results saved to:\n%s\n', summary_out_dir);

%% Local functions
function values = applyStringListOverride(values, env_name)
    override = getenv(env_name);
    if isempty(override)
        return;
    end
    parts = string(strsplit(override, ';'));
    parts = strtrim(parts(:));
    values = parts(parts ~= "");
end

function datasets = buildDatasetList(root_dirs, replicate_names, analysis_dirs, analysis_relative_dir)
    datasets = struct('DatasetID', {}, 'DatasetName', {}, 'DatasetRoot', {}, 'AnalysisDir', {});

    if ~isempty(analysis_dirs)
        analysis_dirs = string(analysis_dirs(:));
        for i = 1:numel(analysis_dirs)
        if ~isfolder(analysis_dirs(i))
            warning('Skipping missing analysis folder: %s', analysis_dirs(i));
            continue;
        end
            [dataset_root, default_name] = inferDatasetRootAndNameFromAnalysisDir(analysis_dirs(i));
            datasets(end+1) = makeDataset(numel(datasets) + 1, default_name, ...
                string(dataset_root), analysis_dirs(i)); %#ok<AGROW>
        end
        return;
    end

    root_dirs = string(root_dirs(:));
    root_dirs = root_dirs(root_dirs ~= "");
    for i = 1:numel(root_dirs)
        root_dir = root_dirs(i);
        analysis_dir = string(fullfile(char(root_dir), analysis_relative_dir));
        if ~isfolder(analysis_dir)
            warning('Skipping root_dir because its AfterSegAnalysis folder is missing: %s', root_dir);
            continue;
        end

        if ~isempty(replicate_names) && numel(replicate_names) >= i && replicate_names(i) ~= ""
            dataset_name = replicate_names(i);
        else
            [~, dataset_name_char] = fileparts(char(root_dir));
            dataset_name = string(dataset_name_char);
        end
        datasets(end+1) = makeDataset(numel(datasets) + 1, dataset_name, ...
            root_dir, analysis_dir); %#ok<AGROW>
    end
end

function [dataset_root, dataset_name] = inferDatasetRootAndNameFromAnalysisDir(analysis_dir)
    analysis_parent = fileparts(char(analysis_dir));
    [candidate_root, parent_name] = fileparts(analysis_parent);
    if strcmpi(parent_name, 'cropped_aligned_imgs')
        dataset_root = candidate_root;
        [~, dataset_name] = fileparts(dataset_root);
    else
        dataset_root = analysis_parent;
        dataset_name = parent_name;
    end
end

function dataset = makeDataset(dataset_id, dataset_name, dataset_root, analysis_dir)
    dataset = struct();
    dataset.DatasetID = dataset_id;
    dataset.DatasetName = string(dataset_name);
    dataset.DatasetRoot = string(dataset_root);
    dataset.AnalysisDir = string(analysis_dir);
end

function merged_table = appendCsvWithDataset(merged_table, dataset, csv_name, warn_if_missing)
    csv_path = fullfile(char(dataset.AnalysisDir), csv_name);
    if ~isfile(csv_path)
        if warn_if_missing
            warning('Missing expected CSV for %s: %s', dataset.DatasetName, csv_name);
        end
        return;
    end

    tbl = readCsvTable(csv_path);
    if height(tbl) == 0
        return;
    end
    tbl = addDatasetColumns(tbl, dataset);
    merged_table = appendTableUnion(merged_table, tbl);
end

function tbl = readCsvTable(csv_path)
    try
        tbl = readtable(csv_path, 'TextType', 'string');
    catch
        tbl = readtable(csv_path);
        tbl = convertTextColumnsToString(tbl);
    end
end

function tbl = convertTextColumnsToString(tbl)
    for i_var = 1:width(tbl)
        values = tbl.(i_var);
        if iscellstr(values) || ischar(values) %#ok<ISCLSTR>
            tbl.(i_var) = string(values);
        end
    end
end

function tbl = addDatasetColumns(tbl, dataset)
    n = height(tbl);
    tbl.DatasetID = repmat(dataset.DatasetID, n, 1);
    tbl.DatasetName = repmat(dataset.DatasetName, n, 1);
    tbl.DatasetRoot = repmat(dataset.DatasetRoot, n, 1);
    tbl.AnalysisDir = repmat(dataset.AnalysisDir, n, 1);

    if ismember('CellID', tbl.Properties.VariableNames)
        tbl.CellID = string(tbl.CellID);
        tbl.DatasetCellID = dataset.DatasetName + "::" + tbl.CellID;
        tbl = movevars(tbl, {'DatasetID', 'DatasetName', 'DatasetRoot', ...
            'AnalysisDir', 'DatasetCellID'}, 'Before', 'CellID');
    else
        tbl = movevars(tbl, {'DatasetID', 'DatasetName', 'DatasetRoot', ...
            'AnalysisDir'}, 'Before', 1);
    end
end

function summary_table = summarizeNumericByGroup(tbl, group_vars, value_vars)
    summary_table = table();
    if height(tbl) == 0
        return;
    end

    group_vars = group_vars(ismember(group_vars, string(tbl.Properties.VariableNames)));
    value_vars = value_vars(ismember(value_vars, string(tbl.Properties.VariableNames)));
    if isempty(group_vars) || isempty(value_vars)
        return;
    end

    group_key = makeGroupKey(tbl, group_vars);
    unique_keys = unique(group_key, 'stable');

    for i_group = 1:numel(unique_keys)
        idx = group_key == unique_keys(i_group);
        first_row = find(idx, 1, 'first');
        row = table();

        for i_var = 1:numel(group_vars)
            var_name = char(group_vars(i_var));
            row.(var_name) = tbl.(var_name)(first_row);
        end
        row.NCells = sum(idx);

        for i_var = 1:numel(value_vars)
            var_name = char(value_vars(i_var));
            vals = double(tbl.(var_name)(idx));
            row.("Mean_" + var_name) = nanMean(vals);
            row.("Median_" + var_name) = nanMedian(vals);
            row.("Std_" + var_name) = nanStd(vals);
        end

        summary_table = appendTableUnion(summary_table, row);
    end
end

function key = makeGroupKey(tbl, group_vars)
    key = strings(height(tbl), 1);
    separator = char(30);
    for i_var = 1:numel(group_vars)
        values = string(tbl.(char(group_vars(i_var))));
        key = key + separator + values;
    end
end

function writeIfNotEmpty(tbl, out_dir, file_name)
    if height(tbl) == 0
        return;
    end
    writetable(tbl, fullfile(out_dir, file_name));
end

function makeMergedSummaryPlots(nn_summary, nn_rows, overlap_summary, ...
        cluster_metrics, cluster_summary, nn_by_rep, overlap_by_rep, ...
        cluster_by_rep, out_dir, save_plot_formats)
    if height(nn_summary) > 0
        plotNearestNeighborCellBars(nn_summary, out_dir, save_plot_formats);
    end
    if height(nn_by_rep) > 0
        plotNearestNeighborReplicateBars(nn_by_rep, out_dir, save_plot_formats);
    end
    if height(nn_rows) > 0
        plotNearestNeighborHistograms(nn_rows, out_dir, save_plot_formats);
    end
    if height(overlap_summary) > 0
        plotOverlapCellSummary(overlap_summary, out_dir, save_plot_formats);
    end
    if height(overlap_by_rep) > 0
        plotOverlapReplicateSummary(overlap_by_rep, out_dir, save_plot_formats);
    end
    if height(cluster_metrics) > 0
        plotClusterMetrics(cluster_metrics, cluster_summary, out_dir, save_plot_formats);
    end
    if height(cluster_by_rep) > 0
        plotClusterMetricsByReplicate(cluster_by_rep, out_dir, save_plot_formats);
    end
end

function plotNearestNeighborCellBars(nn_summary, out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 900 430]);
    directions = unique(nn_summary.Direction, 'stable');
    group_labels = strings(1, numel(directions));

    hold on;
    for i_dir = 1:numel(directions)
        idx = nn_summary.Direction == directions(i_dir);
        actual_vals = nn_summary.ActualMedianDistance_um(idx);
        random_vals = nn_summary.RandomMedianDistanceMean_um(idx);
        x_actual = (i_dir - 1) * 3 + 1;
        x_random = (i_dir - 1) * 3 + 2;

        bar(x_actual, nanMean(actual_vals), 0.8, 'FaceColor', [0.20 0.40 0.80], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.75);
        bar(x_random, nanMean(random_vals), 0.8, 'FaceColor', [0.75 0.75 0.75], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.9);
        scatterWithJitter(x_actual, actual_vals, [0.02 0.12 0.35]);
        scatterWithJitter(x_random, random_vals, [0.20 0.20 0.20]);
        group_labels(i_dir) = directions(i_dir);
    end
    ylabel('Cell median nearest-neighbor distance (um)');
    set(gca, 'XTick', (0:numel(directions)-1) * 3 + 1.5, 'XTickLabel', group_labels, ...
        'TickLabelInterpreter', 'none', 'Box', 'off', 'LineWidth', 1);
    legend({'Actual', 'Randomized'}, 'Location', 'best', 'Box', 'off');
    title('Merged weighted-centroid nearest-neighbor distances');
    saveFigure(fig, out_dir, 'merged_nearest_neighbor_actual_vs_random_by_cell', save_plot_formats);
    close(fig);
end

function plotNearestNeighborReplicateBars(nn_by_rep, out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 900 430]);
    directions = unique(nn_by_rep.Direction, 'stable');
    group_labels = strings(1, numel(directions));

    hold on;
    for i_dir = 1:numel(directions)
        idx = nn_by_rep.Direction == directions(i_dir);
        actual_vals = nn_by_rep.Mean_ActualMedianDistance_um(idx);
        random_vals = nn_by_rep.Mean_RandomMedianDistanceMean_um(idx);
        x_actual = (i_dir - 1) * 3 + 1;
        x_random = (i_dir - 1) * 3 + 2;

        bar(x_actual, nanMean(actual_vals), 0.8, 'FaceColor', [0.20 0.40 0.80], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.75);
        bar(x_random, nanMean(random_vals), 0.8, 'FaceColor', [0.75 0.75 0.75], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.9);
        plotPairedValues(x_actual, x_random, actual_vals, random_vals);
        scatterWithJitter(x_actual, actual_vals, [0.02 0.12 0.35]);
        scatterWithJitter(x_random, random_vals, [0.20 0.20 0.20]);
        group_labels(i_dir) = directions(i_dir);
    end
    ylabel('Replicate mean of cell median distance (um)');
    set(gca, 'XTick', (0:numel(directions)-1) * 3 + 1.5, 'XTickLabel', group_labels, ...
        'TickLabelInterpreter', 'none', 'Box', 'off', 'LineWidth', 1);
    legend({'Actual', 'Randomized'}, 'Location', 'best', 'Box', 'off');
    title('Nearest-neighbor distances by biological replicate');
    saveFigure(fig, out_dir, 'merged_nearest_neighbor_actual_vs_random_by_replicate', save_plot_formats);
    close(fig);
end

function plotNearestNeighborHistograms(nn_rows, out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 900 430]);
    directions = unique(nn_rows.Direction, 'stable');
    tiledlayout(1, numel(directions), 'TileSpacing', 'compact', 'Padding', 'compact');
    for i_dir = 1:numel(directions)
        nexttile;
        idx = nn_rows.Direction == directions(i_dir);
        vals = nn_rows.NearestDistance_um(idx);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            histogram(vals, 'Normalization', 'probability', ...
                'FaceColor', [0.20 0.40 0.80], 'EdgeColor', 'none');
        end
        xlabel('Nearest distance (um)');
        ylabel('Probability');
        title(directions(i_dir), 'Interpreter', 'none');
        set(gca, 'Box', 'off', 'LineWidth', 1);
    end
    saveFigure(fig, out_dir, 'merged_nearest_neighbor_distance_histograms', save_plot_formats);
    close(fig);
end

function plotOverlapCellSummary(overlap_summary, out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 1050 380]);
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plotActualRandomPair(overlap_summary.JaccardIndex, ...
        overlap_summary.RandomJaccardIndexMean, 'Jaccard index');
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

    saveFigure(fig, out_dir, 'merged_voxel_overlap_actual_vs_random_by_cell', save_plot_formats);
    close(fig);
end

function plotOverlapReplicateSummary(overlap_by_rep, out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 1050 380]);
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plotActualRandomPair(overlap_by_rep.Mean_JaccardIndex, ...
        overlap_by_rep.Mean_RandomJaccardIndexMean, 'Replicate mean Jaccard index');
    title('Voxel overlap');

    nexttile;
    plotActualRandomPair(overlap_by_rep.Mean_OverlapFractionOfCh1, ...
        overlap_by_rep.Mean_RandomOverlapFractionOfCh1Mean, ...
        sprintf('Replicate mean fraction of %s', overlap_by_rep.Channel1(1)));
    title('Overlap over channel 1');

    nexttile;
    plotActualRandomPair(overlap_by_rep.Mean_OverlapFractionOfCh2, ...
        overlap_by_rep.Mean_RandomOverlapFractionOfCh2Mean, ...
        sprintf('Replicate mean fraction of %s', overlap_by_rep.Channel2(1)));
    title('Overlap over channel 2');

    saveFigure(fig, out_dir, 'merged_voxel_overlap_actual_vs_random_by_replicate', save_plot_formats);
    close(fig);
end

function plotClusterMetrics(cluster_metrics, cluster_summary, out_dir, save_plot_formats)
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
    plotMetricHistograms(cluster_metrics, channels, 'Volume_um3', ...
        'Cluster volume (um^3)', 'Cluster size');

    nexttile;
    vals_var = 'IntegratedIntensity_VolumeTimesMean';
    hold on;
    for i_ch = 1:numel(channels)
        vals = cluster_metrics.(vals_var)(cluster_metrics.Channel == channels(i_ch));
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
    plotMetricHistograms(cluster_metrics, channels, 'EquiRadius_um', ...
        'Equivalent radius (um)', 'Equivalent radius');

    saveFigure(fig, out_dir, 'merged_cluster_metrics_summary_by_cell', save_plot_formats);
    close(fig);
end

function plotClusterMetricsByReplicate(cluster_by_rep, out_dir, save_plot_formats)
    fig = figure('Color', 'w', 'Position', [100 100 1050 420]);
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    channels = unique(cluster_by_rep.Channel, 'stable');

    nexttile;
    plotReplicateChannelBars(cluster_by_rep, channels, 'Mean_ClusterCount', ...
        'Replicate mean cluster count');
    title('Cluster number');

    nexttile;
    plotReplicateChannelBars(cluster_by_rep, channels, 'Mean_MeanVolume_um3', ...
        'Replicate mean cluster volume (um^3)');
    title('Cluster size');

    nexttile;
    plotReplicateChannelBars(cluster_by_rep, channels, 'Mean_MeanIntegratedIntensity', ...
        'Replicate mean integrated intensity');
    title('Cluster intensity');

    saveFigure(fig, out_dir, 'merged_cluster_metrics_summary_by_replicate', save_plot_formats);
    close(fig);
end

function plotMetricHistograms(tbl, channels, var_name, x_label_text, title_text)
    hold on;
    for i_ch = 1:numel(channels)
        vals = tbl.(var_name)(tbl.Channel == channels(i_ch));
        vals = vals(isfinite(vals) & vals > 0);
        if ~isempty(vals)
            histogram(vals, 'Normalization', 'probability', 'DisplayStyle', 'stairs', ...
                'LineWidth', 1.5, 'EdgeColor', channelColor(i_ch));
        end
    end
    xlabel(x_label_text);
    ylabel('Probability');
    title(title_text);
    legend(channels, 'Interpreter', 'none', 'Box', 'off');
    set(gca, 'Box', 'off', 'LineWidth', 1);
end

function plotReplicateChannelBars(tbl, channels, var_name, y_label_text)
    hold on;
    for i_ch = 1:numel(channels)
        idx = tbl.Channel == channels(i_ch);
        vals = tbl.(var_name)(idx);
        bar(i_ch, nanMean(vals), 0.7, 'FaceColor', channelColor(i_ch), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.75);
        scatterWithJitter(i_ch, vals, [0.05 0.05 0.05]);
    end
    set(gca, 'XTick', 1:numel(channels), 'XTickLabel', channels, ...
        'TickLabelInterpreter', 'none', 'Box', 'off', 'LineWidth', 1);
    ylabel(y_label_text);
end

function plotActualRandomPair(actual_vals, random_vals, y_label_text)
    hold on;
    actual_vals = actual_vals(:);
    random_vals = random_vals(:);
    bar(1, nanMean(actual_vals), 0.7, 'FaceColor', [0.20 0.40 0.80], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.75);
    bar(2, nanMean(random_vals), 0.7, 'FaceColor', [0.75 0.75 0.75], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9);
    plotPairedValues(1, 2, actual_vals, random_vals);
    scatterWithJitter(1, actual_vals, [0.02 0.12 0.35]);
    scatterWithJitter(2, random_vals, [0.20 0.20 0.20]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Actual', 'Randomized'}, ...
        'Box', 'off', 'LineWidth', 1);
    ylabel(y_label_text);
end

function plotPairedValues(x1, x2, vals1, vals2)
    for i = 1:min(numel(vals1), numel(vals2))
        if isfinite(vals1(i)) && isfinite(vals2(i))
            plot([x1 x2], [vals1(i) vals2(i)], '-', 'Color', [0.65 0.65 0.65]);
        end
    end
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

function saveFigure(fig, out_dir, base_name, save_plot_formats)
    for i_fmt = 1:numel(save_plot_formats)
        fmt = lower(save_plot_formats(i_fmt));
        out_path = fullfile(out_dir, base_name + "." + fmt);
        if fmt == "fig"
            savefig(fig, char(out_path));
        elseif fmt == "pdf"
            exportgraphics(fig, char(out_path), 'ContentType', 'vector');
        else
            exportgraphics(fig, char(out_path), 'Resolution', 300);
        end
    end
end

function out = appendTableUnion(base_table, add_table)
    if width(base_table) == 0
        out = add_table;
        return;
    end
    if width(add_table) == 0
        out = base_table;
        return;
    end

    base_vars = string(base_table.Properties.VariableNames);
    add_vars = string(add_table.Properties.VariableNames);
    all_vars = unique([base_vars, add_vars], 'stable');

    base_table = addMissingVariables(base_table, all_vars);
    add_table = addMissingVariables(add_table, all_vars);
    add_table = add_table(:, base_table.Properties.VariableNames);
    out = [base_table; add_table];
end

function tbl = addMissingVariables(tbl, all_vars)
    existing_vars = string(tbl.Properties.VariableNames);
    for i_var = 1:numel(all_vars)
        var_name = char(all_vars(i_var));
        if ~ismember(all_vars(i_var), existing_vars)
            tbl.(var_name) = missingColumn(height(tbl));
        end
    end
    tbl = tbl(:, cellstr(all_vars));
end

function values = missingColumn(n)
    values = strings(n, 1);
    values(:) = missing;
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
