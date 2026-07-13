% Scripts for visualization of 3D segmentation results
clc;clear;close all;

%% Cluster number
root_dir = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20251106_SIM_live_JF549_HaloRPB1';
sub_dir = [
    "Exp1_001_20251106_141854",...
    "Exp1_002_20251106_143014",...
    "Exp1_003_20251106_143342",...
    "Exp1_004_20251106_143649",...
    "Exp1_005_20251106_144119",...
    "Exp1_006_20251106_144811",...
    ];

% root_dir = 'D:\Denglab-DataCenter\Qinhua_WangBo_enhan_promo_project\HIS-SIM\20251211_qinhua_OCT4_BRD4_BD555\processed\seg_result';
% sub_dir = [
%     "Sparse_22_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_23_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_23_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     "Sparse_24_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_24_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     "Sparse_25_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_26_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_26_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     "Sparse_27_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_27_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     "Sparse_28_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_28_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     "Sparse_29_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_29_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     "Sparse_30_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_30_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     "Sparse_31_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_31_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     "Sparse_33_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_33_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     "Sparse_34_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
%     "Sparse_34_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
%     ];

% merge_tab_stats3D = table();
% for idir = 1:length(sub_dir)
%     seg_result_dir = fullfile(root_dir,sub_dir(idir));
%     seg_result_files = dir(fullfile(seg_result_dir,'tiff_new\seg_result','*_segResults.mat')); 
%     for ii = 1:length(seg_result_files)
%         t = 1;
%         load(fullfile(seg_result_files(ii).folder,seg_result_files(ii).name),'stats3D_table');
%     end
% end

name_ch = ["SIM561"]; iCh = 1;
merge_tab_stats3D = table();   % initialize empty table
cell_counter = 0;              % counter for cell_idx

for idir = 1:length(sub_dir)
    seg_result_dir = fullfile(root_dir, sub_dir(idir));
    seg_result_files = dir(fullfile(seg_result_dir, 'tiff_new', 'seg_result', '*_segResults.mat')); 
    % seg_result_files = dir(fullfile(seg_result_dir, '*_segResults.mat')); 
    temp_Filenames = {seg_result_files.name}; idx = contains(temp_Filenames,name_ch(iCh)); seg_result_files = seg_result_files(idx); 
    
    for ii = 1:length(seg_result_files)
        % load stats3D_table from file
        load(fullfile(seg_result_files(ii).folder, seg_result_files(ii).name), 'stats3D_table');
        
        % increment cell counter
        cell_counter = cell_counter + 1;
        
        % add new column cell_idx
        stats3D_table.cell_idx = repmat(cell_counter, height(stats3D_table), 1);
        
        % append to merged table
        merge_tab_stats3D = [merge_tab_stats3D; stats3D_table];
    end
end



fig3 = figure;
tiledlayout(3,1);
ax = nexttile;
histogram(merge_tab_stats3D.Volume,'Normalization','probability');
xlabel('Voxels'); ylabel('Counts');ax.Box = 'off';
nexttile;
histogram(merge_tab_stats3D.VolumePhysical,'Normalization','probability');
xlabel('\mum^3'); ylabel('Counts');
merge_tab_stats3D.equiRadius = ((3*merge_tab_stats3D.VolumePhysical)/(4*pi)).^(1/3);
nexttile;
histogram(merge_tab_stats3D.equiRadius,'Normalization','probability');
xlabel('\mum'); ylabel('Counts');

counts = groupcounts(merge_tab_stats3D.cell_idx);
mean_y = mean(counts);
std_y = std(counts);
x = categorical(name_ch(iCh));

% ================== Cluster Number Counts ================== %
fig = figure('Position', [711 415 560 368]); hold all;
bar(x, mean_y, 'FaceColor','none','EdgeColor','k','LineWidth',2);  
ylabel('Cluster # per cell', 'FontSize', 20); 

% Overlay with error bars
errorbar(x, mean_y, std_y,  '.', 'color', 'k', 'linewidth', 1); 

% Overlay data points
hold on;
spread = 0.5; % Spread factor for jitter
allData = counts; % Extract data points for current condition
plot(rand(size(allData)) * spread - (spread / 2) + 1, allData, 'k.', 'MarkerSize', 5);

% adjust format
ax = gca; ax.FontSize = 12; ax.LineWidth = 0.5; ax.XAxis.TickLabelRotation = 45;  ax.TickLabelInterpreter = 'none';  ax.XColor = 'k'; ax.YColor = 'k';   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
ax.YLim = [0, 3500]; 
ax.PlotBoxAspectRatio = [1 2.5 1];

% save_path = 'D:\Denglab-DataCenter\ChenFei_INTAC_temp\20251106_SIM_live_JF549_HaloRPB1\figures';
% save figure 
% save_path = fullfile(root_dir,'figures'); mkdir(save_path);
% exportgraphics(fig, fullfile(save_path,sprintf('Cluster_numbers_perCell_Ch%s.pdf',name_ch(iCh))));

%% Radius
clc;clear;close all;
root_dir = 'D:\Denglab-DataCenter\Qinhua_WangBo_enhan_promo_project\HIS-SIM\20251211_qinhua_OCT4_BRD4_BD555\processed\seg_result';
sub_dir = [
    "Sparse_22_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_23_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_23_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_24_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_24_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_25_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_26_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_26_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_27_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_27_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_28_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_28_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_29_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_29_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_30_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_30_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_31_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_31_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_33_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_33_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_34_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_34_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    ];

merge_tab_stats3D = {};   % initialize empty table

name_ch = ["C488", "C638"];
pixel_size    = 0.065/2; % pixel size, unit : um

for iCh = 1:length(name_ch)
    
    for idir = 1:length(sub_dir)
        fprintf("Processing Ch %d, File %d ... \n",iCh, idir);

        seg_result_dir = fullfile(root_dir, sub_dir(idir));
        img_list = dir(fullfile(seg_result_dir, '*_seglabels.tif')); 
        temp_Filenames = {img_list.name}; idx = contains(temp_Filenames,name_ch(iCh)); img_list = img_list(idx);           
        
        if length(img_list) > 1; error('img_list must contain only one element.'); end
        V = tiffreadVolume(fullfile(img_list.folder, img_list.name));
        
        % cluster radius using Maximum Intensity Projection (MIP) along Z-axis
        % 1. Identify all unique cluster labels (excluding background 0)
        clusterIDs = unique(V);
        clusterIDs(clusterIDs == 0) = [];
        numClusters = numel(clusterIDs);
        
        % Initialize results array
        mipRadii = zeros(numClusters, 1);
        
        for i = 1:numClusters
            id = clusterIDs(i);
            
            % 2. Isolate the specific cluster
            tempMask = (V == id);
            
            % 3. Calculate Maximum Intensity Projection (MIP) along Z-axis
            % For a binary mask, this is equivalent to 'any' voxel existing at (x,y)
            projection = max(tempMask, [], 3);
            
            % 4. Measure properties of the 2D projection
            stats = regionprops(projection, 'Area');
        
            % % Mean radius = (Major + Minor) / 4
            % mipRadii(i) = mean([stats.MajorAxisLength, stats.MinorAxisLength]) / 2;
            
            % 5. Calculate Equivalent Radius (R = sqrt(Area/pi))
            % If a cluster has multiple disconnected parts in projection, sum their areas
            totalArea = sum([stats.Area]);
            mipRadii(i) = sqrt(totalArea / pi);
        end
               
        mipRadii_real = mipRadii*pixel_size; % unit: um
        merge_tab_stats3D{iCh,idir} = mipRadii_real;
    end
end

save(fullfile(root_dir,"result_mipRadii.mat"),"merge_tab_stats3D","sub_dir","root_dir",'-mat');

%
figure;tiledlayout(2,1);
ax1 = nexttile;
temp = vertcat(merge_tab_stats3D{1,:});
histogram(temp,'Normalization','probability','NumBins',30,'FaceColor','none');
xlabel('r_{mip} (um)');ylabel('Probability');title('OCT4 Alexa488');ax1.XLim = [0 0.25];
ax2 = nexttile;
temp = vertcat(merge_tab_stats3D{2,:});
histogram(temp,'Normalization','probability','NumBins',30,'FaceColor','none');
xlabel('r_{mip} (um)');ylabel('Probability');title('BRD4 BD555');ax2.XLim = [0 0.25];

fig = gcf;
save_path = fullfile(root_dir,'figures'); mkdir(save_path);
exportgraphics(fig, fullfile(save_path,sprintf('Cluster_radius.pdf')));

%% Colocalization
clc;clear;close all;
root_dir = 'D:\Denglab-DataCenter\Qinhua_WangBo_enhan_promo_project\HIS-SIM\20251211_qinhua_OCT4_BRD4_BD555\processed\seg_result';
sub_dir = [
    "Sparse_22_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_23_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_23_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_24_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_24_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_25_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_26_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_26_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_27_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_27_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_28_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_28_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_29_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_29_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_30_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_30_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_31_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_31_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_33_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_33_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    "Sparse_34_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell01_aligned",...
    "Sparse_34_Merge_488_Em525_561_Em609_Wiener_3DSIM-3_Oct4_488_Brd4_561_Cell02_aligned",...
    ];

%% NiCheng OCT4-Lexy

%clc;clear;close all;

root_dir = '/dataA/nicheng/NC-imaging/2026/SIM';

sub_dir = [
    "Fixed-NT",...
    "Fixed-BL",...
    ];

merge_tab_stats2D = table();   % initialize empty table
merge_tab_stats3D = table();   % initialize empty table
cell_counter = 0;              % counter for cell_idx

for idir = 1:length(sub_dir)
    seg_result_dir = fullfile(root_dir, sub_dir(idir),'processed','seg_result_directOtsuLabel');    
    seg_result_files = dir(fullfile(seg_result_dir, '*_segResults.mat')); 
    temp_Filenames = {seg_result_files.name}; % idx = contains(temp_Filenames,name_ch(iCh)); seg_result_files = seg_result_files(idx); 
    
    for ii = 1:length(seg_result_files)
        % load stats3D_table from file
        load(fullfile(seg_result_files(ii).folder, seg_result_files(ii).name), 'stats2D_table','stats3D_table');
        
        % increment cell counter
        cell_counter = cell_counter + 1;

         % add new column cell_idx
        stats2D_table.cell_idx = repmat(cell_counter, height(stats2D_table), 1);
        stats2D_table.condition = repmat(categorical(sub_dir(idir)), height(stats2D_table), 1);
        
        % add new column cell_idx
        stats3D_table.cell_idx = repmat(cell_counter, height(stats3D_table), 1);
        stats3D_table.condition = repmat(categorical(sub_dir(idir)), height(stats3D_table), 1);
       
        % append to merged table
        merge_tab_stats2D = [merge_tab_stats2D; stats2D_table];
        merge_tab_stats3D = [merge_tab_stats3D; stats3D_table];
    end
end

% ================== Cluster Number Counts by Condition ================== %
% Count clusters per cell, grouped by condition
% counts_per_cell = groupsummary(merge_tab_stats3D, {'cell_idx','condition'}, 'numel', 'cell_idx');
% Correct usage: count rows per cell_idx per condition
counts_per_cell = groupsummary(merge_tab_stats3D, {'cell_idx','condition'}, @numel);


% Extract condition labels
conditions = categories(merge_tab_stats3D.condition);

% Preallocate
mean_y = zeros(numel(conditions),1);
std_y  = zeros(numel(conditions),1);
allData = cell(numel(conditions),1);

for ic = 1:numel(conditions)
    % Select rows for this condition
    idx = counts_per_cell.condition == conditions{ic};
    data = counts_per_cell.GroupCount(idx);   % cluster counts per cell
    
    % Store stats
    mean_y(ic) = mean(data);
    std_y(ic)  = std(data);
    allData{ic} = data;
end

% Plot
fig = figure('Position', [711 415 560 368]); hold on;
x = categorical(conditions);
x = reordercats(x, conditions); % enforce your desired order

bar(x, mean_y, 'FaceColor','none','EdgeColor','k','LineWidth',2);
ylabel('Cluster # per cell', 'FontSize', 20);

% Overlay error bars
errorbar(x, mean_y, std_y, '.', 'color', 'k', 'linewidth', 1);

% Overlay jittered datapoints for each condition
spread = 0.3;
for ic = 1:numel(conditions)
    jitterX = rand(size(allData{ic})) * spread - (spread/2) + ic;
    plot(jitterX, allData{ic}, 'k.', 'MarkerSize', 5);
end

% Adjust format
ax = gca;
ax.FontSize = 12; ax.LineWidth = 0.5;
ax.XAxis.TickLabelRotation = 45;
ax.TickLabelInterpreter = 'none';
ax.XColor = 'k'; ax.YColor = 'k';
ax.PlotBoxAspectRatio = [1 2.5 1];

% save_path = 'D:\Denglab-DataCenter\ChenFei_INTAC_temp\20251106_SIM_live_JF549_HaloRPB1\figures';
% % save figure 
% exportgraphics(fig, fullfile(save_path,'Cluster_numbers_perCell.pdf'));

% ================== Cluster EquivRadius_um Counts by Condition ================== %
% Count clusters per cell, grouped by condition
% counts_per_cell = groupsummary(merge_tab_stats3D, {'cell_idx','condition'}, 'numel', 'cell_idx');
% Correct usage: count rows per cell_idx per condition
radius_per_cell = groupsummary(merge_tab_stats2D, {'cell_idx','condition'}, @mean);


% Extract condition labels
conditions = categories(merge_tab_stats2D.condition);

% Preallocate
mean_y = zeros(numel(conditions),1);
std_y  = zeros(numel(conditions),1);
allData = cell(numel(conditions),1);

for ic = 1:numel(conditions)
    % Select rows for this condition
    idx = radius_per_cell.condition == conditions{ic};
    data = radius_per_cell.fun1_EquivRadius_um(idx);   % cluster counts per cell
    
    % Store stats
    mean_y(ic) = mean(data);
    std_y(ic)  = std(data);
    allData{ic} = data;
end

% Plot
fig = figure('Position', [711 415 560 368]); hold on;
x = categorical(conditions);
x = reordercats(x, conditions); % enforce your desired order

bar(x, mean_y, 'FaceColor','none','EdgeColor','k','LineWidth',2);
ylabel('MIP Radius per cell', 'FontSize', 20);

% Overlay error bars
errorbar(x, mean_y, std_y, '.', 'color', 'k', 'linewidth', 1);

% Overlay jittered datapoints for each condition
spread = 0.3;
for ic = 1:numel(conditions)
    jitterX = rand(size(allData{ic})) * spread - (spread/2) + ic;
    plot(jitterX, allData{ic}, 'k.', 'MarkerSize', 5);
end

% Adjust format
ax = gca;
ax.FontSize = 12; ax.LineWidth = 0.5;
ax.XAxis.TickLabelRotation = 45;
ax.TickLabelInterpreter = 'none';
ax.XColor = 'k'; ax.YColor = 'k';
ax.PlotBoxAspectRatio = [1 2.5 1];
ax.YLim = [0.1 0.16];

% save_path = 'D:\Denglab-DataCenter\ChenFei_INTAC_temp\20251106_SIM_live_JF549_HaloRPB1\figures';
% % save figure 
% exportgraphics(fig, fullfile(save_path,'Cluster_numbers_perCell.pdf'));

%% Colocalization
coloc_list_name = {};
tManders_A = [];
tManders_B = [];

name_ch = ["C488", "C638"];
pixel_size    = 0.065/2; % pixel size, unit : um


    
for idir = 1:length(sub_dir)
    fprintf("Processing File %d ... \n", idir);

    seg_result_dir = fullfile(root_dir, sub_dir(idir));
    coloc_list = dir(fullfile(seg_result_dir, '*.csv')); 
    temp_Filenames = {coloc_list.name}; idx = contains(temp_Filenames,'Coloc_thres7_roi1'); coloc_list = coloc_list(idx);           
    
    if length(coloc_list) > 1; error('coloc_list must contain only one element.'); end
    coloc_tab = readtable(fullfile(coloc_list.folder, coloc_list.name));
    
    tManders_A = [tManders_A; coloc_tab{18,2}];
    tManders_B = [tManders_B; coloc_tab{19,2}];
    coloc_list_name{idir} = coloc_list.name;
end

merge_coloc_tab = table(coloc_list_name',tManders_A,tManders_B);

x = categorical(["tM_A","tM_B"]);
x = reordercats(x);
mean_y = mean(merge_coloc_tab{:,2:3},1);
std_y = std(merge_coloc_tab{:,2:3},1);
% ================== Manders Coefficient ================== %
fig = figure('Position', [711 415 560 368]); hold all;
bar(x, mean_y, 'FaceColor','none','EdgeColor','k','LineWidth',2);  
ylabel('Thresholded Manders', 'FontSize', 20); 

% Overlay with error bars
errorbar(x, mean_y, std_y,  '.', 'color', 'k', 'linewidth', 1); 

% Overlay data points
hold on;
spread = 0.5; % Spread factor for jitter
tmp = merge_coloc_tab{:,2:3};
for i = 1:length(x)    
    allData = tmp(:,i);
    plot(rand(size(allData)) * spread - (spread / 2) + i, allData, 'k.', 'MarkerSize', 5);
end

% adjust format
ax = gca; ax.FontSize = 12; ax.LineWidth = 0.5; ax.XAxis.TickLabelRotation = 45;  ax.TickLabelInterpreter = 'none';  ax.XColor = 'k'; ax.YColor = 'k';   %Relative length of each axis, specified as a three-element vector of the form [px py pz] 
ax.YLim = [0.05, 0.15]; 
ax.PlotBoxAspectRatio = [1 2.5 1];

% save figure 
save_path = fullfile(root_dir,'figures'); mkdir(save_path);
exportgraphics(fig, fullfile(save_path,sprintf('Coloc_Manders.pdf')));
