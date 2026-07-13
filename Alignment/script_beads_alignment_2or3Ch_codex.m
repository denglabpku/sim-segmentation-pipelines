% Script align beads in 3D for SIM/Confocal 3D image stacks

% This alignment is compatible to do 2 or 3 channel alignment.
% Same core as script_beads_alignment_3Ch_codex.m

% Zuhui Wang, 2025/12/04
% Optimized by Codex, 2026/06/02

% Always verify the alignment results using a independent beads image
% Best result is obtain when always fix z of beads_roi_origin_um and
% process_roi_origin_um to 0; and fit_transform_3D=true

%% ESTIMATE TRANSFORMATION MATRIX
close all; clear; clc;
addpath(genpath('C:\Users\zuhui\OneDrive - Peking University\Documents\MATLAB\SIM_microscopy\Segmentation_3DSIM'))
addpath(genpath('C:\Users\zuhui\OneDrive - Peking University\Documents\MATLAB\bfmatlab'),'-end')

clc;clear;close all;
% root_path = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned\20260415_8N_JF549-Flag647_beads';
root_path = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned\20260429_8N_JF549-Flag647_beads';

% For multiple bead image pairs (use cell arrays):
% C1 488 , C2 640, C3 560 (fixed)
% C3 channel as fixed reference
C3_beads_csv = {... % Replace with your actual file names
    % 'C2_localizations\beads_002_20260415_160214_C561_bgsub.csv',...
    % 'C2_localizations\beads_003_20260415_160526_C561_bgsub.csv',...
    % 'C2_localizations\beads_004_20260415_160828_C561_bgsub.csv',...
    % 'C2_localizations\beads_005_20260415_161113_C561_bgsub.csv',...  

    'C2_localizations\beads_010_20260429_134425_C561_bgsub.csv',...
    'C2_localizations\beads_011_20260429_134623_C561_bgsub.csv',...
    'C2_localizations\beads_012_20260429_134737_C561_bgsub.csv',...
    'C2_localizations\beads_013_20260429_134905_C561_bgsub.csv',...
    };

C1_beads_csv = {... % Replace with your actual file names
    % 'C1_localizations\beads_002_20260415_160214_C488_bgsub.csv',...
    % 'C1_localizations\beads_003_20260415_160526_C488_bgsub.csv',...
    % 'C1_localizations\beads_004_20260415_160828_C488_bgsub.csv',...
    % 'C1_localizations\beads_005_20260415_161113_C488_bgsub.csv',...  

    'C1_localizations\beads_010_20260429_134425_C488_bgsub.csv',...
    'C1_localizations\beads_011_20260429_134623_C488_bgsub.csv',...
    'C1_localizations\beads_012_20260429_134737_C488_bgsub.csv',...
    'C1_localizations\beads_013_20260429_134905_C488_bgsub.csv',...
    };

% C2 channel to align to C3
C2_beads_csv = {... % Replace with your actual file names    
    % 'C3_localizations\beads_002_20260415_160214_C640_bgsub.csv',...
    % 'C3_localizations\beads_003_20260415_160526_C640_bgsub.csv',...
    % 'C3_localizations\beads_004_20260415_160828_C640_bgsub.csv',...
    % 'C3_localizations\beads_005_20260415_161113_C640_bgsub.csv',...

    'C3_localizations\beads_010_20260429_134425_C640_bgsub.csv',...
    'C3_localizations\beads_011_20260429_134623_C640_bgsub.csv',...
    'C3_localizations\beads_012_20260429_134737_C640_bgsub.csv',...
    'C3_localizations\beads_013_20260429_134905_C640_bgsub.csv',...
    };

px = 0.065; % image pixel size, unit um
py = px;
pz = 0.135; % z step, unit um

% Check if C2 bead files are provided (optional)
if exist('C2_beads_csv', 'var') && ~isempty(C2_beads_csv)
    do_C2_alignment = true;
    % Convert to cell array if char
    if ischar(C2_beads_csv)
        C2_beads_csv = {C2_beads_csv};
    end
else
    do_C2_alignment = false;
    fprintf('Note: C2 bead files not provided, skipping C2 alignment.\n');
end

dist_thresh = 20*px;
use3D = false; % Pairing uses XY only; the transform fit can still use real Z below.
fit_transform_3D = true; % Uses real TrackMate Z values. Falls back to XY affine if Z is underdetermined.
tform_type = "affine";

%% Compute transformation C1 -> C3
fprintf('=== Computing transformation C1 -> C3 ===\n');
all_paired_C3 = [];
all_paired_C1 = [];

% Check if inputs are cell arrays or character arrays (for backward compatibility)
if ischar(C3_beads_csv)
    C3_beads_csv = {C3_beads_csv};
    C1_beads_csv = {C1_beads_csv};
end
if ischar(C1_beads_csv)
    C1_beads_csv = {C1_beads_csv};
end
if numel(C1_beads_csv) ~= numel(C3_beads_csv)
    error('C1_beads_csv and C3_beads_csv must contain the same number of files.');
end
validateBeadCsvFiles(root_path, C3_beads_csv, 'C3');
validateBeadCsvFiles(root_path, C1_beads_csv, 'C1');
if do_C2_alignment
    if numel(C2_beads_csv) ~= numel(C3_beads_csv)
        error('C2_beads_csv and C3_beads_csv must contain the same number of files.');
    end
    validateBeadCsvFiles(root_path, C2_beads_csv, 'C2');
end

% Loop through all bead image pairs for C1->C3
for pair_idx = 1:length(C3_beads_csv)
    fprintf('Processing C1->C3 bead pair %d/%d: %s and %s\n', pair_idx, length(C3_beads_csv), ...
        C3_beads_csv{pair_idx}, C1_beads_csv{pair_idx});

    C3_beads_tab = readtable(fullfile(root_path, C3_beads_csv{pair_idx}));
    C1_beads_tab = readtable(fullfile(root_path, C1_beads_csv{pair_idx}));

    C3_beads_pos = table2array(C3_beads_tab(:,4:6));
    C1_beads_pos = table2array(C1_beads_tab(:,4:6));

    [paired_C3_beads_pos, paired_C1_beads_pos] = position_pairing(C3_beads_pos, C1_beads_pos, dist_thresh, use3D); % position_pairing(fix_beads_pos,moving_beads_pos,dist_thresh,use3D)

    if isempty(paired_C3_beads_pos)
        warning('No paired beads found for C1->C3 pair %d. Skipping.', pair_idx);
        continue;
    end

    % Collect paired positions
    all_paired_C3 = [all_paired_C3; paired_C3_beads_pos];
    all_paired_C1 = [all_paired_C1; paired_C1_beads_pos];
end

% Check if we have any paired positions for C1->C3
if isempty(all_paired_C3)
    error('No paired beads found for C1->C3 alignment.');
end

fprintf('C1->C3: Total paired beads collected: %d from %d bead image pairs\n', size(all_paired_C3,1), length(C3_beads_csv));

% Create summary figure for C1->C3
figure;
hold on;
scatter3(all_paired_C3(:,1), all_paired_C3(:,2), all_paired_C3(:,3), 'r', 'SizeData', 60);
scatter3(all_paired_C1(:,1), all_paired_C1(:,2), all_paired_C1(:,3), 'g', 'SizeData', 60);
for i = 1:size(all_paired_C3, 1)
    line([all_paired_C3(i,1), all_paired_C1(i,1)], ...
         [all_paired_C3(i,2), all_paired_C1(i,2)], ...
         [all_paired_C3(i,3), all_paired_C1(i,3)], ...
         'Color', 'b', 'LineWidth', 2);
end
hold off;
title(sprintf('C1->C3: All pooled paired beads (%d image pairs, %d beads)', length(C3_beads_csv), size(all_paired_C3,1)), 'Interpreter', 'none');
zlim([-10 10]); view(3); grid on;

% Compute transformation C1 -> C3
[tform_C1_to_C3, fit_mode_C1] = fitBeadTransform(all_paired_C1, all_paired_C3, fit_transform_3D, tform_type, 'C1->C3');
reportTransformResiduals(tform_C1_to_C3, all_paired_C1, all_paired_C3, 'C1->C3', fit_mode_C1);

%% Compute transformation C2 -> C3
if do_C2_alignment
    fprintf('\n=== Computing transformation C2 -> C3 ===\n');
    all_paired_C3 = [];
    all_paired_C2 = [];

    % Check if inputs are cell arrays or character arrays
    if ischar(C3_beads_csv)  % Already converted above, but keep for consistency
        C3_beads_csv = {C3_beads_csv};
    end
    if ischar(C2_beads_csv)
        C2_beads_csv = {C2_beads_csv};
    end

    % Loop through all bead image pairs for C2->C3
    for pair_idx = 1:length(C3_beads_csv)
        fprintf('Processing C2->C3 bead pair %d/%d: %s and %s\n', pair_idx, length(C3_beads_csv), ...
            C3_beads_csv{pair_idx}, C2_beads_csv{pair_idx});

        C3_beads_tab = readtable(fullfile(root_path, C3_beads_csv{pair_idx}));
        C2_beads_tab = readtable(fullfile(root_path, C2_beads_csv{pair_idx}));

        C3_beads_pos = table2array(C3_beads_tab(:,4:6));
        C2_beads_pos = table2array(C2_beads_tab(:,4:6));

        [paired_C3_beads_pos, paired_C2_beads_pos] = position_pairing(C3_beads_pos, C2_beads_pos, dist_thresh, use3D);

        if isempty(paired_C3_beads_pos)
            warning('No paired beads found for C2->C3 pair %d. Skipping.', pair_idx);
            continue;
        end

        % Collect paired positions
        all_paired_C3 = [all_paired_C3; paired_C3_beads_pos];
        all_paired_C2 = [all_paired_C2; paired_C2_beads_pos];
    end

    % Check if we have any paired positions for C2->C3
    if isempty(all_paired_C3)
        error('No paired beads found for C2->C3 alignment.');
    end

    fprintf('C2->C3: Total paired beads collected: %d from %d bead image pairs\n', size(all_paired_C3,1), length(C3_beads_csv));

    % Create summary figure for C2->C3
    figure;
    hold on;
    scatter3(all_paired_C3(:,1), all_paired_C3(:,2), all_paired_C3(:,3), 'r', 'SizeData', 60);
    scatter3(all_paired_C2(:,1), all_paired_C2(:,2), all_paired_C2(:,3), 'm', 'SizeData', 60); % magenta for C2
    for i = 1:size(all_paired_C3, 1)
        line([all_paired_C3(i,1), all_paired_C2(i,1)], ...
             [all_paired_C3(i,2), all_paired_C2(i,2)], ...
             [all_paired_C3(i,3), all_paired_C2(i,3)], ...
             'Color', 'b', 'LineWidth', 2);
    end
    hold off;
    title(sprintf('C2->C3: All pooled paired beads (%d image pairs, %d beads)', length(C3_beads_csv), size(all_paired_C3,1)), 'Interpreter', 'none');
    zlim([-10 10]); view(3); grid on;

    % Compute transformation C2 -> C3
    [tform_C2_to_C3, fit_mode_C2] = fitBeadTransform(all_paired_C2, all_paired_C3, fit_transform_3D, tform_type, 'C2->C3');
    reportTransformResiduals(tform_C2_to_C3, all_paired_C2, all_paired_C3, 'C2->C3', fit_mode_C2);
else
    tform_C2_to_C3 = [];
    fprintf('\n=== Skipping C2 alignment (no bead files provided) ===\n');
end

%% Display transformation matrices
fprintf('\n=== Transformation matrices ===\n');
fprintf('tform_C1_to_C3.A =\n');
disp(tform_C1_to_C3.A);
if ~isempty(tform_C2_to_C3)
    fprintf('tform_C2_to_C3.A =\n');
    disp(tform_C2_to_C3.A);
else
    fprintf('tform_C2_to_C3: [] (no C2 alignment performed)\n');
end

%% IMAGE REGISTRATION OPTIONS
run_nd2_registration = false;
run_tif_registration = false;
run_WF_registration = true;
run_mrc_denoise_registration = true;
run_mono_tif_beads_registration = false;
% If beads and processed images use the same scan ROI/crop origin, keep both as [0 0 0].
% If they differ, enter each ROI origin in the same physical coordinate system as TrackMate coordinates (um).
beads_roi_origin_um = [0 0 0];
process_roi_origin_um = [0 0 0];
tform_C1_to_C3_apply = shiftTransformOrigin(tform_C1_to_C3, beads_roi_origin_um, process_roi_origin_um);
if ~isempty(tform_C2_to_C3)
    tform_C2_to_C3_apply = shiftTransformOrigin(tform_C2_to_C3, beads_roi_origin_um, process_roi_origin_um);
else
    tform_C2_to_C3_apply = [];
end

%% REGISTRATION IMAGE USING SEPARATE CHANNEL IMAGES (MRC SIM data)
% all aligns to 560 channel
if run_mrc_denoise_registration
% Debug, using beads image as test sample
root_path = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260612_DLD1_RPB1-JF549_s5p-mintbody488';

img_subdirs = dir(fullfile(root_path));
temp_img_subdirs = {img_subdirs.name};
name_pat1 = "DLD1haloPEIs5p"; 
idx = contains(temp_img_subdirs,name_pat1);
img_subdirs = img_subdirs(idx);

% Pixel sizes for image registration (convert transformations from um to pixels)
px = 0.0325; % image pixel size, unit um; note after Wiener pixel get half !!!
py = px;
pz = 0.135; % z step, unit um

align_mode = 2;
include_beads_channel = [true, false, true]; % C3 is the fixed channel. C1->C3 (488->560), C2-C3 (640->560)
channel_suffix = ["SIM488", "SIM561"];

save_path = fullfile(root_path,'uncropped_aligned_imgs');
ensureDir(save_path);

for idx_dir =  1:length(img_subdirs)
    fprintf('Converting image: %s ... \n', img_subdirs(idx_dir).name);
    
    Vol = [];
    for iCh = 1:length(channel_suffix)   
        fprintf('Converting channel: %s ... \n', channel_suffix{iCh});
        img_list = dir(fullfile(img_subdirs(idx_dir).folder, img_subdirs(idx_dir).name, '*.mrc'));
        name_pat1 = string(channel_suffix{iCh}); name_pat2 = "SIrecon"; name_pat3 = "Denoise";
        temp_Filenames = {img_list.name}; 
        idx = contains(temp_Filenames,name_pat1) & contains(temp_Filenames,name_pat2) & contains(temp_Filenames,name_pat3);
        img_list = img_list(idx);
        
        % 使用 ReadMRC 读取图像数据。假设 now_img 是一个 3D 矩阵 (rows, cols, frames)
        now_img = ReadMRC(fullfile(img_list.folder, img_list.name));   
        now_img = rot90(now_img); % make direction same as in Fiji ImageJ
        
        Vol = cat(4,Vol,now_img);
    end    

    % Determine number of channels and split volume
    % Try to auto-detect: check if divisible by 3 (3 channels) or 2 (2 channels)    
    if align_mode == 3
        % 3 channels assumed: [C1=488, C2=640, C3=561 (fixed)] 
        C1_vol = Vol(:,:,:,1);
        C2_vol = Vol(:,:,:,3);
        C3_vol = Vol(:,:,:,2);
        fprintf('Align mode: 3 channels assumed: [C1=488, C2=640, C3=561 (fixed)]\n');
    elseif align_mode == 2
        % 2 channels: [C1=488, C2=NULL, C3=561 (fixed)]     
        C1_vol = Vol(:,:,:,1);       
        C2_vol = []; % No C2 channel
        C3_vol = Vol(:,:,:,2);
        fprintf('Align mode: 2 channels: [C1=488, C2=NULL, C3=561 (fixed)]\n');   
    elseif align_mode == 1
        % 2 channels: [C1=NULL, C2=640, C3=561 (fixed)]     
        C1_vol = [];        
        C2_vol = Vol(:,:,:,2);
        C3_vol = Vol(:,:,:,1);
        fprintf('Align mode: 2 channels: [C1=NULL, C2=640, C3=561 (fixed)]\n');      
    end

    S = diag([1/px, 1/py, 1/pz, 1]); % um -> pixel
    S_inv = diag([px, py, pz, 1]); % pixel -> um

    % Prepare output reference based on C3 channel size
    sizeC3 = size(C3_vol);
    Rfixed = imref3d(sizeC3);

    if ~include_beads_channel(3)
        error('C3 channel is not using as fixed channel, can not continue');
    end

    % Register C1 to C3
    if include_beads_channel(1)
        T_px_C1 = S * tform_C1_to_C3_apply.A * S_inv;
        tform_px_C1 = affinetform3d(T_px_C1);
        C1_reg = imwarp(C1_vol, tform_px_C1, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);
    else
        C1_reg = C1_vol;
    end

    % Register C2 to C3
    if include_beads_channel(2)
        if isempty(tform_C2_to_C3_apply)
            error('include_channel(2) is true, but no C2->C3 transform is available.');
        end
        T_px_C2 = S * tform_C2_to_C3_apply.A * S_inv;
        tform_px_C2 = affinetform3d(T_px_C2);
        C2_reg = imwarp(C2_vol, tform_px_C2, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);    
    else
        C2_reg = C2_vol;
    end

    % Combine channels into multi-channel volume
    % Order: C3 (reference), C1_aligned, C2_aligned
    % Use cat(4,...) to create XYZCT format where T is channel dimension
    if ~isempty(C1_reg) && ~isempty(C2_reg) && ~isempty(C3_vol)
        % 3-channel output: C3 (reference), C1_aligned, C2_aligned        
        merged_imgstack_final = cat(4, C1_reg, C3_vol, C2_reg);
        fprintf('  Saved as 3-channel image (C1_reg, C3_vol, C2_reg)\n');
    elseif ~isempty(C1_reg) && isempty(C2_reg) && ~isempty(C3_vol)
        merged_imgstack_final = cat(4, C1_reg, C3_vol);
        fprintf('  Saved as 2-channel image (C1_reg, C3_vol)\n');
    elseif isempty(C1_reg) && ~isempty(C2_reg) && ~isempty(C3_vol)
        merged_imgstack_final = cat(4, C3_vol,C2_reg);
        fprintf('  Saved as 2-channel image (C3_vol,C2_reg)\n');
    else
        error('Missing channels for combining');
    end

    % export the registered image, must in OME.TIF
    merged_imgstack_final(merged_imgstack_final<0) = 0;
    bfsave(merged_imgstack_final, fullfile(save_path, sprintf('%s_aligned.ome.tif', img_subdirs(idx_dir).name)), 'dimensionOrder', 'XYZCT');
end    
disp('Convertion finished!')
end


%% REGISTRATION IMAGE using ND2 stacks, save each channel, also pad black image at begin and end frame, do not touch 405
if run_nd2_registration
root_path = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260525_共定位统计\共定位统计\SSB1_JF549';
% Pixel sizes for image registration (convert transformations from um to pixels)
px = 0.06; % image pixel size, unit um; note after Wiener pixel get half
py = px;
pz = 0.2; % z step, unit um
channel_number = 3;
channel_surfix = ["C405", "C488", "C561"];
start_slice = 0; end_hold_slice = 0; % 488/561 settings
include_channel = [false, true, true];
% input_path = fullfile(root_path,'input_imgs');
input_path = fullfile(root_path);
save_path = fullfile(root_path,'split_aligned_imgs_codex');
other_path1 = fullfile(root_path,'aligned_cropped_imgs');
other_path2 = fullfile(root_path,'coloc_results_Coloc2_substractBackground25px');
ensureDir(save_path);
ensureDir(other_path1);
ensureDir(other_path2);

if ~exist(input_path, "dir")
    error('Input path does not exist: %s', input_path);
end
zstack_img_files = dir(fullfile(input_path,'*.nd2'));
if isempty(zstack_img_files)
    warning('No ND2 files found in %s', input_path);
end


for i = 1:length(zstack_img_files)
    zstack_img_file = zstack_img_files(i).name;
    fprintf('Converting image: %s ... \n', zstack_img_file);
    % Vol = tiffreadVolume(fullfile(input_path,zstack_img_file));
    Vol = MemoryEfficientND2reader(fullfile(input_path,zstack_img_file));
    z_len = size(Vol,3);

    % Determine number of channels and split volume
    % Try to auto-detect: check if divisible by 3 (3 channels) or 2 (2 channels)
    if mod(z_len, channel_number) ~= 0
        error('Cannot split %s: z_len=%d is not divisible by channel_number=%d.', zstack_img_file, z_len, channel_number);
    end
    if channel_number == 3
        % 3 channels assumed: [C1, C2, C3] in Z order
        ch_len = z_len / 3;
        C1_vol = Vol(:,:,1:3:end);
        C2_vol = Vol(:,:,2:3:end);
        C3_vol = Vol(:,:,3:3:end);
        fprintf('  Detected 3 channels, each with %d z-slices\n', ch_len);
    elseif channel_number == 2
        % 2 channels (backward compatible): [C3, C1] in Z order
        % ch_len = z_len / 2;
        % C1_vol = Vol(:,:,1:2:end);
        % C3_vol = Vol(:,:,2:2:end);
        % C2_vol = []; % No C2 channel
        % fprintf('  Detected 2 channels, each with %d z-slices (C3, C1)\n', ch_len);
        error('Not implemented yet!');
    else
        error('Cannot determine channel splitting: z_len=%d not divisible by 2 or 3', z_len);
    end
    if start_slice < 1 || start_slice > ch_len
        error('start_slice=%d is outside the valid z-slice range 1:%d for %s.', start_slice, ch_len, zstack_img_file);
    end    

    S = diag([1/px, 1/py, 1/pz, 1]); % um -> pixel
    S_inv = diag([px, py, pz, 1]); % pixel -> um

    % Prepare output reference based on C3 channel size
    sizeC3 = size(C3_vol);
    Rfixed = imref3d(sizeC3);

    if ~include_channel(3)
        error('C3 channel is not using as fixed channel, can not continue');
    end

    % Register C1 to C3
    if include_channel(1)
        T_px_C1 = S * tform_C1_to_C3_apply.A * S_inv;
        tform_px_C1 = affinetform3d(T_px_C1);
        C1_reg = imwarp(C1_vol, tform_px_C1, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);
    else
        C1_reg = C1_vol;
    end

    % Register C2 to C3
    if include_channel(2)
        if isempty(tform_C2_to_C3_apply)
            error('include_channel(2) is true, but no C2->C3 transform is available.');
        end
        T_px_C2 = S * tform_C2_to_C3_apply.A * S_inv;
        tform_px_C2 = affinetform3d(T_px_C2);
        C2_reg = imwarp(C2_vol, tform_px_C2, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);    
    else
        C2_reg = C2_vol;
    end

    % Combine channels into multi-channel volume
    % Order: C3 (reference), C1_aligned, C2_aligned
    % Use cat(4,...) to create XYZCT format where T is channel dimension
    if ~isempty(C1_reg) && ~isempty(C2_reg) && ~isempty(C3_vol)
        % 3-channel output: C3 (reference), C1_aligned, C2_aligned
        zero_frame = zeros(size(C1_reg,1), size(C1_reg,2), 'like', C1_reg); % 构造零帧 (与单帧大小一致)
        C1_reg = cat(3, zero_frame, C1_reg(:,:,start_slice:end-end_hold_slice), zero_frame); % 在前后各加一帧零图像
        C2_reg = cat(3, zero_frame, C2_reg(:,:,start_slice:end-end_hold_slice), zero_frame); % 在前后各加一帧零图像
        C3_vol = cat(3, zero_frame, C3_vol(:,:,start_slice:end-end_hold_slice), zero_frame); % 在前后各加一帧零图像        
        fprintf('Saved as mono-channel image for each channel :C1_aligned, C2_aligned, C3_fixed\n');

        % export the registered image, must in OME.TIF
        bfsave(uint16(C1_reg), fullfile(save_path, sprintf('%s_%s_aligned.ome.tif', channel_surfix(1), zstack_img_file(1:end-4))), 'dimensionOrder', 'XYZCT');
        bfsave(uint16(C2_reg), fullfile(save_path, sprintf('%s_%s_aligned.ome.tif', channel_surfix(2), zstack_img_file(1:end-4))), 'dimensionOrder', 'XYZCT');
        bfsave(uint16(C3_vol), fullfile(save_path, sprintf('%s_%s_aligned.ome.tif', channel_surfix(3), zstack_img_file(1:end-4))), 'dimensionOrder', 'XYZCT');

    elseif ~isempty(C1_reg) && isempty(C2_reg) && ~isempty(C3_vol)
        % 2-channel output (backward compatible): C3, C1_aligned
        % merged_imgstack_final = cat(4, C1_reg, C3_vol);
        error('Not implemented yet!');
        % fprintf('  Saved as 2-channel image (C1_aligned, C3_fixed)\n');        
    else
        error('Missing channels for combining');
    end    
end
disp('Convertion finished!')
end

%% REGISTRATION IMAGE using tif stacks
if run_tif_registration
root_path = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260525_共定位统计\共定位统计\8N_8wg16_JF549\';
% Pixel sizes for image registration (convert transformations from um to pixels)
px = 0.06; % image pixel size, unit um; note after Wiener pixel get half
py = px;
pz = 0.2; % z step, unit um
channel_number = 3;
input_path = fullfile(root_path,'split_imgs');
save_path = fullfile(root_path,'aligned_imgs');
ensureDir(save_path);

% zstack_img_file = {...
%     'Sparse_1_Merge_405_Em450_488_Em525_638_Em667_Wiener_3DSIM-3_OCT4_ASHL2_Cell01.ome.tif',...    
%     };
if ~exist(input_path, "dir")
    error('Input path does not exist: %s', input_path);
end
zstack_img_files = dir(fullfile(input_path,'*.tif'));
if isempty(zstack_img_files)
    warning('No TIF files found in %s', input_path);
end


for i = 1:length(zstack_img_files)
    zstack_img_file = zstack_img_files(i).name;
    fprintf('Converting image: %s ... \n', zstack_img_file);
    Vol = tiffreadVolume(fullfile(input_path,zstack_img_file));    
    z_len = size(Vol,3);

    % Determine number of channels and split volume
    % Try to auto-detect: check if divisible by 3 (3 channels) or 2 (2 channels)
    if mod(z_len, channel_number) ~= 0
        error('Cannot split %s: z_len=%d is not divisible by channel_number=%d.', zstack_img_file, z_len, channel_number);
    end
    if channel_number == 3
        % 3 channels assumed: [C1, C2, C3] in Z order
        ch_len = z_len / 3;
        C1_vol = Vol(:,:,1:ch_len);
        C2_vol = Vol(:,:,ch_len+1:2*ch_len);
        C3_vol = Vol(:,:,2*ch_len+1:end);
        fprintf('  Detected 3 channels, each with %d z-slices\n', ch_len);
    elseif channel_number == 2
        % 2 channels (backward compatible): [C3, C1] in Z order
        ch_len = z_len / 2;
        C1_vol = Vol(:,:,1:ch_len);
        C3_vol = Vol(:,:,ch_len+1:end);
        C2_vol = []; % No C2 channel
        fprintf('  Detected 2 channels, each with %d z-slices (C3, C1)\n', ch_len);
    else
        error('Cannot determine channel splitting: z_len=%d not divisible by 2 or 3', z_len);
    end    

    S = diag([1/px, 1/py, 1/pz, 1]); % um -> pixel
    S_inv = diag([px, py, pz, 1]); % pixel -> um

    % Prepare output reference based on C3 channel size
    sizeC3 = size(C3_vol);
    Rfixed = imref3d(sizeC3);

    % Register C1 to C3
    if ~isempty(C1_vol)
        T_px_C1 = S * tform_C1_to_C3_apply.A * S_inv;
        tform_px_C1 = affinetform3d(T_px_C1);
        C1_reg = imwarp(C1_vol, tform_px_C1, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);
    else
        C1_reg = [];
    end

    % Register C2 to C3
    if ~isempty(C2_vol) && ~isempty(tform_C2_to_C3_apply)
        T_px_C2 = S * tform_C2_to_C3_apply.A * S_inv;
        tform_px_C2 = affinetform3d(T_px_C2);
        C2_reg = imwarp(C2_vol, tform_px_C2, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);
    elseif ~isempty(C2_vol) && isempty(tform_C2_to_C3_apply)
        warning('C2 channel present but no transformation available. Skipping C2 registration.');
        C2_reg = [];
    else
        C2_reg = [];
    end

    % Combine channels into multi-channel volume
    % Order: C3 (reference), C1_aligned, C2_aligned
    % Use cat(4,...) to create XYZCT format where T is channel dimension
    if ~isempty(C1_reg) && ~isempty(C2_reg) && ~isempty(C3_vol)
        % 3-channel output: C3 (reference), C1_aligned, C2_aligned
        merged_imgstack_final = cat(4, C1_reg, C2_reg, C3_vol);
        fprintf('  Saved as 3-channel image (C1_aligned, C2_aligned, C3_fixed)\n');
    elseif ~isempty(C1_reg) && isempty(C2_reg) && ~isempty(C3_vol)
        % 2-channel output (backward compatible): C3, C1_aligned
        merged_imgstack_final = cat(4, C1_reg, C3_vol);
        fprintf('  Saved as 2-channel image (C1_aligned, C3_fixed)\n');
    else
        error('Missing channels for combining');
    end

    % export the registered image, must in OME.TIF
    bfsave(uint16(merged_imgstack_final), fullfile(save_path, sprintf('%s_aligned.ome.tif', stripTiffSuffix(zstack_img_file))), 'dimensionOrder', 'XYZCT');    
end
disp('Convertion finished!')
end




%% BEADS/WF REGISTRATION IMAGE from MRC (CHECK ONLY)
% all aligns to 560 channel
if run_WF_registration
% Debug, using beads image as test sample
root_path = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260612_DLD1_RPB1-JF549_s5p-mintbody488';
save_path = fullfile(root_path,'uncropped_aligned_WF_imgs');
ensureDir(save_path);

img_subdirs = dir(fullfile(root_path));
temp_img_subdirs = {img_subdirs.name};
name_pat1 = "DLD1haloPEIs5p"; 
idx = contains(temp_img_subdirs,name_pat1);
img_subdirs = img_subdirs(idx);

% Pixel sizes for image registration (convert transformations from um to pixels)
px = 0.065; % image pixel size, unit um
py = px;
pz = 0.135; % z step, unit um

align_mode = 2;
include_beads_channel = [true, false, true]; % C3 is the fixed channel. C1->C3 (488->560), C2-C3 (640->560)
channel_suffix = ["SIM488", "SIM561"];

for idx_dir =  1:length(img_subdirs)
    fprintf('Converting image: %s ... \n', img_subdirs(idx_dir).name);    
    Vol = [];
    for iCh = 1:length(channel_suffix)   
        fprintf('Converting channel: %s ... \n', channel_suffix{iCh});
        img_list = dir(fullfile(img_subdirs(idx_dir).folder, img_subdirs(idx_dir).name, '*.mrc'));
        name_pat1 = string(channel_suffix{iCh}); name_pat2 = "SIrecon"; name_pat3 = "Denoise";
        temp_Filenames = {img_list.name}; 
        idx = contains(temp_Filenames,name_pat1) & ~contains(temp_Filenames,name_pat2) & ~contains(temp_Filenames,name_pat3);
        img_list = img_list(idx);
        
        % 使用 ReadMRC 读取图像数据。假设 now_img 是一个 3D 矩阵 (rows, cols, frames)
        now_img = ReadMRC(fullfile(img_list.folder, img_list.name));   
        now_img = rot90(now_img); % make direction same as in Fiji ImageJ

        % 确保数据类型是 single (32位浮点数)
        now_img = single(now_img);        
        [ImHeight,ImWidth,~] = size(now_img);

        % 假设 now_img 的维度为 [row, height, slice]
        [row, col, slice] = size(now_img);
        
        % 计算新的 slice 维度（每15帧求和）
        new_slice = floor(slice / 15);
        
        % 将 now_img 重塑为 [row, height, 15, new_slice]
        reshaped_img = reshape(now_img, [row, col, 15, new_slice]);
        
        % 沿第3维（即每15帧）求和
        now_summed_img = sum(reshaped_img, 3); now_summed_img = squeeze(now_summed_img);
        imgstack_final = now_summed_img;
        
        Vol = cat(4,Vol,imgstack_final);
    end    

    % Determine number of channels and split volume
    % Try to auto-detect: check if divisible by 3 (3 channels) or 2 (2 channels)    
    if align_mode == 3
        % 3 channels assumed: [C1=488, C2=640, C3=561 (fixed)] 
        C1_vol = Vol(:,:,:,1);
        C2_vol = Vol(:,:,:,3);
        C3_vol = Vol(:,:,:,2);
        fprintf('Align mode: 3 channels assumed: [C1=488, C2=640, C3=561 (fixed)]\n');
    elseif align_mode == 2
        % 2 channels: [C1=488, C2=NULL, C3=561 (fixed)]     
        C1_vol = Vol(:,:,:,1);       
        C2_vol = []; % No C2 channel
        C3_vol = Vol(:,:,:,2);
        fprintf('Align mode: 2 channels: [C1=488, C2=NULL, C3=561 (fixed)]\n');   
    elseif align_mode == 1
        % 2 channels: [C1=NULL, C2=640, C3=561 (fixed)]     
        C1_vol = [];        
        C2_vol = Vol(:,:,:,2);
        C3_vol = Vol(:,:,:,1);
        fprintf('Align mode: 2 channels: [C1=NULL, C2=640, C3=561 (fixed)]\n');      
    end    

    S = diag([1/px, 1/py, 1/pz, 1]); % um -> pixel
    S_inv = diag([px, py, pz, 1]); % pixel -> um

    % Prepare output reference based on C3 channel size
    sizeC3 = size(C3_vol);
    Rfixed = imref3d(sizeC3);

    if ~include_beads_channel(3)
        error('C3 channel is not using as fixed channel, can not continue');
    end

    % Register C1 to C3
    if include_beads_channel(1)
        T_px_C1 = S * tform_C1_to_C3_apply.A * S_inv;
        tform_px_C1 = affinetform3d(T_px_C1);
        C1_reg = imwarp(C1_vol, tform_px_C1, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);
    else
        C1_reg = C1_vol;
    end

    % Register C2 to C3
    if include_beads_channel(2)
        if isempty(tform_C2_to_C3_apply)
            error('include_channel(2) is true, but no C2->C3 transform is available.');
        end
        T_px_C2 = S * tform_C2_to_C3_apply.A * S_inv;
        tform_px_C2 = affinetform3d(T_px_C2);
        C2_reg = imwarp(C2_vol, tform_px_C2, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);    
    else
        C2_reg = C2_vol;
    end

    % Combine channels into multi-channel volume
    % Order: C3 (reference), C1_aligned, C2_aligned
    % Use cat(4,...) to create XYZCT format where T is channel dimension
    if ~isempty(C1_reg) && ~isempty(C2_reg) && ~isempty(C3_vol)
        % 3-channel output: C3 (reference), C1_aligned, C2_aligned        
        merged_imgstack_final = cat(4, C1_reg, C3_vol, C2_reg);
        fprintf('  Saved as 3-channel image (C1_reg, C3_vol, C2_reg)\n');
    elseif ~isempty(C1_reg) && isempty(C2_reg) && ~isempty(C3_vol)
        merged_imgstack_final = cat(4, C1_reg, C3_vol);
        fprintf('  Saved as 2-channel image (C1_reg, C3_vol)\n');
    elseif isempty(C1_reg) && ~isempty(C2_reg) && ~isempty(C3_vol)
        merged_imgstack_final = cat(4, C3_vol,C2_reg);
        fprintf('  Saved as 2-channel image (C3_vol,C2_reg)\n');
    else
        error('Missing channels for combining');
    end
    
    merged_imgstack_final(merged_imgstack_final<0) = 0;
    % export the registered image, must in OME.TIF
    bfsave(merged_imgstack_final, fullfile(save_path, sprintf('%s_aligned.ome.tif', img_subdirs(idx_dir).name)), 'dimensionOrder', 'XYZCT');
end
disp('Convertion finished!')
end


%% BEADS/WF REGISTRATION IMAGE from mono channel TIF (CHECK ONLY)
if run_mono_tif_beads_registration
% Debug, using beads image as test sample
root_path = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned\Archive\DONOTUSE_20260415_8N_JF549-Flag647_beads\';
input_path = fullfile(root_path,'converted_WF_tif');
save_path = fullfile(root_path,'aligned_imgs');
ensureDir(save_path);
zstack_img_file = {'beads_002_20260415_160214'};

% Pixel sizes for image registration (convert transformations from um to pixels)
px = 0.065; % image pixel size, unit um
py = px;
pz = 0.135; % z step, unit um

align_mode = 3;
include_beads_channel = [true, true, true]; % C3 is the fixed channel. C1->C3 (488->560), C2-C3 (640->560)
channel_suffix = ["SIM488", "SIM561", "SIM640"];

for idx_dir = 1:length(zstack_img_file)
    fprintf('Converting image: %s ... \n', zstack_img_file{idx_dir});
    
    img_list = dir(fullfile(input_path,'*.tif'));
    name_pat1 = string(zstack_img_file{idx_dir}); name_pat2 = "bgsub";
    temp_Filenames = {img_list.name}; 
    idx = contains(temp_Filenames,name_pat1) & contains(temp_Filenames,name_pat2);
    img_list = img_list(idx);
    Vol = [];
    for ii = 1:length(img_list)
        Vol_tmp = tiffreadVolume(fullfile(img_list(ii).folder,img_list(ii).name));
        Vol = cat(4,Vol,Vol_tmp);
    end    

    % Determine number of channels and split volume
    % Try to auto-detect: check if divisible by 3 (3 channels) or 2 (2 channels)    
    if align_mode == 3
        % 3 channels assumed: [C1=488, C2=640, C3=561 (fixed)] 
        C1_vol = Vol(:,:,:,1);
        C2_vol = Vol(:,:,:,3);
        C3_vol = Vol(:,:,:,2);
        fprintf('Align mode: 3 channels assumed: [C1=488, C2=640, C3=561 (fixed)]\n');
    elseif align_mode == 2
        % 2 channels: [C1=488, C2=NULL, C3=561 (fixed)]     
        C1_vol = Vol(:,:,:,1);       
        C2_vol = []; % No C2 channel
        C3_vol = Vol(:,:,:,2);
        fprintf('Align mode: 2 channels: [C1=488, C2=NULL, C3=561 (fixed)]\n');   
    elseif align_mode == 1
        % 2 channels: [C1=NULL, C2=640, C3=561 (fixed)]     
        C1_vol = [];        
        C2_vol = Vol(:,:,:,2);
        C3_vol = Vol(:,:,:,1);
        fprintf('Align mode: 2 channels: [C1=NULL, C2=640, C3=561 (fixed)]\n');      
    end    

    S = diag([1/px, 1/py, 1/pz, 1]); % um -> pixel
    S_inv = diag([px, py, pz, 1]); % pixel -> um

    % Prepare output reference based on C3 channel size
    sizeC3 = size(C3_vol);
    Rfixed = imref3d(sizeC3);

    if ~include_beads_channel(3)
        error('C3 channel is not using as fixed channel, can not continue');
    end

    % Register C1 to C3
    if include_beads_channel(1)
        T_px_C1 = S * tform_C1_to_C3_apply.A * S_inv;
        tform_px_C1 = affinetform3d(T_px_C1);
        C1_reg = imwarp(C1_vol, tform_px_C1, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);
    else
        C1_reg = C1_vol;
    end

    % Register C2 to C3
    if include_beads_channel(2)
        if isempty(tform_C2_to_C3_apply)
            error('include_channel(2) is true, but no C2->C3 transform is available.');
        end
        T_px_C2 = S * tform_C2_to_C3_apply.A * S_inv;
        tform_px_C2 = affinetform3d(T_px_C2);
        C2_reg = imwarp(C2_vol, tform_px_C2, 'OutputView', Rfixed, ...
            'InterpolationMethod', 'cubic', 'FillValues', 0);    
    else
        C2_reg = C2_vol;
    end

    % Combine channels into multi-channel volume
    % Order: C3 (reference), C1_aligned, C2_aligned
    % Use cat(4,...) to create XYZCT format where T is channel dimension
    if ~isempty(C1_reg) && ~isempty(C2_reg) && ~isempty(C3_vol)
        % 3-channel output: C3 (reference), C1_aligned, C2_aligned        
        merged_imgstack_final = cat(4, C1_reg, C3_vol, C2_reg);
        fprintf('  Saved as 3-channel image (C1_reg, C3_vol, C2_reg)\n');
    elseif ~isempty(C1_reg) && isempty(C2_reg) && ~isempty(C3_vol)
        merged_imgstack_final = cat(4, C1_reg, C3_vol);
        fprintf('  Saved as 2-channel image (C1_reg, C3_vol)\n');
    elseif isempty(C1_reg) && ~isempty(C2_reg) && ~isempty(C3_vol)
        merged_imgstack_final = cat(4, C3_vol,C2_reg);
        fprintf('  Saved as 2-channel image (C3_vol,C2_reg)\n');
    else
        error('Missing channels for combining');
    end

    merged_imgstack_final(merged_imgstack_final<0) = 0;
    % export the registered image, must in OME.TIF
    bfsave(merged_imgstack_final, fullfile(save_path, sprintf('3Dalign_%s_aligned.ome.tif', zstack_img_file{idx_dir})), 'dimensionOrder', 'XYZCT');
end
disp('Convertion finished!')

end


%% AUXILIARY FUNCTIONS
function validateBeadCsvFiles(root_path, fileList, channelLabel)
for k = 1:numel(fileList)
    filePath = fullfile(root_path, fileList{k});
    if ~exist(filePath, 'file')
        error('%s bead CSV does not exist: %s', channelLabel, filePath);
    end
end
end

function [tform, fitMode] = fitBeadTransform(movingPoints, fixedPoints, fit3D, tformType, label)
movingPoints = double(movingPoints);
fixedPoints = double(fixedPoints);
if size(movingPoints, 2) ~= 3 || size(fixedPoints, 2) ~= 3
    error('%s: movingPoints and fixedPoints must be N-by-3 arrays.', label);
end
if size(movingPoints, 1) ~= size(fixedPoints, 1)
    error('%s: movingPoints and fixedPoints must contain the same number of paired beads.', label);
end
if any(~isfinite(movingPoints(:))) || any(~isfinite(fixedPoints(:)))
    error('%s: paired bead coordinates contain NaN or Inf.', label);
end

min2D = minControlPoints(tformType, 2);
min3D = minControlPoints(tformType, 3);
if fit3D && size(movingPoints, 1) >= min3D && hasFull3DSpread(movingPoints) && hasFull3DSpread(fixedPoints)
    try
        tform = fitgeotform3d(movingPoints, fixedPoints, char(tformType));
        fitMode = sprintf('3D %s', char(tformType));
        return
    catch ME
        warning('%s: 3D %s fit failed: %s. Falling back to XY %s.', label, char(tformType), ME.message, char(tformType));
    end
elseif fit3D
    warning('%s: bead coordinates do not have enough stable 3D spread for a 3D %s fit. Falling back to XY %s and preserving Z.', label, char(tformType), char(tformType));
end

if size(movingPoints, 1) < min2D
    error('%s: need at least %d paired beads for an XY %s fit; got %d.', label, min2D, char(tformType), size(movingPoints, 1));
end
tform2D = fitgeotform2d(movingPoints(:,1:2), fixedPoints(:,1:2), char(tformType));
tform = convert2DTformTo3D(tform2D);
fitMode = sprintf('XY %s, Z unchanged', char(tformType));
end

function minPts = minControlPoints(tformType, dimensionality)
switch char(tformType)
    case 'translation'
        minPts = 1;
    case 'rigid'
        if dimensionality == 3
            minPts = 3;
        else
            minPts = 2;
        end
    case 'affine'
        minPts = dimensionality + 1;
    otherwise
        error('Unsupported transform type: %s. Use translation, rigid, or affine.', char(tformType));
end
end

function tf = hasFull3DSpread(points)
centered = points - mean(points, 1);
scale = max(max(points, [], 1) - min(points, [], 1));
if scale == 0
    tf = false;
    return
end
tf = rank(centered, scale * 1e-6) >= 3;
end

function tform3D = convert2DTformTo3D(tform2D)
A2 = tform2D.A;
A3 = eye(4);
A3(1:2, 1:2) = A2(1:2, 1:2);
A3(1:2, 4) = A2(1:2, 3);
tform3D = affinetform3d(A3);
end

function shiftedTform = shiftTransformOrigin(tform, sourceOriginUm, targetOriginUm)
sourceOriginUm = sourceOriginUm(:);
targetOriginUm = targetOriginUm(:);
if numel(sourceOriginUm) ~= 3 || numel(targetOriginUm) ~= 3
    error('ROI origins must be 1-by-3 or 3-by-1 vectors: [x y z] in um.');
end
A = tform.A;
originDelta = targetOriginUm - sourceOriginUm;
A(1:3,4) = A(1:3,4) + (A(1:3,1:3) - eye(3)) * originDelta;
shiftedTform = affinetform3d(A);
end

function reportTransformResiduals(tform, movingPoints, fixedPoints, label, fitMode)
registeredPoints = transformPointsForward(tform, movingPoints);
delta = registeredPoints - fixedPoints;
errXY = hypot(delta(:,1), delta(:,2));
err3D = sqrt(sum(delta.^2, 2));
fprintf('%s residuals (%s): XY median %.4f um, XY p95 %.4f um, 3D median %.4f um, 3D p95 %.4f um, max 3D %.4f um\n', ...
    label, fitMode, median(errXY), percentileValue(errXY, 95), median(err3D), percentileValue(err3D, 95), max(err3D));
end

function p = percentileValue(values, percentile)
values = sort(values(:));
values = values(isfinite(values));
if isempty(values)
    p = NaN;
    return
end
idx = 1 + (numel(values) - 1) * percentile / 100;
lo = floor(idx);
hi = ceil(idx);
if lo == hi
    p = values(lo);
else
    p = values(lo) + (idx - lo) * (values(hi) - values(lo));
end
end

function ensureDir(pathStr)
if ~exist(pathStr, "dir")
    mkdir(pathStr);
end
end

function baseName = stripTiffSuffix(fileName)
[~, name, ext] = fileparts(fileName);
baseName = name;
if strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')
    if numel(baseName) > 4 && strcmpi(baseName(end-3:end), '.ome')
        baseName = baseName(1:end-4);
    end
end
end

