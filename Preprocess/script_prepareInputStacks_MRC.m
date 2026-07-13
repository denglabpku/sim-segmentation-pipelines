%% Prepare mrc SIM image
% Note: This code assumes you have the ReadMRC function available and in
% your MATLAB path. Fred Sigworth (2026). Imagic, MRC, DM and STAR file i/o
% (https://www.mathworks.com/matlabcentral/fileexchange/27021-imagic-mrc-dm-and-star-file-i-o),
% MATLAB Central File Exchange. Retrieved April 17, 2026.

% This script for MRC data is newer than script_prepareInputStacks.m.

% Zuhui Wang 2026/06/08

%% Prepare multi-color mrc WF image
close all; clear; clc;
% Channel configuration
% channel_number = 1;  % Number of channels
channel_suffix = {'488', '561', '640'};  % Suffix for each channel
start_frame = 1; % remove high background far defocus frames to avoid error in localizations, also combine with rolling ball background substraction

raw_img_dir = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned\20260429_8N_JF549-Flag647\beads';
save_img_dir = fullfile(raw_img_dir,'converted_WF_tif');

if ~exist(save_img_dir,"dir")
    mkdir(save_img_dir)
end

img_subdirs = dir(fullfile(raw_img_dir));
temp_img_subdirs = {img_subdirs.name};
name_pat1 = "beads"; 
idx = contains(temp_img_subdirs,name_pat1);
img_subdirs = img_subdirs(idx);

for idx_dir =  1:length(img_subdirs)
        
    for iCh = 1:length(channel_suffix)   
        img_list = dir(fullfile(img_subdirs(idx_dir).folder, img_subdirs(idx_dir).name, '*.mrc'));
        name_pat1 = string(channel_suffix{iCh}); name_pat2 = "SIrecon"; name_pat3 = "Denoise";
        temp_Filenames = {img_list.name}; 
        idx = contains(temp_Filenames,name_pat1) & ~contains(temp_Filenames,name_pat2) & ~contains(temp_Filenames,name_pat3);
        img_list = img_list(idx);
        
        % 使用 ReadMRC 读取图像数据。假设 now_img 是一个 3D 矩阵 (rows, cols, frames)
        now_img = ReadMRC(fullfile(img_list.folder, img_list.name));   
        now_img = rot90(now_img); % make direction same as in Fiji ImageJ
        
        fprintf('Converting image: %s ... Current Channel: %d/%d \n', img_list.name, iCh, length(img_list));
    
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
        imgstack_final = now_summed_img(:,:,start_frame:end);
        
        % Save with channel suffix
        img_name = sprintf('%s_C%s.tif', img_subdirs(idx_dir).name, channel_suffix{iCh});
        full_save_path = fullfile(save_img_dir, img_name);
        
        % 使用 Tiff 类保存为 32-bit float
        t = Tiff(full_save_path, 'w');
        tagstruct.ImageLength = size(imgstack_final,1);
        tagstruct.ImageWidth = size(imgstack_final,2);
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; % 浮点数
        tagstruct.BitsPerSample = 32;                      % 32位
        tagstruct.SamplesPerPixel = 1;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        
        % 如果是多帧，逐帧写入
        for f = 1:size(imgstack_final,3)
            setTag(t, tagstruct);
            write(t, imgstack_final(:,:,f));
            if f < size(imgstack_final,3)
                writeDirectory(t);
            end
        end

        close(t);

    end
end

%% Prepare dual-color tif SIM image (usually after alignment)
% I highly recommend using the fixed 561 channel to draw cell ROI, which is
% easier to reuse ROI if redo alignment.
clc; close all; clear;
addpath(genpath('C:\Users\Zuhui\OneDrive - Peking University\Documents\MATLAB\bfmatlab'),'-end');
raw_img_dir = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260612_DLD1_RPB1-JF549_s5p-mintbody488';
input_img_dir = fullfile(raw_img_dir,'uncropped_aligned_imgs');
save_img_dir = fullfile(raw_img_dir,'cropped_aligned_imgs');

if ~exist(save_img_dir,"dir")
    mkdir(save_img_dir)
end

img_list = dir(fullfile(input_img_dir,'*.tif'));
name_pat1 = "aligned"; name_pat2 = "ome";
temp_Filenames = {img_list.name}; 
idx = contains(temp_Filenames,name_pat1) & contains(temp_Filenames,name_pat2) & ~startsWith(temp_Filenames, ".");
img_list = img_list(idx);
channel_number = 2;

for iImage = 1:length(img_list)  
    now_img = tiffreadVolume(fullfile(img_list(iImage).folder,img_list(iImage).name));   
    z_len = size(now_img,3)/channel_number;

    now_img_chx = {now_img(:,:,1:z_len), now_img(:,:,z_len+1:end)};        
    
    fprintf('Converting image: %s ... %d/%d \n', img_list(iImage).name, iImage, length(img_list));
      
    [ImHeight,ImWidth,~] = size(now_img_chx{1});    

    midFrameIdx=inputdlg("Middle frame index:","Require inputs",[1 50]);
    midFrameIdx = str2double(midFrameIdx{1});
       
    % Draw nuclear masks   
    while (1)   
        close all;
        % Show middle frame to input nucleus number
        midFrame = now_img_chx{1}(:,:,midFrameIdx);        

        fig_ROI = figure('position',[362 137 600 600]);        
        imshow(midFrame,[]);
        hIc = imcontrast; uiwait(hIc);
        title('Middle frame', 'Interpreter','none');

        % Input nuclear number
        seg_prompt = {'Enter number of nuclei:'};
        seg_prompt_name = 'Input';
        dims = [1 35];
        seg_answer = inputdlg(seg_prompt,seg_prompt_name,dims,{'1'});
        Nuc_number = str2double(seg_answer{1});  % number of nucleus mask 
        hROI = cell(Nuc_number,1); 
        Nuc_individual_masks = zeros(ImHeight,ImWidth,Nuc_number);
        roi_info_nuc = cell(Nuc_number,1);         
        
        for iNuc = 1:Nuc_number
            if iNuc == 1                
                ax_ROI = fig_ROI.CurrentAxes;
                hROI{iNuc} = drawpolygon(ax_ROI,'LabelVisible','on',...
                    'FaceAlpha', 0, 'Label', num2str(iNuc),...
                    'MarkerSize', 5,'LineWidth', 2,'Color','green');  
                wait(hROI{iNuc});
                position = hROI{iNuc}.Position; % Each row represents the [x y] coordinates of a ROI vertex
                roi_info_nuc(iNuc) = {position};
                Nuc_individual_masks(:,:,iNuc) = createMask(hROI{iNuc});
    
            else % Nuc_number > 1 && iNuc < Nuc_number
                message = sprintf('Outline another ROI');
                uiwait(msgbox(message));
    
                ax_ROI = fig_ROI.CurrentAxes;
                hROI{iNuc} = drawpolygon(ax_ROI,'LabelVisible','on',...
                    'FaceAlpha', 0, 'Label', num2str(iNuc),...
                    'MarkerSize', 5,'LineWidth', 2,'Color','green');
    
                wait(hROI{iNuc});
                position = hROI{iNuc}.Position; % Each row represents the [x y] coordinates of a ROI vertex
                roi_info_nuc(iNuc) = {position};
                Nuc_individual_masks(:,:,iNuc) = createMask(hROI{iNuc});
    
            end
        end % for
        if iNuc == Nuc_number     
            promptMessage = sprintf('Are you happy with the ROI, yes to continue, no to re-select ROI?');
            button = questdlg(promptMessage, 'Is ROI good?', 'Yes', 'No', 'Yes');  % quest, title, btn1, btn2, defbtn
            if strcmpi(button, 'Yes')  %compare string
                % save ROI                                
                save(fullfile(save_img_dir, sprintf("%s_NuclearMask.mat",img_list(iImage).name(1:end-4))),'roi_info_nuc','ImHeight','ImWidth','midFrameIdx');
                close all;
                break; 
            end
        end
    end % while
    
    
    for iCell = 1:length(roi_info_nuc)

        % --- Step 1: Get polygon vertices for this nucleus ---
        polyVerts = roi_info_nuc{iCell};   % N×2 double [x,y]
        
        % --- Step 2: Convert polygon to bounding rectangle ---
        xMin = floor(min(polyVerts(:,1)));
        xMax = ceil(max(polyVerts(:,1)));
        yMin = floor(min(polyVerts(:,2)));
        yMax = ceil(max(polyVerts(:,2)));

        merged_imgstack_final = [];

        for ich = 1:channel_number        
            % --- Step 3: Crop the rectangle ROI from full image stack ---
            cropped_img = now_img_chx{ich}(yMin:yMax, xMin:xMax, :);
            
            % --- Step 4: Construct zero frames (same size as cropped ROI) --- 
            zero_frame = zeros(size(cropped_img,1), size(cropped_img,2), 'single'); % 构造零帧 (与单帧大小一致)
            imgstack_final = cat(3, zero_frame, cropped_img, zero_frame); % 在前后各加一帧零图像
            
            % --- Step 5: Convert polygon to mask, cropped to ROI ---
            nuclearMaskFull = poly2mask(polyVerts(:,1), polyVerts(:,2), size(now_img,1), size(now_img,2));
            nuclearMaskCrop = nuclearMaskFull(yMin:yMax, xMin:xMax);
            
            % --- Step 6: Apply mask to cropped stack ---
            for f = 1:size(imgstack_final,3)
                imgstack_final(:,:,f) = imgstack_final(:,:,f) .* single(nuclearMaskCrop);
            end    
            merged_imgstack_final(:,:,:,ich) = imgstack_final;
                    
        end

        img_name = sprintf('%s_Cell%02d.ome.tif', img_list(iImage).name(1:end-4),iCell);
        full_save_path = fullfile(save_img_dir, img_name);
    
        bfsave(merged_imgstack_final,full_save_path,'dimensionOrder', 'XYZCT');
    end

end

%% Prepare dual-color tif SIM image (reuse ROI)

clc; close all; clear;
addpath(genpath('C:\Users\Zuhui\OneDrive - Peking University\Documents\MATLAB\bfmatlab'),'-end');
raw_img_dir = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260605_SIM_unaligned\20260417_11_JF549-Flag647';
input_img_dir = fullfile(raw_img_dir,'uncropped_aligned_imgs');
save_img_dir = fullfile(raw_img_dir,'cropped_aligned_imgs');

if ~exist(save_img_dir,"dir")
    mkdir(save_img_dir)
end

img_list = dir(fullfile(input_img_dir,'*.tif'));
name_pat1 = "aligned"; name_pat2 = "ome";
temp_Filenames = {img_list.name}; 
idx = contains(temp_Filenames,name_pat1) & contains(temp_Filenames,name_pat2) & ~startsWith(temp_Filenames, ".");
img_list = img_list(idx);
channel_number = 2;

for iImage = 1:length(img_list)  
    now_img = tiffreadVolume(fullfile(img_list(iImage).folder,img_list(iImage).name));   
    z_len = size(now_img,3)/channel_number;

    now_img_chx = {now_img(:,:,1:z_len), now_img(:,:,z_len+1:end)};        
    
    fprintf('Converting image: %s ... %d/%d \n', img_list(iImage).name, iImage, length(img_list));
      
    [ImHeight,ImWidth,~] = size(now_img_chx{1});        

    load(fullfile(save_img_dir, sprintf("%s_NuclearMask.mat",img_list(iImage).name(1:end-4))),'roi_info_nuc');
    fprintf('Previous defined cell ROI:%s loaded!  \n', img_list(iImage).name(1:end-4));
    
    
    for iCell = 1:length(roi_info_nuc)

        % --- Step 1: Get polygon vertices for this nucleus ---
        polyVerts = roi_info_nuc{iCell};   % N×2 double [x,y]
        
        % --- Step 2: Convert polygon to bounding rectangle ---
        xMin = floor(min(polyVerts(:,1)));
        xMax = ceil(max(polyVerts(:,1)));
        yMin = floor(min(polyVerts(:,2)));
        yMax = ceil(max(polyVerts(:,2)));

        merged_imgstack_final = [];

        for ich = 1:channel_number        
            % --- Step 3: Crop the rectangle ROI from full image stack ---
            cropped_img = now_img_chx{ich}(yMin:yMax, xMin:xMax, :);
            
            % --- Step 4: Construct zero frames (same size as cropped ROI) --- 
            zero_frame = zeros(size(cropped_img,1), size(cropped_img,2), 'single'); % 构造零帧 (与单帧大小一致)
            imgstack_final = cat(3, zero_frame, cropped_img, zero_frame); % 在前后各加一帧零图像
            
            % --- Step 5: Convert polygon to mask, cropped to ROI ---
            nuclearMaskFull = poly2mask(polyVerts(:,1), polyVerts(:,2), size(now_img,1), size(now_img,2));
            nuclearMaskCrop = nuclearMaskFull(yMin:yMax, xMin:xMax);
            
            % --- Step 6: Apply mask to cropped stack ---
            for f = 1:size(imgstack_final,3)
                imgstack_final(:,:,f) = imgstack_final(:,:,f) .* single(nuclearMaskCrop);
            end    
            merged_imgstack_final(:,:,:,ich) = imgstack_final;
                    
        end

        img_name = sprintf('%s_Cell%02d.ome.tif', img_list(iImage).name(1:end-4),iCell);
        full_save_path = fullfile(save_img_dir, img_name);
    
        bfsave(merged_imgstack_final,full_save_path,'dimensionOrder', 'XYZCT');
    end
    clearvars roi_info_nuc

end

disp("Imaging ROI cropping finished!")
