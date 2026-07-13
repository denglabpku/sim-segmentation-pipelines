% Script to segment 3D clusters from Multi-SIM of Naxi's Tech.

% In this script, I segments using the raw V or Vnorm transformed Vn for globalOtsu; Vnorm (w/w.o tophat) for LoG,
% V for quantification.

% When import spot coordinate into Imaris, use the exact image dimension as
% voxelSize defined!

% I choose optimal z slice range using visual check of FFT z-stack of raw
% image, such that artifacts are not obvious (no abnormal radial spikes) or 
% FFT high-score higher than 0.7.

% This is the latest segmentation pipeline for 3D SIM for live cell SIM.

% 1. Generate BW from Otsu inside ROI.
% 2. Label BW components.
% 3. Detect LoG peaks.
% 4. Keep only LoG peaks inside BW.
% 5. Remove BW components with zero LoG peaks.
% 6. Use all remaining LoG peaks as watershed markers.
% 7. Run one marker-controlled watershed inside the kept BW mask.
% 8. For each watershed segment, calculate:
%    - parent BW component ID
%    - volume
%    - weighted centroid in um
%    - mean intensity
% 9. Within each parent BW component, fuse watershed segments whose weighted-centroid distance is below the SIM resolution threshold.
% 10. Relabel continuously.
% 11. Remove too-small objects.
% 12. Export only final `stats3D_table.csv`.

%% Imaging setting
close all; clc; clear;
addpath(genpath('C:\Users\zuhui\OneDrive - Peking University\Documents\MATLAB\bfmatlab'),'-end');
addpath(genpath('C:\Users\zuhui\OneDrive - Peking University\Documents\MATLAB\SIM_microscopy\Segmentation_3DSIM'),'-end');

xy_resolution = 0.1;   % um (lateral resolution)
z_resolution  = 0.3;   % um (axial resolution)
pixel_size    = 0.065/2; % um (actual pixel size in XY)
z_step        = 0.135; % um (slice spacing in Z)
Ch_num        = 2; % channel number
Ch_string = ["C488", "C561"];

% Minimal voxel count at resolution limit
minVoxels_at_resolution = (xy_resolution * xy_resolution * z_resolution) / ...
                          (pixel_size * pixel_size * z_step);

% ----- Parameters (imaging) -----
% Store in params.imaging
params.imaging = struct();
params.imaging.xy_resolution          = xy_resolution;
params.imaging.z_resolution           = z_resolution;
params.imaging.pixel_size             = pixel_size;
params.imaging.z_step                 = z_step;
params.imaging.minVoxels_at_resolution = minVoxels_at_resolution;
params.imaging.Ch_num = Ch_num;


%% Load and preprocess
% ----- Paths -----
root_dir = 'E:\Denglab-DataCenter\ChenFei_INTAC_temp\20260612_DLD1_RPB1-JF549_s5p-mintbody488';
in_tif_dir = fullfile(root_dir,'cropped_aligned_imgs');
in_tif_files = dir(fullfile(in_tif_dir,'*.ome.tif')); % choose color-drift aligned image stacks

for iFile = 1:length(in_tif_files)
    fprintf('Start processing file %s: %d/%d\n',in_tif_files(iFile).name,iFile,length(in_tif_files));
    close all;
    
    in_tif = fullfile(in_tif_files(iFile).folder,in_tif_files(iFile).name);
    [filepath, filename, ext] = fileparts(in_tif);filename = filename(1:end-4); % remove .ome
    out_dir = fullfile(filepath, 'seg_result',filename);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end 
    
    % Initialize 
    merged_imgstack_final = [];
    merged_normed_imgstack_final = [];
    merged_normedTophat_imgstack_final = [];
    hfScore_raw = [];

    % Load volume
    Vol = tiffreadVolume(in_tif);               % returns numeric array    
    z_len = size(Vol,3);
    Vol_chx = {Vol(:,:,1:z_len/2),Vol(:,:,(z_len/2+1):end)}; % make sure use same channel as fixed channel of beads image     

    % Plot average frequency and median intensity across the frames
    hfScore_raw = plotZSliceQuality(Vol_chx, Ch_string);    
    
    for iCh = 1:Ch_num
        % --- Paths ---
        name_raw_tif  = fullfile(out_dir,sprintf('%s_good_raw.tif',Ch_string(iCh)));
        name_norm_tif  = fullfile(out_dir,sprintf('%s_good_norm.tif',Ch_string(iCh)));        
        label_tif = fullfile(out_dir,sprintf('%s_good_seglabels.tif',Ch_string(iCh)));

        V = Vol_chx{iCh};  

        if iCh == 1           
            prompt = {'Enter start slice:', 'Enter end slice:'};
            dlgtitle = 'Good Slice Range';
            dims = [1 35]; definput = {'8','15'};opts.WindowStyle = 'normal';
            % Non-modal dialog (does not block figure interaction)
            answer = inputdlg(prompt, dlgtitle, dims, definput, opts);
            
            % Execution still waits until you press OK/Cancel,
            % but you can interact with fig1 while the dialog is open.
            if ~isempty(answer)
                good_slice_range = str2double(answer);
            else
                error('No input provided.');
            end  
        end
                                    
        % ----- Parameters (processing) -----
        params.preprocessing = struct('good_slice_range', good_slice_range,...
                                      'doMedFilt', false, 'medFiltSize', [3 3 3], ...
                                      'doGauss', false, 'gaussSigma', 0.7,...
                                      'doTopHat',true,'topHat_r', 12,...% use TopHat when background is high such as in live-cell SIM, radius should be larger than your interested cluster size 
                                      'sliceNormType','Mean'); % 'Mean'|'Median'; Mean gets better results when correction photobleach induced signal variation across z slices, which get better Otsu after norm across slices.
        params.thresholding = struct('mode', 'globalOtsu_norm', ... % 'globalOtsu_raw' |'globalOtsu_norm' | 'fixed' | 'sliceAdaptive'
                                     'fixedLevel', 500, ...    % used when mode == 'fixed'
                                     'adaptiveWindow', [25 25]); % used when mode == 'sliceAdaptive'
        params.morphology   = struct('minVoxels', minVoxels_at_resolution, 'closeRadius', 0, 'doFill', false); % [50  2  false]
        params.components   = struct('connectivity', 26);
        params.LoG          = struct('sigmaX', 1.3, 'sigmaY', 1.3, 'sigmaZ', 0.9, 'threshold', 2, 'h_value', 0); % guided by real data results, choose the size close to resolution. most cases I fix h_value to 0.
        params.segmentFusion = struct('xyResolution', xy_resolution, 'zResolution', z_resolution, ...
                                      'distanceThreshold', 2.0); % distanceThreshold=1 means 1 unit of xyResolution, I choose 3 such that large cluster are not segmented while small does. visually check.
        params.measure      = struct('voxelSize', pixel_size*pixel_size*z_step,... % [dy dx dz] in pixels
                                      'target', 'norm'); % raw | norm
        params.save         = struct('out_dir', out_dir);
        params.timing       = struct();

        % ----------- loading ROI mask before thresholding ----------- %
        % filename is like: originalName.ome_Cell01
        cellToken = regexp(filename, '_Cell(\d+)$', 'tokens', 'once');
        cellIdx = str2double(cellToken{1});
        baseName = regexprep(filename, '_Cell\d+$', '');
        
        maskFile = fullfile(filepath, sprintf('%s_NuclearMask.mat', baseName));
        S = load(maskFile, 'roi_info_nuc', 'ImHeight', 'ImWidth');
        
        polyVerts = S.roi_info_nuc{cellIdx};
        
        xMin = floor(min(polyVerts(:,1)));
        xMax = ceil(max(polyVerts(:,1)));
        yMin = floor(min(polyVerts(:,2)));
        yMax = ceil(max(polyVerts(:,2)));
        
        nuclearMaskFull = poly2mask(polyVerts(:,1), polyVerts(:,2), S.ImHeight, S.ImWidth);
        roiMask2D = nuclearMaskFull(yMin:yMax, xMin:xMax);
        
        % ----- Check background -----
        % High background leads to Wiener reconstruct artifact 
        hfScore_raw = plotZSliceQuality(Vol_chx, Ch_string, params.preprocessing.good_slice_range);

        % remove coverslip dirt signal        
        V = V(:,:,params.preprocessing.good_slice_range(1):params.preprocessing.good_slice_range(2));
        merged_imgstack_final(:,:,:,iCh) = V;

        % only using ROI region for segmentation
        roiMask3D = repmat(roiMask2D, 1, 1, size(V,3));
        bk_mean = zeros(size(V,3),1);        
        for z = 1:size(V,3)
            sliceVals = V(:,:,z);
            sliceVals = sliceVals(roiMask2D);   % ROI only, no outside zeros
            switch params.preprocessing.sliceNormType
                case 'Mean'
                    bk_mean(z) = mean(sliceVals);
                case 'Median'
                    bk_mean(z) = median(sliceVals);
            end
        end        
        
        % normalize intensity of each slice to reduce illumination variance
        % BUT only use Vnorm to do segments, do not use it for biological
        % intensity quantification!
        % 3D结构光本质上是宽场照明，物镜探测平面和结构光固定不动，依赖样品的上下移动来完成结构光z轴相位的变化。一般来讲信号量是不会变化的，但是当样品过厚时，宽场光进去样品不同深度会有不同程度的衰减，才会出现这样的现象
        Vnorm = zeros(size(V), 'like', V);
        for z = 1:size(V,3)
            if bk_mean(z) > 0
                Vnorm(:,:,z) = V(:,:,z) ./ bk_mean(z);
            else
                Vnorm(:,:,z) = V(:,:,z); % safeguard if background is zero
            end
        end        
        Vnorm(~roiMask3D) = 0;  
        merged_normed_imgstack_final(:,:,:,iCh) = Vnorm; 
        
        % Preprocessing
        if params.preprocessing.doMedFilt
            Vnorm = medfilt3(Vnorm, params.preprocessing.medFiltSize);
        end
        if params.preprocessing.doGauss
            Vnorm = imgaussfilt3(Vnorm, params.preprocessing.gaussSigma);
        end
        if params.preprocessing.doTopHat
            se2d = strel('disk', params.preprocessing.topHat_r);  % tune: maybe 4-10 px
            for z = 1:size(Vnorm,3)
                Vnorm(:,:,z) = imtophat(Vnorm(:,:,z), se2d);
            end   
            merged_normedTophat_imgstack_final(:,:,:,iCh) = Vnorm;
        end

        %% check slice again
        V_test = {Vnorm};
        hfScore_raw = plotZSliceQuality(V_test, Ch_string(iCh));
        %% Thresholding to binary mask
        switch params.thresholding.mode
            % Only GlobalOstu is optimized, other options may should also
            % use Vn instead of V as input.
            case 'globalOtsu_raw'                        
                % Otsu on ROI masked, percentile cleaned pixels, more robust as suggested by codex  
                % V looks better than Vnorm, more close to human inspect
                validVals = V(roiMask3D);
                
                % This is perc
                lo = prctile(validVals, 0.1);
                hi = prctile(validVals, 99.9);
                
                Vn = (V - lo) ./ max(hi - lo, eps);
                Vn = min(max(Vn, 0), 1);
                
                lvl = graythresh(Vn(roiMask3D));
                
                BW = Vn > lvl;
                BW(~roiMask3D) = false;   

            case 'globalOtsu_norm'                        
                % Otsu on ROI masked, percentile cleaned pixels, more robust as suggested by codex  
                % V looks better than Vnorm, more close to human inspect
                validVals = Vnorm(roiMask3D);
                
                % This is perc
                lo = prctile(validVals, 0.1);
                hi = prctile(validVals, 99.9);
                
                Vn = (Vnorm - lo) ./ max(hi - lo, eps);
                Vn = min(max(Vn, 0), 1);
                
                lvl = graythresh(Vn(roiMask3D));
                
                BW = Vn > lvl;
                BW(~roiMask3D) = false;
        
            case 'fixed'
                % Use absolute threshold on original intensities
                BW = Vnorm > params.thresholding.fixedLevel;
        
            case 'sliceAdaptive'
                % Adaptive per-slice using 2D adaptthresh + imbinarize
                BW = false(size(Vnorm), 'logical');
                win = params.thresholding.adaptiveWindow;
                reV = rescale(Vnorm);
                for z = 1:size(Vnorm,3)
                    % robust normalization for each slice
                    % s = V(:,:,z);
                    % lo = prctile(s(:), 1); hi = prctile(s(:), 99);
                    % sn = min(max((s - lo) / max(hi - lo, eps), 0), 1);
                    sn = reV(:,:,z);
                    local_sensitivity = 0.01;
                    T = adaptthresh(sn,local_sensitivity, 'NeighborhoodSize', 101, 'ForegroundPolarity', 'bright');
                    BW(:,:,z) = imbinarize(sn, T);
                end
        
            otherwise
                error('Unknown thresholding mode.');
        end

        
        
        %% Morphological cleanup and component labeling
        % Morphological smoothing/cleanup
        if params.morphology.minVoxels > 0
            % removes all connected components (objects) that have fewer than P pixels from the binary image
            BW = bwareaopen(BW, ceil(params.morphology.minVoxels), params.components.connectivity);
        end
        if params.morphology.closeRadius > 0
            SE = strel('sphere', params.morphology.closeRadius);
            BW = imclose(BW, SE);
        end
        if params.morphology.doFill
            % 3D hole filling is expensive; often unnecessary for clusters
            % If needed, approximate by slice-wise fill:
            for z = 1:size(BW,3)
                BW(:,:,z) = imfill(BW(:,:,z), 'holes');
            end
        end
        
        %% Not use anymore
        % D = -bwdist(~BW); % Distance transform of binary image
        % Lw = watershed(D,8); Lw(~BW) = 0; BWw = BW; BWw(Lw == 0) = false; % refine mask
        % CC = bwconncomp(BWw, 26); L = labelmatrix(CC); % Connected components
        
        % % without watershed
        % CC = bwconncomp(BW, 26); L = labelmatrix(CC); % Connected components
        
        % % Measurements
        % % voxel spacing lets you convert voxel counts to physical volumes
        % stats = regionprops3(CC, V, 'Volume', 'Centroid', 'BoundingBox', 'PrincipalAxisLength');
        % % Scale voxel to match real voxel spacing
        % vs = params.measure.voxelSpacing;  % [dy dx dz]
        % voxelVol = vs(1) * vs(2) * vs(3);
        % stats.Volume = stats.Volume * voxelVol; % Scale voxel
        % stats.VolumePhysical = stats.Volume * params.measure.voxelSize; % Convert Volume (voxels) to um^3
        
        % % Sort by volume (optional)
        % stats = sortrows(stats, 'Volume', 'descend');
    
        %% Localize LoG peaks used as watershed markers        
        % only use Vnorm to do segments, do not use it for biological
        % intensity quantification!
        [spots, quality, I_log] = log_detector_3d_Sage(Vnorm, params.LoG.sigmaX,params.LoG.sigmaY,params.LoG.sigmaZ,params.LoG.threshold, params.LoG.h_value); % current best by trail-and-error
        if isempty(spots)
            spots = zeros(0,3);
            quality = zeros(0,1);
        end
        quality = quality(:);

        fig2 = figure;histogram(quality);xlabel('LoG quality after threshold');
    
        % Convert spot coordinates to linear indices.
        % Note: MATLAB uses [row, col, slice] which corresponds to [y, x, z].
        int_spots = round(spots);
        inBounds = int_spots(:,1) >= 1 & int_spots(:,1) <= size(BW,2) & ...
                   int_spots(:,2) >= 1 & int_spots(:,2) <= size(BW,1) & ...
                   int_spots(:,3) >= 1 & int_spots(:,3) <= size(BW,3);
        spots = spots(inBounds,:);
        int_spots = int_spots(inBounds,:);
        quality = quality(inBounds);
        linearIndices = sub2ind(size(BW), int_spots(:,2), int_spots(:,1), int_spots(:,3));

        % Keep only LoG peaks inside BW. If multiple peaks round to the same
        % voxel, keep the highest-quality peak as the watershed marker.
        spotInMask = BW(linearIndices);
        spots_filtered = spots(spotInMask, :);
        spotLinearIndices = linearIndices(spotInMask);
        quality_filtered = quality(spotInMask);
        if ~isempty(spotLinearIndices)
            [uniqueSpotLinearIndices,~,uniqueGroup] = unique(spotLinearIndices);
            if numel(uniqueSpotLinearIndices) < numel(spotLinearIndices)
                keepUnique = false(numel(spotLinearIndices),1);
                for iUnique = 1:numel(uniqueSpotLinearIndices)
                    members = find(uniqueGroup == iUnique);
                    [~, bestMemberIdx] = max(quality_filtered(members));
                    keepUnique(members(bestMemberIdx)) = true;
                end
                spots_filtered = spots_filtered(keepUnique,:);
                spotLinearIndices = spotLinearIndices(keepUnique);
                quality_filtered = quality_filtered(keepUnique);
            end
        end
    
        %% Watershed from LoG markers, then fuse close watershed segments
        tWatershedFusion = tic;
        CC_all = bwconncomp(BW, params.components.connectivity);
        BW_component_labels_all = labelmatrix(CC_all);
        spotComponentLabels = BW_component_labels_all(spotLinearIndices);

        % Keep only BW components containing at least one in-mask LoG peak.
        componentsWithSpots = unique(spotComponentLabels);
        componentsWithSpots(componentsWithSpots == 0) = [];
        BW_seeded = ismember(BW_component_labels_all, componentsWithSpots);
        CC = bwconncomp(BW_seeded, params.components.connectivity);
        BW_component_labels = labelmatrix(CC);
        spotComponentLabels = BW_component_labels(spotLinearIndices);

        nMarkers = numel(spotLinearIndices);
        markers = zeros(size(BW), 'uint32');
        markers(spotLinearIndices) = uint32(1:nMarkers);
        spotTable = table((1:nMarkers)', spotComponentLabels(:), spotLinearIndices(:), ...
            spots_filtered(:,1), spots_filtered(:,2), spots_filtered(:,3), quality_filtered(:), ...
            'VariableNames', {'MarkerID','ComponentID','LinearIndex','X','Y','Z','Quality'});

        L_watershed = zeros(size(BW), 'uint32');
        if nMarkers > 0            
            topo = -I_log;
            topo(~BW_seeded) = Inf;
            topo_marked = imimposemin(topo, markers > 0);
            Lsplit = watershed(topo_marked);
            Lsplit(~BW_seeded) = 0;

            for iObj = 1:CC.NumObjects
                pixIdx = CC.PixelIdxList{iObj};
                markerRows = find(spotComponentLabels == iObj);
                markerBasins = Lsplit(spotLinearIndices(markerRows));
                assignedBasins = zeros(numel(markerRows), 1);
                nAssignedBasins = 0;

                for iMarker = 1:numel(markerRows)
                    basinLabel = markerBasins(iMarker);
                    if basinLabel > 0 && ~ismember(basinLabel, assignedBasins(1:nAssignedBasins))
                        basinMask = Lsplit(pixIdx) == basinLabel;
                        L_watershed(pixIdx(basinMask)) = uint32(markerRows(iMarker));
                        nAssignedBasins = nAssignedBasins + 1;
                        assignedBasins(nAssignedBasins) = basinLabel;
                    end
                end
                
                % Must include this, otherwise, later fusing will not work
                % properly
                orphanIdx = pixIdx(L_watershed(pixIdx) == 0);
                if ~isempty(orphanIdx)
                    nearestMarkerRows = nearestMarkersByResolutionDistance( ...
                        orphanIdx, spots_filtered, markerRows, size(BW), ...
                        pixel_size, z_step, ...
                        params.segmentFusion.xyResolution, params.segmentFusion.zResolution);
                    L_watershed(orphanIdx) = uint32(nearestMarkerRows);
                end
            end
        end

        rawSegmentLabels = unique(L_watershed(:));
        rawSegmentLabels(rawSegmentLabels == 0) = [];
        if isempty(rawSegmentLabels)
            L_watershed = zeros(size(BW), 'uint32');
            segmentParentComponentIDs = zeros(0,1);
        else
            [~, loc] = ismember(L_watershed, rawSegmentLabels);
            L_watershed = uint32(loc);
            segmentParentComponentIDs = spotComponentLabels(double(rawSegmentLabels));
        end

        [L, segmentFusionTable, dResolution_segment] = fuseSegmentsByWeightedCentroid( ...
            L_watershed, V, segmentParentComponentIDs, ...
            pixel_size, z_step, ...
            params.segmentFusion.xyResolution, params.segmentFusion.zResolution, ...
            params.segmentFusion.distanceThreshold);

        figure;
        histogram(dResolution_segment);
        xlabel('Resolution-normalized weighted-centroid distance');
        ylabel('Frequency');
        title('Histogram of Pairwise Segment Centroid Distances inside one BW label');

        params.timing.watershedAndCentroidFusion = toc(tWatershedFusion);
        fprintf('%s watershed + centroid fusion: %.2f s\n', Ch_string(iCh), params.timing.watershedAndCentroidFusion);
    
        % Step 5: Measurements the each labeled objects
        % voxel spacing lets you convert voxel counts to physical volumes
        tInitialStats = tic;
        switch params.measure.target
            case 'raw'
                measure_target = V;
            case 'norm'
                measure_target = Vnorm;
        end
        stats = regionprops3(L, measure_target, 'Volume', 'Centroid', 'BoundingBox', 'PrincipalAxisLength');
        params.timing.initialRegionprops3 = toc(tInitialStats);
        fprintf('%s initial regionprops3: %.2f s\n', Ch_string(iCh), params.timing.initialRegionprops3);
        % Scale voxel to match real voxel spacing
        % vs = params.measure.voxelSpacing;  % [dy dx dz]
        % voxelVol = vs(1) * vs(2) * vs(3);
        % stats.Volume = stats.Volume * voxelVol; % Scale voxel
        stats.VolumePhysical = stats.Volume * params.measure.voxelSize; % Convert Volume (voxels) to um^3
           
        %% Remove clusters that are not good
    
        % Define minimum physical volume threshold and remove z-edge objects.
        minVoxelSize = ceil(params.morphology.minVoxels);   % example in µm^3, adjust to your dataset

        % -----------------------------------------------------------
        % Remove edge-touching objects, currently not used since will
        % remove some obvious large clusters.
        % -----------------------------------------------------------
        % edgeLabels = unique([L(:,:,1); L(:,:,end)]);
        % edgeLabels(edgeLabels == 0) = [];         
        % numLabels = height(stats);
        % edgeIdx = false(numLabels,1);
        % edgeLabels = double(edgeLabels(edgeLabels <= numLabels));
        % edgeIdx(edgeLabels) = true;        
        % Find indices of components that meet the size criterion and are not edge-truncated.
        % validIdx = stats.Volume >= minVoxelSize & ~edgeIdx;   

        validIdx = stats.Volume >= minVoxelSize;
    
        validLabels = find(validIdx);
        [~, loc] = ismember(L, validLabels);
        L_filtered = uint32(loc);
        
        % -----------------------------------------------------------
        % Remove edge-touching objects, currently not used since will
        % remove some obvious large clusters.
        % -----------------------------------------------------------
        % assert(isempty(setdiff(unique([L_filtered(:,:,1); L_filtered(:,:,end)]), 0)), ...
        %     'Edge-truncated labels remain after filtering.');
        finalLabels = unique(L_filtered(:));
        finalLabels(finalLabels == 0) = [];
        if isempty(finalLabels)
            assert(max(L_filtered(:)) == 0, 'Empty final labels should contain only background.');
        else
            assert(isequal(finalLabels, uint32((1:max(L_filtered(:)))')), ...
                'L_filtered labels must be continuous from 1.');
        end
    
        % Measurements
        % voxel spacing lets you convert voxel counts to physical volumes
        tFinalStats = tic;        
        stats3D_table = regionprops3(L_filtered, measure_target, 'Volume', 'Centroid', 'BoundingBox', 'PrincipalAxisLength','WeightedCentroid','MeanIntensity');
        params.timing.finalRegionprops3 = toc(tFinalStats);
        fprintf('%s final regionprops3: %.2f s\n', Ch_string(iCh), params.timing.finalRegionprops3);
        % Scale voxel to match real voxel spacing
        stats3D_table.VolumePhysical = stats3D_table.Volume * params.measure.voxelSize; % Convert Volume (voxels) to um^3 
        stats3D_table.ClusterID = (1:height(stats3D_table))';
        stats3D_table = movevars(stats3D_table, 'ClusterID', 'Before', 1);
        if height(stats3D_table) > 0
            stats3D_table.WeightedCentroid_x_um = (stats3D_table.WeightedCentroid(:,1)-0.5).*pixel_size;
            stats3D_table.WeightedCentroid_y_um = (stats3D_table.WeightedCentroid(:,2)-0.5).*pixel_size;
            stats3D_table.WeightedCentroid_z_um = (stats3D_table.WeightedCentroid(:,3)-0.5).*z_step;
        else
            stats3D_table.WeightedCentroid_x_um = zeros(0,1);
            stats3D_table.WeightedCentroid_y_um = zeros(0,1);
            stats3D_table.WeightedCentroid_z_um = zeros(0,1);
        end

        numClusters = double(max(L_filtered(:)));
        assert(height(stats3D_table) == numClusters, 'stats3D_table row count must match final labels.');
        assert(isequal(stats3D_table.ClusterID(:), (1:numClusters)'), ...
            'stats3D_table.ClusterID must match final label IDs.');
        stats3D_table.equiRadius = ((3*stats3D_table.VolumePhysical)/(4*pi)).^(1/3);

        
        
        % % Remove Rod‑like segments artifacts
        % axisLens = stats.PrincipalAxisLength; % Nx3 matrix [major, intermediate, minor]
        % elongationRatio = axisLens(:,1) ./ axisLens(:,3); % major/minor
        % stats.elongationRatio = elongationRatio;
        % elongationCutoff = 5; % discard if major/minor > 5
        % 
        % % Remove clusters below resolution limit
        % minVoxelSize =  ceil(params.morphology.minVoxels);
        % 
        % validIdx = stats.Volume >= minVoxelSize & ...
        %            elongationRatio <= elongationCutoff & ...
        %            stats.elongVersusVol >= 0.005;
        % 
        % CC_valid = CC;
        % CC_valid.PixelIdxList = CC.PixelIdxList(validIdx);
        % CC_valid.NumObjects   = sum(validIdx);
        % 
        % L_filtered = labelmatrix(CC_valid);
        % stats_filtered = stats(validIdx,:);
    
        %% 2D max projection each clusters to avoid z axis elongation
        % Assume:
        % L_filter : 3D label matrix (filtered, sequential labels)
        % voxelSize: [dy dx dz] in µm
        
        numClusters = double(max(L_filtered(:)));
        if numClusters == 0
            stats2D_table = table(zeros(0,1), zeros(0,1), ...
                'VariableNames', {'ClusterID','EquivRadius_um'});
        else
            stats2D = struct([]);
                    
            for k = 1:numClusters
                % Step 1: Extract cluster mask
                mask3D = (L_filtered == k);
            
                if ~any(mask3D(:))
                    continue; % skip empty labels
                end
            
                % Step 2: Max projection along Z
                mask2D = max(mask3D, [], 3);
            
                % Step 3: Measure 2D properties
                area_2d = sum(mask2D,'all');        
            
                % Step 4: Scale to physical units
                stats2D(k).ClusterID          = k;
                stats2D(k).EquivRadius_um   = sqrt(area_2d * params.imaging.pixel_size^2 / pi);        
            end
            
            % Convert to table for convenience
            stats2D_table = struct2table(stats2D);
        end
    
        %% Visualization, QC, and saving
        % --- Paths ---
        % name_raw_tif  = fullfile(out_dir,sprintf('%s_good_raw.tif',Ch_string(iCh)));
        % name_norm_tif  = fullfile(out_dir,sprintf('%s_good_norm.tif',Ch_string(iCh)));
        % label_tif = fullfile(out_dir,sprintf('%s_good_seglabels.tif',Ch_string(iCh)));        
        
        % --- Save bkCor_raw intensity (float32) ---
        % writeFloat32Tiff(single(V), name_raw_tif);
        % writeFloat32Tiff(single(Vnorm), name_norm_tif);
        % merged_imgstack_final(:,:,:,iCh) = V;
        % merged_normed_imgstack_final(:,:,:,iCh) = Vnorm;        
        
        % --- Save labels (uint32) ---
        writeUint32Tiff(uint32(L_filtered), label_tif);
        
        save(fullfile(out_dir,sprintf("%s_segResults.mat",Ch_string(iCh))), ...
            "params","stats3D_table","stats2D_table","spotTable","segmentFusionTable");
    
        writetable(stats3D_table,fullfile(out_dir,sprintf('%s_stats3D_table.csv',Ch_string(iCh))));
        
        %% Quantification
        fig3 = figure;
        tiledlayout(3,1);
        nexttile;
        histogram(stats3D_table.Volume);
        xlabel('Voxels'); ylabel('Counts');
        nexttile;
        histogram(stats3D_table.VolumePhysical);
        xlabel('\mum^3'); ylabel('Counts');
        nexttile;
        histogram(stats3D_table.equiRadius);
        xlabel('\mum'); ylabel('Counts');

    end

    % Export dual-channel images        
    bfsave(single(merged_imgstack_final),fullfile(out_dir, sprintf('%s_good_raw.ome.tif',"Merge")),'dimensionOrder', 'XYZCT');
    bfsave(single(merged_normed_imgstack_final),fullfile(out_dir,sprintf('%s_good_norm.ome.tif',"Merge")),'dimensionOrder', 'XYZCT'); 
    if params.preprocessing.doTopHat
        bfsave(single(merged_normedTophat_imgstack_final),fullfile(out_dir,sprintf('%s_good_normTophat.ome.tif',"Merge")),'dimensionOrder', 'XYZCT'); 
    end

    f = msgbox(["Check result, if everything okay, continue..."]);
    uiwait(f);

end
disp('All files have been processed!')

% %% Visualization, QC, and saving
% % QC: show binary and labels
% figure('Name','Binary mid-slices');
% slice(double(BW), size(BW,2)/2, size(BW,1)/2, size(BW,3)/2);
% colormap(gray); shading interp; axis image ij;
% 
% % 3D rendering (binary)
% figure('Name','3D isosurface (binary)');
% p = patch(isosurface(BW, 0.5)); p.FaceColor = [0.9 0.1 0.1]; p.EdgeColor = 'none';
% camlight; lighting gouraud; axis equal tight off;
% 
% % 3D rendering (intensity masked by BW)
% figure('Name','Masked intensity isosurface');
% Vm = V; Vm(~BW) = 0;
% isoLevel = prctile(Vm(Vm>0), 70); % adjust as needed
% p2 = patch(isosurface(Vm, isoLevel)); p2.FaceColor = [0.2 0.6 1]; p2.EdgeColor = 'none';
% camlight; lighting gouraud; axis equal tight off;
% 
% % Volumetric label visualization (requires Image Processing Toolbox R2019b+)
% try
%     volshow(L); title('Labels');
% catch
%     % Fallback: slice montage
%     figure('Name','Label montage'); montage(uint16(L), 'Indices', round(linspace(1,size(L,3),16)));
% end
% 
% % Save outputs
% [~, base] = fileparts(in_tif);
% if params.save.saveBinary
%     bw_path = fullfile(out_dir, [base '_seg_binary.tif']);
%     writeFloat32Tiff(single(BW), bw_path); % preserve float32 for compatibility; BW is logical
% end
% if params.save.saveLabel
%     lbl_path = fullfile(out_dir, [base '_seg_labels.tif']);
%     writeUint32Tiff(uint32(L), lbl_path);  % labels often exceed uint16 in big volumes
% end
% if params.save.saveStats
%     stats_path = fullfile(out_dir, [base '_seg_stats.csv']);
%     writetable(stats, stats_path);
% end
% 
% 
% sx = 1;
% sy= 1;
% sz = 135/65/2;
% A = [sx 0 0 0; 0 sy 0 0; 0 0 sz 0; 0 0 0 1];
% tform = affinetform3d(A);
% 
% vol = volshow(L,Transformation=tform);

% ----------------------- Use to check LoG peaks  ----------------------- %
% % --- Configuration ---
% filename = fullfile(out_dir,'zstacked_filtered_output0.tif');
% h = 0; % Define your 'h' value for the title string
% [img_height, img_width, z_max] = size(V); % Get original image dimensions
% 
% % Delete the file if it already exists to prevent adding to old runs
% if exist(filename, 'file')
%     delete(filename);
% end
% 
% % Create a hidden figure with a fixed pixel size matching the image
% fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, img_width, img_height]); 
% 
% % Create an axes that fills 100% of the figure with no default boundaries
% ax = axes(fig, 'Position', [0 0 1 1]); 
% 
% % --- Main Export Loop (Increasing Z) ---
% for z = 1:z_max
%     % 1. Render the custom frame directly into the targeted axes
%     imshow(V(:,:,z), [0 6500], 'Parent', ax); 
%     hold(ax, 'on');
% 
%     filt_spots = spots_filtered(spots_filtered(:,3) >= z & spots_filtered(:,3) < z+1, :);
%     if ~isempty(filt_spots)
%         % Plot the scatter points onto the matching axes
%         scatter(ax, filt_spots(:,1), filt_spots(:,2), 8, 'red', 'filled');
%     end
% 
%     % Ensure data limits perfectly track the pixel boundaries
%     xlim(ax, [0.5, img_width + 0.5]);
%     ylim(ax, [0.5, img_height + 0.5]);
%     hold(ax, 'off');
% 
%     % Force MATLAB to completely render the frame layout
%     drawnow; 
% 
%     % 2. Capture only the axes contents (bypasses window frames and gray margins)
%     frame = getframe(ax);
%     im = frame.cdata;
% 
%     % Force dimensions: Handles high-DPI scaling (Mac Retina/4K) if present
%     if size(im, 1) ~= img_height || size(im, 2) ~= img_width
%         im = imresize(im, [img_height, img_width]);
%     end
% 
%     % 3. Write/Append the frame into the Z-stack TIFF
%     if z == 1
%         imwrite(im, filename, 'tiff');
%     else
%         imwrite(im, filename, 'tiff', 'WriteMode', 'append');
%     end
% end
% 
% % Clean up by closing the hidden figure workspace
% close(fig);
% disp('Z-stacked TIFF export complete.');


%% Auxiliary function

function hfScore_raw = plotZSliceQuality(Vol_chx, Ch_string, varargin)

    if nargin == 3
        x_line_range = varargin{1};        
    end
    % PLOTZSLICEQUALITY Plots background noise vs FFT quality metrics across Z-slices
    Ch_num = length(Vol_chx);    
    
    figure;
    % Get the default MATLAB color palette to ensure each channel gets a unique color
    colors = get(gca, 'ColorOrder');
    
    for iCh = 1:Ch_num  
        V = Vol_chx{iCh}; 
        
        % Pick a color for this channel (loop around if Ch_num > 7)
        colorIdx = mod(iCh-1, size(colors,1)) + 1;
        current_color = colors(colorIdx, :);
        
        % ----- 1. Left Axis: Background Intensity -----       
        bk_mean = median(V,[1 2]);
        bk_mean = bk_mean(:);
        
        yyaxis left
        % Store the handle, apply the channel color, use '--o' for left data
        p_left(iCh) = plot(1:length(bk_mean), bk_mean/max(bk_mean), ...
                           'LineStyle', '--', 'Marker', 'o', 'Color', current_color, 'LineWidth', 1.5);   
        hold on; % Hold on must be applied to the active axis
    
        % ----- 2. Right Axis: High-Frequency Fourier Score -----
        hfScore = zeros(size(V,3), 1); % Fixed: preallocate 1 column per channel loop
        for z = 1:size(V,3)
            slice = single(V(:,:,z));
            % 1. Create the 2D Hanning Window matching your image size,
            % this is to remove cruciform artifacts or FFT leakage at image
            % edge
            [rows, cols] = size(slice);
            win_x = hann(cols); % Horizontal window vector
            win_y = hann(rows); % Vertical window vector
            win_2d = win_y * win_x'; % Full 2D window grid
            
            % 2. Apply the window to the image
            slice_windowed = slice .* win_2d;

            F = fftshift(fft2(slice_windowed));
            P = abs(F).^2;
            P = log1p(P);
            hfScore(z) = mean(P,"all");
        end
        hfScore_raw(:,iCh) = hfScore;

        yyaxis right
        % Store the handle, apply the same channel color, use '-s' for right data
        p_right(iCh) = plot(1:length(bk_mean), hfScore/max(hfScore), ...
                            'LineStyle', '-', 'Marker', 's', 'Color', current_color, 'LineWidth', 1.5);  
        hold on;
    end
    
    % ----- 3. Formatting & Axis-Specific Labels -----
    yyaxis left
    ylabel('Normalized Median Intensity (Dashed, --o)');
    % Set left y-axis spine to a neutral dark gray color instead of matching the last line
    ax = gca; ax.YAxis(1).Color = [0.15 0.15 0.15]; 
    
    yyaxis right
    ylabel('Normalized FFT High-Freq Score (Solid, -s)');
    ax.YAxis(2).Color = [0.15 0.15 0.15];
    
    xlabel('Z slices');     
    grid on;
    
    % ----- 4. Intelligent Dual-Layer Legend -----
    % Create distinct legend entries combining Channel Name + Metric Type
    legend_labels = cell(1, Ch_num * 2);
    for iCh = 1:Ch_num
        legend_labels{iCh} = [Ch_string{iCh} ' (Bk Mean)'];
        legend_labels{Ch_num + iCh} = [Ch_string{iCh} ' (FFT Score)'];
    end
    
    if nargin == 3
        xline(x_line_range(1),'r--', 'LineWidth', 2);   % red dashed line
        xline(x_line_range(2),'r--', 'LineWidth', 2);   % red dashed line
        title(sprintf('Z slices outside two lines will be discarded.'));
    end
    
    % Feed all line handles and strings into a single consolidated legend
    legend([p_left, p_right], legend_labels, 'Location', 'best'); 
end

function [L_fused, segmentFusionTable, dResolution_log] = fuseSegmentsByWeightedCentroid( ...
    L_watershed, V, parentComponentIDs, pixelSize, zStep, xyResolution, zResolution, distanceThreshold)

    dResolution_log = [];
    nSegments = double(max(L_watershed(:)));
    if nSegments == 0
        L_fused = zeros(size(L_watershed), 'uint32');
        segmentFusionTable = table();
        return;
    end

    parentComponentIDs = parentComponentIDs(:);
    statsPre = regionprops3(L_watershed, V, 'Volume', 'Centroid', 'WeightedCentroid');
    weightedCentroids = statsPre.WeightedCentroid;

    % WeightedCentroid can be undefined if a segment has zero total intensity.
    fallbackIdx = any(isnan(weightedCentroids), 2);
    weightedCentroids(fallbackIdx,:) = statsPre.Centroid(fallbackIdx,:);

    parent = 1:nSegments;
    componentIDs = unique(parentComponentIDs);
    componentIDs(componentIDs == 0) = [];

    for iComponent = 1:numel(componentIDs)
        segmentIDs = find(parentComponentIDs == componentIDs(iComponent));
        if numel(segmentIDs) < 2
            continue;
        end

        for iSegment = 1:numel(segmentIDs)-1
            segI = segmentIDs(iSegment);
            for jSegment = iSegment+1:numel(segmentIDs)
                segJ = segmentIDs(jSegment);

                dx = (weightedCentroids(segI,1) - weightedCentroids(segJ,1)) * pixelSize / xyResolution;
                dy = (weightedCentroids(segI,2) - weightedCentroids(segJ,2)) * pixelSize / xyResolution;
                dz = (weightedCentroids(segI,3) - weightedCentroids(segJ,3)) * zStep / zResolution;
                dResolution = sqrt(dx^2 + dy^2 + dz^2);
                dResolution_log(end+1,1) = dResolution; %#ok<AGROW>

                if dResolution <= distanceThreshold
                    [rootI, parent] = findRoot(parent, segI);
                    [rootJ, parent] = findRoot(parent, segJ);
                    if rootI ~= rootJ
                        parent(rootJ) = rootI;
                    end
                end
            end
        end
    end

    roots = zeros(nSegments, 1);
    for iSegment = 1:nSegments
        [roots(iSegment), parent] = findRoot(parent, iSegment);
    end

    groupRoots = unique(roots, 'stable');
    segmentToFused = zeros(nSegments, 1);
    for iGroup = 1:numel(groupRoots)
        segmentToFused(roots == groupRoots(iGroup)) = iGroup;
    end

    labelMap = zeros(nSegments+1, 1, 'uint32');
    labelMap(2:end) = uint32(segmentToFused);
    L_fused = labelMap(double(L_watershed) + 1);

    segmentFusionTable = table((1:nSegments)', parentComponentIDs, statsPre.Volume, ...
        weightedCentroids(:,1), weightedCentroids(:,2), weightedCentroids(:,3), segmentToFused, ...
        'VariableNames', {'SegmentID','ParentBWComponentID','Volume','WeightedCentroidX', ...
        'WeightedCentroidY','WeightedCentroidZ','FusedLabelID'});
end

function [root, parent] = findRoot(parent, idx)
    root = idx;
    while parent(root) ~= root
        root = parent(root);
    end

    while parent(idx) ~= idx
        nextIdx = parent(idx);
        parent(idx) = root;
        idx = nextIdx;
    end
end

function nearestMarkerRows = nearestMarkersByResolutionDistance(voxelIdx, markerSpots, markerRows, volSize, pixelSize, zStep, xyResolution, zResolution)
    if isempty(voxelIdx)
        nearestMarkerRows = zeros(0,1);
        return;
    end

    [vy, vx, vz] = ind2sub(volSize, voxelIdx);
    vx = double(vx);
    vy = double(vy);
    vz = double(vz);

    nearestD2 = inf(numel(voxelIdx), 1);
    nearestMarkerRows = zeros(numel(voxelIdx), 1);

    for iMarker = 1:numel(markerRows)
        markerRow = markerRows(iMarker);
        dx = (vx - markerSpots(markerRow,1)) .* pixelSize ./ xyResolution;
        dy = (vy - markerSpots(markerRow,2)) .* pixelSize ./ xyResolution;
        dz = (vz - markerSpots(markerRow,3)) .* zStep ./ zResolution;
        d2 = dx.^2 + dy.^2 + dz.^2;

        updateIdx = d2 < nearestD2;
        nearestD2(updateIdx) = d2(updateIdx);
        nearestMarkerRows(updateIdx) = markerRow;
    end
end
