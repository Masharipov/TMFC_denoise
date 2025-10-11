% Example of TMFC_denoise usage #2: Event-related design
%
% Face repetition fMRI dataset from the SPM website: 
% http://www.fil.ion.ucl.ac.uk/spm/data/face_rep/
%==========================================================================

clear

%% Directory containing the facerep data
%--------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end
data_path = fullfile(data_path,'facerep');

%% Select SPM.mat files
%--------------------------------------------------------------------------
% Option 1 - Command line:
nSub = 1; % This example uses data from a single subject
for iSub = 1:nSub
    SPM_paths{iSub,1} = fullfile(data_path,['sub-' sprintf('%02d', iSub)],'GLM_categorical','SPM.mat');  % Paths to SPM.mat files
    subject_paths{iSub,1} = fullfile(data_path,['sub-' sprintf('%02d', iSub)]); % Paths to subject folders containing GLMs
end

% Option 2 - GUI:
% Select the sub-01 folder and the SPM.mat file inside the GLM subfolder
% [SPM_paths,subject_paths] = tmfc_select_subjects_GUI(0);

%% Define denoising options
%--------------------------------------------------------------------------
% Option 1 - Command line:
% Default parameters have been changed for this dataset
options.motion = '6HMP';            % Alternatives: 12HMP, 24HMP (Default)
options.translation_idx = [1 2 3];  % For FSL/AFNI use [4 5 6]
options.rotation_idx = [4 5 6];     % For FSL/AFNI use [1 2 3]
options.rotation_unit = 'rad';      % For HCP/AFNI use 'deg'
options.head_radius = 50;           % Human head radius in [mm]
options.DVARS = 1;                  % DVARS calculation: 1 = enabled (default), 0 = disabled
options.aCompCor = [0.5 0.5];       % [0.5 0.5] for aCompcor50;
                                    % [5 5] for fixed PCs aCompCor (Default: five WM PCs, five CSF PCs);
                                    % [0 0] do not calculate aCompCor
options.aCompCor_ort = 1;           % Disable pre-orthogonalization = 0
options.rWLS = 0;                   % rWLS estimation: 1 = enabled, 0 = disabled
options.spikereg = 0;               % Spike regression: 1 = enabled, 0 = disabled
options.spikeregFDthr = 0.5;        % Select FD threshold for spike regression
options.WM_CSF = 'none';            % Alternatives: 2Phys, 4Phys, 8Phys
options.GSR = 'none';               % Alternatives: GSR, 2GSR, 4GSR 
options.parallel = 0;               % Parallel computations: 1 = enabled, 0 = disabled

% Option 2 - GUI:
% options = tmfc_denoise_options_GUI;

%% Select structural T1 images in native space
%--------------------------------------------------------------------------
% Option 1 - Command line:
anat_subfolder = 'Structural';
for iSub = 1:length(SPM_paths)
    anat_file = spm_select('ExtList', fullfile(data_path,['sub-' sprintf('%02d', iSub)],anat_subfolder), '^sM.*\.img$');
    anat_paths{iSub,1} = fullfile(data_path,['sub-' sprintf('%02d', iSub)], ...
        anat_subfolder,anat_file);  % Paths to structural images
end

% Option 2 - GUI:
% Select the 'Structural' folder and apply the 'sM*.img' text filter
% anat_paths = tmfc_select_anat_GUI(subject_paths);

%% Select realigned and unsmoothed functional images in MNI space
%--------------------------------------------------------------------------
% Option 1 - Command line:
func_subfolder = 'RawEPI';
for iSub = 1:length(SPM_paths)
    func_files = spm_select('ExtList', fullfile(data_path,['sub-' sprintf('%02d', iSub)],func_subfolder), '^wars.*\.img$');
    for iScan = 1:length(func_files)
        func_paths.fname{iScan,1} = fullfile(data_path,['sub-' sprintf('%02d', iSub)], ...
            func_subfolder,deblank(func_files(iScan,:)));  % Paths to functional images
    end
end

% Option 2 - GUI:
% Select the 'RawEPI' folder and apply the '^wars*.img' text filter
% func_paths = tmfc_select_func_GUI(SPM_paths,subject_paths);

%% Create TMFC_denoise subfolders
%--------------------------------------------------------------------------
% Delete previously created TMFC_denoise files?
clear_all = 0; % Do not delete

for iSub = 1:length(SPM_paths)
    GLM_subfolder = fileparts(SPM_paths{iSub});
    if clear_all == 1 && exist(fullfile(GLM_subfolder,'TMFC_denoise'),'dir')
        rmdir(fullfile(GLM_subfolder,'TMFC_denoise'),'s');
    end
    if ~exist(fullfile(GLM_subfolder,'TMFC_denoise'),'dir')
        mkdir(fullfile(GLM_subfolder,'TMFC_denoise'));
    end
    clear GLM_subfolder
end

%% Calculate head motion parameters (HMP) and framewise displacement (FD)
%--------------------------------------------------------------------------
FD = tmfc_head_motion(SPM_paths,subject_paths,options);

%% Plot FD time series
%--------------------------------------------------------------------------
% Option 1: Allows saving group FD statistics only (using the 'Save' button)
% tmfc_plot_FD(FD);

% Option 2: Allows saving group FD statistics and TMFC_denoise settings (using the 'Save' button)
tmfc_plot_FD(FD,options,SPM_paths,subject_paths,anat_paths,func_paths);

%% Create spike regressors
%--------------------------------------------------------------------------
% Select the FD threshold for spike regression
% Option 1 - Command line:
% options.spikeregFDthr = FDthr;
% Option 2 - GUI:
% FDthr = tmfc_plot_FD(FD);

% Create regressors using the selected FD threshold:
if options.spikereg == 1
    tmfc_spikereg(SPM_paths,options);
end

%% Create GM/WM/CSF and whole-brain (WB) masks
%--------------------------------------------------------------------------
% Masking parameters:
% (1) Probability thresholds for GM, WM, and CSF maps.
% (2) Number of dilation/erosion cycles for GM, WM, and CSF masks.

% Option 1 - Command line:
% Default parameters have been changed for this dataset
options.GMmask.prob = 0.99;   % Default: 0.95
options.WMmask.prob = 0.95;   % Default: 0.99
options.CSFmask.prob = 0.95;  % Default: 0.99
options.GMmask.dilate = 1;    % Default: 2
options.WMmask.erode = 2;     % Default: 3
options.CSFmask.erode = 2;

% Option 2 - GUI:
% [options.GMmask.prob, options.WMmask.prob, options.CSFmask.prob, ...
%  options.GMmask.dilate, options.WMmask.erode, options.CSFmask.erode] = tmfc_masks_GUI();

% Create masks:
if sum(options.aCompCor)~=0 || ~strcmpi(options.WM_CSF,'none') || ~strcmpi(options.GSR,'none')  || options.DVARS == 1
    masks = tmfc_create_masks(SPM_paths,anat_paths,func_paths,options);
end

%% Calculate physiological regressors
%--------------------------------------------------------------------------
if sum(options.aCompCor)~=0 || ~strcmpi(options.WM_CSF,'none') || ~strcmpi(options.GSR,'none')
    tmfc_physioreg(SPM_paths,subject_paths,func_paths,masks,options);
end

%% Estimate updated GLMs with noise regressors
%--------------------------------------------------------------------------
output_paths = tmfc_estimate_updated_GLMs(SPM_paths,masks,options);

%% Calculate and plot DVARS
%--------------------------------------------------------------------------
% Calculate DVARS
if options.DVARS == 1
    [preDVARS,postDVARS] = tmfc_calculate_DVARS(FD,SPM_paths,options,masks,output_paths);

    % Plot DVARS
    % Option 1: Allows saving group FD-DVARS statistics only (using the 'Save' button)
    % tmfc_plot_DVARS(preDVARS,postDVARS,FD);
    
    % Option 2: Allows saving group FD-DVARS statistics and TMFC_denoise settings (using the 'Save' button)
    tmfc_plot_DVARS(preDVARS,postDVARS,FD,options,SPM_paths,subject_paths,anat_paths,func_paths,masks);
end
