% Prepare face repetition fMRI dataset from the SPM website: 
% http://www.fil.ion.ucl.ac.uk/spm/data/face_rep/
%
% The dataset is described in the SPM manual:
% http://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf#Chap:data:faces
%
% Reference:
% Henson, R.N.A., Shallice, T., Gorno-Tempini, M.-L. and Dolan, R.J. (2002)
% Face repetition effects in implicit and explicit memory tests as 
% measured by fMRI. Cerebral Cortex, 12, 178-186.
%==========================================================================

clear

%==========================================================================
% Data preparation steps are adapted from the face_rep_spm12_batch.m sript:
% 
% Directory containing the facerep data
%--------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end
data_path = fullfile(data_path,'facerep');
if ~exist(data_path,'dir'); mkdir(data_path); end
fprintf('%-40s:', 'Downloading face repetition dataset from the SPM website (58 MB) ...'); 
urlwrite('http://www.fil.ion.ucl.ac.uk/spm/download/data/face_rep/face_rep.zip',fullfile(data_path,'facerep.zip'));
unzip(fullfile(data_path,'facerep.zip'),data_path);
movefile(fullfile(data_path,'face_rep'),fullfile(data_path,'sub-01'));
data_path = fullfile(data_path,'sub-01');
fprintf(' %20s\n', '...done');

% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPATIAL PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

% Select functional and structural scans
%--------------------------------------------------------------------------
f = spm_select('FPList', fullfile(data_path,'RawEPI'), '^sM.*\.img$');
a = spm_select('FPList', fullfile(data_path,'Structural'), '^sM.*\.img$');

% Realign
%--------------------------------------------------------------------------
matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = cellstr(f);

% Slice Timing Correction
%--------------------------------------------------------------------------
matlabbatch{2}.spm.temporal.st.scans{1} = cellstr(spm_file(f,'prefix','r'));
matlabbatch{2}.spm.temporal.st.nslices = 24;
matlabbatch{2}.spm.temporal.st.tr = 2;
matlabbatch{2}.spm.temporal.st.ta = 2-2/24;
matlabbatch{2}.spm.temporal.st.so = 24:-1:1;
matlabbatch{2}.spm.temporal.st.refslice = 12;

% Coregister
%--------------------------------------------------------------------------
matlabbatch{3}.spm.spatial.coreg.estimate.ref    = cellstr(spm_file(f(1,:),'prefix','mean'));
matlabbatch{3}.spm.spatial.coreg.estimate.source = cellstr(a);

% Segment
%--------------------------------------------------------------------------
matlabbatch{4}.spm.spatial.preproc.channel.vols  = cellstr(a);
matlabbatch{4}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{4}.spm.spatial.preproc.warp.write    = [0 1];

% Normalise: Write
%--------------------------------------------------------------------------
matlabbatch{5}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample = cellstr(char(spm_file(f,'prefix','ar'),spm_file(f(1,:),'prefix','mean')));
matlabbatch{5}.spm.spatial.normalise.write.woptions.vox  = [3 3 3];

matlabbatch{6}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
matlabbatch{6}.spm.spatial.normalise.write.subj.resample = cellstr(spm_file(a,'prefix','m','ext','nii'));
matlabbatch{6}.spm.spatial.normalise.write.woptions.vox  = [1 1 1.5];

% Smooth
%--------------------------------------------------------------------------
matlabbatch{7}.spm.spatial.smooth.data = cellstr(spm_file(f,'prefix','war'));
matlabbatch{7}.spm.spatial.smooth.fwhm = [8 8 8];

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASSICAL STATISTICAL ANALYSIS (CATEGORICAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

% Load onsets
%--------------------------------------------------------------------------
onsets    = load(fullfile(data_path,'sots.mat'));
condnames = {'N1' 'N2' 'F1' 'F2'};

% Output Directory
%--------------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM_categorical';

% Model Specification
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM_categorical'));
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t = 24;
matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t0 = 12;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans = cellstr(spm_file(f,'prefix','swar'));
for i=1:numel(condnames)
    matlabbatch{2}.spm.stats.fmri_spec.sess.cond(i).name = condnames{i};
    matlabbatch{2}.spm.stats.fmri_spec.sess.cond(i).onset = onsets.sot{i};
    matlabbatch{2}.spm.stats.fmri_spec.sess.cond(i).duration = 0;
end
matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg   = cellstr(spm_file(f(1,:),'prefix','rp_','ext','.txt'));
matlabbatch{2}.spm.stats.fmri_spec.fact(1).name     = 'Fame';
matlabbatch{2}.spm.stats.fmri_spec.fact(1).levels   = 2;
matlabbatch{2}.spm.stats.fmri_spec.fact(2).name     = 'Rep';
matlabbatch{2}.spm.stats.fmri_spec.fact(2).levels   = 2;
matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];

% Model Estimation
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM_categorical','SPM.mat'));


spm_jobman('run',matlabbatch);
close all


