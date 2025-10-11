% Prepare auditory fMRI dataset from the SPM website: 
% http://www.fil.ion.ucl.ac.uk/spm/data/auditory/
%
% The dataset is described in the SPM manual:
% http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#Chap:data:auditory
%
% This experiment was conducted by Geraint Rees under the direction of 
% Karl Friston and the FIL methods group. The purpose was to explore
% equipment and techniques in the early days of our fMRI experience.
% As such it has not been formally written up, and is freely available
% for personal education and evaluation purposes.
%==========================================================================

clear

%==========================================================================
% Data preparation steps are adapted from the auditory_spm12_batch.m sript:
% 
% Directory containing the auditory data
%--------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end
data_path = fullfile(data_path,'auditory','sub-01');
if ~exist(data_path,'dir'); mkdir(data_path); end
fprintf('%-40s:', 'Downloading auditory dataset from the SPM website (32.6 MB) ...'); 
urlwrite('http://www.fil.ion.ucl.ac.uk/spm/download/data/MoAEpilot/MoAEpilot.zip',fullfile(data_path,'auditory.zip'));
unzip(fullfile(data_path,'auditory.zip'),data_path);
fprintf(' %20s\n', '...done');


% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREAMBLE: DUMMY SCANS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = spm_select('FPList', fullfile(data_path,'fM00223'), '^f.*\.img$') ;

clear matlabbatch

matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'dummy';

matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(f(1:12,:));
matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = cellstr(fullfile(data_path,'dummy'));

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SPATIAL PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = spm_select('FPList', fullfile(data_path,'fM00223'), '^f.*\.img$');
a = spm_select('FPList', fullfile(data_path,'sM00223'), '^s.*\.img$');

clear matlabbatch

% Realign
%--------------------------------------------------------------------------
matlabbatch{1}.spm.spatial.realign.estwrite.data = {cellstr(f)};
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];

% Coregister
%--------------------------------------------------------------------------
matlabbatch{2}.spm.spatial.coreg.estimate.ref    = cellstr(spm_file(f(1,:),'prefix','mean'));
matlabbatch{2}.spm.spatial.coreg.estimate.source = cellstr(a);

% Segment
%--------------------------------------------------------------------------
matlabbatch{3}.spm.spatial.preproc.channel.vols  = cellstr(a);
matlabbatch{3}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{3}.spm.spatial.preproc.warp.write    = [0 1];

% Normalise: Write
%--------------------------------------------------------------------------
matlabbatch{4}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample = cellstr(f);
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox  = [3 3 3];

matlabbatch{5}.spm.spatial.normalise.write.subj.def      = cellstr(spm_file(a,'prefix','y_','ext','nii'));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample = cellstr(spm_file(a,'prefix','m','ext','nii'));
matlabbatch{5}.spm.spatial.normalise.write.woptions.vox  = [1 1 3];

% Smooth
%--------------------------------------------------------------------------
matlabbatch{6}.spm.spatial.smooth.data = cellstr(spm_file(f,'prefix','w'));
matlabbatch{6}.spm.spatial.smooth.fwhm = [6 6 6];

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLM SPECIFICATION, ESTIMATION, INFERENCE, RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = spm_select('FPList', fullfile(data_path,'fM00223'), '^swf.*\.img$');

clear matlabbatch

% Output Directory
%--------------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM';

% Model Specification
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM'));
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 7;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans = cellstr(f);
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.name = 'active';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.onset = 6:12:84;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond.duration = 6;
matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = cellstr(spm_file(f(1,:),'prefix','rp_','ext','.txt')); 

% Model Estimation
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM','SPM.mat'));

spm_jobman('run',matlabbatch);
close all