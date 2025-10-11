Spike Regression
================

The function ``tmfc_spikereg`` calculates spike regressors. 
The number of spike regressors equals the number of flagged time points. 
It is called automatically by the main function TMFC_denoise if the user has selected the SpikeReg option, or it can be run manually::

    tmfc_spikereg(SPM_paths,options);

Outputs are saved in the ``TMFC_denoise`` subfolder within each subject’s first-level GLM directory, 
with filenames of the form ``SpikeReg_[FDthr_0.50mm].mat``, where ``FDthr`` denotes the selected FD threshold. 

**Note:** Adding spike regressors can introduce negative FD-DVARS correlations. 
Future work should be devoted to understanding this cautionary behavior.


