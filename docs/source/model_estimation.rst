Model estimation
================

The ``tmfc_estimate_updated_GLMs`` function estimates updated GLMs with nuisance regressors. 
It is called automatically by the main function TMFC_denoise if the user has selected the corresponding options, or it can be run manually::

    output_paths = tmfc_estimate_updated_GLMs(SPM_paths,masks,options);

The outputs are saved in the ``TMFC_denoise/[WM*e*]_[CSF*e*]_[GM*d*]/GLM_*`` subfolders, where folder names encode the selected denoising and masking parameters (see :ref:`masks`). 

 - ``WM*e*`` — probability threshold and number of erosion cycles.
 - ``CSF*e*`` — threshold and number of erosion cycles. 
 - ``GM*d*`` — threshold and number of dilation cycles.

Each selected denoising option appends a corresponding suffix to the updated GLM subfolder (see :ref:`options`). 
For example, ``GLM_[24HMP]_[aCompCor50]_[rWLS]`` indicates that the updated GLM includes 24 head-motion regressors, a variable number of aCompCor regressors explaining 50% of WM and CSF variance, and was estimated using rWLS.

The updated GLM subfolders contain the standard outputs from SPM model estimation, as well as ``GLM_batch.m`` files, which store ``matlabbatch`` structures that can be reopened in the SPM batch system. 

The ``SPM.mat`` files in these subfolders can be used as input to the ``TMFC toolbox``, which implements **gPPI and BSC-LSS methods** with or without FIR task regression (Masharipov et al., 2024). 
gPPI and LSS models automatically include nuisance regressors and, optionally, FIR regressors, along with high-pass filter regressors. 
Therefore, noise regression, FIR task co-activation regression (optional), and high-pass filtering are performed in a single step, which avoids reintroducing signal related to nuisance covariates (Lindquist et al., 2019).


