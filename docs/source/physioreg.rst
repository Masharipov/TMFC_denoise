Tissue-based nuisance regressors
================================

The ``tmfc_physioreg`` function calculates tissue-based nuisance regressors. 
It is called automatically by the main function TMFC_denoise if the user has selected the corresponding options, or it can be run manually::

    tmfc_physioreg(SPM_paths,subject_paths,func_paths,masks,options);

The outputs are saved in the ``TMFC_denoise/[WM*e*]_[CSF*e*]_[GM*d*]`` subfolders, where folder names encode the selected mask parameters (see :ref:`masks`): 

 - ``WM*e*`` — probability threshold and number of erosion cycles.
 - ``CSF*e*`` — threshold and number of erosion cycles. 
 - ``GM*d*`` — threshold and number of dilation cycles.

Generated files include (depending on user-selected options, see :ref:`options`):

 - ``2Phys.mat``, ``4Phys.mat``, ``8Phys.mat`` — WM/CSF signals.
 - ``GSR.mat``, ``2GSR.mat``, ``4GSR.mat`` — whole-brain signals.
 - ``[aCompCor_*WM_*CSF_Ort].mat`` — a fixed number of WM and CSF principal components (PCs). These files also include information on the variance explained by WM/CSF PCs per session and the mean variance explained across sessions. 
 - ``[aCompCor50_Ort].mat`` — a variable number of PCs explaining 50% of WM/CSF signals variance. This file also reports the mean and total number of PCs for WM and CSF across sessions.

**Note:** When pre-orthogonalization of WM and CSF signals is enabled before PC extraction,
the suffix ``_Ort`` is appended to all aCompCor output files.

