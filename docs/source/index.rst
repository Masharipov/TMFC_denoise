TMFC_denoise Documentation
=====================================

**TMFC_denoise** is a MATLAB toolbox for SPM12/SPM25 that performs GLM-based denoising (**noise regression**).
 
This toolbox allows you to **add noise regressors** to the original general linear model (GLM), calculate **framewise displacement (FD)**, 
Derivative of root mean square VARiance over voxelS **(DVARS)**, 
and **FD-DVARS correlation** before and after denoising.

The updated GLMs can be used for **task-based activation analysis** or for **task-modulated functional connectivity (TMFC) analysis**.

.. toctree::
   :maxdepth: 2
   :caption: Introduction

   install
   overview
   prepare_data

.. toctree::
   :maxdepth: 2
   :caption: Usage Guide

   select_subjects
   denoising_options
   select_anat
   select_func
   FD_plot
   spikereg
   mask_generation
   physioreg
   model_estimation
   DVARS

.. toctree::
   :caption: Command-line usage
   :maxdepth: 1

   command_line_usage

.. toctree::
   :caption: FAQ
   :maxdepth: 1

   FAQ

.. toctree::
   :caption: References
   :maxdepth: 1

   references
   
  
   


   
   



