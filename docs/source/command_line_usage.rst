.. _commandline:

Command Line Usage
==================

All functions in **TMFC_denoise** can be executed directly from the MATLAB command line or used in custom scripts without the GUI.

Example scripts demonstrating typical command-line usage are available in the 
`examples <https://github.com/IHB-IBR-department/TMFC_denoise/tree/main/examples>`_ folder of the TMFC_denoise GitHub repository:

- ``prepare_auditory_dataset.m`` — prepares the **auditory fMRI dataset** from the SPM website:  
  http://www.fil.ion.ucl.ac.uk/spm/data/auditory/
- ``prepare_facerep_dataset.m`` — prepares the **face repetition fMRI dataset** from the SPM website:  
  http://www.fil.ion.ucl.ac.uk/spm/data/face_rep/
- ``TMFC_denoise_Example_01_Auditory_dataset.m`` — example of TMFC_denoise usage (**block design**, auditory dataset).
- ``TMFC_denoise_Example_02_Facerep_dataset.m`` — example of TMFC_denoise usage (**event-related design**, face repetition dataset).

These scripts illustrate the full processing workflow, including dataset preparation, mask generation, denoising, and quality-control steps.  
They can be used as templates for adapting **TMFC_denoise** to other fMRI studies.


