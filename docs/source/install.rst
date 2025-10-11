Installation
============

You can install **TMFC_denoise** either:

1. As a **standalone toolbox**, or  
2. As part of the **TMFC toolbox**
   (https://github.com/IHB-IBR-department/TMFC_toolbox).

Dependencies
------------

Before using **TMFC_denoise**, make sure the following software is installed:

- **MATLAB** R2021b or newer
- **SPM12** or **SPM25**
- (Optional) **Parallel Computing Toolbox** – for parallel computations

Download
--------

You can obtain the toolbox from any of the following sources:

- **GitHub repository:**  
  https://github.com/IHB-IBR-department/TMFC_denoise
- **Zenodo DOI:**  
  https://doi.org/10.5281/zenodo.17176264
- **Included in the TMFC toolbox**:  
  https://github.com/IHB-IBR-department/TMFC_toolbox

Installation Steps
------------------

1. **Download and unzip** the toolbox archive from GitHub or Zenodo.
2. In MATLAB, open :menuselection:`Home --> Set Path --> Add with Subfolders`.
3. Select the ``TMFC_denoise`` folder  
   (or the ``TMFC_toolbox`` folder, if using the integrated version).
4. Click **Save** and **Close**.
5. Test the installation::

    TMFC_denoise

   This command should open the TMFC_denoise GUI.

Alternatively, if you installed the full TMFC toolbox::

    TMFC

Then, in the main TMFC GUI, choose :menuselection:`Tools --> Denoise`.

Command-Line Usage
------------------

The **TMFC_denoise** toolbox can also be executed without using the GUI::


    output_paths = TMFC_denoise(SPM_paths,subject_paths,options)
 
    output_paths = TMFC_denoise(SPM_paths,subject_paths,options,anat_paths,func_paths)

    output_paths = TMFC_denoise(SPM_paths,subject_paths,options,anat_paths,func_paths,display_FD,estimate_GLMs,clear_all)


For command-line examples, see the ``examples`` folder in the GitHub repository:
`examples <https://github.com/IHB-IBR-department/TMFC_denoise/tree/main/examples>`_ folder of the TMFC_denoise GitHub repository


