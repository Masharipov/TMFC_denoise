FAQ
===

If your question is not covered here, please refer to the  
`GitHub Issues page <https://github.com/IHB-IBR-department/TMFC_denoise/issues>`_  
or contact the developer directly at: masharipov@ihb.spb.ru

---

**Q1. I only want to calculate motion metrics (FD and DVARS). Is that possible?**  
Yes. In *Denoising options*, set all denoising options to ``none`` and select ``6HMP`` for head motion.  
This will compute FD and FD-DVARS correlations (if DVARS calculation is enabled) without adding regressors to the model.

---

**Q2. Can TMFC_denoise be used for resting-state fMRI?**  
Yes. While optimized for task-based analyses, TMFC_denoise can also denoise resting-state data organized in SPM format.  
Simply skip task regressors in your GLM.

---

**Q3. What if I have data preprocessed with fMRIPrep, HCP, or FSL?**  
TMFC_denoise is compatible with GLMs created in SPM using externally preprocessed data.  
Adjust the order of motion regressors (``translation_idx`` and ``rotation_idx``) and rotation units (``deg`` or ``rad``) in the *Denoising options* window accordingly (see :ref:`options`).

---

**Q4. Does TMFC_denoise work only with the TMFC toolbox?**  
No. TMFC_denoise can also be used for task-activation analyses or other purposes independently of the TMFC toolbox (see :ref:`overview`).

---

**Q5. How can I resume a partially finished run?**  
TMFC_denoise saves intermediate outputs (e.g., masks, FD, aCompCor) in subject subfolders.  
If processing was interrupted, simply re-run the same command — previously completed subjects will be skipped automatically.

---

**Q6. Do I need to recompute masks or regressors each time I change options?**  
No. Segmentation and already created masks with the same parameters are not recomputed.  
Similarly, previously generated HMP expansions and tissue-based regressors for a given mask configuration are reused automatically.

---

**Q7. My GLM already contains nuisance regressors (except motion parameters).**  
That’s fine. Check the indices of motion regressors in ``SPM.Sess.C`` and specify them in *Denoising options*.  
All existing regressors from the original model will be preserved in the updated GLM.

---

**Q8. My model includes time/dispersion derivatives and parametric modulators.**  
No problem — these regressors will remain in the updated model after denoising.

---

**Q9. I don’t have unsmoothed preprocessed functional images. Can I still run TMFC_denoise?**  
Yes, although using unsmoothed images for tissue-based regressors is recommended.  
Smoothed data can still be used, but gray matter contamination may slightly reduce the specificity of WM/CSF regressors.  
This effect is minimized by erosion of WM and CSF masks.

---

**Q10. My sessions were concatenated using ``spm_fmri_concatenate``. Can I use TMFC_denoise?**  
Yes. TMFC_denoise fully supports concatenated multi-session models created with ``spm_fmri_concatenate``.

---

**Q11. My model was already estimated with rWLS. Should I select the rWLS option again?**  
You should use models estimated without rWLS (``none``, ``AR(1)``, or ``FAST``).  
If your model was previously estimated with rWLS but you do *not* select the rWLS option, the updated model will default to AR(1).  
To apply rWLS again, enable the rWLS option explicitly.

---

**Q12. My eroded WM/CSF mask contains only a few voxels. Is that acceptable?**  
Consider using a more liberal probability threshold and reducing the number of erosion cycles.  
You may also decrease the number of GM dilation cycles, since the GM mask is subtracted from the CSF mask.

---

**Q13. Can I use functional images in native space?**  
For tissue-based regressors and DVARS, no — they must be in MNI space.  
However, you can still perform HMP expansions, rWLS estimation, and FD calculation on native-space data.

---

**Q14. What are the best denoising options for gPPI and BSC analyses?**  
There is currently no universally optimal denoising strategy for gPPI or BSC.  
Large-scale benchmarking is still needed. Select options empirically for your dataset and inspect FD–DVARS correlations —  
successful denoising should reduce the group mean correlation toward zero.

---

**Q15. Can I use noise regressors from TMFC_denoise in other software?**  
Yes. All generated nuisance regressors are saved in the subject’s ``TMFC_denoise`` subfolder and can be used in other analysis frameworks.

---

