# GRACE-Data-Processing-Forward-Modeling
# GRACE Data Processing & Forward Modeling

This repository contains MATLAB scripts for processing GRACE Level-2 (Spherical Harmonic Coefficients) data and improving its spatial resolution using a forward modeling approach with WATERGAP Total Water Storage (TWS) data.  

## 1. Overview
GRACE (Gravity Recovery and Climate Experiment) mission provides time-variable gravity field data, which can be used to estimate surface mass changes. This project processes GRACE Level-2 data by:  
- Converting geoid/potential data into surface mass changes  
- Applying necessary corrections (Degree-1 correction, C20/C30 replacement, GIA correction)  
- Filtering to reduce noise and improve spatial localization
- Using WATERGAP TWS data to refine GRACE-derived mass changes via forward modeling  

## 2. Data Sources
### GRACE Level-2 Data  
- File: `GRACE_200204_201706.mat`  
- Source: [CSR Release-06 GRACE Level-2 Data](https://www2.csr.utexas.edu/grace/RL06.html)  

### WATERGAP TWS Data  
- File: `Watergap_TWS_sh_2003.mat`  
- Source: [PANGAEA Dataset](https://doi.pangaea.de/10.1594/PANGAEA.918447)  

## 3. Code Structure
### P01. GRACE Data Processing  
The first step processes GRACE Level-2 data and converts it into surface mass change.  

#### Key Steps:  
1. Convert geoid/potential to surface mass (water equivalent)  
   - Apply Degree-1 correction using TN-13
     (Source: https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_CSR_RL0603.txt)
   - Replace C20/C30 coefficients using SLR solutions  
     (Source: https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt)
   - Apply Glacial Isostatic Adjustment (GIA) correction 
     (Peltier et al., 2018)

2. Filtering
   - De-correlation filtering (PM filtering):
     Reduces correlated errors  
   - Slepian localization: Enhances spatial resolution  
   - Gaussian smoothing: Reduces high-frequency noise  

3. Save Processed Data
   - Output: `./output/P01_Slepian_GRACE_GSM_surfacemass_P4M8_DL60_cg_G400_mmH20.mat`  
### P02. Forward Modeling
This step refines the GRACE-derived mass changes using WATERGAP TWS data.  

#### Key Steps:
1. Load GRACE and WATERGAP Data 
   - Load the processed GRACE GSM data from Step 1  
   - Load WATERGAP TWS data(initial estimates)  

2. Forward Modeling Process
   - Compute the difference between GRACE and WATERGAP data
   - Convert the difference into spatial grid format  
   - Update WATERGAP TWS using the difference correction  
   - Iterate until convergence (RMS difference < 0.1)  

3. Save Forward Modeling Results
   - Output: `./output/P02_FM_CSR_GSM_RL06_sp60_cg_P4M8_G400_deg1_c20_c30_slr_PE_ICE6GD_watergap.mat`  

## 4. How to Run  
### Requirements 
- MATLAB  
- Data files placed in `./Data/` directory  
- Supporting Functions in current directory
  ⚠️ **Due to file size limitations, these functions are not included in this repository.**  
  :You can download them from: https://drive.google.com/file/d/1pUBvslOX3uTJpjSi5cAkmDR9UOp_G1NE/view?usp=drivesdk
  :Alternatively, please contact me if you need access to these functions.  

### Execution
1. Run `P01_GRACE_Processing.m` to process GRACE data  
2. Run `P02_Forward_Modeling.m` to perform forward modeling  

## 5. References 
- Loomis, B. D., et al. (2019). "Replacing GRACE/GRACE-FO C20 with SLR-derived estimates." *Geophysical Research Letters*. [DOI:10.1029/2019GL082929](https://doi.org/10.1029/2019GL082929)  
- Peltier, W. R., et al. (2018). "ICE6G-D Glacial Isostatic Adjustment Model."  
- CSR GRACE/GRACE-FO Release-06: [https://www2.csr.utexas.edu/grace/RL06.html](https://www2.csr.utexas.edu/grace/RL06.html)  

