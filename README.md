# CFPS-Multivariate-Statistical-Data-Analysis
Multivariate Statistical Data Analysis of a Cell-Free Protein synthesis: This work provides necessary data, scripts and functions useful to provide important insights that can help support the operation and control of CFPS processes. This work was used to provide the results provided in the article "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control.", AIChE Journal. The prediction accuracy resuslts of the different algorithms change significantly from the published work due to the normalization of yield, aggregate% and monomer% done for confidentiality protection reasons. 
DOI: 10.5281/zenodo.4469875

Data-sets: D1-D5.xlsx


Scripts:

  -PCA_all.m: PCA analysis of all the data-sets (initial and final values)

  -PCA_D1D2.m: PCA analysis of D1 and D2

  -PCA_D5.m: PCA analysis of D5

  -PLS_D2D3D5.m: PLS regression and confidence limits calculation of D1 & D2, D5

  -MPLS_D3D4.m: Multiway-PLS regression and confidence limits calculation  of D3 & D4

  -Pred_accuracy.m: Compares the prediction accuracy of different algorithms (initial and final values)

  -Pred_accuracy_bytime.m: Compares the prediction accuracy of different algorithms using different time slots (Multiway).

  -Uneven_vectors_accuracy.m: Compares the prediction accuracy of different algorithms using uneven vector techniques.
  
  -Monitoring_T2_SPE.m: Shows MPLS monitoring techniques


Functions:

  -predm.m: Function to do leave-one-out predictions using PLS regression.

  -predpoly.m: Function to do leave-one-out predictions using Polynomial PLS regression.

  -predlinear.m: Function to do leave-one-out predictions using OLS regression.

  -clim.m: Function to do PLS regression and confidence limits using bootstrap calculation 

  -polypls.m: Function to do PLS regression with polynomial inner-relation.

  -polypred.m: Function to predict using polynomial PLS.

  -plsnip.m: Function to calculate PLS regressionfor only one LV.

  -plscv.m: K-fold Cross-validation for PLS regression

  -plsnipals.m: NIPALS algorithm for PLS regression


This work was supported by the UK Engineering & Physical Sciences Research Council (EPSRC) [EP/P006485/1] and a consortium of industrial users and sector organizations in the Future Targeted Healthcare Manufacturing Hub hosted by UCL Biochemical Engineering in collaboration with UK universities. CFPS is envisioned as a potential solution for the simplified, robust, flexible and local production of stratified or personalised biotherapeutic medicines, as might be necessary in a future pharmacy or hospital setting.
