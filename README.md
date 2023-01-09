# lehner19natcc
Collection of scripts to produce figures in Lehner et al., 2019, Nature Climate Change (https://doi.org/10.1038/s41558-019-0639-x)

Author: Flavio Lehner, 2019, flehner@ucar.edu or flavio.lehner@cornell.edu

Requirements: 
- Matlab and NCL. Custom-built Matlab functions are provided.

Main scripts:
- runoff_efficiency_paper_tas_pr_ro_re_maps_for_sharing.ncl: NCL script to create Fig. 1
- runoff_efficiency_paper_global_regression_analysis_map_for_sharing.ncl: NCL script to create Fig. S2
- runoff_efficiency_cmip5_no_recon_for_sharing.m: Matlab script to create Figs. 2-4 and Figs. S1, S3, and S4. Includes additional unpublished figures that support the analysis and conclusions of the main paper.

Data needs: 
- CMIP5: Monthly mean 'tas', 'pr', 'mrros' from CMIP5 (not provided here, but available in its raw form at https://esgf-node.llnl.gov/search/cmip5/)
- Observations: forthcoming
- Shapefiles: forthcoming
