close all
clear all

% ==============================================================================
%
% runoff_efficiency_cmip5_no_recon_for_sharing.m
% author: flehner@ucar.edu
%
% Scripts to produce Figs. 2-4 and Figs. S1, S3, and S4 in Lehner et al. 2019.
% Includes additional unpublished figures that support the analysis and
% conclusions of the main paper.
%
% Full citation:
% Lehner, F., A. W. Wood, J. A. Vano, D. M. Lawrence, M. P. Clark, J. S. Mankin (2019):
% The potential to reduce uncertainty in regional runoff projections from climate models
% Nature Climate Change, DOI: 10.1038/s41558-019-0639-x
% https://www.nature.com/articles/s41558-019-0639-x
%
% Required functions:
% -
%
% Notes:
% - This is a terribly long script with controls at the beginning that
% allow to plot individual figures. Ideally, this will be converted to a more
% user-friendly format one day, but it will likely/hopefully not be
% in Matlab anymore.
% - Postprocessing of GCMs is not included here and involves croping and averaging
% the relevant quantities for the watersheds of interest using HUC shapefiles.
%
% ==============================================================================

% -- models --
models_piControl = {'bcc-csm1-1' 'bcc-csm1-1-m' 'BNU-ESM' 'CanESM2' 'CCSM4' 'CESM1-BGC' 'CESM1-CAM5' 'CMCC-CESM' 'CMCC-CM' 'CMCC-CMS' 'CNRM-CM5' 'CNRM-CM5-2' 'CSIRO-Mk3-6-0' 'FGOALS-g2' 'FIO-ESM' 'GFDL-CM2p1' 'GFDL-CM3' 'GFDL-ESM2G' 'GFDL-ESM2M' 'GISS-E2-H' 'GISS-E2-H-CC' 'GISS-E2-R' 'GISS-E2-R-CC' 'HadCM3' 'HadGEM2-ES' 'inmcm4' 'IPSL-CM5A-LR' 'IPSL-CM5A-MR' 'IPSL-CM5B-LR' 'MIROC4h' 'MIROC5' 'MIROC-ESM' 'MIROC-ESM-CHEM' 'MPI-ESM-LR' 'MPI-ESM-MR' 'MPI-ESM-P' 'MRI-CGCM3' 'MRI-ESM1' 'NorESM1-M' 'NorESM1-ME'};
models_1pctCO2 = {'ACCESS1-0' 'ACCESS1-3' 'bcc-csm1-1' 'bcc-csm1-1-m' 'BNU-ESM' 'CanESM2' 'CCSM4' 'CESM1-BGC' 'CESM1-CAM5' 'CNRM-CM5' 'CNRM-CM5-2' 'CSIRO-Mk3-6-0' 'CSIRO-Mk3L-1-2' 'GFDL-CM3' 'GFDL-ESM2G' 'GFDL-ESM2M' 'GISS-E2-H' 'GISS-E2-R' 'inmcm4' 'IPSL-CM5A-LR' 'IPSL-CM5A-MR' 'IPSL-CM5B-LR' 'MIROC5' 'MIROC-ESM' 'MPI-ESM-LR' 'MPI-ESM-MR' 'MPI-ESM-P' 'MRI-CGCM3' 'NorESM1-M' 'NorESM1-ME'};
models_common = intersect(models_1pctCO2,models_piControl); % 'FGOALS-g2' excluded on purpose
% -- create a selection of models with balanced P-E=R in piControl (within 3%):
models_common = {'CCSM4' 'CESM1-BGC' 'CESM1-CAM5' 'CNRM-CM5' 'CanESM2' 'GFDL-ESM2G' 'GFDL-ESM2M' 'GISS-E2-H' 'IPSL-CM5A-MR' 'MIROC-ESM' 'MIROC5' 'MPI-ESM-LR' 'MPI-ESM-MR' 'MRI-CGCM3' 'NorESM1-M' 'NorESM1-ME' 'bcc-csm1-1' 'bcc-csm1-1-m' 'inmcm4'};
% -- carbon cycle models:
models_cc = {'CanESM2' 'CESM1-BGC' 'GFDL-ESM2M' 'IPSL-CM5A-LR' 'MPI-ESM-LR' 'NorESM1-ME'};
% -- models that have CLM (3 or 4) as their land component:
models_clm = {'CCSM4' 'CESM1-BGC' 'CESM1-CAM5' 'NorESM1-M' 'NorESM1-ME' 'bcc-csm1-1' 'bcc-csm1-1-m'};

scen        = {'piControl' '1pctCO2' 'historical'};
scen2       = {'piControl','1pctCO2','esmFixClim1','esmFdbk1'};
scen2_names = {'piControl','1pctCO2','CO_{2phys}','CO_{2rad}'};

% -- figure pathout --
figpath = '~/Dropbox/publication/lehner18_runoff_efficiency_constraint/fig/';

% -- which figure to plot --
% -- main paper figures (Fig. 2-4):
fig2      = 1; % precip vs runoff scatter plots
fig3      = 0; % runoff sensitivities
fig4      = 0; % less complicated future plot (one future temperature)

% -- supplementary figures (only Fig. S1, S3, and S4 were actually used in Lehner et al. 2019):
figS1     = 0; % methodological uncertainties in sensitivities (only set up for 'UC' basin)
figS3     = 0; % CO2 effects (Abby Swann's simulations)
figS4     = 0; % runoff vs streamflow, does it behave the same?
figS5     = 0; % moving window sensitivities
figS6     = 0; % comparing Reclamation and GRDC runoff (only for UC and CO_DALLES)
figS8     = 0; % check the agreement of different sensitivity calculations (2-term OLS vs 3-term OLS with storage coefficient)
fig5      = 0; % less complicated future plot (several future temperatures)

% -- which basin to plot --
vari = 'UC' % UC CO_DALLES NS
% UC = Upper Colorado
% CO_DALLES = Columbia at Dalles
% NS = Northern Sierras

% -- observations
obs_data = 3 ; % 1=PRISM/PRISM, 2=GPCC/BEST, 3=Livneh

% -- reference period
refstart  = 1950; % 1929 1950
refende   = 2008;

tf = 1.23348; % transfer factor from 10^6 acre feet to km^3

anoml = 0; % express everything as anomalies to mean? 1=yes, 0=no

reg_or_med = 1; % sensitivities as 1=regression or 2=epsilon (Sankar and Vogel)

% -- carry over term or not
co = 0; % 1=yes, 0=no -- previous-year carry over term (proxy for storage, see Milly et al 2018 WRR)
if co == 1
  co_tag = '_with_carryover';
else
  co_tag = '';
end

% -- regress out precip from temp time series? 0=no, 1=yes
regout0 = 0;
% -- regress out temp from precip time series? 0=no, 1=yes
regout1 = 0;
% -- regress out precip from runoff time series? 0=no, 1=yes
regout2 = 0;
% -- regress out temp from runoff time series? 0=no, 1=yes
regout3 = 0;

if strcmp(vari,'UC')==1
  vari_name = 'Upper Colorado';
elseif strcmp(vari,'CO_DALLES')==1
  vari_name = 'Columbia';
else
  vari_name = 'Northern Sierras';
end


% === CMIP5 data ===============================================================
% -- mrro   = total runoff (incl. subsurface drainage)
% -- mrros  = surface runoff
for s = 1:3
  pathin = ['~/Dropbox/work/cmip5_ncar/runoff_efficiency_ts_collection/cmip5_historical_cut/'];
  % -- for now limited to historical 1900-2005 (106 years), maybe develop script to account for longer simulations like piControl
  min_length_start = 1;
  min_length = 106;
  if s == 1 % ======> piControl
    period = '';
  elseif s == 2 % ======> 1pctCO2
    period = '';
  elseif s == 3 % ======> historical
    period = '_190001-200512';
  end
  for m = 1:length(models_common)
    % -- runoff
    tmp1 = ncread([pathin 'mrro_' scen{s} '_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
    ro(m,s,:) = tmp1((min_length_start-1)*12+1:min_length*12);
    for i = 1:(length(tmp1)/12)-1
      tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
    end
    % -- snippet to determine min length (not needed after it has been determined and hard-coded the first time)
    % if (length(tmp1)/12)-1 < min_length
    %   min_length = (length(tmp1)/12)-1;
    % end
    % ------------------------------------------------------
    ro_wy(m,s,:)      = tmp2(min_length_start:min_length-1);
    ro_wy_mean(m,s)   = nanmean(ro_wy(m,s,:));
    ro_wy_median(m,s) = nanmedian(ro_wy(m,s,:));
    ro_wy_percent_mean(m,s,:)   = (ro_wy(m,s,:)./ro_wy_mean(m,s))*100;
    ro_wy_percent_median(m,s,:) = (ro_wy(m,s,:)./ro_wy_median(m,s))*100;
    clear('tmp1','tmp2')
    % % -- runoff surface
    % tmp1 = ncread([pathin 'mrros_' scen{s} '_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
    % ros(m,s,:) = tmp1((min_length_start-1)*12+1:min_length*12);
    % for i = 1:(length(tmp1)/12)-1
    %   tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
    % end
    % ros_wy(m,s,:)      = tmp2(min_length_start:min_length-1);
    % ros_wy_mean(m,s)   = nanmean(ros_wy(m,s,:));
    % ros_wy_median(m,s) = nanmedian(ros_wy(m,s,:));
    % ros_wy_percent_mean(m,s,:)   = (ros_wy(m,s,:)./ros_wy_mean(m,s))*100;
    % ros_wy_percent_median(m,s,:) = (ros_wy(m,s,:)./ros_wy_median(m,s))*100;
    % clear('tmp1','tmp2')
    % -- precip
    tmp1 = ncread([pathin 'pr_' scen{s} '_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
    pr(m,s,:) = tmp1((min_length_start-1)*12+1:min_length*12);
    for i = 1:(length(tmp1)/12)-1
      tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
      tmp3(i) = sum(tmp1((i-1)*12+10:i*12+4));
    end
    pr_wy(m,s,:)          = tmp2(min_length_start:min_length-1);
    pr_oct_apr(m,s,:)     = tmp3(min_length_start:min_length-1);
    pr_wy_mean(m,s)       = nanmean(pr_wy(m,s,:));
    pr_oct_apr_mean(m,s)  = nanmean(pr_oct_apr(m,s,:));
    pr_wy_median(m,s)       = nanmedian(pr_wy(m,s,:));
    pr_oct_apr_median(m,s)  = nanmedian(pr_oct_apr(m,s,:));
    pr_wy_percent_mean(m,s,:) = (pr_wy(m,s,:)./pr_wy_mean(m,s))*100;
    pr_wy_percent_median(m,s,:) = (pr_wy(m,s,:)./pr_wy_median(m,s))*100;
    clear('tmp1','tmp2','tmp3')
    % -- evapotransp (et)
    % tmp1 = ncread([pathin 'evspsbl_' scen{s} '_' models_common{m} '_remapcon2_1x1_ts.nc'],[vari '_TS']);
    tmp1 = ncread([pathin 'hfls_' scen{s} '_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
    et(m,s,:) = tmp1((min_length_start-1)*12+1:min_length*12);
    for i = 1:(length(tmp1)/12)-1
      tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
      tmp3(i) = sum(tmp1((i-1)*12+10:i*12+4));
    end
    et_wy(m,s,:)          = tmp2(min_length_start:min_length-1);
    et_oct_apr(m,s,:)     = tmp3(min_length_start:min_length-1);
    et_wy_mean(m,s)       = nanmean(et_wy(m,s,:));
    et_oct_apr_mean(m,s)  = nanmean(et_oct_apr(m,s,:));
    et_wy_median(m,s)       = nanmedian(et_wy(m,s,:));
    et_oct_apr_median(m,s)  = nanmedian(et_oct_apr(m,s,:));
    clear('tmp1','tmp2','tmp3')
    % -- temperature
    tmp1 = ncread([pathin 'tas_' scen{s} '_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS'])-273.16;
    tas(m,s,:) = tmp1((min_length_start-1)*12+1:min_length*12);
    for i = 1:(length(tmp1)/12)-1
      tmp2(i) = nanmean(tmp1((i-1)*12+10:i*12+9));
      tmp3(i) = nanmean(tmp1(i*12+5:i*12+6));
    end
    tas_wy(m,s,:)         = tmp2(min_length_start:min_length-1);
    tas_may_jun(m,s,:)    = tmp3(min_length_start:min_length-1);
    tas_wy_mean(m,s)      = nanmean(tas_wy(m,s,:));
    tas_wy_anom(m,s,:)      = tas_wy(m,s,:)-tas_wy_mean(m,s);
    tas_may_jun_mean(m,s) = nanmean(tas_may_jun(m,s,:));
    tas_wy_median(m,s)      = nanmedian(tas_wy(m,s,:));
    tas_may_jun_median(m,s) = nanmedian(tas_may_jun(m,s,:));
    clear('tmp1','tmp2','tmp3')
  end
end

if regout0 == 1
  % -- regress out precip from temp
  for s = 1:2
    for m = 1:length(models_common)
      b = regress(squeeze(tas_wy(m,s,:)),[ones(size(squeeze(pr_wy(m,s,:)))) squeeze(pr_wy(m,s,:))]);
      tas_wy_pr_regressed_out(m,s,:) = squeeze(tas_wy(m,s,:)) - (b(1)+b(2)*squeeze(pr_wy(m,s,:)));
    end
  end
else
  tas_wy_pr_regressed_out = tas_wy;
end
if regout1 == 1
  % -- regress out temp from precip
  for s = 1:2
    for m = 1:length(models_common)
      % b = regress(squeeze(pr_wy_percent_mean(m,s,:)),[ones(size(squeeze(tas_wy(m,s,:)))) squeeze(tas_wy(m,s,:))]);
      % pr_wy_percent_mean_resid(m,s,:) = squeeze(pr_wy_percent_mean(m,s,:)) - (b(1)+b(2)*squeeze(tas_wy(m,s,:)));
      b = regress(squeeze(pr_wy(m,s,:)),[ones(size(squeeze(tas_wy(m,s,:)))) squeeze(tas_wy(m,s,:))]);
      pr_wy_tas_regressed_out(m,s,:) = squeeze(pr_wy(m,s,:)) - (b(1)+b(2)*squeeze(tas_wy(m,s,:))) + nanmean(squeeze(pr_wy(m,s,:)));
      pr_wy_tas_regressed_out_percent_mean(m,s,:) = (pr_wy_tas_regressed_out(m,s,:)./mean(pr_wy_tas_regressed_out(m,s,:)))*100;
    end
  end
else
  pr_wy_tas_regressed_out = pr_wy;
  pr_wy_tas_regressed_out_percent_mean = pr_wy_percent_mean;
end
if regout2 == 1
  % -- regress out precip from runoff
  for s = 1:2
    for m = 1:length(models_common)
      % b = regress(squeeze(ro_wy_percent_mean(m,s,:)),[ones(size(squeeze(pr_wy_percent_mean(m,s,:)))) squeeze(pr_wy_percent_mean(m,s,:))]);
      % ro_wy_percent_mean_resid_pr(m,s,:) = squeeze(ro_wy_percent_mean(m,s,:)) - (b(1)+b(2)*squeeze(pr_wy_percent_mean(m,s,:))) + nanmean(squeeze(ro_wy_percent_mean(m,s,:)));
      b = regress(squeeze(ro_wy(m,s,:)),[ones(size(squeeze(pr_wy(m,s,:)))) squeeze(pr_wy(m,s,:))]);
      ro_wy_pr_regressed_out(m,s,:) = squeeze(ro_wy(m,s,:)) - (b(1)+b(2)*squeeze(pr_wy(m,s,:))) + nanmean(squeeze(ro_wy(m,s,:)));
      ro_wy_pr_regressed_out_percent_mean(m,s,:) = (ro_wy_pr_regressed_out(m,s,:)./mean(ro_wy_pr_regressed_out(m,s,:)))*100;
    end
  end
else
  ro_wy_pr_regressed_out = ro_wy;
  ro_wy_pr_regressed_out_percent_mean = ro_wy_percent_mean;
end
if regout3 == 1
  % -- regress out temp from runoff
  for s = 1:2
    for m = 1:length(models_common)
      b = regress(squeeze(ro_wy(m,s,:)),[ones(size(squeeze(tas_wy(m,s,:)))) squeeze(tas_wy(m,s,:))]);
      ro_wy_tas_regressed_out(m,s,:) = squeeze(ro_wy(m,s,:)) - (b(1)+b(2)*squeeze(tas_wy(m,s,:))) + nanmean(squeeze(ro_wy(m,s,:)));
      ro_wy_tas_regressed_out_percent_mean(m,s,:) = (ro_wy_tas_regressed_out(m,s,:)./mean(ro_wy_tas_regressed_out(m,s,:)))*100;
    end
  end
else
  ro_wy_tas_regressed_out = ro_wy;
  ro_wy_tas_regressed_out_percent_mean = ro_wy_percent_mean;
end
% -- anomalies
for s = 1:2
  for m = 1:length(models_common)
    tas_wy_pr_regressed_out_anom(m,s,:)      = tas_wy_pr_regressed_out(m,s,:)-nanmean(tas_wy_pr_regressed_out(m,s,:));
  end
end

% -- calculate P-E
pme    = pr-et;
pme_wy = pr_wy-et_wy;

% -- calculate runoff efficiency
re_wy = (ro_wy./pr_wy)*100;
re_wy_median = nanmedian(re_wy,3);



% === CMIP5 data historical2rcp85 ===============================================================
pathin    = ['~/Dropbox/work/cmip5_ncar/runoff_efficiency_ts_collection/'];
period    = '_190001-200512';
start     = 1901;
ende      = 2099;
time_h2r  = start:ende;
for m = 1:length(models_common)
  % -- runoff
  tmp0 = ncread([pathin '/cmip5_historical_cut/mrro_historical_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
  tmp1 = ncread([pathin '/cmip5_rcp85_cut/mrro_rcp85_' models_common{m} '_remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
  xx0 = size(tmp0);
  xx1 = size(tmp1);
  if xx0(1)>xx0(2)
    tmp0 = tmp0';
  end
  if xx1(1)>xx1(2)
    tmp1 = tmp1';
  end
  tmpc = [tmp0 tmp1(1:1128)];
  ro_h2r(m,:) = tmpc;
  for i = 1:(length(tmpc)/12)-1
    tmp2(i) = sum(tmpc((i-1)*12+10:i*12+9));
  end
  ro_h2r_wy(m,:)        = tmp2;
  ro_h2r_wy_mean(m)   = nanmean(ro_h2r_wy(m,refstart-start+1:refende-start+1));
  ro_h2r_wy_median(m) = nanmedian(ro_h2r_wy(m,refstart-start+1:refende-start+1));
  ro_h2r_wy_percent_mean(m,:)   = (ro_h2r_wy(m,:)./ro_h2r_wy_mean(m))*100;
  ro_h2r_wy_percent_median(m,:) = (ro_h2r_wy(m,:)./ro_h2r_wy_median(m))*100;
  clear('tmp0','tmp1','tmp2')
  % -- runoff surface
  tmp0 = ncread([pathin '/cmip5_historical_cut/mrros_historical_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
  tmp1 = ncread([pathin '/cmip5_rcp85_cut/mrros_rcp85_' models_common{m} '_remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
  xx0 = size(tmp0);
  xx1 = size(tmp1);
  if xx0(1)>xx0(2)
    tmp0 = tmp0';
  end
  if xx1(1)>xx1(2)
    tmp1 = tmp1';
  end
  tmpc = [tmp0 tmp1(1:1128)];
  ros_h2r(m,:) = tmpc;
  for i = 1:(length(tmpc)/12)-1
    tmp2(i) = sum(tmpc((i-1)*12+10:i*12+9));
  end
  ros_h2r_wy(m,:)        = tmp2;
  ros_h2r_wy_mean(m)   = nanmean(ros_h2r_wy(m,refstart-start+1:refende-start+1));
  ros_h2r_wy_median(m) = nanmedian(ros_h2r_wy(m,refstart-start+1:refende-start+1));
  ros_h2r_wy_percent_mean(m,:)   = (ros_h2r_wy(m,:)./ros_h2r_wy_mean(m))*100;
  ros_h2r_wy_percent_median(m,:) = (ros_h2r_wy(m,:)./ros_h2r_wy_median(m))*100;
  clear('tmp0','tmp1','tmp2')
  % -- precip
  tmp0 = ncread([pathin '/cmip5_historical_cut/pr_historical_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
  tmp1 = ncread([pathin '/cmip5_rcp85_cut/pr_rcp85_' models_common{m} '_remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
  xx0 = size(tmp0);
  xx1 = size(tmp1);
  if xx0(1)>xx0(2)
    tmp0 = tmp0';
  end
  if xx1(1)>xx1(2)
    tmp1 = tmp1';
  end
  tmpc = [tmp0 tmp1(1:1128)];
  pr_h2r(m,:) = tmpc;
  for i = 1:(length(tmpc)/12)-1
    tmp2(i) = sum(tmpc((i-1)*12+10:i*12+9));
  end
  pr_h2r_wy(m,:)        = tmp2;
  pr_h2r_wy_mean(m)   = nanmean(pr_h2r_wy(m,refstart-start+1:refende-start+1));
  pr_h2r_wy_median(m) = nanmedian(pr_h2r_wy(m,refstart-start+1:refende-start+1));
  pr_h2r_wy_percent_mean(m,:)   = (pr_h2r_wy(m,:)./pr_h2r_wy_mean(m))*100;
  pr_h2r_wy_percent_median(m,:) = (pr_h2r_wy(m,:)./pr_h2r_wy_median(m))*100;
  clear('tmp0','tmp1','tmp2')
  % -- evapotransp (et)
  tmp0 = ncread([pathin '/cmip5_historical_cut/hfls_historical_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
  tmp1 = ncread([pathin '/cmip5_rcp85_cut/hfls_rcp85_' models_common{m} '_remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
  xx0 = size(tmp0);
  xx1 = size(tmp1);
  if xx0(1)>xx0(2)
    tmp0 = tmp0';
  end
  if xx1(1)>xx1(2)
    tmp1 = tmp1';
  end
  tmpc = [tmp0 tmp1(1:1128)];
  et_h2r(m,:) = tmpc;
  for i = 1:(length(tmpc)/12)-1
    tmp2(i) = sum(tmpc((i-1)*12+10:i*12+9));
  end
  et_h2r_wy(m,:)        = tmp2;
  et_h2r_wy_mean(m)   = nanmean(et_h2r_wy(m,refstart-start+1:refende-start+1));
  et_h2r_wy_median(m) = nanmedian(et_h2r_wy(m,refstart-start+1:refende-start+1));
  et_h2r_wy_percent_mean(m,:)   = (et_h2r_wy(m,:)./et_h2r_wy_mean(m))*100;
  et_h2r_wy_percent_median(m,:) = (et_h2r_wy(m,:)./et_h2r_wy_median(m))*100;
  clear('tmp0','tmp1','tmp2')
  % -- temperature
  tmp0 = ncread([pathin '/cmip5_historical_cut/tas_historical_' models_common{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS'])-273.16;
  tmp1 = ncread([pathin '/cmip5_rcp85_cut/tas_rcp85_' models_common{m} '_remapcon2_1x1_ts.nc'],[vari '_TS'])-273.16;
  xx0 = size(tmp0);
  xx1 = size(tmp1);
  if xx0(1)>xx0(2)
    tmp0 = tmp0';
  end
  if xx1(1)>xx1(2)
    tmp1 = tmp1';
  end
  tmpc = [tmp0 tmp1(1:1128)];
  tas_h2r(m,:) = tmpc;
  for i = 1:(length(tmpc)/12)-1
    tmp2(i) = nanmean(tmpc((i-1)*12+10:i*12+9));
  end
  tas_h2r_wy(m,:)        = tmp2;
  tas_h2r_wy_mean(m)     = nanmean(tas_h2r_wy(m,refstart-start+1:refende-start+1));
  tas_h2r_wy_median(m)   = nanmedian(tas_h2r_wy(m,refstart-start+1:refende-start+1));
  tas_h2r_wy_anom(m,:)  = tas_h2r_wy(m,:)-tas_h2r_wy_mean(m);
  clear('tmp0','tmp1','tmp2')
end

% -- calculate P-E
pme_h2r    = pr_h2r-et_h2r;
pme_h2r_wy = pr_h2r_wy-et_h2r_wy;

% -- calculate runoff efficiency
re_h2r_wy = (ro_h2r_wy./pr_h2r_wy)*100;
re_h2r_wy_median = nanmedian(re_h2r_wy,3);


if regout0 == 1
  % -- regress out precip from temp
  for m = 1:length(models_common)
    b = regress(tas_h2r_wy(m,:)',[ones(size(pr_h2r_wy(m,:)))' pr_h2r_wy(m,:)']);
    tas_h2r_wy_pr_regressed_out(m,:) = squeeze(tas_h2r_wy(m,:)) - (b(1)+b(2)*squeeze(pr_h2r_wy(m,:)));
  end
else
  tas_h2r_wy_pr_regressed_out = tas_h2r_wy;
end
if regout1 == 1
  % -- regress out temp from precip
  for m = 1:length(models_common)
    b = regress(pr_h2r_wy(m,:)',[ones(size(tas_h2r_wy(m,:)))' tas_h2r_wy(m,:)']);
    pr_h2r_wy_tas_regressed_out(m,:) = squeeze(pr_h2r_wy(m,:)) - (b(1)+b(2)*squeeze(tas_h2r_wy(m,:)));
    pr_h2r_wy_tas_regressed_out_percent_mean(m,:) = (pr_h2r_wy_tas_regressed_out(m,:)./mean(pr_h2r_wy_tas_regressed_out(m,:)))*100;
  end
else
  pr_h2r_wy_tas_regressed_out = pr_h2r_wy;
  pr_h2r_wy_tas_regressed_out_percent_mean = pr_h2r_wy_percent_mean;
end
if regout2 == 1
  % -- regress out precip from runoff
  for m = 1:length(models_common)
    b = regress(ro_h2r_wy(m,:)',[ones(size(pr_h2r_wy(m,:)))' pr_h2r_wy(m,:)']);
    ro_h2r_wy_pr_regressed_out(m,:) = squeeze(ro_h2r_wy(m,:)) - (b(1)+b(2)*squeeze(pr_h2r_wy(m,:))) + nanmean(ro_h2r_wy(m,:));
    ro_h2r_wy_pr_regressed_out_percent_mean(m,:) = (ro_h2r_wy_pr_regressed_out(m,:)./mean(ro_h2r_wy_pr_regressed_out(m,:)))*100;
  end
else
  ro_h2r_wy_pr_regressed_out = ro_h2r_wy;
  ro_h2r_wy_pr_regressed_out_percent_mean = ro_h2r_wy_percent_mean;
end
if regout3 == 1
  % -- regress out temp from runoff
  for m = 1:length(models_common)
    b = regress(ro_h2r_wy(m,:)',[ones(size(tas_h2r_wy(m,:)))' tas_h2r_wy(m,:)']);
    ro_h2r_wy_tas_regressed_out(m,:) = ro_h2r_wy(m,:) - (b(1)+b(2)*squeeze(tas_h2r_wy(m,:))) + nanmean(ro_h2r_wy(m,:));
    ro_h2r_wy_tas_regressed_out_percent_mean(m,:) = (ro_h2r_wy_tas_regressed_out(m,:)./mean(ro_h2r_wy_tas_regressed_out(m,:)))*100;
  end
else
  ro_h2r_wy_tas_regressed_out = ro_h2r_wy;
  ro_h2r_wy_tas_regressed_out_percent_mean = ro_h2r_wy_percent_mean;
end
% -- anomalies
for m = 1:length(models_common)
  tas_h2r_wy_pr_regressed_out_anom(m,:) = tas_h2r_wy_pr_regressed_out(m,:)-nanmean(tas_h2r_wy_pr_regressed_out(m,refstart-start+1:refende-start+1));
end

corr(ro_h2r_wy_tas_regressed_out(1,:)',pr_h2r_wy_percent_mean(1,:)')
nn=80;
for m = 1:19
  corr(ro_h2r_wy_tas_regressed_out(m,1:nn)',pr_h2r_wy_percent_mean(m,1:nn)')
end


% === CMIP5 ESM data ===============================================================
% -- minimum length of simulation in years:
min_length_cc = 140;
period = '';
refstart_cc = 1;
refende_cc  = 50;
for s = 2:4 % -- the different experimental setups (2=1pctCO2, 3=esmFixClim1, 4=esmFdbk1)
  pathin = ['~/Dropbox/work/cmip5_ncar/runoff_efficiency_ts_collection/cmip5_esm/'];
  for m = 1:length(models_cc)
    % -- runoff
    tmp1 = ncread([pathin 'mrro_' models_cc{m} '_' scen2{s} '_remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
    ro_cc(m,s,:) = tmp1((min_length_start-1)*12+1:min_length_cc*12);
    for i = 1:(length(tmp1)/12)-1
      tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
    end
    ro_cc_wy(m,s,:)      = tmp2(min_length_start:min_length_cc-1);
    ro_cc_wy_mean(m,s)   = nanmean(ro_cc_wy(m,s,refstart_cc:refende_cc));
    ro_cc_wy_median(m,s) = nanmedian(ro_cc_wy(m,s,refstart_cc:refende_cc));
    ro_cc_wy_percent_median(m,s,:) = (ro_cc_wy(m,s,:)./ro_cc_wy_median(m,s))*100;
    ro_cc_wy_percent_mean(m,s,:) = (ro_cc_wy(m,s,:)./ro_cc_wy_mean(m,s))*100;
    clear('tmp1','tmp2')
    % -- precip
    tmp1 = ncread([pathin 'pr_' models_cc{m} '_' scen2{s} '_remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
    pr_cc(m,s,:) = tmp1((min_length_start-1)*12+1:min_length_cc*12);
    for i = 1:(length(tmp1)/12)-1
      tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
      tmp3(i) = sum(tmp1((i-1)*12+10:i*12+4));
    end
    pr_cc_wy(m,s,:)          = tmp2(min_length_start:min_length_cc-1);
    pr_cc_wy_mean(m,s)       = nanmean(pr_cc_wy(m,s,refstart_cc:refende_cc));
    pr_cc_wy_median(m,s)       = nanmedian(pr_cc_wy(m,s,refstart_cc:refende_cc));
    pr_cc_wy_percent_median(m,s,:) = (pr_cc_wy(m,s,:)./pr_cc_wy_median(m,s))*100;
    pr_cc_wy_percent_mean(m,s,:) = (pr_cc_wy(m,s,:)./pr_cc_wy_mean(m,s))*100;
    rt_cc_wy(m,s,:)          = (ro_cc_wy(m,s,:)-ro_cc_wy_mean(m,s))-(pr_cc_wy(m,s,:)-pr_cc_wy_mean(m,s));
    clear('tmp1','tmp2','tmp3')
    % -- evapotransp (et)
    tmp1 = ncread([pathin 'hfls_' models_cc{m} '_' scen2{s} '_remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
    et_cc(m,s,:) = tmp1((min_length_start-1)*12+1:min_length_cc*12);
    for i = 1:(length(tmp1)/12)-1
      tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
      tmp3(i) = sum(tmp1((i-1)*12+10:i*12+4));
    end
    et_cc_wy(m,s,:)          = tmp2(min_length_start:min_length_cc-1);
    et_cc_wy_mean(m,s)       = nanmean(et_cc_wy(m,s,refstart_cc:refende_cc));
    et_cc_wy_median(m,s)       = nanmedian(et_cc_wy(m,s,refstart_cc:refende_cc));
    et_cc_wy_percent_median(m,s,:) = (et_cc_wy(m,s,:)./et_cc_wy_median(m,s))*100;
    et_cc_wy_percent_mean(m,s,:) = (et_cc_wy(m,s,:)./et_cc_wy_mean(m,s))*100;
    clear('tmp1','tmp2','tmp3')
    % -- temperature
    tmp1 = ncread([pathin 'tas_' models_cc{m} '_' scen2{s} '_remapcon2_1x1_ts.nc'],[vari '_TS'])-273.16;
    tas_cc(m,s,:) = tmp1((min_length_start-1)*12+1:min_length_cc*12);
    for i = 1:(length(tmp1)/12)-1
      tmp2(i) = nanmean(tmp1((i-1)*12+10:i*12+9));
      tmp3(i) = nanmean(tmp1(i*12+5:i*12+6));
    end
    tas_cc_wy(m,s,:)         = tmp2(min_length_start:min_length_cc-1);
    tas_cc_wy_mean(m,s)      = nanmean(tas_cc_wy(m,s,refstart_cc:refende_cc));
    tas_cc_wy_median(m,s)    = nanmedian(tas_cc_wy(m,s,refstart_cc:refende_cc));
    tas_cc_wy_anom(m,s,:)    = tas_cc_wy(m,s,:)-tas_cc_wy_mean(m,s);
    clear('tmp1','tmp2','tmp3')
  end
end
% -- fill in respective piControl into same matrix --
s = 1;
pathin = ['~/Dropbox/work/cmip5_ncar/runoff_efficiency_ts_collection/cmip5_historical_cut/'];
% -- for now limited to historical 1900-2005 (106 years), maybe develop script to account for longer simulations like piControl
min_length_start = 1;
min_length = min_length_cc; %100;
for m = 1:length(models_cc)
  % -- runoff
  tmp1 = ncread([pathin 'mrro_' scen2{s} '_' models_cc{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
  ro_cc(m,s,(min_length_start-1)*12+1:min_length*12) = tmp1((min_length_start-1)*12+1:min_length*12);
  for i = 1:(length(tmp1)/12)-1
    tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
  end
  ro_cc_wy(m,s,min_length_start:min_length-1)      = tmp2(min_length_start:min_length-1);
  ro_cc_wy(m,s,ro_cc_wy(m,s,:)==0) = NaN;
  ro_cc_wy_mean(m,s)   = nanmean(ro_cc_wy(m,s,:));
  ro_cc_wy_median(m,s) = nanmedian(ro_cc_wy(m,s,:));
  ro_cc_wy_percent_mean(m,s,:)   = (ro_cc_wy(m,s,:)./ro_cc_wy_mean(m,s))*100;
  ro_cc_wy_percent_median(m,s,:) = (ro_cc_wy(m,s,:)./ro_cc_wy_median(m,s))*100;
  clear('tmp1','tmp2')
  % -- precip
  tmp1 = ncread([pathin 'pr_' scen{s} '_' models_cc{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
  pr_cc(m,s,(min_length_start-1)*12+1:min_length*12) = tmp1((min_length_start-1)*12+1:min_length*12);
  for i = 1:(length(tmp1)/12)-1
    tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
    tmp3(i) = sum(tmp1((i-1)*12+10:i*12+4));
  end
  pr_cc_wy(m,s,min_length_start:min_length-1)          = tmp2(min_length_start:min_length-1);
  pr_cc_wy(m,s,pr_cc_wy(m,s,:)==0) = NaN;
  pr_cc_wy_mean(m,s)       = nanmean(pr_cc_wy(m,s,:));;
  pr_cc_wy_median(m,s)       = nanmedian(pr_cc_wy(m,s,:));
  pr_cc_wy_percent_mean(m,s,:) = (pr_cc_wy(m,s,:)./pr_cc_wy_mean(m,s))*100;
  pr_cc_wy_percent_median(m,s,:) = (pr_cc_wy(m,s,:)./pr_cc_wy_median(m,s))*100;
  clear('tmp1','tmp2','tmp3')
  % -- evapotransp (et)
  tmp1 = ncread([pathin 'hfls_' scen{s} '_' models_cc{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS']) * tf;
  et_cc(m,s,(min_length_start-1)*12+1:min_length*12) = tmp1((min_length_start-1)*12+1:min_length*12);
  for i = 1:(length(tmp1)/12)-1
    tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
    tmp3(i) = sum(tmp1((i-1)*12+10:i*12+4));
  end
  et_cc_wy(m,s,min_length_start:min_length-1)          = tmp2(min_length_start:min_length-1);
  et_cc_wy(m,s,et_cc_wy(m,s,:)==0) = NaN;
  et_cc_wy_mean(m,s)       = nanmean(et_cc_wy(m,s,:));
  et_cc_wy_median(m,s)       = nanmedian(et_cc_wy(m,s,:));
  et_cc_wy_percent_mean(m,s,:) = (et_cc_wy(m,s,:)./et_cc_wy_mean(m,s))*100;
  et_cc_wy_percent_median(m,s,:) = (et_cc_wy(m,s,:)./et_cc_wy_median(m,s))*100;
  clear('tmp1','tmp2','tmp3')
  % -- temperature
  tmp1 = ncread([pathin 'tas_' scen{s} '_' models_cc{m} '_remapcon2_1x1_ts' period '.nc'],[vari '_TS'])-273.16;
  tas_cc(m,s,(min_length_start-1)*12+1:min_length*12) = tmp1((min_length_start-1)*12+1:min_length*12);
  for i = 1:(length(tmp1)/12)-1
    tmp2(i) = nanmean(tmp1((i-1)*12+10:i*12+9));
    tmp3(i) = nanmean(tmp1(i*12+5:i*12+6));
  end
  tas_cc_wy(m,s,min_length_start:min_length-1)         = tmp2(min_length_start:min_length-1);
  tas_cc_wy(m,s,tas_cc_wy(m,s,:)==0) = NaN;
  tas_cc_wy_mean(m,s)      = nanmean(tas_cc_wy(m,s,:));
  tas_cc_wy_median(m,s)      = nanmedian(tas_cc_wy(m,s,:));
  tas_cc_wy_anom(m,s,:)      = tas_cc_wy(m,s,:)-tas_cc_wy_mean(m,s);
  clear('tmp1','tmp2','tmp3')
end

if regout0 == 1
  % -- regress out precip from temp
  for s = 1:4
    for m = 1:length(models_cc)
      b = regress(squeeze(tas_cc_wy(m,s,:)),[ones(size(squeeze(pr_cc_wy(m,s,:)))) squeeze(pr_cc_wy(m,s,:))]);
      tas_cc_wy_pr_regressed_out(m,s,:) = squeeze(tas_cc_wy(m,s,:)) - (b(1)+b(2)*squeeze(pr_cc_wy(m,s,:)));
    end
  end
else
  tas_cc_wy_pr_regressed_out = tas_cc_wy;
end
if regout1 == 1
  % -- regress out temp from precip
  for s = 1:4
    for m = 1:length(models_cc)
      b = regress(squeeze(pr_cc_wy(m,s,:)),[ones(size(squeeze(tas_cc_wy(m,s,:)))) squeeze(tas_cc_wy(m,s,:))]);
      pr_cc_wy_tas_regressed_out(m,s,:) = squeeze(pr_cc_wy(m,s,:)) - (b(1)+b(2)*squeeze(tas_cc_wy(m,s,:))) + nanmean(squeeze(pr_cc_wy(m,s,:)));
      pr_cc_wy_tas_regressed_out_percent_mean(m,s,:) = (pr_cc_wy_tas_regressed_out(m,s,:)./mean(pr_cc_wy_tas_regressed_out(m,s,:)))*100;
    end
  end
else
  pr_cc_wy_tas_regressed_out = pr_cc_wy;
  pr_cc_wy_tas_regressed_out_percent_mean = pr_cc_wy_percent_mean;
end
% -- calculate P-E
pme_cc    = pr_cc-et_cc;
pme_cc_wy = pr_cc_wy-et_cc_wy;
% -- calculate runoff efficiency
re_cc_wy = (ro_cc_wy./pr_cc_wy)*100;
re_cc_wy_median = nanmedian(re_cc_wy,3);
re_cc_wy_mean = nanmean(re_cc_wy,3);
% -- anomalies
for s = 1:4
  for m = 1:length(models_cc)
    tas_cc_wy_pr_regressed_out(m,s,:)   = tas_cc_wy_pr_regressed_out(m,s,:)-nanmean(tas_cc_wy_pr_regressed_out(m,s,:));
    tas_cc_wy_anom(m,s,:)         = tas_cc_wy(m,s,:)-nanmean(tas_cc_wy_mean(m,s));
    pme_cc_wy_percent_mean(m,s,:) = ((pme_cc_wy(m,s,:)-nanmean(pme_cc_wy(m,s,refstart_cc:refende_cc)))/nanmean(pme_cc_wy(m,s,refstart_cc:refende_cc)))*100;
    re_cc_wy_anom(m,s,:)          = re_cc_wy(m,s,:)-nanmean(re_cc_wy_mean(m,s));
  end
end

% ------------------------------------------------------------------------------




% -- find 50-year average futures that are 2, 3, 4 degC warmer than piControl:
wl = 50;
for m = 1:length(models_common)
  tas_wy_50rm(m,:) = rm(squeeze(tas_wy(m,2,:)),wl)-nanmean(squeeze(tas_wy(m,1,:)));
  if min(abs(1-tas_wy_50rm(m,:))) > 0.2
    idx_1C(m,:) = NaN;
  else
    idx_1C(m,:) = find(abs(1-tas_wy_50rm(m,:)) == min(abs(1-tas_wy_50rm(m,:))));
  end
  if min(abs(2-tas_wy_50rm(m,:))) > 0.2
    idx_2C(m,:) = NaN;
  else
    idx_2C(m,:) = find(abs(2-tas_wy_50rm(m,:)) == min(abs(2-tas_wy_50rm(m,:))));
  end
  if min(abs(3-tas_wy_50rm(m,:))) > 0.2
    idx_3C(m,:) = NaN;
  else
    idx_3C(m,:) = find(abs(3-tas_wy_50rm(m,:)) == min(abs(3-tas_wy_50rm(m,:))));
  end
  if min(abs(4-tas_wy_50rm(m,:))) > 0.2
    idx_4C(m,:) = NaN;
  else
    idx_4C(m,:) = find(abs(4-tas_wy_50rm(m,:)) == min(abs(4-tas_wy_50rm(m,:))));
  end
  if min(abs(5-tas_wy_50rm(m,:))) > 0.2
    idx_5C(m,:) = NaN;
  else
    idx_5C(m,:) = find(abs(5-tas_wy_50rm(m,:)) == min(abs(5-tas_wy_50rm(m,:))));
  end
end

% -- calculate future change
ro_wy_delta_absolute  = nanmean(ro_wy(:,2,end-49:end),3) - nanmean(ro_wy(:,1,:),3); % absolute
ro_wy_delta           = ((mean(ro_wy(:,2,end-49:end),3) - nanmean(ro_wy(:,1,:),3)) ./ nanmean(ro_wy(:,1,:),3))*100; % percent
re_wy_delta           = nanmean(re_wy(:,2,end-49:end),3) - nanmean(re_wy(:,1,:),3);
pr_wy_delta_absolute  = nanmean(pr_wy(:,2,end-49:end),3) - nanmean(pr_wy(:,1,:),3); % absolute
pr_wy_delta           = ((mean(pr_wy(:,2,end-49:end),3) - nanmean(pr_wy(:,1,:),3)) ./ nanmean(pr_wy(:,1,:),3))*100; % percent
et_wy_delta           = nanmean(et_wy(:,2,end-49:end),3) - nanmean(et_wy(:,1,:),3);
pme_wy_delta          = nanmean(pme_wy(:,2,end-49:end),3) - nanmean(pme_wy(:,1,:),3);
tas_wy_delta          = nanmean(tas_wy(:,2,end-49:end),3) - nanmean(tas_wy(:,1,:),3);




% === CESM LE data ===============================================================
pathin = ['~/Dropbox/work/cmip5_ncar/runoff_efficiency_ts_collection/cesm_le/'];
start     = 1920;
ende      = 2100;
time_cesm = start:ende;
for e = 1:40
  % -- runoff
  tmp1 = ncread([pathin 'QRUNOFF.r' num2str(e) 'i1p1.192001-210012.remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
  ro_cesm(e,:) = tmp1;
  for i = 1:(length(tmp1)/12)-1
    tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
  end
  ro_cesm_wy(e,:)      = tmp2;
  ro_cesm_am(e,:)      = am(tmp1)*12;
  ro_cesm_wy_mean(e)   = nanmean(ro_cesm_wy(e,refstart-start+1:refende-start+1));
  ro_cesm_am_mean(e)   = nanmean(ro_cesm_am(e,refstart-start+1:refende-start+1));
  ro_cesm_wy_median(e) = nanmedian(ro_cesm_wy(e,refstart-start+1:refende-start+1));
  ro_cesm_wy_percent_mean(e,:) = (ro_cesm_wy(e,:)./ro_cesm_wy_mean(e))*100;
  ro_cesm_am_percent_mean(e,:) = (ro_cesm_am(e,:)./ro_cesm_am_mean(e))*100;
  ro_cesm_wy_percent_median(e,:) = (ro_cesm_wy(e,:)./ro_cesm_wy_median(e))*100;
  clear('tmp1','tmp2')
  % -- precip
  tmp1 = ncread([pathin 'PRECT.r' num2str(e) 'i1p1.192001-210012.remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
  pr_cesm(e,:) = tmp1;
  for i = 1:(length(tmp1)/12)-1
    tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
    tmp3(i) = sum(tmp1((i-1)*12+10:i*12+4));
  end
  pr_cesm_wy(e,:)          = tmp2;
  pr_cesm_am(e,:)          = am(tmp1)*12;
  pr_cesm_oct_apr(e,:)     = tmp3;
  pr_cesm_wy_mean(e)       = nanmean(pr_cesm_wy(e,refstart-start+1:refende-start+1));
  pr_cesm_am_mean(e)       = nanmean(pr_cesm_am(e,refstart-start+1:refende-start+1));
  pr_cesm_wy_median(e)       = nanmedian(pr_cesm_wy(e,refstart-start+1:refende-start+1));
  pr_cesm_wy_percent_mean(e,:) = (pr_cesm_wy(e,:)./pr_cesm_wy_mean(e))*100;
  pr_cesm_am_percent_mean(e,:) = (pr_cesm_am(e,:)./pr_cesm_am_mean(e))*100;
  pr_cesm_wy_percent_median(e,:) = (pr_cesm_wy(e,:)./pr_cesm_wy_median(e))*100;
  clear('tmp1','tmp2','tmp3')
  % -- temperature
  tmp1 = ncread([pathin 'TREFHT.r' num2str(e) 'i1p1.192001-210012.remapcon2_1x1_ts.nc'],[vari '_TS'])-273.16;
  tas_cesm(e,:) = tmp1;
  for i = 1:(length(tmp1)/12)-1
    tmp2(i) = nanmean(tmp1((i-1)*12+10:i*12+9));
    tmp3(i) = nanmean(tmp1(i*12+5:i*12+6));
  end
  tas_cesm_wy(e,:)         = tmp2;
  tas_cesm_am(e,:)         = am(tmp1);
  tas_cesm_may_jun(e,:)    = tmp3;
  tas_cesm_wy_mean(e)      = nanmean(tas_cesm_wy(e,refstart-start+1:refende-start+1));
  tas_cesm_am_mean(e)      = nanmean(tas_cesm_am(e,refstart-start+1:refende-start+1));
  tas_cesm_wy_median(e)    = nanmedian(tas_cesm_wy(e,refstart-start+1:refende-start+1));
  tas_cesm_wy_anom(e,:)    = tas_cesm_wy(e,:)-tas_cesm_wy_mean(e);
  tas_cesm_am_anom(e,:)    = tas_cesm_am(e,:)-tas_cesm_am_mean(e);
  clear('tmp1','tmp2','tmp3')
  % -- streamflow (QCHANR from River Transport Model) --
  tmp1 = ncread([pathin 'QCHANR.r' num2str(e) 'i1p1.192001-210012.lees_ferry.nc'],['QCHANR']) * 86400 *30.4 * 1e-9; % convert from m3/s to km3/year
  flow_cesm(e,:) = tmp1;
  for i = 1:(length(tmp1)/12)-1
    tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
    tmp3(i) = sum(tmp1(i*12+5:i*12+6));
  end
  flow_cesm_wy(e,:)         = tmp2;
  flow_cesm_am(e,:)         = am(squeeze(tmp1))*12;
  flow_cesm_may_jun(e,:)    = tmp3;
  flow_cesm_wy_mean(e)      = nanmean(flow_cesm_wy(e,refstart-start+1:refende-start+1));
  flow_cesm_am_mean(e)      = nanmean(flow_cesm_am(e,refstart-start+1:refende-start+1));
  flow_cesm_wy_median(e)    = nanmedian(flow_cesm_wy(e,refstart-start+1:refende-start+1));
  flow_cesm_wy_percent_mean(e,:) = (flow_cesm_wy(e,:)./flow_cesm_wy_mean(e))*100;
  flow_cesm_am_percent_mean(e,:) = (flow_cesm_am(e,:)./flow_cesm_am_mean(e))*100;
  flow_cesm_wy_anom(e,:)    = flow_cesm_wy(e,:)-flow_cesm_wy_mean(e);
  clear('tmp1','tmp2','tmp3')
end

if regout0 == 1
  % -- regress out precip from temp
  for e = 1:40
    b = regress(tas_cesm_wy(e,:)',[ones(size(pr_cesm_wy(e,:)))' pr_cesm_wy(e,:)']);
    tas_cesm_wy_pr_regressed_out(e,:) = squeeze(tas_cesm_wy(e,:)) - (b(1)+b(2)*squeeze(pr_cesm_wy(e,:)));
  end
else
  tas_cesm_wy_pr_regressed_out = tas_cesm_wy;
end
if regout1 == 1
  % -- regress out temp from precip
  for e = 1:40
    b = regress(pr_cesm_wy(e,:)',[ones(size(tas_cesm_wy(e,:)))' tas_cesm_wy(e,:)']);
    pr_cesm_wy_tas_regressed_out(e,:) = squeeze(pr_cesm_wy(e,:)) - (b(1)+b(2)*squeeze(tas_cesm_wy(e,:))) + nanmean(squeeze(pr_cesm_wy(e,:)));
    pr_cesm_wy_tas_regressed_out_percent_mean(e,:) = (pr_cesm_wy_tas_regressed_out(e,:)./mean(pr_cesm_wy_tas_regressed_out(e,:)))*100;
  end
else
  pr_cesm_wy_tas_regressed_out = pr_cesm_wy;
  pr_cesm_wy_tas_regressed_out_percent_mean = pr_cesm_wy_percent_mean;
end
if regout2 == 1
  % -- regress out precip from runoff
  for e = 1:40
    b = regress(ro_cesm_wy(e,:)',[ones(size(pr_cesm_wy(e,:)))' pr_cesm_wy(e,:)']);
    ro_cesm_wy_pr_regressed_out(e,:) = squeeze(pr_cc_wy(e,:)) - (b(1)+b(2)*squeeze(tas_cc_wy(e,:))) + nanmean(squeeze(pr_cc_wy(e,:)));
    ro_cesm_wy_pr_regressed_out_percent_mean(m,:) = (ro_cesm_wy_pr_regressed_out(e,:)./mean(ro_cesm_wy_pr_regressed_out(e,:)))*100;
  end
else
  for e = 1:40
    ro_cesm_wy_pr_regressed_out(e,:) = ro_cesm_wy(e,:);
    ro_cesm_wy_pr_regressed_out_percent_mean(e,:) = ro_cesm_wy_percent_mean(e,:);
  end
end
if regout3 == 1
  % -- regress out temp from runoff
  for e = 1:40
    b = regress(ro_cesm_wy(e,:)',[ones(size(tas_cesm_wy(e,:)))' tas_cesm_wy(e,:)']);
    ro_cesm_wy_tas_regressed_out(e,:) = squeeze(ro_cesm_wy(e,:)) - (b(1)+b(2)*squeeze(tas_cesm_wy(e,:))) + nanmean(ro_cesm_wy(e,:));
    ro_cesm_wy_tas_regressed_out_percent_mean(m,:) = (ro_cesm_wy_tas_regressed_out(e,:)./mean(ro_cesm_wy_tas_regressed_out(e,:)))*100;
  end
else
  for e = 1:40
    ro_cesm_wy_tas_regressed_out(e,:) = ro_cesm_wy(e,:);
    ro_cesm_wy_tas_regressed_out_percent_mean(e,:) = ro_cesm_wy_percent_mean(e,:);
  end
end
for e = 1:40
  tas_cesm_wy_pr_regressed_out(e,:) = tas_cesm_wy_pr_regressed_out(e,:)-nanmean(tas_cesm_wy_pr_regressed_out(e,1:86));
end

% -- calculate runoff efficiency
re_cesm_wy = (ro_cesm_wy./pr_cesm_wy)*100;




% === CESM LE control data ===============================================================
pathin = ['~/Dropbox/work/LE_control/'];
start     = 0400;
ende      = 2200;
time_cesm_control = start:ende;
% -- runoff
tmp1 = ncread([pathin 'lnd/QRUNOFF/QRUNOFF.040001-220012.remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
ro_cesm_control = tmp1;
for i = 1:(length(tmp1)/12)-1
  tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
end
ro_cesm_control_wy      = tmp2;
ro_cesm_control_wy_mean   = nanmean(ro_cesm_control_wy); %(refstart-start+1:refende-start+1));
ro_cesm_control_wy_median = nanmedian(ro_cesm_control_wy); %(refstart-start+1:refende-start+1));
ro_cesm_control_wy_percent_mean = (ro_cesm_control_wy./ro_cesm_control_wy_mean)*100;
ro_cesm_control_wy_percent_median = (ro_cesm_control_wy./ro_cesm_control_wy_median)*100;
clear('tmp1','tmp2')
% -- precip
tmp1 = ncread([pathin 'atm/PRECT/PRECT.040001-220012.remapcon2_1x1_ts.nc'],[vari '_TS']) * tf;
pr_cesm = tmp1;
for i = 1:(length(tmp1)/12)-1
  tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
  tmp3(i) = sum(tmp1((i-1)*12+10:i*12+4));
end
pr_cesm_control_wy          = tmp2;
pr_cesm_control_oct_apr     = tmp3;
pr_cesm_control_wy_mean       = nanmean(pr_cesm_control_wy); %(refstart-start+1:refende-start+1));
pr_cesm_control_wy_median       = nanmedian(pr_cesm_control_wy); %(refstart-start+1:refende-start+1));
pr_cesm_control_wy_percent_mean = (pr_cesm_control_wy./pr_cesm_control_wy_mean)*100;
pr_cesm_control_wy_percent_median = (pr_cesm_control_wy./pr_cesm_control_wy_median)*100;
clear('tmp1','tmp2','tmp3')
% -- temperature
tmp1 = ncread([pathin 'atm/TREFHT/TREFHT.040001-220012.remapcon2_1x1_ts.nc'],[vari '_TS'])-273.16;
tas_cesm = tmp1;
for i = 1:(length(tmp1)/12)-1
  tmp2(i) = nanmean(tmp1((i-1)*12+10:i*12+9));
  tmp3(i) = nanmean(tmp1(i*12+5:i*12+6));
end
tas_cesm_control_wy         = tmp2;
tas_cesm_control_may_jun    = tmp3;
tas_cesm_control_wy_mean      = nanmean(tas_cesm_control_wy); %(refstart-start+1:refende-start+1));
tas_cesm_control_wy_median      = nanmedian(tas_cesm_control_wy); %(refstart-start+1:refende-start+1));
tas_cesm_control_wy_anom  = tas_cesm_control_wy-tas_cesm_control_wy_mean;
clear('tmp1','tmp2','tmp3')
% -- streamflow (QCHANR from River Transport Model)
tmp1 = ncread([pathin 'rof/QCHANR/lees_ferry/QCHANR.040001-220012.lees_ferry.nc'],['QCHANR']) * 86400 * 30.4 * 1e-9; % convert from m3/s to km3/year
flow_cesm = tmp1;
for i = 1:(length(tmp1)/12)-1
  tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
  tmp3(i) = sum(tmp1(i*12+5:i*12+6));
end
flow_cesm_control_wy         = tmp2;
flow_cesm_control_may_jun    = tmp3;
flow_cesm_control_wy_mean    = nanmean(flow_cesm_control_wy); %(refstart-start+1:refende-start+1));
flow_cesm_control_wy_median  = nanmedian(flow_cesm_control_wy); %(refstart-start+1:refende-start+1));
flow_cesm_control_wy_anom    = flow_cesm_control_wy-flow_cesm_control_wy_mean;
clear('tmp1','tmp2','tmp3')
% -- water storage (WT from CLM)
tmp1 = ncread([pathin 'lnd/WT/WT.040001-220012.remapcon2_1x1_ts.nc'],[vari '_TS']) * tf * 1e-6; % acre feet --> km3/year
wt_cesm = tmp1;
for i = 1:(length(tmp1)/12)-1
  tmp2(i) = sum(tmp1((i-1)*12+10:i*12+9));
end
wt_cesm_control_wy         = tmp2;
wt_cesm_control_wy_mean    = nanmean(wt_cesm_control_wy); %(refstart-start+1:refende-start+1));
wt_cesm_control_wy_median  = nanmedian(wt_cesm_control_wy); %(refstart-start+1:refende-start+1));
wt_cesm_control_wy_anom    = wt_cesm_control_wy-wt_cesm_control_wy_mean;
clear('tmp1','tmp2','tmp3')

if regout0 == 1
  % -- regress out precip from temp
  b = regress(tas_cesm_control_wy',[ones(size(pr_cesm_control_wy))' pr_cesm_control_wy']);
  tas_cesm_control_wy_pr_regressed_out = squeeze(tas_cesm_control_wy) - (b(1)+b(2)*squeeze(pr_cesm_control_wy));
else
  tas_cesm_control_wy_pr_regressed_out = tas_cesm_control_wy;
end
if regout1 == 1
  % -- regress out temp from precip
  b = regress(pr_cesm_control_wy',[ones(size(tas_cesm_control_wy))' tas_cesm_control_wy']);
  pr_cesm_control_wy_tas_regressed_out = pr_cesm_control_wy - (b(1)+b(2)*tas_cesm_control_wy) + nanmean(pr_cesm_control_wy);
  pr_cesm_control_wy_tas_regressed_out_percent_mean = (pr_cesm_control_wy_tas_regressed_out/mean(pr_cesm_control_wy_tas_regressed_out))*100;
else
  pr_cesm_control_wy_tas_regressed_out = pr_cesm_control_wy;
  pr_cesm_control_wy_tas_regressed_out_percent_mean = pr_cesm_control_wy_percent_mean;
end
if regout2 == 1
  % -- regress out precip from runoff
  b = regress(ro_cesm_control_wy',[ones(size(pr_cesm_control_wy))' pr_cesm_control_wy']);
  ro_cesm_control_wy_pr_regressed_out = ro_cesm_control_wy - (b(1)+b(2)*squeeze(pr_cesm_control_wy)) + nanmean(ro_cesm_control_wy);
  ro_cesm_control_wy_pr_regressed_out_percent_mean = (ro_cesm_control_wy_pr_regressed_out/mean(ro_cesm_control_wy_pr_regressed_out))*100;
end
if regout3 == 1
  % -- regress out temp from runoff
  b = regress(ro_cesm_control_wy_percent_mean',[ones(size(tas_cesm_control_wy))' tas_cesm_control_wy']);
  ro_cesm_control_wy_tas_regressed_out = ro_cesm_control_wy - (b(1)+b(2)*tas_cesm_control_wy) + nanmean(ro_cesm_control_wy);
  ro_cesm_control_wy_tas_regressed_out_percent_mean = (ro_cesm_control_wy_tas_regressed_out/mean(ro_cesm_control_wy_tas_regressed_out))*100;
end
tas_cesm_control_wy_pr_regressed_out_anom = tas_cesm_control_wy_pr_regressed_out-nanmean(tas_cesm_control_wy_pr_regressed_out); %(1:86));

% -- calculate runoff efficiency
re_cesm_control_wy = (ro_cesm_control_wy./pr_cesm_control_wy)*100;





% === OBSERVATIONS =============================================================
pathin = '~/Dropbox/work/';

time_obs_common = 1950:2008; %1929:2008; %1906:2012;

% -- Flow: from USBR naturalized flow (in acre feet/month)
if strcmp(vari,'UC')==1
  % -- begin Upper Colorado data --
  % -- Reclamation data:
  data              = load([pathin '/streamflow_data/UC_NaturalFlows.csv']);
  timem_obs         = 1905+1/12:1/12:2016;
  ro_obs            = NaN(length(timem_obs),1);
  % ro_obs(10:end)    = data(:,3)*1e-6; % convert to million acre feet (MAF)
  ro_obs(10:end)    = data(:,3)*1e-6 * tf; % in mm
  time_ro_obs       = 1906:2015;
  % -- UC River Commission data:
  data2             = load([pathin '/streamflow_data/upper_colorado_river_commission_virgin_flow.edit.txt']);
  % --> http://www.ucrcommission.com/RepDoc/UCRCAnnualReports/68_UCRC_Annual_Report.pdf
  % ro_obs2_wy        = data2(:,3); % already in million acre feet (MAF)
  ro_obs2_wy        = data2(:,3) * tf; % in mm
  time_ro_obs2      = 1896:2016;
  clear('data2')
  % -- Global River Data Centre (GRDC) data (in m3/s):
  grdc_id           = 4152450;
  data              = read_mixed_csv([pathin '/streamflow_data/GRDC_data/Export/' num2str(grdc_id) '_Q_Month.txt'],';');
  for i = 39:length(data(:,2))
    ro_obs3(i-38)         = str2num(cell2mat(data(i,4))) * 86400*30.4*1e-9; % convert to km3/month
  end
  [y,m,d]           = datevec(data(39:end,1));
  timem_ro_obs3     = y(1)+(m(1)/12):1/12:y(end)+(m(end)/12);
  time_ro_obs3      = y(1)+1:y(end)-1;
elseif strcmp(vari,'CO_DALLES')==1
  data    = read_mixed_csv([pathin '/streamflow_data/NRNI_Flows_The_Dalles.csv'],',');
  tmp     = cell2mat(data(8:end,1));
  for i = 1:length(tmp)
    yyyy(i)    = str2double(tmp(i,1:4));
    mm(i)      = str2double(tmp(i,6:7));
    dd(i)      = str2double(tmp(i,9:10));
  end
  yyyy(yyyy==2028) = 1928;
  yyyy(yyyy==2029) = 1929;
  flow    = str2double(data(8:end,2));
  years   = unique(yyyy);
  timem   = years(1)+1/12:1/12:years(end)+1;
  for y = 1:length(years)
    yidx = find(yyyy==years(y));
    for m = 1:12
      midx = find(mm==m);
      idx = intersect(yidx,midx);
      flow2(y,m,1:length(idx)) = flow(idx);
    end
  end
  ro_obs = (sum(flow2,3))';
  ro_obs = ro_obs(:) * 86400 * 2.29569e-5 * 1e-6 *tf; % in mm per month
  ro_obs(ro_obs==0) = NaN;
  time_ro_obs = 1929:2008;
  % -- Global River Data Centre (GRDC) data (in m3/s):
  grdc_id           = 4115200;
  data              = read_mixed_csv([pathin '/streamflow_data/GRDC_data/Export/' num2str(grdc_id) '_Q_Month.txt'],';');
  for i = 46:length(data(:,2))
    ro_obs3(i-45)       = str2num(cell2mat(data(i,4))) * 86400*30.4*1e-9; % convert to km3/month
  end
  [y,m,d]           = datevec(data(46:end,1));
  timem_ro_obs3     = y(1)+(m(1)/12):1/12:y(end)+(m(end)/12);
  time_ro_obs3      = y(1)+1:y(end)-1;
elseif strcmp(vari,'NS')==1
  data              = load([pathin '/streamflow_data/cali_flows/northern_sierras_four_index_flow_190510-201805.edit.csv']);
  data              = data(4:end-5,:);
  timem_obs         = 1906+1/12:1/12:2018;
  ro_obs            = data(:,end)*1e-6 * tf; % mm per month
  time_ro_obs       = 1907:2017;
end
% -- make water year totals:
for i = 2:length(ro_obs)/12
  ro_obs_wy(i-1)  = nansum(ro_obs((i-2)*12+10:(i-1)*12+9));
end
if strcmp(vari,'NS')==1
  'no third party obs'
else
  for i = 2:length(ro_obs3)/12
    ro_obs3_wy(i-1)  = nansum(ro_obs3((i-2)*12+10:(i-1)*12+9));
  end
  ro_obs3_wy(ro_obs3_wy<=0)= NaN;
end


% -- Precipitation (obs) --
if obs_data == 1
  obs_data_name = 'PRISM';
  pr_obs 	= ncread(['' pathin 'prism/precip/prism_precip_189501-201712_uc_ts.nc'],[vari '_TS']) * tf;
  for i = 2:length(pr_obs)/12
    pr_obs_wy(i-1) = sum(pr_obs((i-2)*12+10:(i-1)*12+9));
  end
  time_pr_obs	= 1896:2017;
  % -- alternative data
  obs2_data_name = 'GPCC';
  pr_obs2 	= ncread(['' pathin 'observations/precipitation/gpcc/gpcc_precip_colo_colu_ts.nc'],[vari '_TS']) * tf;
  for i = 2:length(pr_obs2)/12
    pr_obs2_wy(i-1) = sum(pr_obs2((i-2)*12+10:(i-1)*12+9));
  end
  time_pr_obs2	= 1902:2016;
elseif obs_data == 2
  obs_data_name = 'GPCC';
  pr_obs 	= ncread(['' pathin 'observations/precipitation/gpcc/gpcc_precip_colo_colu_ts.nc'],[vari '_TS']) * tf;
  for i = 2:length(pr_obs)/12
    pr_obs_wy(i-1) = sum(pr_obs((i-2)*12+10:(i-1)*12+9));
  end
  time_pr_obs	= 1902:2016;
  % -- alternative data
  if strcmp(vari,'UC')==1
    obs2_data_name = 'PRISM';
    pr_obs2 	= ncread(['' pathin 'prism/precip/prism_precip_189501-201712_uc_ts.nc'],[vari '_TS']) * tf;
    for i = 2:length(pr_obs2)/12
      pr_obs2_wy(i-1) = sum(pr_obs2((i-2)*12+10:(i-1)*12+9));
    end
    time_pr_obs2	= 1896:2017;
  else
    obs2_data_name = 'GPCC';
    pr_obs2_wy    = pr_obs_wy;
    time_pr_obs2  = time_pr_obs;
  end
elseif obs_data == 3
  obs_data_name = 'Livneh';
  pr_obs 	= ncread(['' pathin 'livneh/livneh_precip_colo_colu_ts.nc'],[vari '_TS']) * tf;
  for i = 2:length(pr_obs)/12
    pr_obs_wy(i-1) = sum(pr_obs((i-2)*12+10:(i-1)*12+9));
  end
  time_pr_obs	= 1916:2011;
  obs2_data_name = 'Livneh';
  pr_obs2_wy    = pr_obs_wy;
  time_pr_obs2  = time_pr_obs;
end

% -- cut to common length
pr_obs_wy      = pr_obs_wy(find(time_pr_obs==time_obs_common(1)):find(time_pr_obs==time_obs_common(end)));
ro_obs_wy      = ro_obs_wy(find(time_ro_obs==time_obs_common(1)):find(time_ro_obs==time_obs_common(end)));
ro_obs_wy_percent_mean = (ro_obs_wy./mean(ro_obs_wy))*100;
pr_obs_wy_percent_mean = (pr_obs_wy./mean(pr_obs_wy))*100;
ro_obs_wy_percent_median = (ro_obs_wy./median(ro_obs_wy))*100;
pr_obs_wy_percent_median = (pr_obs_wy./median(pr_obs_wy))*100;
pr_obs2_wy     = pr_obs2_wy(find(time_pr_obs2==time_obs_common(1)):find(time_pr_obs2==time_obs_common(end)));
pr_obs2_wy_percent_mean = (pr_obs2_wy./mean(pr_obs2_wy))*100;

% -- UC temp (obs) --
if obs_data == 1
  tas_obs 	= ncread(['' pathin 'prism/tmean/prism_tmean_189501-201712_uc_ts.nc'],[vari '_TS']);
  tas_obs_am = am(tas_obs);
  tas_obs_am = tas_obs_am(2:end); % adjust to same length as WY
  for i = 2:length(tas_obs)/12
    tas_obs_wy(i-1) = nanmean(tas_obs((i-2)*12+10:(i-1)*12+9));
  end
  time_tas_obs	= 1896:2017;
elseif obs_data == 2
  tas_obs 	= ncread(['' pathin 'observations/temperature/best/best_temperature_185001-201712_colo_colu_ts.nc'],[vari '_TS']);
  tas_obs_am = am(tas_obs);
  tas_obs_am = tas_obs_am(2:end); % adjust to same length as WY
  for i = 2:length(tas_obs)/12
    tas_obs_wy(i-1) = nanmean(tas_obs((i-2)*12+10:(i-1)*12+9));
  end
  time_tas_obs	= 1851:2017;
elseif obs_data == 3
  obs_data_name = 'Livneh';
  tas_obs 	= ncread(['' pathin 'livneh/livneh_tmean_colo_colu_ts.nc'],[vari '_TS']);
  for i = 2:length(pr_obs)/12
    tas_obs_wy(i-1) = mean(tas_obs((i-2)*12+10:(i-1)*12+9));
  end
  time_tas_obs	= 1916:2011;
end
tas_obs_wy = tas_obs_wy-mean(tas_obs_wy(refstart-time_tas_obs(1)+1:refende-time_tas_obs(1)));

% -- cut to common length
tas_obs_wy      = tas_obs_wy(find(time_tas_obs==time_obs_common(1)):find(time_tas_obs==time_obs_common(end)));
tas_obs_wy_anom = tas_obs_wy;

if regout0 == 1
  % -- regress out precip from temp
  b = regress(tas_obs_wy',[ones(size(pr_obs_wy')) pr_obs_wy']);
  tas_obs_wy_pr_regressed_out = tas_obs_wy - (b(1)+b(2)*pr_obs_wy);
else
  tas_obs_wy_pr_regressed_out = tas_obs_wy;
end
if regout1 == 1
  % -- regress out temp from precip
  b = regress(pr_obs_wy',[ones(size(tas_obs_wy')) tas_obs_wy']);
  pr_obs_wy_tas_regressed_out = pr_obs_wy - (b(1)+b(2)*tas_obs_wy) + nanmean(pr_obs_wy);
  pr_obs_wy_tas_regressed_out_percent_mean = (pr_obs_wy_tas_regressed_out./mean(pr_obs_wy_tas_regressed_out))*100;
else
  pr_obs_wy_tas_regressed_out = pr_obs_wy;
  pr_obs_wy_tas_regressed_out_percent_mean = pr_obs_wy_percent_mean;
end
if regout2 == 1
  % -- regress out precip from runoff
  b = regress(ro_obs_wy',[ones(size(pr_obs_wy')) pr_obs_wy']);
  ro_obs_wy_pr_regressed_out = ro_obs_wy - (b(1)+b(2)*pr_obs_wy) + nanmean(ro_obs_wy);
  ro_obs_wy_pr_regressed_out_percent_mean = (ro_obs_wy_pr_regressed_out./mean(ro_obs_wy_pr_regressed_out))*100;
else
  ro_obs_wy_pr_regressed_out = ro_obs_wy;
  ro_obs_wy_pr_regressed_out_percent_mean = ro_obs_wy_percent_mean;
end
if regout3 == 1
  % -- regress out temp from runoff
  b = regress(ro_obs_wy',[ones(size(tas_obs_wy')) tas_obs_wy']);
  ro_obs_wy_tas_regressed_out = ro_obs_wy - (b(1)+b(2)*tas_obs_wy) + nanmean(ro_obs_wy);
  ro_obs_wy_tas_regressed_out_percent_mean = (ro_obs_wy_tas_regressed_out./mean(ro_obs_wy_tas_regressed_out))*100;
else
  ro_obs_wy_tas_regressed_out = ro_obs_wy;
  ro_obs_wy_tas_regressed_out_percent_mean = ro_obs_wy_percent_mean;
end

% -- UC RE --
re_obs_wy       = (ro_obs_wy./pr_obs_wy)*100;

% -- make sure matrix size is correct:
re_obs_wy = re_obs_wy';
ro_obs_wy = ro_obs_wy';


if anoml == 1
  % -- make everything as anomalies to their long-term median:
  for m = 1:length(models_common)
    for s = 1:length(scen)
      tas_wy(m,s,:) = tas_wy(m,s,:)-tas_wy_median(m,1,:);
      pr_wy(m,s,:)  = pr_wy(m,s,:)-pr_wy_median(m,1,:);
      et_wy(m,s,:)  = et_wy(m,s,:)-et_wy_median(m,1,:);
      ro_wy(m,s,:)  = ro_wy(m,s,:)-ro_wy_median(m,1,:);
      re_wy(m,s,:)  = re_wy(m,s,:)-re_wy_median(m,1,:);
    end
  end
  % -- make everything as anomalies to the first 30 years of record:
  for m = 1:length(models_cc)
    for s = 1:length(scen2)
      % -- relative anomalies (for pr, et, pme, ro)
      pr_cc_wy_percent(m,s,:)  = (pr_cc_wy(m,s,:)-mean(pr_cc_wy(m,s,1:30),3)) ./ nanmean(pr_cc_wy(m,s,1:30),3);
      et_cc_wy_percent(m,s,:)  = (et_cc_wy(m,s,:)-mean(et_cc_wy(m,s,1:30),3)) ./ nanmean(et_cc_wy(m,s,1:30),3);
      pme_cc_wy_percent(m,s,:) = (pme_cc_wy(m,s,:)-mean(pme_cc_wy(m,s,1:30),3)) ./ nanmean(pme_cc_wy(m,s,1:30),3);
      ro_cc_wy_percent(m,s,:)  = (ro_cc_wy(m,s,:)-mean(ro_cc_wy(m,s,1:30),3)) ./ nanmean(ro_cc_wy(m,s,1:30),3);
    end
  end
  for m = 1:length(models_cc)
    for s = 1:length(scen2)
      % -- absolute anomalies
      tas_cc_wy(m,s,:) = tas_cc_wy(m,s,:)-mean(tas_cc_wy(m,s,1:30),3);
      pr_cc_wy(m,s,:)  = pr_cc_wy(m,s,:)-mean(pr_cc_wy(m,s,1:30),3);
      et_cc_wy(m,s,:)  = et_cc_wy(m,s,:)-mean(et_cc_wy(m,s,1:30),3);
      pme_cc_wy(m,s,:) = pme_cc_wy(m,s,:)-mean(pme_cc_wy(m,s,1:30),3);
      ro_cc_wy(m,s,:)  = ro_cc_wy(m,s,:)-mean(ro_cc_wy(m,s,1:30),3);
      re_cc_wy(m,s,:)  = re_cc_wy(m,s,:)-mean(re_cc_wy(m,s,1:30),3);
    end
  end
  tas_obs_wy = tas_obs_wy-median(tas_obs_wy);
  pr_obs_wy = pr_obs_wy-median(pr_obs_wy);
  ro_obs_wy = ro_obs_wy-median(ro_obs_wy);
  re_obs_wy   = re_obs_wy-median(re_obs_wy);
end

% -- CO2 vmr from CESM1 rcp85:
co2 = am(ncread('~/Dropbox/work/LE/co2vmr.hist2rcp85.185001-208012.nc','CO2VMR'))*1e6; % in ppm
time_co2 = 1850:2080;

% -- find models that get wetter and drier
idx_dry_delta = find(pr_wy_delta<0);
idx_wet_delta = find(pr_wy_delta>0);




% === MAKE COMMON TIME PERIOD, at least start year ===
ro_h2r_wy   = ro_h2r_wy(:,find(time_h2r==time_obs_common(1)):end);
pme_h2r_wy  = pme_h2r_wy(:,find(time_h2r==time_obs_common(1)):end);
pr_h2r_wy   = pr_h2r_wy(:,find(time_h2r==time_obs_common(1)):end);
re_h2r_wy   = re_h2r_wy(:,find(time_h2r==time_obs_common(1)):end);
tas_h2r_wy  = tas_h2r_wy(:,find(time_h2r==time_obs_common(1)):end);
tas_h2r_wy_pr_regressed_out_anom  = tas_h2r_wy_pr_regressed_out_anom(:,find(time_h2r==time_obs_common(1)):end);
tas_h2r_wy_anom  = tas_h2r_wy_anom(:,find(time_h2r==time_obs_common(1)):end);
co2_h2r_anom  = co2(find(time_co2==time_obs_common(1)):end)-300;

ro_h2r_wy_percent_mean            = ro_h2r_wy_percent_mean(:,find(time_h2r==time_obs_common(1)):end);
ro_h2r_wy_pr_regressed_out_percent_mean   = ro_h2r_wy_pr_regressed_out_percent_mean(:,find(time_h2r==time_obs_common(1)):end);
ro_h2r_wy_tas_regressed_out_percent_mean  = ro_h2r_wy_tas_regressed_out_percent_mean(:,find(time_h2r==time_obs_common(1)):end);

pr_h2r_wy_percent_mean            = pr_h2r_wy_percent_mean(:,find(time_h2r==time_obs_common(1)):end);
pr_h2r_wy_tas_regressed_out_percent_mean      = pr_h2r_wy_tas_regressed_out_percent_mean(:,find(time_h2r==time_obs_common(1)):end);

ro_h2r_wy_percent_median          = ro_h2r_wy_percent_median(:,find(time_h2r==time_obs_common(1)):end);
pr_h2r_wy_percent_median          = pr_h2r_wy_percent_median(:,find(time_h2r==time_obs_common(1)):end);



% --------------------------------------------------------------------------------
% === calculate runoff sensitivities ===

% -- CMIP5 historical -------
nn = 80;
for m = 1:length(models_common)
  if co == 0
    [pr_ro_reg_h2r_wy(m,1) pr_ro_reg_h2r_wy(m,2:3) tas_ro_reg_h2r_wy(m,1) tas_ro_reg_h2r_wy(m,2:3) int_ro_reg_h2r_wy(m,1) int_ro_reg_h2r_wy(m,2:3)]...
                           = runoff_sens_reg(ro_h2r_wy_percent_mean(m,1:nn)-100,pr_h2r_wy_percent_mean(m,1:nn)-100,tas_h2r_wy_anom(m,1:nn));
    [pr_ro_reg_h2r_wy_dt(m,1) pr_ro_reg_h2r_wy_dt(m,2:3) tas_ro_reg_h2r_wy_dt(m,1) tas_ro_reg_h2r_wy_dt(m,2:3) int_ro_reg_h2r_wy_dt(m,1) int_ro_reg_h2r_wy_dt(m,2:3)]...
                          = runoff_sens_reg(detrend(ro_h2r_wy_percent_mean(m,1:nn)-100),detrend(pr_h2r_wy_percent_mean(m,1:nn)-100),detrend(tas_h2r_wy_anom(m,1:nn)));
  else
    [pr_ro_reg_h2r_wy(m,1) pr_ro_reg_h2r_wy(m,2:3) tas_ro_reg_h2r_wy(m,1) tas_ro_reg_h2r_wy(m,2:3) int_ro_reg_h2r_wy(m,1) int_ro_reg_h2r_wy(m,2:3)]...
                          = runoff_sens_reg_with_carryover(ro_h2r_wy_percent_mean(m,1:nn)-100,pr_h2r_wy_percent_mean(m,1:nn)-100,tas_h2r_wy_anom(m,1:nn));
  end

  [pr_ro_reg_h2r_wy_co2(m,1) pr_ro_reg_h2r_wy_co2(m,2:3) tas_ro_reg_h2r_wy_co2(m,1) tas_ro_reg_h2r_wy_co2(m,2:3) co2_ro_reg_h2r_wy_co2(m,1) co2_ro_reg_h2r_wy_co2(m,2:3) int_ro_reg_h2r_wy_co2(m,1) int_ro_reg_h2r_wy_co2(m,2:3)]...
                          = runoff_sens_reg_with_co2(ro_h2r_wy_percent_mean(m,1:nn)-100,pr_h2r_wy_percent_mean(m,1:nn)-100,tas_h2r_wy_anom(m,1:nn),co2_h2r_anom(1:nn));
  pr_ro_reg2_h2r_wy(m,1)   = elast_sankar(ro_h2r_wy_percent_mean(m,1:nn),pr_h2r_wy_percent_mean(m,1:nn));
  pr_ro_reg2_h2r_wy(m,2:3) = elast_sankar(ro_h2r_wy_percent_mean(m,1:nn),pr_h2r_wy_percent_mean(m,1:nn));
end


% -- CMIP5 piControl and 1pctCO2 -------
for m = 1:length(models_common)
  % -- ro vs pr
  for s = 1:2 % -- piControl and 1pctCO2
    if co == 0
      [pr_ro_reg_wy(m,s,1) pr_ro_reg_wy(m,s,2:3) tas_ro_reg_wy(m,s,1) tas_ro_reg_wy(m,s,2:3) int_ro_reg_wy(m,s,1) int_ro_reg_wy(m,s,2:3)]...
       = runoff_sens_reg(squeeze(ro_wy_percent_mean(m,s,1:nn))-100,squeeze(pr_wy_percent_mean(m,s,1:nn))-100,squeeze(tas_wy_anom(m,s,1:nn)));
    else
      [pr_ro_reg_wy(m,s,1) pr_ro_reg_wy(m,s,2:3) tas_ro_reg_wy(m,s,1) tas_ro_reg_wy(m,s,2:3) int_ro_reg_wy(m,s,1) int_ro_reg_wy(m,s,2:3)]...
       = runoff_sens_reg_with_carryover(squeeze(ro_wy_percent_mean(m,s,1:nn))-100,squeeze(pr_wy_percent_mean(m,s,1:nn))-100,squeeze(tas_wy_anom(m,s,1:nn)));
    end
    pr_ro_reg2_wy(m,s,1) = elast_sankar(squeeze(ro_wy_percent_mean(m,s,1:nn)),squeeze(pr_wy_percent_mean(m,s,1:nn)));
    pr_ro_reg2_wy(m,s,2:3) = elast_sankar(squeeze(ro_wy_percent_mean(m,s,1:nn)),squeeze(pr_wy_percent_mean(m,s,1:nn)));
  end
end


% -- CMIP5 carbon cycle models -------
nn = 100; % number of years in piControl
s = 2;
for m = 1:length(models_cc)
  if co == 0
    [pr_ro_reg_cc_wy(m,s,1) pr_ro_reg_cc_wy(m,s,2:3) tas_ro_reg_cc_wy(m,s,1) tas_ro_reg_cc_wy(m,s,2:3) int_ro_reg_cc_wy(m,s,1) int_ro_reg_cc_wy(m,s,2:3)]...
     = runoff_sens_reg(squeeze(ro_cc_wy_percent_mean(m,s,1:nn))-100,squeeze(pr_cc_wy_percent_mean(m,s,1:nn))-100,squeeze(tas_cc_wy_anom(m,s,1:nn)));
  else
    [pr_ro_reg_cc_wy(m,s,1) pr_ro_reg_cc_wy(m,s,2:3) tas_ro_reg_cc_wy(m,s,1) tas_ro_reg_cc_wy(m,s,2:3) int_ro_reg_cc_wy(m,s,1) int_ro_reg_cc_wy(m,s,2:3)]...
     = runoff_sens_reg_with_carryover(squeeze(ro_cc_wy_percent_mean(m,s,1:nn))-100,squeeze(pr_cc_wy_percent_mean(m,s,1:nn))-100,squeeze(tas_cc_wy_anom(m,s,1:nn)));
  end
  pr_ro_reg2_cc_wy(m,s,1) = elast_sankar(ro_cc_wy_percent_mean(m,s,1:nn),pr_cc_wy_percent_mean(m,s,1:nn));
  pr_ro_reg2_cc_wy(m,s,2:3) = elast_sankar(ro_cc_wy_percent_mean(m,s,1:nn),pr_cc_wy_percent_mean(m,s,1:nn));
end

% -- Obs -------
if co == 0
  % -- multi-var regression (using Livneh and Livneh)
  [pr_ro_reg_obs_wy(1) pr_ro_reg_obs_wy(2:3) tas_ro_reg_obs_wy(1) tas_ro_reg_obs_wy(2:3) int_ro_reg_obs_wy(1) int_ro_reg_obs_wy(2:3)]...
   = runoff_sens_reg(ro_obs_wy_percent_mean-100,pr_obs_wy_percent_mean-100,tas_obs_wy_anom);
  [pr_ro_reg_obs_wy_dt(1) pr_ro_reg_obs_wy_dt(2:3) tas_ro_reg_obs_wy_dt(1) tas_ro_reg_obs_wy_dt(2:3) int_ro_reg_obs_wy_dt(1) int_ro_reg_obs_wy_dt(2:3)]...
   = runoff_sens_reg(detrend(ro_obs_wy_percent_mean-100),detrend(pr_obs_wy_percent_mean-100),detrend(tas_obs_wy_anom));
else
  % -- multi-var regression with carry-over term (Milly et al. 2018, WRR, Eq. (18)):
  [pr_ro_reg_obs_wy(1) pr_ro_reg_obs_wy(2:3) tas_ro_reg_obs_wy(1) tas_ro_reg_obs_wy(2:3) co_term_obs_wy(1) co_term_obs_wy(2:3) int_ro_reg_obs_wy(1) int_ro_reg_obs_wy(2:3)]...
   = runoff_sens_reg_with_carryover(ro_obs_wy_percent_mean-100,pr_obs_wy_percent_mean-100,tas_obs_wy_anom);
end
% % -- multi-var regression after applying a running mean to all time series
% rx = 5;
% [pr_ro_reg_obs_wy(1) pr_ro_reg_obs_wy(2:3) tas_ro_reg_obs_wy(1) tas_ro_reg_obs_wy(2:3) int_ro_reg_obs_wy(1) int_ro_reg_obs_wy(2:3)]...
%  = runoff_sens_reg(rm(ro_obs_wy_percent_mean-100,rx),rm(pr_obs_wy_percent_mean-100,rx),rm(tas_obs_wy_anom,rx))
 % -- added term for CO2
[pr_ro_reg_obs_wy_co2(1) pr_ro_reg_obs_wy_co2(2:3) tas_ro_reg_obs_wy_co2(1) tas_ro_reg_obs_wy_co2(2:3) co2_ro_reg_obs_wy_co2(1) co2_ro_reg_obs_wy_co2(2:3) int_ro_reg_obs_wy_co2(1) int_ro_reg_obs_wy_co2(2:3)]...
 = runoff_sens_reg_with_co2(ro_obs_wy_percent_mean-100,pr_obs_wy_percent_mean-100,tas_obs_wy_anom,co2_h2r_anom(1:59));

% -- elasticity from Sankar and Vogel 2002
pr_ro_reg2_obs_wy(1)    = elast_sankar(ro_obs_wy_percent_mean,pr_obs_wy_percent_mean);
pr_ro_reg2_obs_wy(2:3)  = elast_sankar(ro_obs_wy_percent_mean,pr_obs_wy_percent_mean);
pr_ro_reg2_obs2_wy(1)   = elast_sankar(ro_obs_wy_percent_mean,pr_obs2_wy_percent_mean);
pr_ro_reg2_obs2_wy(2:3) = elast_sankar(ro_obs_wy_percent_mean,pr_obs2_wy_percent_mean);
% -- create random but plausible combination of P and T sensitivities for uncertainty estimate later:
ni = 1000;
nj = 1000;
tmp1 = normrnd(tas_ro_reg_obs_wy(1),abs(tas_ro_reg_obs_wy(3)-tas_ro_reg_obs_wy(2))/4,ni,1);
tmp2 = normrnd(pr_ro_reg_obs_wy(1),abs(pr_ro_reg_obs_wy(3)-pr_ro_reg_obs_wy(2))/4,nj,1);

pred_rmse = NaN(ni,nj,length(ro_obs_wy_percent_mean));
for i = 1:ni
  pred_tas = tmp1(i)*tas_obs_wy_anom;
  for j = 1:nj
    pred_pr         = tmp2(j)*(pr_obs_wy_percent_mean-100);
    pred            = pred_tas+pred_pr;
    pred_rmse(i,j)  = rmse((ro_obs_wy_percent_mean-100),pred);
  end
end
for i = 1:ni
   [m,idx] = min(pred_rmse(i,:));
   pr_index_given_tas(i) = idx;
end
tas_ro_reg_obs_wy_mc  = tmp1;
pr_ro_reg_obs_wy_mc   = tmp2;%(pr_index_given_tas);


% -- CESM LE -------
for e = 1:40
  pr_ro_reg2_cesm_wy(e,1) = elast_sankar(ro_cesm_wy_percent_mean(e,1:nn),pr_cesm_wy_percent_mean(e,1:nn));
  pr_ro_reg2_cesm_wy(e,2:3) = elast_sankar(ro_cesm_wy_percent_mean(e,1:nn),pr_cesm_wy_percent_mean(e,1:nn));
  if co == 0
    % -- multi-var regression:
    [pr_ro_reg_cesm_wy(e,1) pr_ro_reg_cesm_wy(e,2:3) tas_ro_reg_cesm_wy(e,1) tas_ro_reg_cesm_wy(e,2:3) int_ro_reg_cesm_wy(e,1) int_ro_reg_cesm_wy(e,2:3)]...
     = runoff_sens_reg(ro_cesm_wy_percent_mean(e,1:nn)-100,pr_cesm_wy_percent_mean(e,1:nn)-100,tas_cesm_wy_anom(e,1:nn));
    % -- multi-var regression but using flow instead of runoff:
    [pr_ro_reg_cesm_wy_with_flow(e,1) pr_ro_reg_cesm_wy_with_flow(e,2:3) tas_ro_reg_cesm_wy_with_flow(e,1) tas_ro_reg_cesm_wy_with_flow(e,2:3) int_ro_reg_cesm_wy_with_flow(e,1) int_ro_reg_cesm_wy_with_flow(e,2:3)]...
     = runoff_sens_reg(flow_cesm_wy_percent_mean(e,1:nn)-100,pr_cesm_wy_percent_mean(e,1:nn)-100,tas_cesm_wy_anom(e,1:nn));
  else
    % -- multi-var regression with carry-over term (Milly et al. 2018, WRR, Eq. (18)):
    % [pr_ro_reg_cesm_wy_with_ro_and_carryover(e,1) pr_ro_reg_cesm_wy_with_ro_and_carryover(e,2:3) tas_ro_reg_cesm_wy_with_ro_and_carryover(e,1) tas_ro_reg_cesm_wy_with_ro_and_carryover(e,2:3) co_term_cesm_wy(e,1) co_term_cesm_wy(e,2:3) int_ro_reg_cesm_wy_with_ro_and_carryover(e,1) int_ro_reg_cesm_wy_with_ro_and_carryover(e,2:3)]...
    % = runoff_sens_reg_with_carryover(ro_cesm_wy_percent_mean(e,1:nn)-100,pr_cesm_wy_percent_mean(e,1:nn)-100,tas_cesm_wy_anom(e,1:nn));
    [pr_ro_reg_cesm_wy(e,1) pr_ro_reg_cesm_wy(e,2:3) tas_ro_reg_cesm_wy(e,1) tas_ro_reg_cesm_wy(e,2:3) co_term_cesm_wy(e,1) co_term_cesm_wy(e,2:3) int_ro_reg_cesm_wy(e,1) int_ro_reg_cesm_wy(e,2:3)]...
    = runoff_sens_reg_with_carryover(ro_cesm_wy_percent_mean(e,1:nn)-100,pr_cesm_wy_percent_mean(e,1:nn)-100,tas_cesm_wy_anom(e,1:nn));
    % -- multi-var regression with carry-over term (Milly et al. 2018, WRR, Eq. (18)), but using flow instead of runoff:
    [pr_ro_reg_cesm_wy_with_flow_and_carryover(e,1) pr_ro_reg_cesm_wy_with_flow_and_carryover(e,2:3) tas_ro_reg_cesm_wy_with_flow_and_carryover(e,1) tas_ro_reg_cesm_wy_with_flow_and_carryover(e,2:3) co_term_cesm_wy(e,1) co_term_cesm_wy(e,2:3) int_ro_reg_cesm_wy_with_flow_and_carryover(e,1) int_ro_reg_cesm_wy_with_flow_and_carryover(e,2:3)]...
    = runoff_sens_reg_with_carryover(flow_cesm_wy_percent_mean(e,1:nn)-100,pr_cesm_wy_percent_mean(e,1:nn)-100,tas_cesm_wy_anom(e,1:nn));
  end
end


% --------------------------------------------------------------------------------



% -- calc future change and spread for basins ----
ys = 2021;
ye = 2050;
ave = nanmean(pr_h2r_wy_percent_mean(:,ys-time_obs_common(1)+1:ye-time_obs_common(1)+1),2);
max(ave)
min(ave)
prctile(ave,95)
prctile(ave,5)
plot(ave)
hline(100,'k')
% return
% ------------------------------------------------



% -- list sensitivities for text part in paper --
pr_ro_reg_obs_wy = pr_ro_reg_obs_wy
tas_ro_reg_obs_wy = tas_ro_reg_obs_wy
pr_ro_reg_h2r_wy_min = min(pr_ro_reg_h2r_wy(:,1))
pr_ro_reg_h2r_wy_max = max(pr_ro_reg_h2r_wy(:,1))
tas_ro_reg_h2r_wy_min = min(tas_ro_reg_h2r_wy(:,1))
tas_ro_reg_h2r_wy_max = max(tas_ro_reg_h2r_wy(:,1))

% return











%% ==========================================================================================================================
%% === PLOTTING =============================================================================================================
close all




if figS6 == 1
% -----------------------------------------------------------------------------
% Figure S6?
% ----
  % -- comparing Reclamation and GRDC runoff (only for UC and CO_DALLES)
  close all
  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 20 10])
  hold on
  h1 = plot(time_ro_obs,ro_obs_wy)
  h2 = plot(time_ro_obs3,ro_obs3_wy,'r')
  box on
  title(['Water year flow ' vari_name ])
  xlabel('Time (Year)')
  ylabel('Flow (km^3/year)')
  legend([h1 h2],'Official record (naturalized)','Global Runoff Data Centre data','Location','NorthWest')
  % legend boxoff

  set(gcf,'PaperPositionMode','auto');
  fileo = [figpath 'Fig_S6_' vari];
  print('-r300','-loose', '-depsc', [fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  return
end







if figS8 == 1
% -----------------------------------------------------------------------------
% Figure S8
% ----
close all

  t = 2; % tas degC anomaly relative to reference period

  % futstart  = 2021; % 2001 2050
  % futende   = 2070; % 2050 2099
  % or...
  futtas    = t; % tas degC anomaly relative to reference period
  futwl     = 40; % window length over which 'futtas' is defined (in years)

  magnif1   = 50; % T sensitivity marker size
  magnif2   = 5;  % amplification of delta R PDFs on y-axis

  % -- find future period that is 2C warmer than reference period
  for e = 1:40
   tmp     = rm(tas_cesm_wy_anom(e,:),futwl);
   idx(e)  = find(abs(futtas-tmp)==min(abs(futtas-tmp)))
  end
  % -- calculate simulated P changes
  for e = 1:40
   pr_cesm_change(e)              = nanmean(pr_cesm_wy_percent_mean(e,idx(e)-(futwl/2):idx(e)+(futwl/2)-1),2)-100;
  end

  % -- calculate actually simulated and predicted R changes
  for e = 1:40
    ro_cesm_change(e)                  = nanmean(ro_cesm_wy_percent_mean(e,idx(e)-(futwl/2):idx(e)+(futwl/2)-1),2)-100;
    ro_cesm_constr_change(e)           = tas_ro_reg_cesm_wy(e,1)*futtas + pr_ro_reg_cesm_wy(e,1)*(mean(pr_cesm_wy_percent_mean(e,idx(e)-(futwl/2):idx(e)+(futwl/2)-1),2)-100);
    ro_cesm_constr_change2(e)          = tas_ro_reg_cesm_wy_with_flow_and_carryover(e,1)*futtas + pr_ro_reg_cesm_wy_with_flow_and_carryover(e,1)*(mean(pr_cesm_wy_percent_mean(e,idx(e)-(futwl/2):idx(e)+(futwl/2)-1),2)-100);
  end

  close all
  hold on
  plot(ro_cesm_change,ro_cesm_constr_change,'bo')
  plot(ro_cesm_change,ro_cesm_constr_change2,'ro')

  rmse(ro_cesm_change',ro_cesm_constr_change')
  rmse(ro_cesm_change',ro_cesm_constr_change2')



  return
  set(gcf,'PaperPositionMode','auto');
  fileo = [figpath 'Fig_S8_' vari];
  print('-r300','-loose', '-depsc', [fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  return
end









if figS4 == 1
  % -----------------------------------------------------------------------------
  % Figure S4?
  % ----

  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 20 20])

  subplot(2,1,1)
  title('(a) Upper Colorado River Basin flow volume')
  hold on
  x = [ro_obs_wy; ro_cesm_control_wy'; flow_cesm_control_wy'];
  g = [zeros(length(ro_obs_wy), 1); ones(length(ro_cesm_control_wy), 1); 2*ones(length(flow_cesm_control_wy), 1)];
  boxplot(x,g,'labels',{'Observed streamflow','CESM control runoff','CESM control streamflow'})
  text(.7,100,['r(CESM_{runoff},CESM_{streamflow}) = ' num2str(corr(ro_cesm_control_wy',flow_cesm_control_wy')) ''])
  box on
  ylabel('Runoff or streamflow (km^3/year)')

  subplot(2,1,2)
  title('(b) Upper Colorado River Basin water year autocorrelation')
  hold on
  ii = 30;
  tmp0  = xcov(ro_obs_wy,'coeff');
  tmp00 = xcov(pr_obs_wy,'coeff');
  for i = 1:floor(length(ro_cesm_control_wy)/59)
    tmp1  = xcov(ro_cesm_control_wy((i-1)*59+1:i*59),'coeff');
    tmp2  = xcov(flow_cesm_control_wy((i-1)*59+1:i*59),'coeff');
    h1 = plot(0:ii,tmp1(length(tmp1)/2:length(tmp1)/2+ii),'Color',[.5 .5 1],'MarkerSize',12)
    h2 = plot(0:ii,tmp2(length(tmp2)/2:length(tmp2)/2+ii),'Color',[1 .5 .5],'MarkerSize',12)
  end
  h0  = plot(0:ii,tmp0(length(tmp0)/2:length(tmp0)/2+ii),'k.-','MarkerSize',12)
  box on
  hline(0,'k')
  set(gca,'Layer','top')
  ylabel('Corraltion coefficient')
  xlabel('Time (Years lagged)')
  legend([h0(1) h1(1) h2(1)],'Observed streamflow','CESM control runoff','CESM control streamflow','Location','NorthEast')

  set(gcf,'PaperPositionMode','auto');
  fileo = [figpath 'Fig_S4'];
  print('-r300','-loose', '-depsc', [fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  return
end








if fig4 == 1
% ------------------------------------------------------------------------
% Figure 4? in paper (simpler version of old Fig 4)
% ---------
close all


t = 2; % tas degC anomaly relative to reference period

% futstart  = 2021; % 2001 2050
% futende   = 2070; % 2050 2099
% or...
futtas    = t; % tas degC anomaly relative to reference period
futwl     = 40; % window length over which 'futtas' is defined (in years)

magnif1   = 50; % T sensitivity marker size
magnif2   = 5;  % amplification of delta R PDFs on y-axis

% -- find future period that is 2C warmer than reference period
for m = 1:length(models_common)
  tmp     = rm(tas_h2r_wy_anom(m,:),futwl);
  idx(m)  = find(abs(futtas-tmp)==min(abs(futtas-tmp)))
end
idx'+1850
% -- also find corresponding CO2 concentration (from CESM rcp85 co2vmr file):
for m = 1:length(models_common)
  co2_values(m) = round(co2(time_co2==idx(m)+1850));
end
co2_values'
% -- calculate simulated P changes
for m = 1:length(models_common)
  pr_h2r_change(m)              = nanmean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100;
end
% -- calculate actually simulated and predicted R changes
for m = 1:length(models_common)
  ro_h2r_change(m)                  = nanmean(ro_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100;
  ro_h2r_constr_change(m)           = tas_ro_reg_h2r_wy(m,1)*futtas + pr_ro_reg_h2r_wy(m,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
    % -- incl. interaction term in prediction:
    % ro_h2r_constr_change(m)           = tas_ro_reg_h2r_wy(m,1)*futtas + pr_ro_reg_h2r_wy(m,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100) + ...
    %                                     int_ro_reg_h2r_wy(m,1)*(futtas*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100));
  ro_h2r_constr_change_co2(m)       = tas_ro_reg_h2r_wy_co2(m,1)*futtas + pr_ro_reg_h2r_wy_co2(m,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100) + co2_ro_reg_h2r_wy_co2(m,1)*co2_h2r_anom(idx(m));
  ro_h2r_constr_change_only_pr(m)   = pr_ro_reg_h2r_wy(m,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
  ro_h2r_constr_change_only_tas(m)  = tas_ro_reg_h2r_wy(m,1)*futtas;
  ro_h2r_constr_change_w_uncertainty(m,:) = normrnd(tas_ro_reg_h2r_wy(m,1),abs(tas_ro_reg_h2r_wy(m,2)-tas_ro_reg_h2r_wy(m,1))/1.96,1000,1)*futtas...
   + normrnd(pr_ro_reg_h2r_wy(m,1),abs(pr_ro_reg_h2r_wy(m,2)-pr_ro_reg_h2r_wy(m,1))/1.96,1000,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
  ro_h2r_obs_constr_change(m,1) = tas_ro_reg_obs_wy(1)*futtas + pr_ro_reg_obs_wy(1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
    % % -- using Hoerling et al 2019 temp sens of -2.5%/C for UC:
    % ro_h2r_obs_constr_change(m,1) = -2.5*futtas + 2.0*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
  ro_h2r_obs_constr_change_co2(m,1) = tas_ro_reg_obs_wy_co2(1)*futtas + pr_ro_reg_obs_wy_co2(1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100) + co2_ro_reg_obs_wy_co2(1)*co2_h2r_anom(idx(m));
  % ro_h2r_obs_constr_change_w_uncertainty(m,:) = normrnd(tas_ro_reg_obs_wy(1),abs(tas_ro_reg_obs_wy(3)-tas_ro_reg_obs_wy(2))/3,1000,1)*futtas...
  %  + normrnd(pr_ro_reg_obs_wy(1),abs(pr_ro_reg_obs_wy(3)-pr_ro_reg_obs_wy(2))/3,1000,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
   ro_h2r_obs_constr_change_w_uncertainty(m,:) = tas_ro_reg_obs_wy_mc*futtas + pr_ro_reg_obs_wy_mc*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
end
for m = 1:length(models_clm)
  ii = find(strcmp(models_common, models_clm(m)));
  ro_h2r_change_clm(m) = nanmean(ro_h2r_wy_percent_mean(ii,idx(ii)-(futwl/2):idx(ii)+(futwl/2)-1),2)-100;
end

xlim = [-60 30];

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [10 10 11 20])

subplot(2,1,1)
% title([num2str(t) 'C warming'])
title(vari_name)
hold on
h1 = plot(ro_h2r_change,ro_h2r_constr_change,'ko')%,'MarkerFaceColor','k')
h4 = plot(ro_h2r_change,ro_h2r_constr_change_co2,'mo')%,'MarkerFaceColor','k')
[x,idx] = sort(ro_h2r_change);
plot(ro_h2r_change(idx),polyval(polyfit(ro_h2r_change(idx),ro_h2r_constr_change(idx),1),ro_h2r_change(idx)),'k','LineWidth',1.5)
plot(ro_h2r_change(idx),polyval(polyfit(ro_h2r_change(idx),ro_h2r_constr_change_co2(idx),1),ro_h2r_change(idx)),'m','LineWidth',1.5)
plot([xlim(1) xlim(2)],[xlim(1) xlim(2)],'Color',[.5 .5 .5])
h2 = plot(ro_h2r_change,ro_h2r_constr_change_only_pr,'bo')
h3 = plot(ro_h2r_change,ro_h2r_constr_change_only_tas,'ro')
[x,idx] = sort(ro_h2r_change);
plot(ro_h2r_change(idx),polyval(polyfit(ro_h2r_change(idx),ro_h2r_constr_change_only_pr(idx),1),ro_h2r_change(idx)),'b','LineWidth',1.5)
[x,idx] = sort(ro_h2r_change);
plot(ro_h2r_change(idx),polyval(polyfit(ro_h2r_change(idx),ro_h2r_constr_change_only_tas(idx),1),ro_h2r_change(idx)),'r','LineWidth',1.5)
plot([xlim(1) xlim(2)],[xlim(1) xlim(2)],'Color',[.5 .5 .5])
plot([0 0],[xlim(1) xlim(2)],'Color',[.5 .5 .5])
plot([xlim(1) xlim(2)],[0 0],'Color',[.5 .5 .5])
set(gca,'XLim',xlim,'YLim',xlim)
box on
xlabel('Simulated \DeltaQ (%)')
ylabel('Predicted \DeltaQ (%)')
% -- correlations:
[r,p] = corr(ro_h2r_change',ro_h2r_constr_change');
n = 80;
k = 4;
r2_adj  = 1 - (((1-r^2)*(n-1)) / (n-k-1));
mse     = rmse(ro_h2r_change,ro_h2r_constr_change)^2;
aic     = log(mse)+(k+1)*2/n;
bic     = log(mse)+(k+1)*log(n)/n;
if p < 0.05
  signi = '*';
else
  signi = '';
end
text(5,-22,['R^2'],'Color',[.4 .4 .4])
text(5,-28,[num2str(round(r^2*100)/100) signi ],'Color','k')
% text(5,-28,[num2str(round(r2_adj*100)/100) signi ],'Color','k')
[r,p] = corr(ro_h2r_change',ro_h2r_constr_change_co2');
k = 4;%5;
r2_adj  = 1 - (((1-r^2)*(n-1)) / (n-k-1));
mse     = rmse(ro_h2r_change,ro_h2r_constr_change)^2;
aic     = log(mse)+(k+1)*2/n;
bic     = log(mse)+(k+1)*log(n)/n;
if p < 0.05
  signi = '*';
else
  signi = '';
end
text(5,-33,[num2str(round(r^2*100)/100) signi],'Color','m')
% text(5,-33,[num2str(round(r2_adj*100)/100) signi],'Color','m')
[r,p] = corr(ro_h2r_change',ro_h2r_constr_change_only_pr');
k = 2;
r2_adj  = 1 - (((1-r^2)*(n-1)) / (n-k-1));
if p < 0.05
  signi = '*';
else
  signi = '';
end
text(5,-38,[num2str(round(r^2*100)/100) signi],'Color','b')
% text(5,-38,[num2str(round(r2_adj*100)/100) signi],'Color','b')
[r,p] = corr(ro_h2r_change',ro_h2r_constr_change_only_tas');
k = 2;
r2_adj  = 1 - (((1-r^2)*(n-1)) / (n-k-1));
if p < 0.05
  signi = '*';
else
  signi = '';
end
text(5,-43,[num2str(round(r^2*100)/100) signi],'Color','r')
% text(5,-43,[num2str(round(r2_adj*100)/100) signi],'Color','r')
if strcmp(vari,'UC')==1
  legend([h1 h4 h2 h3],'\DeltaQ','\DeltaQ-[CO2]','P-related \DeltaQ','T-related \DeltaQ','Location','NorthWest')
  legend boxoff
  set(0,'DefaultLegendAutoUpdate','off')
end

% ---
if strcmp(vari,'UC')==1
  offset1 = 0.45+.1-.08;
  offset2 = 0.45-.08;
elseif strcmp(vari,'CO_DALLES')==1
  offset1 = 0.94+.1-.08;
  offset2 = 0.94-.08;
elseif strcmp(vari,'NS')==1
  offset1 = 0.78+.1-.08;
  offset2 = 0.78-.08;
end
ylim = [0 1.1];

subplot(2,1,2)
hold on
% bw = 2;
[f0,x0,u] = ksdensity( ro_h2r_change)%,'width',bw );
h1 = line(x0, f0*magnif2, 'Color','k', 'LineWidth',2);
[f0,x0,u] = ksdensity( ro_h2r_constr_change)%,'width',bw );
h2 = line(x0, f0*magnif2, 'Color','k', 'LineWidth',1,'LineStyle','--');
[f0,x0,u] = ksdensity( ro_h2r_constr_change_only_pr(:))%,'width',bw );
h4 = line(x0, f0*magnif2, 'Color','b', 'LineWidth',1,'LineStyle','-');
[f0,x0,u] = ksdensity( ro_h2r_constr_change_only_tas(:))%,'width',bw );
h5 = line(x0, f0*magnif2, 'Color','r', 'LineWidth',1,'LineStyle','-');
[f0,x0,u] = ksdensity( ro_h2r_obs_constr_change(:,1))%,'width',bw );
h3 = line(x0, f0*magnif2, 'Color',[6 160 39]/255, 'LineWidth',2,'LineStyle','-');

[f0,x0,u] = ksdensity( ro_h2r_obs_constr_change_w_uncertainty(:) );
line(x0, f0*magnif2, 'Color',[6 160 39]/255, 'LineWidth',1.5,'LineStyle','--');
plot([0 0],[ylim(1) ylim(2)],'Color',[.5 .5 .5])
set(gca,'XLim',xlim,'YLim',[ylim(1) ylim(2)],'YTick',[],'YTickLabel',[])
box on
if strcmp(vari,'UC')==1
  legend([h1 h2 h4 h5 h3],...
  ['Simulated \DeltaQ'],...
  ['Predicted \DeltaQ'],...
  ['P-related \DeltaQ'],...
  ['T-related \DeltaQ'],...
  ['Obs-constrained \DeltaQ'],...
  'Location','NorthWest')
  legend boxoff
end
xlabel('\DeltaQ (%)')
ylabel('Density')
plot([min(ro_h2r_change) max(ro_h2r_change)],[offset1 offset1],'Color','k','LineWidth',2)
plot([mean(ro_h2r_change) mean(ro_h2r_change)],[offset1 offset1],'.','Color','k','LineWidth',2,'MarkerSize',18)
% plot([median(ro_h2r_change) median(ro_h2r_change)],[offset1-.02 offset1+.02],'-','Color','k','LineWidth',2)

plot([min(ro_h2r_obs_constr_change) max(ro_h2r_obs_constr_change)],[offset2 offset2],'Color',[6 160 39]/255,'LineWidth',2)
plot([mean(ro_h2r_obs_constr_change) mean(ro_h2r_obs_constr_change)],[offset2 offset2],'.','Color',[6 160 39]/255,'LineWidth',2,'MarkerSize',18)
% plot([median(ro_h2r_obs_constr_change) median(ro_h2r_obs_constr_change)],[offset2-.02 offset2+.02],'-','Color',[6 160 39]/255,'LineWidth',2)

% plot([min(ro_h2r_obs_constr_change_co2) max(ro_h2r_obs_constr_change_co2)],[offset-.05 offset-.05],'Color','m','LineWidth',2)
% plot([mean(ro_h2r_obs_constr_change_co2) mean(ro_h2r_obs_constr_change_co2)],[offset-.05 offset-.05],'.','Color','m','LineWidth',2,'MarkerSize',18)
% plot([median(ro_h2r_obs_constr_change_co2) median(ro_h2r_obs_constr_change_co2)],[offset-.05-.02 offset-.05+.02],'-','Color','m','LineWidth',2)

plot([prctile(ro_h2r_obs_constr_change_w_uncertainty(:),5) prctile(ro_h2r_obs_constr_change_w_uncertainty(:),95)],[offset2 offset2],'-','Color',[6 160 39]/255,'LineWidth',1)
% plot([mean(ro_h2r_obs_constr_change_w_uncertainty(:)) mean(ro_h2r_obs_constr_change_w_uncertainty(:))],[offset2 offset2],'.','Color',[6 160 39]/255,'LineWidth',2,'MarkerSize',18)
% plot([median(ro_h2r_obs_constr_change_w_uncertainty(:)) median(ro_h2r_obs_constr_change_w_uncertainty(:))],[offset2-.02 offset2+.02],'-','Color',[6 160 39]/255,'LineWidth',2)

% -- text labels for range bars:
% -- models
text(min(ro_h2r_change),offset1+.05,[num2str(round(mean(ro_h2r_change)*10)/10) ' ('...
num2str(round(((max(ro_h2r_change)+100)-(min(ro_h2r_change)+100))*10)/10) ')'],'Color','k')
% % -- obs-constraint w/o uncertainty
% text(min(ro_h2r_obs_constr_change),offset2+.05,[num2str(round(mean(ro_h2r_obs_constr_change)*10)/10) '  '...
% num2str(round(((max(ro_h2r_obs_constr_change)+100)-(min(ro_h2r_obs_constr_change)+100))*10)/10)],'Color',[6 160 39]/255)
% -- obs-constraint w/ uncertainty
text(prctile(ro_h2r_obs_constr_change_w_uncertainty(:),5),offset2+.05,[num2str(round(mean(ro_h2r_obs_constr_change(:))*10)/10) ' ('...
num2str(round(((max(ro_h2r_obs_constr_change)+100)-(min(ro_h2r_obs_constr_change)+100))*10)/10) ' | '...
num2str(round(((prctile(ro_h2r_obs_constr_change_w_uncertainty(:),95)+100)-(prctile(ro_h2r_obs_constr_change_w_uncertainty(:),5)+100))*10)/10) ')'],'Color',[6 160 39]/255)


% return
set(gcf,'PaperPositionMode','auto');
fileo = [figpath 'Fig_4_' vari co_tag];
print('-r300','-loose', '-depsc', [fileo '.eps'])
saveas(gcf,fileo,'jpg')
return

end












if fig5 == 1
% ------------------------------------------------------------------------
% Figure 5? in paper (simpler version of Fig 4, but with several temp changes)
% ---------
close all

xlim = [-60 30];

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [10 10 24 35])

for t = 1:2

  futstart  = 2021; % 2001 2050
  futende   = 2070; % 2050 2099
  % or...
  futtas    = t; % tas degC anomaly relative to reference period 1929-2008
  futwl     = 40; % window length over which 'futtas' is defined (in years)

  magnif1   = 50; % T sensitivity marker size
  magnif2   = 5;  % amplification of delta R PDFs on y-axis

  for m = 1:length(models_common)
    tmp     = rm(tas_h2r_wy_anom(m,:),futwl);
    idx(m)  = find(abs(futtas-tmp)==min(abs(futtas-tmp)));
  end
  for m = 1:length(models_common)
    ro_h2r_change(m)              = nanmean(ro_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100;
    ro_h2r_constr_change(m)       = tas_ro_reg_h2r_wy(m,1)*futtas + pr_ro_reg_h2r_wy(m,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
    ro_h2r_constr_change_only_pr(m)   = pr_ro_reg_h2r_wy(m,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
    ro_h2r_constr_change_only_tas(m)  = tas_ro_reg_h2r_wy(m,1)*futtas;
    ro_h2r_constr_change_w_uncertainty(m,:) = normrnd(tas_ro_reg_h2r_wy(m,1),abs(tas_ro_reg_h2r_wy(m,2)-tas_ro_reg_h2r_wy(m,1))/1.96,1000,1)*futtas...
     + normrnd(pr_ro_reg_h2r_wy(m,1),abs(pr_ro_reg_h2r_wy(m,2)-pr_ro_reg_h2r_wy(m,1))/1.96,1000,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
    ro_h2r_obs_constr_change(m,1) = tas_ro_reg_obs_wy(1)*futtas + pr_ro_reg_obs_wy(1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
    ro_h2r_obs_constr_change_w_uncertainty(m,:) = normrnd(tas_ro_reg_obs_wy(1),abs(tas_ro_reg_obs_wy(2)-tas_ro_reg_obs_wy(1))/1.96,1000,1)*futtas...
     + normrnd(pr_ro_reg_obs_wy(1),abs(pr_ro_reg_obs_wy(2)-pr_ro_reg_obs_wy(1))/1.96,1000,1)*(mean(pr_h2r_wy_percent_mean(m,idx(m)-(futwl/2):idx(m)+(futwl/2)-1),2)-100);
  end
  for m = 1:length(models_clm)
    ii = find(strcmp(models_common, models_clm(m)));
    ro_h2r_change_clm(m) = nanmean(ro_h2r_wy_percent_mean(ii,idx(ii)-(futwl/2):idx(ii)+(futwl/2)-1),2)-100;
  end

  subplot(3,2,1+t-1)
  title([num2str(t) 'C warming'])
  hold on
  h1 = plot(ro_h2r_constr_change,ro_h2r_change,'ko')%,'MarkerFaceColor','k')
  [x,idx] = sort(ro_h2r_constr_change);
  plot(ro_h2r_constr_change(idx),polyval(polyfit(ro_h2r_constr_change(idx),ro_h2r_change(idx),1),ro_h2r_constr_change(idx)),'k','LineWidth',1.5)
  plot([xlim(1) xlim(2)],[xlim(1) xlim(2)],'Color',[.5 .5 .5])
  corr(ro_h2r_change',ro_h2r_constr_change')
  plot([0 0],[xlim(1) xlim(2)],'Color',[.5 .5 .5])
  plot([xlim(1) xlim(2)],[0 0],'Color',[.5 .5 .5])
  set(gca,'XLim',xlim,'YLim',xlim)
  box on
  xlabel('Predicted R change (%)')
  ylabel('Simulated R change (%)')
  text(10,-30,['r = ' num2str(round(corr(ro_h2r_change',ro_h2r_constr_change')*100)/100)],'Color','k')
  legend([h1],'R change','Location','NorthWest')
  legend boxoff

  subplot(3,2,3+t-1)
  hold on
  h1 = plot(ro_h2r_constr_change_only_pr,ro_h2r_change,'bo')
  h2 = plot(ro_h2r_constr_change_only_tas,ro_h2r_change,'ro')
  [x,idx] = sort(ro_h2r_constr_change_only_pr);
  plot(ro_h2r_constr_change_only_pr(idx),polyval(polyfit(ro_h2r_constr_change_only_pr(idx),ro_h2r_change(idx),1),ro_h2r_constr_change_only_pr(idx)),'b','LineWidth',1.5)
  [x,idx] = sort(ro_h2r_constr_change_only_tas);
  plot(ro_h2r_constr_change_only_tas(idx),polyval(polyfit(ro_h2r_constr_change_only_tas(idx),ro_h2r_change(idx),1),ro_h2r_constr_change_only_tas(idx)),'r','LineWidth',1.5)
  plot([xlim(1) xlim(2)],[xlim(1) xlim(2)],'Color',[.5 .5 .5])
  plot([0 0],[xlim(1) xlim(2)],'Color',[.5 .5 .5])
  plot([xlim(1) xlim(2)],[0 0],'Color',[.5 .5 .5])
  set(gca,'XLim',xlim,'YLim',xlim)
  box on
  xlabel('Predicted R change (%)')
  ylabel('Simulated R change (%)')
  text(10,-30,['r = ' num2str(round(corr(ro_h2r_change',ro_h2r_constr_change_only_pr')*100)/100)],'Color','b')
  text(10,-35,['r = ' num2str(round(corr(ro_h2r_change',ro_h2r_constr_change_only_tas')*100)/100)],'Color','r')
  legend([h1 h2],'P-related R change','T-related R change','Location','NorthWest')
  legend boxoff

  subplot(3,2,5+t-1)
  hold on
  [f0,x0,u] = ksdensity( ro_h2r_change );
  h1 = line(x0, f0*magnif2, 'Color','k', 'LineWidth',1.5);
  [f0,x0,u] = ksdensity( ro_h2r_constr_change );
  h2 = line(x0, f0*magnif2, 'Color','k', 'LineWidth',1,'LineStyle','--');
  % [f0,x0,u] = ksdensity( ro_h2r_constr_change_w_uncertainty(:) );
  % line(x0, f0*magnif2, 'Color','b', 'LineWidth',1,'LineStyle',':');
  [f0,x0,u] = ksdensity( ro_h2r_constr_change_only_pr(:) );
  h4 = line(x0, f0*magnif2, 'Color','b', 'LineWidth',1,'LineStyle','-');
  [f0,x0,u] = ksdensity( ro_h2r_constr_change_only_tas(:) );
  h5 = line(x0, f0*magnif2, 'Color','r', 'LineWidth',1,'LineStyle','-');
  [f0,x0,u] = ksdensity( ro_h2r_obs_constr_change(:,1) );
  h3 = line(x0, f0*magnif2, 'Color',[6 160 39]/255, 'LineWidth',1.5,'LineStyle','-');
  % [f0,x0,u] = ksdensity( ro_h2r_change_clm(:) );
  % h0 = line(x0, f0*magnif2, 'Color','m', 'LineWidth',1.5,'LineStyle','-');
  % [f0,x0,u] = ksdensity( ro_h2r_obs_constr_change_w_uncertainty(:) );
  % line(x0, f0*magnif2, 'Color','m', 'LineWidth',1,'LineStyle','-');
  plot([0 0],[0 1],'Color',[.5 .5 .5])
  % plot([median(ro_h2r_obs_constr_change) median(ro_h2r_obs_constr_change)],[0 1],'Color',[6 160 39]/255,'LineWidth',1.5)
  % plot([median(ro_h2r_change) median(ro_h2r_change)],[0 1],'k','LineWidth',1.5)
  % plot([median(ro_h2r_constr_change) median(ro_h2r_constr_change)],[0 1],'k--')
  % plot([median(ro_h2r_constr_change_only_pr) median(ro_h2r_constr_change_only_pr)],[0 1],'b')
  % plot([median(ro_h2r_constr_change_only_tas) median(ro_h2r_constr_change_only_tas)],[0 1],'r')
  % text(median(ro_h2r_change),.61,[num2str(round(median(ro_h2r_change)*10)/10)])
  % text(median(ro_h2r_constr_change),.63,[num2str(round(median(ro_h2r_constr_change)*10)/10)])
  % text(median(ro_h2r_obs_constr_change),.65,[num2str(round(median(ro_h2r_obs_constr_change)*10)/10)],'Color',[6 160 39]/255)
  set(gca,'XLim',xlim,'YLim',[0 1])
  box on
  legend([h1 h2 h3 h4 h5],...
  ['Simulated R change (' num2str(round(median(ro_h2r_change)*10)/10) '%)'],...
  ['Predicted R change (' num2str(round(median(ro_h2r_constr_change)*10)/10) '%)'],...
  ['Obs-predicted R change (' num2str(round(median(ro_h2r_obs_constr_change)*10)/10) '%)'],...
  ['P-related R change (' num2str(round(median(ro_h2r_constr_change_only_pr)*10)/10) '%)'],...
  ['T-related R change (' num2str(round(median(ro_h2r_constr_change_only_tas)*10)/10) '%)'],...
  'Location','NorthWest')
  % legend boxoff
  xlabel('Predicted/Simulated R change (%)')
  ylabel('Density')
end


% return
set(gcf,'PaperPositionMode','auto');
fileo = [figpath 'Fig_5_' vari];
print('-r300','-loose', '-depsc', [fileo '.eps'])
saveas(gcf,fileo,'jpg')
% return
end








if figS1 == 1
% -----------------------------------------------------------------------------
% Figure S1?
% ----

% -- methodological sensitivity tests:
% temporal variability of Sankar's epsilon --

% -- Precip data sets --
% -- obs data 1:
obs1_data_name = 'PRISM';
pr_obs1 	= ncread(['' pathin 'prism/precip/prism_precip_189501-201712_uc_ts.nc'],[vari '_TS']); % in MAF
for i = 2:length(pr_obs1)/12
  pr_obs1_wy(i-1) = sum(pr_obs1((i-2)*12+10:(i-1)*12+9));
end
time_pr_obs1	= 1896:2017;
% -- obs data 2:
obs2_data_name = 'GPCC';
pr_obs2 	= ncread(['' pathin 'observations/precipitation/gpcc/gpcc_precip_colo_colu_ts.nc'],[vari '_TS']); % in MAF
for i = 2:length(pr_obs2)/12
  pr_obs2_wy(i-1) = sum(pr_obs2((i-2)*12+10:(i-1)*12+9));
end
time_pr_obs2	= 1902:2016;
% -- obs data 3:
obs3_data_name = 'Livneh';
pr_obs3 	= ncread(['' pathin 'livneh/livneh_precip_colo_colu_ts.nc'],[vari '_TS']); % in MAF
for i = 2:length(pr_obs3)/12
  pr_obs3_wy(i-1) = sum(pr_obs3((i-2)*12+10:(i-1)*12+9));
end
time_pr_obs3	= 1916:2011;
% -- cut to common length
pr_obs1_wy      = pr_obs1_wy(find(time_pr_obs1==time_obs_common(1)):find(time_pr_obs1==time_obs_common(end)));
pr_obs2_wy      = pr_obs2_wy(find(time_pr_obs2==time_obs_common(1)):find(time_pr_obs2==time_obs_common(end)));
pr_obs3_wy      = pr_obs3_wy(find(time_pr_obs3==time_obs_common(1)):find(time_pr_obs3==time_obs_common(end)));
pr_obs1_wy_percent_mean = (pr_obs1_wy./mean(pr_obs1_wy))*100;
pr_obs2_wy_percent_mean = (pr_obs2_wy./mean(pr_obs2_wy))*100;
pr_obs3_wy_percent_mean = (pr_obs3_wy./mean(pr_obs3_wy))*100;

% -- Temp data sets --
% -- obs data 1:
clear('tas_obs1_wy')
obs1_data_name = 'Livneh';
tas_obs1 	= ncread(['' pathin 'livneh/livneh_tmean_colo_colu_ts.nc'],[vari '_TS']);
for i = 2:length(tas_obs1)/12
  tas_obs1_wy(i-1) = mean(tas_obs1((i-2)*12+10:(i-1)*12+9));
end
time_tas_obs1	= 1916:2011;
tas_obs1_wy_anom = tas_obs1_wy-mean(tas_obs1_wy(refstart-time_tas_obs1(1)+1:refende-time_tas_obs1(1)));
% -- obs data 2:
tas_obs2 	= ncread(['' pathin 'observations/temperature/best/best_temperature_185001-201712_colo_colu_ts.nc'],[vari '_TS']);
tas_obs2_am = am(tas_obs2);
tas_obs2_am = tas_obs2_am(2:end); % adjust to same length as WY
for i = 2:length(tas_obs2)/12
  tas_obs2_wy(i-1) = nanmean(tas_obs2((i-2)*12+10:(i-1)*12+9));
end
time_tas_obs2	= 1851:2017;
tas_obs2_wy_anom = tas_obs2_wy-mean(tas_obs2_wy(refstart-time_tas_obs2(1)+1:refende-time_tas_obs2(1)));
% -- obs data 3:
tas_obs3 	= ncread(['' pathin 'prism/tmean/prism_tmean_189501-201712_uc_ts.nc'],[vari '_TS']);
tas_obs3_am = am(tas_obs3);
tas_obs3_am = tas_obs3_am(2:end); % adjust to same length as WY
for i = 2:length(tas_obs3)/12
  tas_obs3_wy(i-1) = nanmean(tas_obs3((i-2)*12+10:(i-1)*12+9));
end
time_tas_obs3	= 1896:2017;
tas_obs3_wy_anom = tas_obs3_wy-mean(tas_obs3_wy(refstart-time_tas_obs3(1)+1:refende-time_tas_obs3(1)));

% -- cut to common length
tas_obs1_wy_anom      = tas_obs1_wy_anom(find(time_tas_obs1==time_obs_common(1)):find(time_tas_obs1==time_obs_common(end)));
tas_obs2_wy_anom      = tas_obs2_wy_anom(find(time_tas_obs2==time_obs_common(1)):find(time_tas_obs2==time_obs_common(end)));
tas_obs3_wy_anom      = tas_obs3_wy_anom(find(time_tas_obs3==time_obs_common(1)):find(time_tas_obs3==time_obs_common(end)));


% -- load streamflow and precip data that Julie Vano used:
data = [1418743.917	1.045659
955518.4167	0.897887
453309.5	0.732259
1260119	1.05651
1505322.417	1.170061
1476323.167	1.208161
748600.3333	1.058718
1421620.583	1.142557
2032799.583	1.242915
2109522.75	1.377675
1754051.833	1.202608
1901946.833	1.273484
1355472.917	1.068602
982124.4167	0.984579
836128.6667	0.881206
790779.8333	0.971848
1006257.167	1.015006
919331.9167	1.086383
1526693	1.318106
885648.75	0.920851
1688850.417	1.294706
1204230.917	1.073559
1801807	1.438768
1401438.167	1.146939
1380431	1.09642
909142.8333	0.958439
910805.4167	0.896307
517506.75	0.803036
876463.1667	0.886962
802812.1667	1.06447];
ro_vano_percent_mean = (data(:,1)./mean(data(:,1)))*100;
pr_vano_percent_mean = data(:,2)*100;
% ro_vano_percent_mean = data(:,1);
% pr_vano_percent_mean = data(:,2);
elast_vano = elast_sankar(ro_vano_percent_mean,pr_vano_percent_mean);
[t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_vano_percent_mean-100,pr_vano_percent_mean-100,tas_obs1_wy_anom(1975-1950+1:2004-1950+1));
reg_vano = t1;
[t1 t2 t3 t4 t5 t6] = runoff_sens_reg_with_carryover(ro_vano_percent_mean-100,pr_vano_percent_mean-100,tas_obs1_wy_anom(1975-1950+1:2004-1950+1));
reg_co_vano = t1;

wl      = 30; % window length
% -- create arrays for different calculation methods (just for obs here)
elast1  = NaN(size(pr_obs1_wy));
elast2  = elast1;
elast3  = elast1;
reg1    = elast1;
reg2    = elast1;
reg3    = elast1;
reg_co1     = elast1;
reg_co2     = elast1;
reg_co3     = elast1;
tas_reg1    = elast1;
tas_reg2    = elast1;
tas_reg3    = elast1;
tas_reg_co1    = elast1;
tas_reg_co2    = elast1;
tas_reg_co3    = elast1;
year1   = 1955;
year2   = 1990;
% -- obs (all different calculation methods) --
for i = 1:length(pr_obs1_wy)-wl
  elast1(i+floor(wl/2)) = elast_sankar(ro_obs_wy(i:i+wl-1),pr_obs1_wy(i:i+wl-1));
  elast2(i+floor(wl/2)) = elast_sankar(ro_obs_wy(i:i+wl-1),pr_obs2_wy(i:i+wl-1));
  elast3(i+floor(wl/2)) = elast_sankar(ro_obs_wy(i:i+wl-1),pr_obs3_wy(i:i+wl-1));
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_obs_wy_percent_mean(i:i+wl-1)-100,pr_obs1_wy_percent_mean(i:i+wl-1)-100,tas_obs1_wy_anom(i:i+wl-1));
  reg1(i+floor(wl/2)) = t1;
  tas_reg1(i+floor(wl/2)) = t3;
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_obs_wy_percent_mean(i:i+wl-1)-100,pr_obs2_wy_percent_mean(i:i+wl-1)-100,tas_obs1_wy_anom(i:i+wl-1));
  reg2(i+floor(wl/2)) = t1;
  tas_reg2(i+floor(wl/2)) = t3;
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_obs_wy_percent_mean(i:i+wl-1)-100,pr_obs3_wy_percent_mean(i:i+wl-1)-100,tas_obs1_wy_anom(i:i+wl-1));
  reg3(i+floor(wl/2)) = t1;
  tas_reg3(i+floor(wl/2)) = t3;
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg_with_carryover(ro_obs_wy_percent_mean(i:i+wl-1)-100,pr_obs1_wy_percent_mean(i:i+wl-1)-100,tas_obs1_wy_anom(i:i+wl-1));
  reg_co1(i+floor(wl/2)) = t1;
  tas_reg_co1(i+floor(wl/2)) = t3;
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg_with_carryover(ro_obs_wy_percent_mean(i:i+wl-1)-100,pr_obs2_wy_percent_mean(i:i+wl-1)-100,tas_obs1_wy_anom(i:i+wl-1));
  reg_co2(i+floor(wl/2)) = t1;
  tas_reg_co2(i+floor(wl/2)) = t3;
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg_with_carryover(ro_obs_wy_percent_mean(i:i+wl-1)-100,pr_obs3_wy_percent_mean(i:i+wl-1)-100,tas_obs1_wy_anom(i:i+wl-1));
  reg_co3(i+floor(wl/2)) = t1;
  tas_reg_co3(i+floor(wl/2)) = t3;
end
% -- bootstrap ("bs") obs:
bl = 2; % block length in years
iterations = 100;
elast1_bs  = NaN(iterations,length(ro_obs_wy_percent_mean));
elast2_bs  = NaN(iterations,length(ro_obs_wy_percent_mean));
reg1_bs    = NaN(iterations,length(ro_obs_wy_percent_mean));
reg2_bs    = NaN(iterations,length(ro_obs_wy_percent_mean));
for iter = 1:iterations
  idx = randi([1 length(ro_obs_wy)-bl+1],floor(length(ro_obs_wy)/bl),1);
  j = 1;
  for i = 1:length(idx)
    ro_tmp1(j:j+bl-1)   = ro_obs_wy(idx(i):idx(i)+bl-1);
    pr_tmp1(j:j+bl-1)   = pr_obs_wy(idx(i):idx(i)+bl-1);
    pr2_tmp1(j:j+bl-1)  = pr_obs2_wy(idx(i):idx(i)+bl-1);
    ro_tmp2(j:j+bl-1)   = ro_obs_wy_percent_mean(idx(i):idx(i)+bl-1)-100;
    pr_tmp2(j:j+bl-1)   = pr_obs_wy_percent_mean(idx(i):idx(i)+bl-1)-100;
    pr2_tmp2(j:j+bl-1)  = pr_obs2_wy_percent_mean(idx(i):idx(i)+bl-1)-100;
    tas_tmp2(j:j+bl-1)  = tas_obs_wy_anom(idx(i):idx(i)+bl-1);
    j = j+bl;
  end
  for i = 1:length(ro_obs_wy)-wl
    elast1_bs(iter,i+floor(wl/2)) = elast_sankar(ro_tmp1(i:i+wl-1),pr_tmp1(i:i+wl-1));
    elast2_bs(iter,i+floor(wl/2)) = elast_sankar(ro_tmp1(i:i+wl-1),pr2_tmp1(i:i+wl-1));
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_tmp2(i:i+wl-1),pr_tmp2(i:i+wl-1),tas_tmp2(i:i+wl-1));
    reg1_bs(iter,i+floor(wl/2)) = t1;
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_tmp2(i:i+wl-1),pr2_tmp2(i:i+wl-1),tas_tmp2(i:i+wl-1));
    reg2_bs(iter,i+floor(wl/2)) = t1;
  end
end
elast1_bs(elast1_bs==0) = NaN;
elast2_bs(elast2_bs==0) = NaN;
elast1_bs_pp = prctile(elast1_bs,[95 5]);
elast2_bs_pp = prctile(elast2_bs,[95 5]);
reg1_bs(reg1_bs==0) = NaN;
reg2_bs(reg2_bs==0) = NaN;
reg1_bs_pp = prctile(reg1_bs,[95 5]);
reg2_bs_pp = prctile(reg2_bs,[95 5]);

% -- cesm piControl --
elast_cesm_control    = NaN(size(ro_cesm_control_wy));
reg_cesm_control      = NaN(size(ro_cesm_control_wy));
tas_reg_cesm_control  = NaN(size(ro_cesm_control_wy));
reg_co_cesm_control   = NaN(size(ro_cesm_control_wy));
tas_reg_co_cesm_control  = NaN(size(ro_cesm_control_wy));
for i = 1:length(ro_cesm_control_wy)-wl
  elast_cesm_control(i+floor(wl/2)) = elast_sankar(ro_cesm_control_wy_percent_mean(i:i+wl-1),pr_cesm_control_wy_percent_mean(i:i+wl-1));
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_cesm_control_wy_percent_mean(i:i+wl-1)-100,pr_cesm_control_wy_percent_mean(i:i+wl-1)-100,tas_cesm_control_wy_anom(i:i+wl-1));
  reg_cesm_control(i+floor(wl/2)) = t1;
  tas_reg_cesm_control(i+floor(wl/2)) = t3;
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg_with_carryover(ro_cesm_control_wy_percent_mean(i:i+wl-1)-100,pr_cesm_control_wy_percent_mean(i:i+wl-1)-100,tas_cesm_control_wy_anom(i:i+wl-1));
  reg_co_cesm_control(i+floor(wl/2)) = t1;
  tas_co_reg_cesm_control(i+floor(wl/2)) = t3;
end
% -- cesm piControl segments --
segments1 = [10 15 20 25 30 35 40 45 50 60 70 80 90 100 150 200 300 400 500 1000];
elast_cesm_control_seg    = NaN(length(segments1),floor(1800/segments1(1)));
reg_cesm_control_seg      = NaN(length(segments1),floor(1800/segments1(1)));
tas_reg_cesm_control_seg  = NaN(length(segments1),floor(1800/segments1(1)));
reg_co_cesm_control_seg      = NaN(length(segments1),floor(1800/segments1(1)));
tas_reg_co_cesm_control_seg  = NaN(length(segments1),floor(1800/segments1(1)));
xx = 0;
for n = 1:length(segments1) % window length nn
  nn = segments1(n);
  xx = xx+1;
  for i = 1:floor(1800/nn) % no of segment of window length nn
    % -- ro vs pr
    elast_cesm_control_seg(xx,i) = elast_sankar(ro_cesm_control_wy_percent_mean((i-1)*nn+1:i*nn),pr_cesm_control_wy_percent_mean((i-1)*nn+1:i*nn));
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_cesm_control_wy_percent_mean((i-1)*nn+1:i*nn)'-100,pr_cesm_control_wy_percent_mean((i-1)*nn+1:i*nn)'-100,tas_cesm_control_wy_anom((i-1)*nn+1:i*nn));
    reg_cesm_control_seg(xx,i) = t1;
    tas_reg_cesm_control_seg(xx,i) = t3;
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg_with_carryover(ro_cesm_control_wy_percent_mean((i-1)*nn+1:i*nn)'-100,pr_cesm_control_wy_percent_mean((i-1)*nn+1:i*nn)'-100,tas_cesm_control_wy_anom((i-1)*nn+1:i*nn));
    reg_co_cesm_control_seg(xx,i) = t1;
    tas_reg_co_cesm_control_seg(xx,i) = t3;
  end
end

% -- cesm le segments --
segments2 = [10 15 20 25 30 35 40 45 50 60 70 80 90 100];
% length of historical2rcp85 to be considered:
ll = 100;
elast_cesm_seg    = NaN(40,floor(ll/segments2(1)),length(segments2));
reg_cesm_seg      = elast_cesm_seg;
tas_reg_cesm_seg  = elast_cesm_seg;
reg_co_cesm_seg      = elast_cesm_seg;
tas_reg_co_cesm_seg  = elast_cesm_seg;
for e = 1:40
  xx = 0;
  for n = 1:length(segments2) % window length nn
    nn = segments2(n);
    xx = xx+1;
    for i = 1:floor(ll/nn) % no of segment of window length nn
      % -- ro vs pr
      elast_cesm_seg(e,i,xx) = elast_sankar(ro_cesm_wy_percent_mean(e,(i-1)*nn+1:i*nn),pr_cesm_wy_percent_mean(e,(i-1)*nn+1:i*nn));
      [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_cesm_wy_percent_mean(e,(i-1)*nn+1:i*nn)'-100,pr_cesm_wy_percent_mean(e,(i-1)*nn+1:i*nn)'-100,tas_cesm_wy_anom(e,(i-1)*nn+1:i*nn));
      reg_cesm_seg(e,i,xx) = t1;
      tas_reg_cesm_seg(e,i,xx) = t3;
      [t1 t2 t3 t4 t5 t6] = runoff_sens_reg_with_carryover(ro_cesm_wy_percent_mean(e,(i-1)*nn+1:i*nn)'-100,pr_cesm_wy_percent_mean(e,(i-1)*nn+1:i*nn)'-100,tas_cesm_wy_anom(e,(i-1)*nn+1:i*nn));
      reg_co_cesm_seg(e,i,xx) = t1;
      tas_reg_co_cesm_seg(e,i,xx) = t3;
    end
  end
end

% -- cesm le --
elast_cesm      = NaN(size(ro_obs_wy_percent_mean));
reg_cesm        = NaN(size(ro_obs_wy_percent_mean));
reg_co_cesm     = NaN(size(ro_obs_wy_percent_mean));
elast_cesm_bs   = NaN(iterations,length(ro_obs_wy_percent_mean));
reg_cesm_bs     = NaN(iterations,length(ro_obs_wy_percent_mean));
reg_co_cesm_bs  = NaN(iterations,length(ro_obs_wy_percent_mean));
for e = 1:40
  for i = 1:length(ro_obs_wy_percent_mean)-wl
    elast_cesm(e,i+floor(wl/2)) = elast_sankar(ro_cesm_wy(e,8+i:8+i+wl-1),pr_cesm_wy(e,8+i:8+i+wl-1));
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_cesm_wy_percent_mean(e,8+i:8+i+wl-1)-100,pr_cesm_wy_percent_mean(e,8+i:8+i+wl-1)-100,tas_cesm_wy_anom(e,8+i:8+i+wl-1));
    reg_cesm(e,i+floor(wl/2)) = t1;
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg_with_carryover(ro_cesm_wy_percent_mean(e,8+i:8+i+wl-1)-100,pr_cesm_wy_percent_mean(e,8+i:8+i+wl-1)-100,tas_cesm_wy_anom(e,8+i:8+i+wl-1));
    reg_co_cesm(e,i+floor(wl/2)) = t1;
  end
end
% -- bootstrap cesm le #e:
e = 2; % ensemble member
for iter = 1:iterations
  idx = randi([1 length(ro_obs_wy)-bl+1],floor(length(ro_obs_wy)/bl),1);
  j = 1;
  for i = 1:length(idx)
    ro_tmp1(j:j+bl-1)   = ro_cesm_wy(e,8+idx(i):8+idx(i)+bl-1);
    pr_tmp1(j:j+bl-1)   = pr_cesm_wy(e,8+idx(i):8+idx(i)+bl-1);
    ro_tmp2(j:j+bl-1)   = ro_cesm_wy_percent_mean(e,8+idx(i):8+idx(i)+bl-1)-100;
    pr_tmp2(j:j+bl-1)   = pr_cesm_wy_percent_mean(e,8+idx(i):8+idx(i)+bl-1)-100;
    tas_tmp2(j:j+bl-1)  = tas_cesm_wy_anom(e,8+idx(i):8+idx(i)+bl-1);
    j = j+bl;
  end
  for i = 1:length(ro_obs_wy)-wl
    elast_cesm_bs(iter,i+floor(wl/2)) = elast_sankar(ro_tmp1(i:i+wl-1),pr_tmp1(i:i+wl-1));
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_tmp2(i:i+wl-1),pr_tmp2(i:i+wl-1),tas_tmp2(i:i+wl-1));
    reg_cesm_bs(iter,i+floor(wl/2)) = t1;
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg_with_carryover(ro_tmp2(i:i+wl-1),pr_tmp2(i:i+wl-1),tas_tmp2(i:i+wl-1));
    reg_co_cesm_bs(iter,i+floor(wl/2)) = t1;
  end
end
elast_cesm_pp = prctile(elast_cesm,[95 5]);
elast_cesm_pp(isnan(elast_cesm_pp)) = 0;
reg_cesm_pp = prctile(reg_cesm,[95 5]);
reg_cesm_pp(isnan(reg_cesm_pp)) = 0;
reg_co_cesm_pp = prctile(reg_co_cesm,[95 5]);
reg_co_cesm_pp(isnan(reg_co_cesm_pp)) = 0;
elast_cesm_bs(elast_cesm_bs==0) = NaN;
elast_cesm_bs_pp = prctile(elast_cesm_bs,[95 5]);
reg_cesm_bs(reg_cesm_bs==0) = NaN;
reg_cesm_bs_pp = prctile(reg_cesm_bs,[95 5]);
reg_co_cesm_bs(reg_co_cesm_bs==0) = NaN;
reg_co_cesm_bs_pp = prctile(reg_co_cesm_bs,[95 5]);



% -- Plotting ------------------
close all

xlim    = [1948 2010];
ylim    = [0 200];

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [5 5 35 30])

subplot(3,2,1)
title('(a) Precipitation and streamflow (Observations)')
% hold on
% jbfill((year1-floor(wl/2):year1+floor(wl/2)-1),zeros(1,wl)+ylim(2),zeros(1,wl)+ylim(1),[.9 .9 .9],'none')
% jbfill((year2-floor(wl/2):year2+floor(wl/2)-1),zeros(1,wl)+ylim(2),zeros(1,wl)+ylim(1),[.9 .9 .9],'none')
hold on
h1 = plot(time_obs_common,pr_obs3_wy_percent_mean,'b')
h2 = plot(time_obs_common,pr_obs1_wy_percent_mean,'b--')
h3 = plot(time_obs_common,pr_obs2_wy_percent_mean,'b:')
h4 = plot(1975:2004,pr_vano_percent_mean,'b-.')
h5 = plot(time_obs_common,ro_obs_wy_percent_mean,'r')
% h6 = plot(1975:2004,ro_vano_percent_mean,'r--')
set(gca,'Layer','top','XLim',xlim,'YLim',[-50 200])
hline(100,'k')
box on
ylabel(['Percent of mean ' num2str(refstart) '-' num2str(refende) ' (%)'])
xlabel('Time (Year)')
legend([h1 h2 h4 h3 h5],'Precipitation (Livneh)','Precipitation (GPCC)','Precipitation (PRISM)','Precipitation (Vano)','Streamflow (Reclamation)','Location','SouthWest')
legend boxoff

subplot(3,2,2)
title('(b) Temperature (Observations)')
hold on
h1 = plot(time_obs_common,tas_obs3_wy_anom,'b')
h2 = plot(time_obs_common,tas_obs1_wy_anom,'b--')
h3 = plot(time_obs_common,tas_obs2_wy_anom,'b:')
set(gca,'Layer','top','XLim',xlim,'YLim',[-2.5 2])
hline(0,'k')
box on
ylabel(['Temperature anomaly relative to ' num2str(refstart) '-' num2str(refende) ' (\circC)'])
xlabel('Time (Year)')
legend([h1 h2 h3],'Temperature (Livneh)','Temperature (BEST)','Temperature (PRISM)','Streamflow (Reclamation)','Location','SouthWest')
legend boxoff

subplot(3,2,3)
title('(c) P sensitivity (Observations)')
hold on
h11 = plot(1950:2008,reg3,'b')
h12 = plot(1950:2008,reg1,'b--')
h13 = plot(1950:2008,reg2,'b:')
h21 = plot(1950:2008,reg_co3,'m')
h22 = plot(1950:2008,reg_co1,'m--')
h23 = plot(1950:2008,reg_co2,'m:')
h31 = plot(1950:2008,elast3,'r')
h32 = plot(1950:2008,elast1,'r--')
h33 = plot(1950:2008,elast2,'r:')
h4 = plot([1990 1990],[reg_vano reg_vano],'bs','MarkerSize',8)
h5 = plot([1990 1990],[reg_co_vano reg_co_vano],'ms','MarkerSize',8)
h6 = plot([1990 1990],[elast_vano elast_vano],'rs','MarkerSize',8)
set(gca,'Layer','top','XLim',[1948 2010],'YLim',[.3 2.35])
% hline(0,'k')
box on
ylabel([num2str(wl) '-year running P sensitivity (%/%)'])
xlabel('Time (Year)')
% legend([h11 h12 h13 h4 h31 h32 h33 h6],'Regression (P Livneh / T Livneh)','Regression (P GPCC / T Livneh)','Regression (P PRISM / T Livneh)','Regression (P Vano / T Livneh)','Elasticity (P Livneh)','Elasticity (P GPCC)','Elasticity (P PRISM)','Elasticity (P Vano)','Location','SouthWest')
legend([h11 h12 h13 h4 h21 h22 h23 h5 h31 h32 h33 h6],'Regression 1 (P Livneh / T Livneh)','Regression 1 (P GPCC / T Livneh)','Regression 1 (P PRISM / T Livneh)','Regression 1 (P Vano / T Livneh)',...
'Regression 2 (P Livneh / T Livneh)','Regression 2 (P GPCC / T Livneh)','Regression 2 (P PRISM / T Livneh)','Regression 2 (P Vano / T Livneh)',...
'Elasticity (P Livneh)','Elasticity (P GPCC)','Elasticity (P PRISM)','Elasticity (P Vano)','Location','SouthWest')
legend boxoff

subplot(3,2,4)
title('(d) T sensitivity (Observations)')
hold on
h1 = plot(1950:2008,tas_reg3,'b')
h2 = plot(1950:2008,tas_reg1,'b--')
h3 = plot(1950:2008,tas_reg2,'b:')
plot(1950:2008,tas_reg_co3,'m')
plot(1950:2008,tas_reg_co1,'m--')
plot(1950:2008,tas_reg_co2,'m:')
set(gca,'Layer','top','XLim',[1948 2010],'YLim',[-24 -5])
% hline(0,'k')
box on
ylabel([num2str(wl) '-year running T sensitivity (%/\circC)'])
xlabel('Time (Year)')
legend([h1 h2 h3],'Regression (P Livneh / T Livneh)','Regression (P Livneh / T BEST)','Regression (P Livneh / T PRISM)','Location','SouthWest')
legend boxoff


subplot(3,2,5)
title('(e) P sensitivity (CESM simulations)')
hold on
% plot(reg_cesm_control)
s = size(reg_cesm_seg);
clear('tmp1','tmp2','tmp3')
for i = 1:s(2)
  tmp1_95(i,:) = prctile(squeeze(elast_cesm_seg(:,i,:)),95);
  tmp2_95(i,:) = prctile(squeeze(reg_cesm_seg(:,i,:)),95);
  tmp3_95(i,:) = prctile(squeeze(tas_reg_cesm_seg(:,i,:)),95);
  tmp4_95(i,:) = prctile(squeeze(tas_reg_co_cesm_seg(:,i,:)),95);
  tmp1_5(i,:) = prctile(squeeze(elast_cesm_seg(:,i,:)),5);
  tmp2_5(i,:) = prctile(squeeze(reg_cesm_seg(:,i,:)),5);
  tmp3_5(i,:) = prctile(squeeze(tas_reg_cesm_seg(:,i,:)),5);
  tmp4_5(i,:) = prctile(squeeze(tas_reg_co_cesm_seg(:,i,:)),5);
end
h3 = jbfill([segments2],squeeze(prctile(tmp2_95,95)),squeeze(prctile(tmp2_5,5)),[1 .7 .7],'none')
h4 = jbfill([segments2],squeeze(prctile(tmp1_95,95)),squeeze(prctile(tmp1_5,5)),[1 .9 .9],'none')
h5 = jbfill([segments2],squeeze(prctile(tmp4_95,95)),squeeze(prctile(tmp4_5,5)),[1 .5 .5],'none')
% h2 = jbfill([segments1],max(elast_cesm_control_seg'),min(elast_cesm_control_seg'),[.9 .9 1],'none')
% h1 = jbfill([segments1],max(reg_cesm_control_seg'),min(reg_cesm_control_seg'),[.7 .7 .7],'none')
h2 = jbfill([segments1],prctile(elast_cesm_control_seg',95),prctile(elast_cesm_control_seg',5),[.9 .9 1],'none')
h1 = jbfill([segments1],prctile(reg_cesm_control_seg',95),prctile(reg_cesm_control_seg',5),[.7 .7 .7],'none')
hold on
% % plot([segments1],min(reg_cesm_control_seg'),'k')
% % plot([segments1],max(reg_cesm_control_seg'),'k')
% % plot([segments1],min(elast_cesm_control_seg'),'b')
% % plot([segments1],max(elast_cesm_control_seg'),'b')
% plot([segments1],prctile(reg_cesm_control_seg',[5 95]),'k')
% plot([segments1],prctile(reg_co_cesm_control_seg',[5 95]),'m')
% plot([segments1],prctile(elast_cesm_control_seg',[5 95]),'b')
% % plot([segments2],squeeze(min(min(elast_cesm_seg))),'m')
% % plot([segments2],squeeze(max(max(elast_cesm_seg))),'m')
% plot([segments2],squeeze(prctile(tmp2_95,95)),'r')
% plot([segments2],squeeze(prctile(tmp2_5,5)),'r')
% plot([segments2],squeeze(prctile(tmp1_95,95)),'m')
% plot([segments2],squeeze(prctile(tmp1_5,5)),'m')
set(gca,'XScale','log','Layer','top','YLim',[0 3.5])
h0 = vline(refende-refstart+1,'r--')
hline(0,'k--')
plot([segments1],nanmedian(reg_cesm_control_seg,2),'k')
plot([segments1],nanmedian(elast_cesm_control_seg,2),'b')
plot([segments2],squeeze(nanmedian(nanmedian(reg_cesm_seg(:,:,:)))),'r')
plot([segments2],squeeze(nanmedian(nanmedian(elast_cesm_seg(:,:,:)))),'m')
% plot([segments1(1) segments1(end)],[reg_cesm_control_seg(end,1) reg_cesm_control_seg(end,1)],'k')
% plot([segments1(1) segments1(end)],[elast_cesm_control_seg(end,1) elast_cesm_control_seg(end,1)],'b')
% plot([segments2(1) segments2(end)],[median(reg_cesm_seg(:,1,end)) median(reg_cesm_seg(:,1,end))],'r')
% plot([segments2(1) segments2(end)],[median(elast_cesm_seg(:,1,end)) median(elast_cesm_seg(:,1,end))],'m')
box on
ylabel(['P sensitivity (%/%)'])
xlabel('Segment length (Years)')
legend([h1 h2 h3 h4 h0],'Regression (CESM piControl)','Elasticity (CESM piControl)','Regression (CESM LE historical)','Elasticity (CESM LE historical)','Length used in paper','Location','SouthEast')
legend boxoff

subplot(3,2,6)
title('(f) T sensitivity (CESM simulations)')
hold on
h2 = jbfill([segments2],squeeze(prctile(tmp3_95,95)),squeeze(prctile(tmp3_5,5)),[1 .7 .7],'none')
% h1 = jbfill([segments1],max(tas_reg_cesm_control_seg'),min(tas_reg_cesm_control_seg'),[.7 .7 .7],'none')
h1 = jbfill([segments1],prctile(tas_reg_cesm_control_seg',95),prctile(tas_reg_cesm_control_seg',5),[.7 .7 .7],'none')
hold on
% plot([segments1],min(tas_reg_cesm_control_seg'),'k')
% plot([segments1],max(tas_reg_cesm_control_seg'),'k')
plot([segments1],prctile(tas_reg_cesm_control_seg',[5 95]),'k')
plot([segments2],squeeze(prctile(tmp3_95,95)),'r')
plot([segments2],squeeze(prctile(tmp3_5,5)),'r')
set(gca,'XScale','log','Layer','top','YLim',[-25 20])
h0 = vline(refende-refstart+1,'r--')
hline(0,'k--')
% plot([segments1(1) segments1(end)],[tas_reg_cesm_control_seg(end,1) tas_reg_cesm_control_seg(end,1)],'k')
plot([segments1],nanmedian(tas_reg_cesm_control_seg,2),'k')
plot([segments2],squeeze(nanmedian(nanmedian(tas_reg_cesm_seg(:,:,:)))),'r')
box on
ylabel(['T sensitivity (%/\circC)'])
xlabel('Segment length (Years)')
legend([h1 h2 h0],'Regression (CESM piControl)','Regression (CESM LE historical)','Length used in paper','Location','SouthEast')
legend boxoff

% return
set(gcf,'PaperPositionMode','auto');
fileo = [figpath 'Fig_S1_' vari];
print('-r300','-loose', '-depsc', [fileo '.eps'])
saveas(gcf,fileo,'jpg')
return
end







if figS5 == 1
% -----------------------------------------------------------------------------
% Figure S5?
% ----

% -- temporal variability of Sankar's epsilon --

% -- load streamflow and precip data that Julie Vano used:
data = [1418743.917	1.045659
955518.4167	0.897887
453309.5	0.732259
1260119	1.05651
1505322.417	1.170061
1476323.167	1.208161
748600.3333	1.058718
1421620.583	1.142557
2032799.583	1.242915
2109522.75	1.377675
1754051.833	1.202608
1901946.833	1.273484
1355472.917	1.068602
982124.4167	0.984579
836128.6667	0.881206
790779.8333	0.971848
1006257.167	1.015006
919331.9167	1.086383
1526693	1.318106
885648.75	0.920851
1688850.417	1.294706
1204230.917	1.073559
1801807	1.438768
1401438.167	1.146939
1380431	1.09642
909142.8333	0.958439
910805.4167	0.896307
517506.75	0.803036
876463.1667	0.886962
802812.1667	1.06447];
ro_vano_percent_mean = (data(:,1)./mean(data(:,1)))*100;
pr_vano_percent_mean = data(:,2)*100;
% ro_vano_percent_mean = data(:,1);
% pr_vano_percent_mean = data(:,2);

wl      = 30; % window length
elast1  = NaN(size(ro_obs_wy_percent_mean));
elast2  = elast1;
reg1    = elast1;
reg2    = elast1;
year1   = 1955;
year2   = 1990;
% -- obs --
for i = 1:length(ro_obs_wy_percent_mean)-wl
  elast1(i+floor(wl/2)) = elast_sankar(ro_obs_wy(i:i+wl-1),pr_obs_wy(i:i+wl-1));
  elast2(i+floor(wl/2)) = elast_sankar(ro_obs_wy(i:i+wl-1),pr_obs2_wy(i:i+wl-1));
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_obs_wy_percent_mean(i:i+wl-1)-100,pr_obs_wy_percent_mean(i:i+wl-1)-100,tas_obs_wy_anom(i:i+wl-1));
  reg1(i+floor(wl/2)) = t1;
  [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_obs_wy_percent_mean(i:i+wl-1)-100,pr_obs2_wy_percent_mean(i:i+wl-1)-100,tas_obs_wy_anom(i:i+wl-1));
  reg2(i+floor(wl/2)) = t1;
end
% -- bootstrap obs:
bl = 2; % block length in years
iterations = 100;
elast1_bs  = NaN(iterations,length(ro_obs_wy_percent_mean));
elast2_bs  = NaN(iterations,length(ro_obs_wy_percent_mean));
reg1_bs    = NaN(iterations,length(ro_obs_wy_percent_mean));
reg2_bs    = NaN(iterations,length(ro_obs_wy_percent_mean));
for iter = 1:iterations
  idx = randi([1 length(ro_obs_wy)-bl+1],floor(length(ro_obs_wy)/bl),1);
  j = 1;
  for i = 1:length(idx)
    ro_tmp1(j:j+bl-1)   = ro_obs_wy(idx(i):idx(i)+bl-1);
    pr_tmp1(j:j+bl-1)   = pr_obs_wy(idx(i):idx(i)+bl-1);
    pr2_tmp1(j:j+bl-1)  = pr_obs2_wy(idx(i):idx(i)+bl-1);
    ro_tmp2(j:j+bl-1)   = ro_obs_wy_percent_mean(idx(i):idx(i)+bl-1)-100;
    pr_tmp2(j:j+bl-1)   = pr_obs_wy_percent_mean(idx(i):idx(i)+bl-1)-100;
    pr2_tmp2(j:j+bl-1)  = pr_obs2_wy_percent_mean(idx(i):idx(i)+bl-1)-100;
    tas_tmp2(j:j+bl-1)  = tas_obs_wy_anom(idx(i):idx(i)+bl-1);
    j = j+bl;
  end
  for i = 1:length(ro_obs_wy)-wl
    elast1_bs(iter,i+floor(wl/2)) = elast_sankar(ro_tmp1(i:i+wl-1),pr_tmp1(i:i+wl-1));
    elast2_bs(iter,i+floor(wl/2)) = elast_sankar(ro_tmp1(i:i+wl-1),pr2_tmp1(i:i+wl-1));
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_tmp2(i:i+wl-1),pr_tmp2(i:i+wl-1),tas_tmp2(i:i+wl-1));
    reg1_bs(iter,i+floor(wl/2)) = t1;
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_tmp2(i:i+wl-1),pr2_tmp2(i:i+wl-1),tas_tmp2(i:i+wl-1));
    reg2_bs(iter,i+floor(wl/2)) = t1;
  end
end
elast1_bs(elast1_bs==0) = NaN;
elast2_bs(elast2_bs==0) = NaN;
elast1_bs_pp = prctile(elast1_bs,[95 5]);
elast2_bs_pp = prctile(elast2_bs,[95 5]);
reg1_bs(reg1_bs==0) = NaN;
reg2_bs(reg2_bs==0) = NaN;
reg1_bs_pp = prctile(reg1_bs,[95 5]);
reg2_bs_pp = prctile(reg2_bs,[95 5]);


% -- cesm le --
elast_cesm    = NaN(size(ro_obs_wy_percent_mean));
reg_cesm      = NaN(size(ro_obs_wy_percent_mean));
elast_cesm_bs = NaN(iterations,length(ro_obs_wy_percent_mean));
reg_cesm_bs   = NaN(iterations,length(ro_obs_wy_percent_mean));
for e = 1:40
  for i = 1:length(ro_obs_wy_percent_mean)-wl
    elast_cesm(e,i+floor(wl/2)) = elast_sankar(ro_cesm_wy(e,8+i:8+i+wl-1),pr_cesm_wy(e,8+i:8+i+wl-1));
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_cesm_wy_percent_mean(e,8+i:8+i+wl-1)-100,pr_cesm_wy_percent_mean(e,8+i:8+i+wl-1)-100,tas_cesm_wy_anom(e,8+i:8+i+wl-1));
    reg_cesm(e,i+floor(wl/2)) = t1;
  end
end
% -- bootstrap cesm le #e:
e = 2; % ensemble member
for iter = 1:iterations
  idx = randi([1 length(ro_obs_wy)-bl+1],floor(length(ro_obs_wy)/bl),1);
  j = 1;
  for i = 1:length(idx)
    ro_tmp1(j:j+bl-1)   = ro_cesm_wy(e,8+idx(i):8+idx(i)+bl-1);
    pr_tmp1(j:j+bl-1)   = pr_cesm_wy(e,8+idx(i):8+idx(i)+bl-1);
    ro_tmp2(j:j+bl-1)   = ro_cesm_wy_percent_mean(e,8+idx(i):8+idx(i)+bl-1)-100;
    pr_tmp2(j:j+bl-1)   = pr_cesm_wy_percent_mean(e,8+idx(i):8+idx(i)+bl-1)-100;
    tas_tmp2(j:j+bl-1)  = tas_cesm_wy_anom(e,8+idx(i):8+idx(i)+bl-1);
    j = j+bl;
  end
  for i = 1:length(ro_obs_wy)-wl
    elast_cesm_bs(iter,i+floor(wl/2)) = elast_sankar(ro_tmp1(i:i+wl-1),pr_tmp1(i:i+wl-1));
    [t1 t2 t3 t4 t5 t6] = runoff_sens_reg(ro_tmp2(i:i+wl-1),pr_tmp2(i:i+wl-1),tas_tmp2(i:i+wl-1));
    reg_cesm_bs(iter,i+floor(wl/2)) = t1;
  end
end
elast_cesm_pp = prctile(elast_cesm,[95 5]);
elast_cesm_pp(isnan(elast_cesm_pp)) = 0;
reg_cesm_pp = prctile(reg_cesm,[95 5]);
reg_cesm_pp(isnan(reg_cesm_pp)) = 0;
elast_cesm_bs(elast_cesm_bs==0) = NaN;
elast_cesm_bs_pp = prctile(elast_cesm_bs,[95 5]);
reg_cesm_bs(reg_cesm_bs==0) = NaN;
reg_cesm_bs_pp = prctile(reg_cesm_bs,[95 5]);


% -- Plotting -------
close all

xlim    = [1925 2010];
ylim    = [-50 200];

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [5 5 35 22])

subplot(2,4,[2 3])
title('Precipitation and streamflow time series')
hold on
jbfill((year1-floor(wl/2):year1+floor(wl/2)-1),zeros(1,wl)+ylim(2),zeros(1,wl)+ylim(1),[.9 .9 .9],'none')
jbfill((year2-floor(wl/2):year2+floor(wl/2)-1),zeros(1,wl)+ylim(2),zeros(1,wl)+ylim(1),[.9 .9 .9],'none')
hold on
h1 = plot(time_obs_common,pr_obs_wy_percent_mean,'k')
h3 = plot(time_obs_common,ro_obs_wy_percent_mean,'r')
% plot(1929:2008,rm(pr_obs_wy_percent_mean,10),'k','LineWidth',1.5)
% plot(1929:2008,rm(ro_obs_wy_percent_mean,10),'r','LineWidth',1.5)
if strcmp(vari,'UC')==1
  h2 = plot(time_obs_common,pr_obs2_wy_percent_mean,'b')
  h4 = plot(1975:2004,pr_vano_percent_mean,'b--')
  h5 = plot(1975:2004,ro_vano_percent_mean,'r--')
end
set(gca,'XLim',xlim,'YLim',ylim,'Layer','top')
hline(100,'k')
box on
ylabel('Percent of mean (%)')
if strcmp(vari,'UC')==1
  legend([h1 h2 h4 h3 h5],['Precipitation (' obs_data_name ')'],['Precipitation (' obs2_data_name ')'],'Precipitation (Vano)','Streamflow (Reclamation)','Streamflow (Vano)','Location','SouthWest')
else
  legend([h1 h3],['Precipitation (' obs_data_name ')'],'Streamflow (Reclamation)','Location','SouthWest')
end

subplot(2,4,[5 6])
title('P elasticity based on Sankar''s epsilon')
% jbfill(1921:2100,max(elast_cesm,1),min(elast_cesm,1),[1 0 0],'none') % min/max of CESM LE
jbfill(1929+(wl/2):2008-(wl/2),elast_cesm_pp(1,(wl/2)+1:end-(wl/2)),elast_cesm_pp(2,(wl/2)+1:end-(wl/2)),[1 .7 .7],'none') % 5-95% of CESM LE
% jbfill(1929+(wl/2):2008-(wl/2)-1,elast_cesm_bs_pp(1,(wl/2)+1:end-(wl/2)-1),elast_cesm_bs_pp(2,(wl/2)+1:end-(wl/2)-1),[1 .9 .9],'none') % 5-95% of bootstrapped CESM LE
jbfill(1929+(wl/2):2008-(wl/2),elast1_bs_pp(1,(wl/2)+1:end-(wl/2)),elast1_bs_pp(2,(wl/2)+1:end-(wl/2)),[.7 .7 .7],'none') % 5-95% of bootstrapped obs
jbfill(1929+(wl/2):2008-(wl/2),elast2_bs_pp(1,(wl/2)+1:end-(wl/2)),elast2_bs_pp(2,(wl/2)+1:end-(wl/2)),[.7 .7 1],'none') % 5-95% of bootstrapped obs
hold on
h4 = plot(1929+(wl/2):2008-(wl/2),nanmean(elast_cesm(:,(wl/2)+1:end-(wl/2)),1),'r')
h1 = plot(1929:2008,elast1,'k')
h5 = plot(1929:2008,elast1*0+pr_ro_reg2_obs_wy(1),'k')
if strcmp(vari,'UC')==1
  h2 = plot(1929:2008,elast2,'b')
  h6 = plot(1929:2008,elast1*0+pr_ro_reg2_obs2_wy(1),'b')
end
plot([year1 year1],[elast1(year1-1929+1) elast1(year1-1929+1)],'kx','MarkerSize',10,'LineWidth',2)
plot([year2 year2],[elast1(year2-1929+1) elast1(year2-1929+1)],'k.','MarkerSize',18)
if strcmp(vari,'UC')==1
  plot([year1 year1],[elast2(year1-1929+1) elast2(year1-1929+1)],'bx','MarkerSize',10,'LineWidth',2)
  plot([year2 year2],[elast2(year2-1929+1) elast2(year2-1929+1)],'b.','MarkerSize',18)
  h3 = plot([year2 year2],[elast_sankar(ro_vano_percent_mean,pr_vano_percent_mean) elast_sankar(ro_vano_percent_mean,pr_vano_percent_mean)],'m.','MarkerSize',18)
end
set(gca,'Layer','top','XLim',xlim,'YLim',[0.7 3.2])
hline(1,'k')
hline(1.5,'k--')
hline(2,'k--')
hline(2.5,'k--')
box on
xlabel('Time (Year)')
ylabel('P elasticity (%/%)')
if strcmp(vari,'UC')==1
  % legend([h1 h2 h3],[ num2str(wl) '-year running window (using ' obs_data_name ')'],[ num2str(wl) '-year running window (using ' obs2_data_name ')'],'Vano','Location','NorthWest')
  legend([h1 h2 h3 h4],[ num2str(wl) '-year running window (using ' obs_data_name ')'],[ num2str(wl) '-year running window (using ' obs2_data_name ')'],'Vano 1975-2004','CESM LE (40 members)','Location','NorthWest')
else
  legend([h1],[ num2str(wl) '-year running window (using ' obs_data_name ')'],'Location','NorthWest')
end


subplot(2,4,[7 8])
title('P elasticity based on regression')
% jbfill(1921:2100,max(reg_cesm,1),min(reg_cesm,1),[1 0 0],'none') % min/max of CESM LE
jbfill(1929+(wl/2):2008-(wl/2),reg_cesm_pp(1,(wl/2)+1:end-(wl/2)),reg_cesm_pp(2,(wl/2)+1:end-(wl/2)),[1 .7 .7],'none') % 5-95% of CESM LE
% jbfill(1929+(wl/2):2008-(wl/2)-1,reg_cesm_bs_pp(1,(wl/2)+1:end-(wl/2)-1),reg_cesm_bs_pp(2,(wl/2)+1:end-(wl/2)-1),[1 .9 .9],'none') % 5-95% of bootstrapped CESM LE
jbfill(1929+(wl/2):2008-(wl/2),reg1_bs_pp(1,(wl/2)+1:end-(wl/2)),reg1_bs_pp(2,(wl/2)+1:end-(wl/2)),[.7 .7 .7],'none') % 5-95% of bootstrapped obs
jbfill(1929+(wl/2):2008-(wl/2),reg2_bs_pp(1,(wl/2)+1:end-(wl/2)),reg2_bs_pp(2,(wl/2)+1:end-(wl/2)),[.7 .7 1],'none') % 5-95% of bootstrapped obs
hold on
h4 = plot(1929+(wl/2):2008-(wl/2),nanmean(reg_cesm(:,(wl/2)+1:end-(wl/2)),1),'r')
h1 = plot(1929:2008,reg1,'k')
h5 = plot(1929:2008,reg1*0+pr_ro_reg_obs_wy(1),'k')
if strcmp(vari,'UC')==1
  h2 = plot(1929:2008,reg2,'b')
  h6 = plot(1929:2008,reg2*0+pr_ro_reg_obs2_wy(1),'b')
end
plot([year1 year1],[reg1(year1-1929+1) reg1(year1-1929+1)],'kx','MarkerSize',10,'LineWidth',2)
plot([year2 year2],[reg1(year2-1929+1) reg1(year2-1929+1)],'k.','MarkerSize',18)
if strcmp(vari,'UC')==1
  plot([year1 year1],[reg2(year1-1929+1) reg2(year1-1929+1)],'bx','MarkerSize',10,'LineWidth',2)
  plot([year2 year2],[reg2(year2-1929+1) reg2(year2-1929+1)],'b.','MarkerSize',18)
  tmp = runoff_sens_reg(ro_vano_percent_mean-100,pr_vano_percent_mean-100,tas_obs_wy_anom(1975-1929+1:2004-1929+1))
  h3 = plot([year2 year2],[tmp tmp],'m.','MarkerSize',18)
end
set(gca,'Layer','top','XLim',xlim,'YLim',[0.7 3.2])
hline(1,'k')
hline(1.5,'k--')
hline(2,'k--')
hline(2.5,'k--')
box on
xlabel('Time (Year)')
ylabel('P elasticity (%/%)')


% return
set(gcf,'PaperPositionMode','auto');
fileo = [figpath 'Fig_S5_' vari];
print('-r300','-loose', '-depsc', [fileo '.eps'])
saveas(gcf,fileo,'jpg')
% return
end










if fig3 == 1
% ------------------------------------------------------------------------
% Figure 3? in paper: Obs and CMIP5 historical2rcp85 sensitivities
% ---------
close all


xx = 6; % x-axis range over which CESM LE members are spread
ylim = [0 length(models_common)+xx+3]; % can be ylim for flipped plots

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [10 10 25 20])

subplot(1,2,1)
% subplot(6,4,[1 2 5 6 9 10 13 14 17 18])
hold on
% jbfill(ylim(1):ylim(2),(ylim(1):ylim(2))*0+pr_ro_reg_obs_wy(3),(ylim(1):ylim(2))*0+pr_ro_reg_obs_wy(2),[1 .7 .7],'none')
jbfill([pr_ro_reg_obs_wy(3) pr_ro_reg_obs_wy(2)], [pr_ro_reg_obs_wy(3) pr_ro_reg_obs_wy(2)]*0+ylim(1), [pr_ro_reg_obs_wy(3) pr_ro_reg_obs_wy(2)]*0+ylim(2),[1 .7 .7],'none')
hold on
% plot([ylim(1) ylim(2)],[pr_ro_reg_obs_wy(1) pr_ro_reg_obs_wy(1)],'r--')
plot([pr_ro_reg_obs_wy(1) pr_ro_reg_obs_wy(1)],[ylim(1) ylim(2)],'r--')
plot([pr_ro_reg2_obs_wy(1) pr_ro_reg2_obs_wy(1)],[ylim(1) ylim(2)],'r-.')
[tmp,idx1] = sort(pr_ro_reg_h2r_wy,1);
for m = 1:length(models_common)
  % h2 = plot([m+.0 m+.0],[pr_ro_reg_wy(idx1(m,1),1,2) pr_ro_reg_wy(idx1(m,1),1,3)],'k')
  %      plot([m+.0 m+.0],[pr_ro_reg_wy(idx1(m,1),1,1) pr_ro_reg_wy(idx1(m,1),1,1)],'k.','MarkerSize',12)
  % h3 = plot([m+.2 m+.2],[pr_ro_reg_wy(idx1(m,1),2,2) pr_ro_reg_wy(idx1(m,1),2,3)],'m')
  %      plot([m+.2 m+.2],[pr_ro_reg_wy(idx1(m,1),2,1) pr_ro_reg_wy(idx1(m,1),2,1)],'m.','MarkerSize',12)
  % h1 = plot([m+.4 m+.4],[pr_ro_reg_h2r_wy(idx1(m,1),2) pr_ro_reg_h2r_wy(idx1(m,1),3)],'b')
  %      plot([m+.4 m+.4],[pr_ro_reg_h2r_wy(idx1(m,1),1) pr_ro_reg_h2r_wy(idx1(m,1),1)],'b.','MarkerSize',12)
  h2 = plot([pr_ro_reg_wy(idx1(m,1),1,2) pr_ro_reg_wy(idx1(m,1),1,3)],[m+.0 m+.0],'k')
       plot([pr_ro_reg_wy(idx1(m,1),1,1) pr_ro_reg_wy(idx1(m,1),1,1)],[m+.0 m+.0],'k.','MarkerSize',12)
  h3 = plot([pr_ro_reg_wy(idx1(m,1),2,2) pr_ro_reg_wy(idx1(m,1),2,3)],[m+.2 m+.2],'m')
       plot([pr_ro_reg_wy(idx1(m,1),2,1) pr_ro_reg_wy(idx1(m,1),2,1)],[m+.2 m+.2],'m.','MarkerSize',12)
  h1 = plot([pr_ro_reg_h2r_wy(idx1(m,1),2) pr_ro_reg_h2r_wy(idx1(m,1),3)],[m+.4 m+.4],'b')
       plot([pr_ro_reg_h2r_wy(idx1(m,1),1) pr_ro_reg_h2r_wy(idx1(m,1),1)],[m+.4 m+.4],'b.','MarkerSize',12)
       plot([pr_ro_reg2_h2r_wy(idx1(m,1),1) pr_ro_reg2_h2r_wy(idx1(m,1),1)],[m+.4 m+.4],'b^','MarkerSize',6)
end
for e = 1:40
  % plot([m+1+e*(xx/40) m+1+e*(xx/40)],[pr_ro_reg_cesm_wy(e,2) pr_ro_reg_cesm_wy(e,3)],'b')
  % plot([m+1+e*(xx/40) m+1+e*(xx/40)],[pr_ro_reg_cesm_wy(e,1) pr_ro_reg_cesm_wy(e,1)],'b.','MarkerSize',12)
  plot([pr_ro_reg_cesm_wy(e,2) pr_ro_reg_cesm_wy(e,3)],[m+1+e*(xx/40) m+1+e*(xx/40)],'b')
  plot([pr_ro_reg_cesm_wy(e,1) pr_ro_reg_cesm_wy(e,1)],[m+1+e*(xx/40) m+1+e*(xx/40)],'b.','MarkerSize',12)
  plot([pr_ro_reg2_cesm_wy(e,1) pr_ro_reg2_cesm_wy(e,1)],[m+1+e*(xx/40) m+1+e*(xx/40)],'b^','MarkerSize',6)
end
% h0 = plot([m+xx+2 m+xx+2],[pr_ro_reg_obs_wy(2) pr_ro_reg_obs_wy(3)],'r')
%      plot([m+xx+2 m+xx+2],[pr_ro_reg_obs_wy(1) pr_ro_reg_obs_wy(1)],'r.','MarkerSize',12)
h0 = plot([pr_ro_reg_obs_wy(2) pr_ro_reg_obs_wy(3)],[m+xx+2 m+xx+2],'r')
plot([pr_ro_reg_obs_wy(1) pr_ro_reg_obs_wy(1)],[m+xx+2 m+xx+2],'r.','MarkerSize',12)
plot([pr_ro_reg2_obs_wy(1) pr_ro_reg2_obs_wy(1)],[m+xx+2 m+xx+2],'r^','MarkerSize',6)
% set(gca,'Layer','top','XLim',ylim,'XTick',[1:length(models_common) length(models_common)+1+(xx/2) ylim(2)-1 ],'XTickLabel',[models_common(idx1(:,1)) '' 'CESM LE' '' '' 'Observations'])
set(gca,'Layer','top','YLim',ylim,'YTick',[1:length(models_common) length(models_common)+1+(xx/2) ylim(2)-1 ],'YTickLabel',[models_common(idx1(:,1)) '' 'CESM LE' '' '' 'Observations'])
% set(gca,'XLim',[0 4],'XTick',[0:1:4])
% xticklabel_rotate([],45,[])%,'Fontsize',14)
% ylabel('dR/dP (%/%)')
xlabel('\DeltaQ/\DeltaP (%/%)')
% xlabel('Models')
ylabel('Models')
% legend([h0(1) h2(1) h3(1) h1(1)],'Observations','CMIP5 piControl','CMIP5 1pctCO2','CMIP5 historical','Location','NorthWest')
legend([h0(1) h2(1) h3(1) h1(1)],'Observations','CMIP5 piControl','CMIP5 1pctCO2','CMIP5 historical','Location','SouthEast')
box on


subplot(1,2,2)
% subplot(6,4,[3 4 7 8 11 12 15 16 19 20])
hold on
% jbfill(ylim(1):ylim(2),(ylim(1):ylim(2))*0+tas_ro_reg_obs_wy(3),(ylim(1):ylim(2))*0+tas_ro_reg_obs_wy(2),[1 .7 .7],'none')
jbfill([tas_ro_reg_obs_wy(3) tas_ro_reg_obs_wy(2)], [tas_ro_reg_obs_wy(3) tas_ro_reg_obs_wy(2)]*0+ylim(1), [tas_ro_reg_obs_wy(3) tas_ro_reg_obs_wy(2)]*0+ylim(2),[1 .7 .7],'none')
hold on
% plot([ylim(1) ylim(2)],[tas_ro_reg_obs_wy(1) tas_ro_reg_obs_wy(1)],'r--')
plot([tas_ro_reg_obs_wy(1) tas_ro_reg_obs_wy(1)],[ylim(1) ylim(2)],'r--')
% [tmp,idx1] = sort(tas_ro_reg_h2r_wy,1); % re-sort according to T sens
for m = 1:length(models_common)
  % h2 = plot([m+.0 m+.0],[tas_ro_reg_wy(idx1(m,1),1,2) tas_ro_reg_wy(idx1(m,1),1,3)],'k')
  %      plot([m+.0 m+.0],[tas_ro_reg_wy(idx1(m,1),1,1) tas_ro_reg_wy(idx1(m,1),1,1)],'k.','MarkerSize',12)
  % h3 = plot([m+.2 m+.2],[tas_ro_reg_wy(idx1(m,1),2,2) tas_ro_reg_wy(idx1(m,1),2,3)],'m')
  %      plot([m+.2 m+.2],[tas_ro_reg_wy(idx1(m,1),2,1) tas_ro_reg_wy(idx1(m,1),2,1)],'m.','MarkerSize',12)
  % h1 = plot([m+.4 m+.4],[tas_ro_reg_h2r_wy(idx1(m,1),2) tas_ro_reg_h2r_wy(idx1(m,1),3)],'b')
  %      plot([m+.4 m+.4],[tas_ro_reg_h2r_wy(idx1(m,1),1) tas_ro_reg_h2r_wy(idx1(m,1),1)],'b.','MarkerSize',12)
  h2 = plot([tas_ro_reg_wy(idx1(m,1),1,2) tas_ro_reg_wy(idx1(m,1),1,3)],[m+.0 m+.0],'k')
       plot([tas_ro_reg_wy(idx1(m,1),1,1) tas_ro_reg_wy(idx1(m,1),1,1)],[m+.0 m+.0],'k.','MarkerSize',12)
  h3 = plot([tas_ro_reg_wy(idx1(m,1),2,2) tas_ro_reg_wy(idx1(m,1),2,3)],[m+.2 m+.2],'m')
       plot([tas_ro_reg_wy(idx1(m,1),2,1) tas_ro_reg_wy(idx1(m,1),2,1)],[m+.2 m+.2],'m.','MarkerSize',12)
  h1 = plot([tas_ro_reg_h2r_wy(idx1(m,1),2) tas_ro_reg_h2r_wy(idx1(m,1),3)],[m+.4 m+.4],'b')
       plot([tas_ro_reg_h2r_wy(idx1(m,1),1) tas_ro_reg_h2r_wy(idx1(m,1),1)],[m+.4 m+.4],'b.','MarkerSize',12)
end
for e = 1:40
  % plot([m+1+e*(xx/40) m+1+e*(xx/40)],[tas_ro_reg_cesm_wy(e,2) tas_ro_reg_cesm_wy(e,3)],'b')
  % plot([m+1+e*(xx/40) m+1+e*(xx/40)],[tas_ro_reg_cesm_wy(e,1) tas_ro_reg_cesm_wy(e,1)],'b.','MarkerSize',12)
  plot([tas_ro_reg_cesm_wy(e,2) tas_ro_reg_cesm_wy(e,3)],[m+1+e*(xx/40) m+1+e*(xx/40)],'b')
  plot([tas_ro_reg_cesm_wy(e,1) tas_ro_reg_cesm_wy(e,1)],[m+1+e*(xx/40) m+1+e*(xx/40)],'b.','MarkerSize',12)
end
% h0 = plot([m+xx+2 m+xx+2],[tas_ro_reg_obs_wy(2) tas_ro_reg_obs_wy(3)],'r')
%      plot([m+xx+2 m+xx+2],[tas_ro_reg_obs_wy(1) tas_ro_reg_obs_wy(1)],'r.','MarkerSize',12)
h0 = plot([tas_ro_reg_obs_wy(2) tas_ro_reg_obs_wy(3)],[m+xx+2 m+xx+2],'r')
     plot([tas_ro_reg_obs_wy(1) tas_ro_reg_obs_wy(1)],[m+xx+2 m+xx+2],'r.','MarkerSize',12)
% set(gca,'Layer','top','XLim',ylim,'XTick',[1:length(models_common) length(models_common)+1+(xx/2) ylim(2)-1 ],'XTickLabel',[models_common(idx1(:,1)) '' 'CESM LE' '' '' 'Observations'])
set(gca,'Layer','top','YLim',ylim,'YTick',[1:length(models_common) length(models_common)+1+(xx/2) ylim(2)-1 ],'YTickLabel',[models_common(idx1(:,1)) '' 'CESM LE' '' '' 'Observations'])
set(gca,'XLim',[-40 30],'XTick',[-40:10:30])
% xticklabel_rotate([],45,[])%,'Fontsize',14)
% ylabel('dR/dT (%/\circC)')
xlabel('\DeltaQ/\DeltaT (%/\circC)')
% xlabel('Models')
ylabel('Models')
vline(0,'k')
% legend([h0(1) h2(1) h3(1) h1(1)],'Observations','CMIP5 piControl','CMIP5 1pctCO2','CMIP5 historical','Location','NorthWest')
% legend([h0(1) h2(1) h3(1) h1(1)],'Observations','CMIP5 piControl','CMIP5 1pctCO2','CMIP5 historical','Location','SouthEast')
box on

% return
set(gcf,'PaperPositionMode','auto');
fileo = [figpath 'Fig_3_' vari '_part1' co_tag];
print('-r300','-loose', '-depsc', [fileo '.eps'])
saveas(gcf,fileo,'jpg')
% return


% --------------

xlim = [0 3.5]; % P sens
ylim = [-30 25]; % T sens

figure2 = figure;
set(figure2, 'units', 'centimeters', 'pos', [10 10 12 10])
hold on
jbfill([pr_ro_reg_obs_wy(3) pr_ro_reg_obs_wy(2)], [pr_ro_reg_obs_wy(3) pr_ro_reg_obs_wy(2)]*0+30, [pr_ro_reg_obs_wy(3) pr_ro_reg_obs_wy(2)]*0-40,[1 .7 .7],'none')
jbfill([0 4],[tas_ro_reg_obs_wy(3) tas_ro_reg_obs_wy(3)], [tas_ro_reg_obs_wy(2) tas_ro_reg_obs_wy(2)],[1 .7 .7],'none')
hold on
% for m = 1:length(models_common)
%   plot([pr_ro_reg_h2r_wy(m,2) pr_ro_reg_h2r_wy(m,3)],[tas_ro_reg_h2r_wy(m,1) tas_ro_reg_h2r_wy(m,1)],'b')
%   plot([pr_ro_reg_h2r_wy(m,1) pr_ro_reg_h2r_wy(m,1)],[tas_ro_reg_h2r_wy(m,2) tas_ro_reg_h2r_wy(m,3)],'b')
% end
% for e = 1:40
%   plot([pr_ro_reg_cesm_wy(e,2) pr_ro_reg_cesm_wy(e,3)],[tas_ro_reg_cesm_wy(e,1) tas_ro_reg_cesm_wy(e,1)],'k')
%   plot([pr_ro_reg_cesm_wy(e,1) pr_ro_reg_cesm_wy(e,1)],[tas_ro_reg_cesm_wy(e,2) tas_ro_reg_cesm_wy(e,3)],'k')
% end
% for e = 1:40
%   h = ellipse((pr_ro_reg_cesm_wy(e,3)-pr_ro_reg_cesm_wy(e,2))/2,((tas_ro_reg_cesm_wy(e,3)+100)-(tas_ro_reg_cesm_wy(e,2)+100))/2,0,pr_ro_reg_cesm_wy(e,1),tas_ro_reg_cesm_wy(e,1))
%   set(h,'Color',[.5 .5 .5])
% end
% --- beginning of CESM part ---
% -- uncertainty from regression and solely due to internal variability using CESM LE:
% h = ellipse((prctile(pr_ro_reg_cesm_wy(:,3),95)-prctile(pr_ro_reg_cesm_wy(:,2),5))/2,(prctile(tas_ro_reg_cesm_wy(:,3)+100,95)-prctile(tas_ro_reg_cesm_wy(:,2)+100,5))/2,0,median(pr_ro_reg_cesm_wy(:,1)),median(tas_ro_reg_cesm_wy(:,1)))
% % set(h,'Color',[0 0 0],'LineStyle','-')
% set(h,'Color',[.7 .7 .7])
% x=get(h,'Xdata');
% y=get(h,'Ydata');
% patch(x,y,[.7 .7 .7],'edgecolor','none');
% h = ellipse((prctile(pr_ro_reg_cesm_wy(:,1),95)-prctile(pr_ro_reg_cesm_wy(:,1),5))/2,(prctile(tas_ro_reg_cesm_wy(:,1)+100,95)-prctile(tas_ro_reg_cesm_wy(:,1)+100,5))/2,0,median(pr_ro_reg_cesm_wy(:,1)),median(tas_ro_reg_cesm_wy(:,1)))
% % set(h,'Color',[0 0 0],'LineStyle','--')
% set(h,'Color',[.5 .5 .5])
% x=get(h,'Xdata');
% y=get(h,'Ydata');
% patch(x,y,[.5 .5 .5],'edgecolor','none');
% plot(median(pr_ro_reg_cesm_wy(:,1)),median(tas_ro_reg_cesm_wy(:,1)),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
% --- end of CESM part ---
for m = 1:length(models_common)
  ellipse((pr_ro_reg_h2r_wy(m,3)-pr_ro_reg_h2r_wy(m,2))/2,((tas_ro_reg_h2r_wy(m,3)+100)-(tas_ro_reg_h2r_wy(m,2)+100))/2,0,pr_ro_reg_h2r_wy(m,1),tas_ro_reg_h2r_wy(m,1),'b')
end
h = ellipse((pr_ro_reg_obs_wy(3)-pr_ro_reg_obs_wy(2))/2,((tas_ro_reg_obs_wy(3)+100)-(tas_ro_reg_obs_wy(2)+100))/2,0,pr_ro_reg_obs_wy(1),tas_ro_reg_obs_wy(1))
set(h,'Color',[1 0 0])
plot(pr_ro_reg_h2r_wy(:,1),tas_ro_reg_h2r_wy(:,1),'bo','MarkerFaceColor','b')
plot(pr_ro_reg_obs_wy(1),tas_ro_reg_obs_wy(1),'ro','MarkerFaceColor','r')
plot([pr_ro_reg_obs_wy(2) pr_ro_reg_obs_wy(3)],[tas_ro_reg_obs_wy(1) tas_ro_reg_obs_wy(1)],'r')
plot([pr_ro_reg_obs_wy(1) pr_ro_reg_obs_wy(1)],[tas_ro_reg_obs_wy(2) tas_ro_reg_obs_wy(3)],'r')
hline(tas_ro_reg_obs_wy(1),'r--')
vline(pr_ro_reg_obs_wy(1),'r--')
for m = 1:length(models_common)
  plot(pr_ro_reg_h2r_wy(m,1), tas_ro_reg_h2r_wy(m,1),'bo','MarkerSize',14,'MarkerFaceColor',[0 0 1])
  text(pr_ro_reg_h2r_wy(m,1), tas_ro_reg_h2r_wy(m,1),num2str(m), 'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',10,'Color',[1 1 1])
end
set(gca,'Layer','top','XLim',xlim,'YLim',ylim)
box on
hline(0,'k')
% vline(1,'k')
ylabel('\DeltaQ/\DeltaT (%/\circC)')
xlabel('\DeltaQ/\DeltaP (%/%)')

% return
set(gcf,'PaperPositionMode','auto');
fileo = [figpath 'Fig_3_' vari '_part2' co_tag];
print('-r300','-loose', '-depsc', [fileo '.eps'])
saveas(gcf,fileo,'jpg')
return
end













% ------------------------------------------------------------------------
tmp1 = pr_obs_wy_percent_mean > 100;
tmp2 = ro_obs_wy_percent_mean > 100;
tmp1c = zeros(size(tmp1));
tmp2c = zeros(size(tmp2));
for i = 1:length(tmp1)-1
  if tmp1(i) == tmp1(i+1)
    tmp1c(i+1) = tmp1c(i)+1;
  end
  if tmp2(i) == tmp2(i+1)
    tmp2c(i+1) = tmp2c(i)+1;
  end
end
close all
figure
hold on
plot(time_obs_common,tmp1c,'b')
plot(time_obs_common,tmp2c,'r--')
xlabel('Time (Year)')
ylabel('Number of consecutive same-sign anomalies (count)')
box on
legend('Precipitation','Streamflow','Location','NorthWest')
% return




if fig2 == 1
% ------------------------------------------------------------------------
% Figure 2? in paper
% ---------
close all

col_scen = [.6 .6 1; 1 .6 .6];
col_model = jet(length(models_common));
% col_model = parula(length(models_common));
% col_model = viridis(length(models_common));

col_model_light = (1-col_model)*.7+col_model;

xvar = 'pr'; % pr tas
yvar = 'ro'; % ro re
clear('var1','var2')

if strcmp(vari,'UC')==1
  magnif1 = 5e2;
  magnif2 = 7e2;
end
if strcmp(vari,'CO_DALLES')==1
  magnif1 = 1e4;
  magnif2 = 5e2;
end
if strcmp(vari,'NS')==1
  magnif1 = 2e2;
  magnif2 = 2e3;
end

ylim12     = [0 12];
xlim1      = [0 1.2];
xlim2      = [0 3.5];
binranges1 = [0:xlim1(2)/10:xlim1(2)];
binranges2 = [0:xlim2(2)/10:xlim2(2)];

figure1 = figure;
% set(figure1, 'units', 'centimeters', 'pos', [5 5 32 11])
set(figure1, 'units', 'centimeters', 'pos', [5 5 22 8])

subplot(1,2,1)
hold on
if strcmp(xvar,'pr')==1
  % xlabel('Precipitation (10^6 acre feet/year)')
  xlabel('Precipitation (km^3/year)')
  if strcmp(vari,'UC')==1
    xlim = [0 300]*tf;
  elseif strcmp(vari,'CO_DALLES')==1
    xlim = [0 900]*tf;
  elseif strcmp(vari,'NS')==1
    xlim = [0 100]*tf;
  end
elseif strcmp(xvar,'tas')==1
  xlabel('Temperature (\circC)')
end
if strcmp(yvar,'ro')==1
  % ylabel('Runoff (10^6 acre feet/year)')
  ylabel('Runoff (km^3/year)')
  if strcmp(vari,'UC')==1
    ylim = [0 120]*tf;
  elseif strcmp(vari,'CO_DALLES')==1
    ylim = [0 550]*tf;
  elseif strcmp(vari,'NS')==1
    ylim = [0 70]*tf;
  end
elseif strcmp(yvar,'re')==1
  ylabel('Runoff efficiency (%)')
end
% plot([xlim(1) xlim(2)],[xlim(1) xlim(2)],'m','LineWidth',2) % 1:1 line
for m = 1:length(models_common)
  for s = 3:3
    if strcmp(xvar,'tas')==1
      % var1_abs(m,:) = eval(['squeeze(' xvar '_wy(m,s,:))-mean(squeeze(' xvar '_wy(m,s,:)))']);
      var1_abs(m,:) = eval([xvar '_h2r_wy(m,1:59)-mean(' xvar '_h2r_wy(m,1:59))']);
    else
      % var1_abs(m,:) = eval(['squeeze(' xvar '_wy(m,s,:))']);
      var1_abs(m,:) = eval([xvar '_h2r_wy(m,1:59)']);
    end
    var2_abs(m,:) = eval([yvar '_h2r_wy(m,1:59)']);
    plot(var1_abs(m,:),var2_abs(m,:),'o','Color',col_model_light(m,:))
  end
end
for m = 1:length(models_common)
  beta_abs_model(m,:) = polyfit(var1_abs(m,:)',var2_abs(m,:)',1);
  h1(m) = plot(sort(var1_abs(m,:)), polyval(polyfit(var1_abs(m,:)',var2_abs(m,:)',1),sort(var1_abs(m,:))), 'Color',col_model(m,:)*0.8,'LineWidth',2);
end
var1_abs_obs = eval([ xvar '_obs_wy']);
var2_abs_obs = eval([ yvar '_obs_wy']);
plot(var1_abs_obs,var2_abs_obs, 'o','Color',[.4 .4 .4])
beta_abs_obs = polyfit(var1_abs_obs,var2_abs_obs',1);
h2 = plot(sort(var1_abs_obs), polyval(polyfit(var1_abs_obs,var2_abs_obs',1),sort(var1_abs_obs)), 'k','LineWidth',3)
box on
title(['(a) Absolute values ' num2str(time_obs_common(1)) '-' num2str(time_obs_common(end))])
set(gca,'Layer','top','XLim',xlim,'YLim',ylim)
if strcmp(vari,'UC')==1
  % legend([h2 h1(1):h1(end)],['Observations' models_common],'Location','NorthWest','FontSize',5)
end
% -- PDFs at side:
for m = 1:length(models_common)
  [f0,x0,u] = ksdensity( var1_abs(m,:) );
  line(x0, f0*magnif1, 'Color',col_model_light(m,:), 'LineWidth',1);
end
[f0,x0,u] = ksdensity( var1_abs_obs );
line(x0, f0*magnif1, 'Color',[.4 .4 .4], 'LineWidth',2);
for m = 1:length(models_common)
  [f0,x0,u] = ksdensity( var2_abs(m,:) );
  line(f0*magnif1, x0, 'Color',col_model_light(m,:), 'LineWidth',1);
end
[f0,x0,u] = ksdensity( var2_abs_obs );
line(f0*magnif1, x0, 'Color',[.4 .4 .4], 'LineWidth',2);
% -- inset histogram:
axes('Position',[.37 .24 .08 .15])
hold on
[n1,x1] = hist(beta_abs_model(:,1),binranges1);
h1 = bar(x1,n1,'hist');
set(h1,'EdgeColor',[.4 .4 .4],'FaceColor','none')
plot([beta_abs_obs(1) beta_abs_obs(1)],ylim12,'k','LineWidth',2)
set(gca,'YLim',ylim12,'XLim',xlim1,'ycolor',[1 1 1],'xtick',[0:.5:3],'FontSize',8)

subplot(1,2,2)
hold on
if strcmp(xvar,'pr')==1
  xlabel('Precipitation (%)')
  if strcmp(vari,'UC')==1
    xlim = [-50 70];
  elseif strcmp(vari,'CO_DALLES')==1
    xlim = [-50 50];
  elseif strcmp(vari,'NS')==1
    xlim = [-100 120];
  end
elseif strcmp(xvar,'tas')==1
  xlabel('Temperature (\circC)')
end
if strcmp(yvar,'ro')==1
  ylabel('Runoff (%)')
  if strcmp(vari,'UC')==1
    ylim = [-130 200];
  elseif strcmp(vari,'CO_DALLES')==1
    ylim = [-80 80];
  elseif strcmp(vari,'NS')==1
    ylim = [-180 300];
  end
elseif strcmp(yvar,'re')==1
  ylabel('Runoff efficiency (%)')
end
% plot([xlim(1) xlim(2)],[xlim(1) xlim(2)],'m','LineWidth',2) % 1:1 line
for m = 1:length(models_common)
  for s = 3:3
    if strcmp(xvar,'tas')==1
      var1_rel(m,:) = eval([xvar '_h2r_wy(m,1:59)-mean(' xvar '_h2r_wy(m,1:59))']);
    else
      var1_rel(m,:) = eval([xvar '_h2r_wy_percent_mean(m,1:59)-100;']);
    end
    if strcmp(yvar,'ro')==1
      var2_rel(m,:) = eval([yvar '_h2r_wy_percent_mean(m,1:59)-100;']);
    elseif strcmp(yvar,'re')==1
      var2_rel(m,:) = eval([yvar '_h2r_wy(m,1:59)) - nanmean(' yvar '_h2r_wy(m,1:59)));']);
    end
    plot(var1_rel(m,:),var2_rel(m,:), 'o','Color',col_model_light(m,:))
  end
end
for m = 1:length(models_common)
  beta_rel_model(m,:) = polyfit(var1_rel(m,:)',var2_rel(m,:)',1);
  plot(sort(var1_rel(m,:)), polyval(polyfit(var1_rel(m,:)',var2_rel(m,:)',1),sort(var1_rel(m,:))), 'Color',col_model(m,:)*0.8,'LineWidth',2)
end
if strcmp(xvar,'tas')==1
  var1_rel_obs = eval([xvar '_obs_wy-mean(' xvar '_obs_wy);']);
else
  var1_rel_obs = eval([xvar '_obs_wy_percent_mean-100;']);
end
if strcmp(yvar,'ro')==1
  var2_rel_obs = eval([yvar '_obs_wy_percent_mean-100;']);
elseif strcmp(yvar,'re')==1
  var2_rel_obs = eval([yvar '_obs_wy-mean(' yvar '_obs_wy);']);
end
plot(var1_rel_obs,var2_rel_obs, 'o','Color',[.4 .4 .4])
beta_rel_obs = polyfit(var1_rel_obs',var2_rel_obs',1);
plot(sort(var1_rel_obs), polyval(polyfit(var1_rel_obs',var2_rel_obs',1),sort(var1_rel_obs)), 'k','LineWidth',3)
box on
title(['(b) Relative anomalies ' num2str(time_obs_common(1)) '-' num2str(time_obs_common(end))])
set(gca,'Layer','top','XLim',xlim,'YLim',ylim)
% -- PDFs at side:
for m = 1:length(models_common)
  [f0,x0,u] = ksdensity( var1_rel(m,:) );
  line(x0, f0*magnif2+ylim(1), 'Color',col_model_light(m,:), 'LineWidth',1);
end
[f0,x0,u] = ksdensity( var1_rel_obs );
line(x0, f0*magnif2+ylim(1), 'Color',[.4 .4 .4], 'LineWidth',2);
for m = 1:length(models_common)
  [f0,x0,u] = ksdensity( var2_rel(m,:) );
  line(f0*magnif2+xlim(1), x0, 'Color',col_model_light(m,:), 'LineWidth',1);
end
[f0,x0,u] = ksdensity( var2_rel_obs );
line(f0*magnif2+xlim(1), x0, 'Color',[.4 .4 .4], 'LineWidth',2);
% -- inset histogram:
axes('Position',[.81 .24 .08 .15])
hold on
[n1,x1] = hist(beta_rel_model(:,1),binranges2);
h1 = bar(x1,n1,'hist');
set(h1,'EdgeColor',[.4 .4 .4],'FaceColor','none')
plot([beta_rel_obs(1) beta_rel_obs(1)],ylim12,'k','LineWidth',2)
set(gca,'YLim',ylim12,'XLim',xlim2,'ycolor',[1 1 1],'xtick',[0:1:3],'FontSize',8)
set(gcf, 'InvertHardcopy', 'off')



% return
set(gcf,'PaperPositionMode','auto');
fileo = [figpath 'Fig2_' xvar '_vs_' yvar '_' vari];
print('-r300','-loose', '-depsc', [fileo '.eps'])
saveas(gcf,fileo,'jpg')
return



close all
figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [5 5 35 15])

subplot(2,4,1)
hold on
plot(mean(tas_h2r_wy(:,1:59),2),pr_ro_reg_h2r_wy(:,1),'o')
hline(pr_ro_reg_obs_wy(1),'r')
vline(mean(tas_obs_wy),'r')
xlabel('T')
ylabel('P sensitivity')
subplot(2,4,2)
hold on
plot(mean(re_h2r_wy(:,1:59),2),pr_ro_reg_h2r_wy(:,1),'o')
hline(pr_ro_reg_obs_wy(1),'r')
vline(mean(re_obs_wy),'r')
xlabel('RE (%)')
ylabel('P sensitivity')
subplot(2,4,3)
hold on
plot(mean(pr_h2r_wy(:,1:59),2),pr_ro_reg_h2r_wy(:,1),'o')
hline(pr_ro_reg_obs_wy(1),'r')
vline(mean(pr_obs_wy),'r')
xlabel('P')
ylabel('P sensitivity')
subplot(2,4,4)
hold on
plot(mean(ro_h2r_wy(:,1:59),2),pr_ro_reg_h2r_wy(:,1),'o')
hline(pr_ro_reg_obs_wy(1),'r')
vline(mean(ro_obs_wy),'r')
xlabel('R')
ylabel('P sensitivity')


subplot(2,4,5)
hold on
plot(mean(tas_h2r_wy(:,1:59),2),tas_ro_reg_h2r_wy(:,1),'o')
hline(tas_ro_reg_obs_wy(1),'r')
vline(mean(tas_obs_wy),'r')
xlabel('T')
ylabel('T sensitivity')
subplot(2,4,6)
hold on
plot(mean(re_h2r_wy(:,1:59),2),tas_ro_reg_h2r_wy(:,1),'o')
hline(tas_ro_reg_obs_wy(1),'r')
vline(mean(re_obs_wy),'r')
xlabel('RE (%)')
ylabel('T sensitivity')
subplot(2,4,7)
hold on
plot(mean(pr_h2r_wy(:,1:59),2),tas_ro_reg_h2r_wy(:,1),'o')
hline(tas_ro_reg_obs_wy(1),'r')
vline(mean(pr_obs_wy),'r')
xlabel('P')
ylabel('T sensitivity')
subplot(2,4,8)
hold on
plot(mean(ro_h2r_wy(:,1:59),2),tas_ro_reg_h2r_wy(:,1),'o')
hline(tas_ro_reg_obs_wy(1),'r')
vline(mean(ro_obs_wy),'r')
xlabel('R')
ylabel('T sensitivity')

['r(pr,ro) = ' num2str(corr(mean(pr_h2r_wy(:,1:59),2),mean(ro_h2r_wy(:,1:59),2)))]
['r(pr,re) = ' num2str(corr(mean(pr_h2r_wy(:,1:59),2),mean(re_h2r_wy(:,1:59),2)))]
['r(ro,re) = ' num2str(corr(mean(ro_h2r_wy(:,1:59),2),mean(re_h2r_wy(:,1:59),2)))]

['r(tas,ro) = ' num2str(corr(mean(tas_h2r_wy(:,1:59),2),mean(ro_h2r_wy(:,1:59),2)))]
['r(tas,pr) = ' num2str(corr(mean(tas_h2r_wy(:,1:59),2),mean(pr_h2r_wy(:,1:59),2)))]
['r(tas,re) = ' num2str(corr(mean(tas_h2r_wy(:,1:59),2),mean(re_h2r_wy(:,1:59),2)))]

['r(ro,pr_ro_reg) = ' num2str(corr(mean(ro_h2r_wy(:,1:59),2),pr_ro_reg_h2r_wy(:,1)))]
['r(pme,pr_ro_reg) = ' num2str(corr(mean(pme_h2r_wy(:,1:59),2),pr_ro_reg_h2r_wy(:,1)))]
['r(pr,pr_ro_reg) = ' num2str(corr(mean(pr_h2r_wy(:,1:59),2),pr_ro_reg_h2r_wy(:,1)))]
['r(tas,pr_ro_reg) = ' num2str(corr(mean(tas_h2r_wy(:,1:59),2),pr_ro_reg_h2r_wy(:,1)))]
['r(re,pr_ro_reg) = ' num2str(corr(mean(re_h2r_wy(:,1:59),2),pr_ro_reg_h2r_wy(:,1)))]
['r(pme,pme_change) = ' num2str(corr(mean(pme_h2r_wy(:,1:59),2), mean(pme_h2r_wy(:,101:end),2)-mean(pme_h2r_wy(:,1:59),2)))]

['r(ro,tas_ro_reg) = ' num2str(corr(mean(ro_h2r_wy(:,1:59),2),tas_ro_reg_h2r_wy(:,1)))]
['r(pme,tas_ro_reg) = ' num2str(corr(mean(pme_h2r_wy(:,1:59),2),tas_ro_reg_h2r_wy(:,1)))]
['r(pr,tas_ro_reg) = ' num2str(corr(mean(pr_h2r_wy(:,1:59),2),tas_ro_reg_h2r_wy(:,1)))]
['r(tas,tas_ro_reg) = ' num2str(corr(mean(tas_h2r_wy(:,1:59),2),tas_ro_reg_h2r_wy(:,1)))]
['r(re,tas_ro_reg) = ' num2str(corr(mean(re_h2r_wy(:,1:59),2),tas_ro_reg_h2r_wy(:,1)))]

['r(ro,tas_ro_reg+pr_ro_reg) = ' num2str(corr(mean(ro_h2r_wy(:,1:59),2), norml(tas_ro_reg_h2r_wy(:,1)).^2+norml(pr_ro_reg_h2r_wy(:,1)).^2))]
['r(pme,tas_ro_reg+pr_ro_reg) = ' num2str(corr(mean(pme_h2r_wy(:,1:59),2), norml(tas_ro_reg_h2r_wy(:,1)).^2+norml(pr_ro_reg_h2r_wy(:,1)).^2))]
['r(pr,tas_ro_reg+pr_ro_reg) = ' num2str(corr(mean(pr_h2r_wy(:,1:59),2), norml(tas_ro_reg_h2r_wy(:,1)).^2+norml(pr_ro_reg_h2r_wy(:,1)).^2))]
['r(tas,tas_ro_reg+pr_ro_reg) = ' num2str(corr(mean(tas_h2r_wy(:,1:59),2), norml(tas_ro_reg_h2r_wy(:,1)).^2+norml(pr_ro_reg_h2r_wy(:,1)).^2))]
['r(re,tas_ro_reg+pr_ro_reg) = ' num2str(corr(mean(re_h2r_wy(:,1:59),2), norml(tas_ro_reg_h2r_wy(:,1)).^2+norml(pr_ro_reg_h2r_wy(:,1)).^2))]

return

end














if figS3 == 1
  % ------------------------------------------------------------------------
  % Figure S3? for paper: Showing with the carbon cycle models that the
  % plant physiology effect doesn't matter much for the runoff sensitivities
  % the different experimental setups: 1=piControl, 2=1pctCO2, 3=esmFixClim1, 4=esmFdbk1
  % ---------
  close all

  cols  = [0 0 0; 0 0 0; 0 1 0; 1 0 0];
  vars  = {'tas','pr','et','ro','pme','re'};
  % units = {'\circC','MAF','MAF','MAF','MAF','%'};
  units = {'\circC','%','%','%','%','%-points'};
  mkrs  = {'o','>','p','d','s','v','*'};
  clear('tmp')

  % -- 56:85 is 30 years centered on CO2 doubling
  start0  = 56; % 56 110
  ende0   = 85; % 85 139


  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 35 35])

  m = 1:length(models_cc);
  for s = 2:4
    subplot(3,6,1+(s-2)*6)
    plot(squeeze(tas_cc_wy_anom(m,s,:))')
    title(vars(1))
    subplot(3,6,2+(s-2)*6)
    plot(squeeze(pr_cc_wy_percent_mean(m,s,:))')
    title(vars(2))
    subplot(3,6,3+(s-2)*6)
    plot(squeeze(et_cc_wy_percent_mean(m,s,:))')
    title(vars(3))
    subplot(3,6,4+(s-2)*6)
    plot(squeeze(ro_cc_wy_percent_mean(m,s,:))')
    title(vars(4))
    subplot(3,6,5+(s-2)*6)
    plot(squeeze(pme_cc_wy_percent_mean(m,s,:))')
    title(vars(5))
    subplot(3,6,6+(s-2)*6)
    plot(squeeze(re_cc_wy(m,s,:))')
    title(vars(6))
  end


  % -- calculate running window runoff sensitivities --
  wl = 50;
  for m = 1:length(models_cc)
    for s = 1:4
      for i = 1:length(ro_cc_wy_percent_mean)-wl+1
        [pr_ro_reg_cc_wy(m,s,i,1) pr_ro_reg_cc_wy(m,s,i,2:3) tas_ro_reg_cc_wy(m,s,i,1) tas_ro_reg_cc_wy(m,s,i,2:3)]...
         = runoff_sens_reg(squeeze(ro_cc_wy_percent_mean(m,s,i:i+wl-1))-100,squeeze(pr_cc_wy_percent_mean(m,s,i:i+wl-1))-100,squeeze(tas_cc_wy_anom(m,s,i:i+wl-1)));
      end
    end
  end

  weak_cols = [170 170 170;...
               136 180 252;...
               78 178 90;...
               242 121 123]/255;
  strong_cols = [0 0 0;...
              0 0 255;...
              31 155 3;...
              255 0 0]/255;

  % -- which models
  m = 1:length(models_cc);

  xlim = [20 120];
  ylim1 = [-0.5 4];
  ylim2 = [-20 20];
  ylim3 = [-20 20];
  ylim4 = [.5 3]; %[-15 20]; % P sens
  ylim5 = [-12 3]; %[-15 20]; % T sens


  % -- line width
  lw = [1 1.5 1 1];

  % -- running means:
  for ii = 1:length(m)
    jj = m(ii);
    for s = 1:4
      tas_cc_wy_anom_rm(ii,s,:)         = rm(squeeze(tas_cc_wy_anom(jj,s,:)),wl);
      ro_cc_wy_rm(ii,s,:)               = rm(squeeze(ro_cc_wy(jj,s,:)),wl);
      pr_cc_wy_percent_mean_rm(ii,s,:)  = rm(squeeze(pr_cc_wy_percent_mean(jj,s,:)),wl);
      ro_cc_wy_percent_mean_rm(ii,s,:)  = rm(squeeze(ro_cc_wy_percent_mean(jj,s,:)),wl);
      rt_cc_wy_rm(ii,s,:)               = rm(squeeze(rt_cc_wy(jj,s,:)),wl);
      pr_cc_wy_rm(ii,s,:)               = rm(squeeze(pr_cc_wy(jj,s,:)),wl);
      pme_cc_wy_percent_mean_rm(ii,s,:) = rm(squeeze(pme_cc_wy_percent_mean(jj,s,:)),wl);
    end
  end



  figure1 = figure;
  set(figure1, 'units', 'centimeters', 'pos', [10 10 32 6])

  idx = wl/2:length(ro_cc_wy_percent_mean)-(wl/2)-1;

  subplot(1,4,1)
  hold on
  x  = (wl/2)+1:length(ro_cc_wy_percent_mean)-(wl/2)+1;
  y1 = prctile(squeeze(tas_cc_wy_anom_rm(:,2,:)),25);
  y2 = prctile(squeeze(tas_cc_wy_anom_rm(:,2,:)),75);
  y1 = y1(~isnan(y1));
  y2 = y2(~isnan(y2));
  patch([x fliplr(x)],[y1 fliplr(y2)],weak_cols(2,:),'EdgeColor','none')
  y1 = prctile(squeeze(tas_cc_wy_anom_rm(:,3,:)),25);
  y2 = prctile(squeeze(tas_cc_wy_anom_rm(:,3,:)),75);
  y1 = y1(~isnan(y1));
  y2 = y2(~isnan(y2));
  patch([x fliplr(x)],[y1 fliplr(y2)],weak_cols(3,:),'EdgeColor','none')
  y1 = prctile(squeeze(tas_cc_wy_anom_rm(:,4,:)),25);
  y2 = prctile(squeeze(tas_cc_wy_anom_rm(:,4,:)),75);
  y1 = y1(~isnan(y1));
  y2 = y2(~isnan(y2));
  patch([x fliplr(x)],[y1 fliplr(y2)],weak_cols(4,:),'EdgeColor','none')
  aa = 1;
  for s = 2:4
    h1(aa) = plot(squeeze(nanmean(tas_cc_wy_anom_rm(:,s,:))),'Color',strong_cols(s,:),'LineWidth',lw(s))
    aa = aa+1;
  end
  legend([h1(1:3)],scen2_names(2:4),'Location','NorthWest')
  legend boxoff
  set(gca,'XLim',xlim,'YLim',ylim1)
  hline(0,'k')
  box on
  ylabel('dT (\circC)')
  xlabel('Time [Simulation year]')

  subplot(1,4,2)
  hold on
  x  = (wl/2)+1:length(ro_cc_wy_percent_mean)-(wl/2)+1;
  y1 = prctile(squeeze(pr_cc_wy_percent_mean_rm(:,2,:))-100,25);
  y2 = prctile(squeeze(pr_cc_wy_percent_mean_rm(:,2,:))-100,75);
  y1 = y1(~isnan(y1));
  y2 = y2(~isnan(y2));
  patch([x fliplr(x)],[y1 fliplr(y2)],weak_cols(2,:),'EdgeColor','none')
  y1 = prctile(squeeze(pr_cc_wy_percent_mean_rm(:,3,:))-100,25);
  y2 = prctile(squeeze(pr_cc_wy_percent_mean_rm(:,3,:))-100,75);
  y1 = y1(~isnan(y1));
  y2 = y2(~isnan(y2));
  patch([x fliplr(x)],[y1 fliplr(y2)],weak_cols(3,:),'EdgeColor','none')
  y1 = prctile(squeeze(pr_cc_wy_percent_mean_rm(:,4,:))-100,25);
  y2 = prctile(squeeze(pr_cc_wy_percent_mean_rm(:,4,:))-100,75);
  y1 = y1(~isnan(y1));
  y2 = y2(~isnan(y2));
  patch([x fliplr(x)],[y1 fliplr(y2)],weak_cols(4,:),'EdgeColor','none')
  for s = 2:4
    plot(squeeze(nanmean(pr_cc_wy_percent_mean_rm(:,s,:)))-100,'Color',strong_cols(s,:),'LineWidth',lw(s))
  end
  set(gca,'XLim',xlim,'YLim',ylim2)
  hline(0,'k')
  box on
  ylabel('dP (%)')
  xlabel('Time [Simulation year]')

  subplot(1,4,3)
  hold on
  for s = 2:4
    x  = (wl/2)+1:length(ro_cc_wy_percent_mean)-(wl/2)+1;
    y1 = prctile(squeeze(ro_cc_wy_percent_mean_rm(:,s,:))-100,20);
    y2 = prctile(squeeze(ro_cc_wy_percent_mean_rm(:,s,:))-100,80);
    idx = find(~isnan(y1));
    patch([x fliplr(x)],[y1(idx) fliplr(y2(idx))],weak_cols(s,:),'EdgeColor','none')
  end
  for s = 2:4
    plot(squeeze(nanmean(ro_cc_wy_percent_mean_rm(:,s,:)))-100,'Color',strong_cols(s,:),'LineWidth',lw(s))
  end
  set(gca,'XLim',xlim,'YLim',ylim3)
  hline(0,'k')
  box on
  ylabel('dR (%)')
  xlabel('Time [Simulation year]')


  % -- calculate predicted runoff changes (transiently) based on runoff sensitivities --
  for mi = 1:length(m)
    for s = 2:4 % scenario
      % -- predicted R change based on period yi:
      yi = [1 10];
      ro_cc_wy_percent_mean_rm_predicted(mi,s,:) = (squeeze(pr_cc_wy_percent_mean_rm(mi,s,idx))-100) .* squeeze(pr_ro_reg_cc_wy(mi,s,:,1)) + (squeeze(tas_cc_wy_anom_rm(mi,s,idx))) .* squeeze(tas_ro_reg_cc_wy(mi,s,:,1));
      ro_cc_wy_percent_mean_rm_predicted_fixed(mi,s,:) = (squeeze(pr_cc_wy_percent_mean_rm(mi,s,idx))-100) .* nanmean(squeeze(pr_ro_reg_cc_wy(mi,s,yi,1))) + squeeze(tas_cc_wy_anom_rm(mi,s,idx)) .* nanmean(squeeze(tas_ro_reg_cc_wy(mi,s,yi,1)));
    end
  end

  subplot(1,4,4)
  hold on
  for s = [3 4 2] % 2:4
    tmp0 = squeeze(nanmean(ro_cc_wy_percent_mean_rm(:,s,:)))-100;
    tmp1 = squeeze(nanmean(ro_cc_wy_percent_mean_rm_predicted(:,s,:),1));
    tmp2 = squeeze(nanmean(ro_cc_wy_percent_mean_rm_predicted_fixed(:,s,:),1));
    y1 = tmp1'*0;
    y2 = (tmp0(idx)-tmp2)';
    y1 = y1(~isnan(y1));
    y2 = y2(~isnan(y2));
    x  = idx(~isnan(y1));
    y3 = tmp1'*0;
    y4 = (tmp0(idx)-tmp1)';
    y3 = y3(~isnan(y3));
    y4 = y4(~isnan(y4));
    h4 = patch([x+1 fliplr(x+1)],[abs(y2)-abs(y4) fliplr(abs(y2)+abs(y4))],weak_cols(s,:),'EdgeColor','none')
  end
  for s = [3 4 2]
    tmp0 = squeeze(nanmean(ro_cc_wy_percent_mean_rm(:,s,:)))-100;
    tmp1 = squeeze(nanmean(ro_cc_wy_percent_mean_rm_predicted(:,s,:),1));
    tmp2 = squeeze(nanmean(ro_cc_wy_percent_mean_rm_predicted_fixed(:,s,:),1));
    y1 = tmp1'*0;
    y2 = (tmp0(idx)-tmp2)';
    y1 = y1(~isnan(y1));
    y2 = y2(~isnan(y2));
    x  = idx(~isnan(y1));
    plot(x+1,abs(y2),'Color',strong_cols(s,:))
  end
  set(gca,'XLim',xlim,'YLim',ylim3)
  hline(0,'k')
  box on
  ylabel('dR (%)')
  xlabel('Time [Simulation year]')


  % return
  set(gcf,'PaperPositionMode','auto');
  fileo = [figpath 'Fig_S3_new_' vari];
  print('-r300','-loose', '-depsc', [fileo '.eps'])
  saveas(gcf,fileo,'jpg')
  return
end
