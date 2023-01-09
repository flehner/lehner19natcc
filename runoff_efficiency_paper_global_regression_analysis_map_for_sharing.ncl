load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; ==============================================================================
;
; runoff_efficiency_paper_global_regression_analysis_map_for_sharing.ncl
; author: flehner@ucar.edu
;
; Script to produce Fig. S2 in Lehner et al. 2019.
;
; Full citation:
; Lehner, F., A. W. Wood, J. A. Vano, D. M. Lawrence, M. P. Clark, J. S. Mankin (2019):
; The potential to reduce uncertainty in regional runoff projections from climate models
; Nature Climate Change, DOI: 10.1038/s41558-019-0639-x
; https://www.nature.com/articles/s41558-019-0639-x
;
;
; Note: if you have obtained the pre-processed data and adjusted the paths to your
; environment, all "recalc..." parameters can be set to "no" for it to
; run through and produce the figure w/o redoing some of the calculations.
;
; ==============================================================================

system("echo  ")

; -- recalc time series
recalc1 = "no"
; -- recalc sensitivities
recalc2 = "no"
; -- recalc prediction
recalc3 = "no"


seas    = "wateryear" ; annual DJF JJA

function wateryear_sum(infield)
begin
  dims  = dimsizes(infield)
  ntime = dims(0)
  nlat  = dims(1)
  nlon  = dims(2)
  tmp1  = new((/(ntime/12)-1,nlat,nlon/),float)
  do i = 0,dimsizes(tmp1(:,0,0))-1
      tmp1(i,:,:) = tofloat(dim_sum_n_Wrap(infield((i*12)+9:(i*12)+21,:,:),0))
  end do
  return(tmp1)
end

function wateryear_avg(infield)
begin
  dims  = dimsizes(infield)
  ntime = dims(0)
  nlat  = dims(1)
  nlon  = dims(2)
  tmp1  = new((/(ntime/12)-1,nlat,nlon/),float)
  do i = 0,dimsizes(tmp1(:,0,0))-1
      tmp1(i,:,:) = tofloat(dim_avg_n_Wrap(infield((i*12)+9:(i*12)+21,:,:),0))
  end do
  return(tmp1)
end

; -- paths
pathin = "/project/cas/flehner/cmip5_ncar/"
pathout = pathin

; -- CMIP5 models
models_common = (/"CCSM4","CESM1-BGC","CESM1-CAM5","CNRM-CM5","CanESM2","GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H","IPSL-CM5A-MR","MIROC-ESM","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME","bcc-csm1-1","bcc-csm1-1-m","inmcm4"/)
nmodels = dimsizes(models_common)

; -- load dummy file
ifile1 = pathin+"historical/Amon/tas/"+models_common(0)+"/r1i1p1/tas_remapcon2_1x1.nc"
a = addfile(ifile1(0),"r")
dummy   = a->tas
lat     = a->lat

; -- river basins
ifile1 = pathin+"major_basins_of_the_world_area_for_cmip5_1x1.nc"
a = addfile(ifile1(0),"r")
tmp     = lonFlip(a->area)
basins  = where(tmp .gt. 0,1,dummy@_FillValue)
nbasins = dimsizes(basins(:,0,0))
delete(tmp)


; -- weights
rad    = 4.0*atan(1.0)/180.0
clat   = cos(lat*rad)
clat!0 = "lat"
clat&lat = lat


; -- load CMIP5 data --
ende_hist   = 200512
start_year  = 1950
ende_year   = 2100
start       = start_year+"01"
ende        = ende_year+"12"
refstart    = 0; 195001
refende     = 58; 200812

; -- CO2 vmr
ifile1 = "/project/cas/flehner/LE/co2vmr.hist2rcp85.185001-208012.nc"
a = addfile(ifile1(0),"r")
tmp         = a->CO2VMR
time_co2    = cd_calendar(tmp&TIME,-1)
ind_start   = ind(time_co2 .eq. start)
co2         = month_to_annual(tmp(ind_start+1::),1)

if recalc1 .eq. "yes"
  ; -- create arrays
  tas_ts  = new((/nmodels,nbasins,ende_year-start_year/),float)
  pr_ts   = new((/nmodels,nbasins,ende_year-start_year/),float)
  mrro_ts = new((/nmodels,nbasins,ende_year-start_year/),float)

  do m = 0,nmodels-1
    print_clock("model "+(m+1)+"/"+nmodels+": "+models_common(m))
    print_clock("  tas"); ------------------
    var = "tas"
    comp = "Amon"
    ; ifile1 = systemfunc("ls -1 "+pathin+"historical/"+comp+"/"+var+"/"+models_common(m)+"/r1i1p1/"+var+"_"+comp+"_"+models_common(m)+"_historical_r1i1p1_*.nc_remapcon2_1x1.nc")
    ifile1 = pathin+"historical/"+comp+"/"+var+"/"+models_common(m)+"/r1i1p1/"+var+"_remapcon2_1x1.nc"
    a = addfile(ifile1(0),"r")
    tmp0    = a->$var$
    time0   = cd_calendar(tmp0&time,-1)
    ind0    = ind(time0 .eq. ende_hist)
    tmp1    = tmp0(0:ind0,:,:)
    delete([/tmp0,time0,ind0/])
    ifile1  = pathin+"rcp85/"+comp+"/"+var+"/"+models_common(m)+"/r1i1p1/"+var+"_remapcon2_1x1.nc"
    a = addfile(ifile1(0),"r")
    tmp2        = a->$var$
    tmp3        = array_append_record(tmp1,tmp2,0)
    time1       = cd_calendar(tmp1&time,-1)
    time2       = cd_calendar(tmp2&time,-1)
    cld         = array_append_record(time1,time2,0)
    ind_start   = ind(cld .eq. start)
    ind_ende    = ind(cld .eq. ende)
    tmp4        = tmp3(ind_start:ind_ende,:,:)
    tas         = wateryear_avg(tmp4)
    copy_VarCoords(tmp4(18::12,:,:),tas)
    ; tas_anom = (tas-conform_dims(dimsizes(tas),dim_avg_n_Wrap(tas,0),(/1,2/)))*100
    ; copy_VarCoords(tas,tas_anom)
    delete([/tmp1,tmp2,tmp3,tmp4,time1,time2,cld/])
    print_clock("  pr"); ------------------
    var = "pr"
    comp = "Amon"
    ifile1 = pathin+"historical/"+comp+"/"+var+"/"+models_common(m)+"/r1i1p1/"+var+"_remapcon2_1x1.nc"
    a = addfile(ifile1(0),"r")
    tmp0    = a->$var$
    time0   = cd_calendar(tmp0&time,-1)
    ind0    = ind(time0 .eq. ende_hist)
    tmp1    = tmp0(0:ind0,:,:)
    delete([/tmp0,time0,ind0/])
    ifile1 = pathin+"rcp85/"+comp+"/"+var+"/"+models_common(m)+"/r1i1p1/"+var+"_remapcon2_1x1.nc"
    a = addfile(ifile1(0),"r")
    tmp2        = a->$var$
    tmp3        = array_append_record(tmp1,tmp2,0)
    time1       = cd_calendar(tmp1&time,-1)
    time2       = cd_calendar(tmp2&time,-1)
    cld         = array_append_record(time1,time2,0)
    ind_start   = ind(cld .eq. start)
    ind_ende    = ind(cld .eq. ende)
    tmp4        = tmp3(ind_start:ind_ende,:,:)
    pr          = wateryear_avg(tmp4)
    copy_VarCoords(tmp4(18::12,:,:),pr)
    ; pr_percent_mean = (pr/conform_dims(dimsizes(pr),dim_avg_n_Wrap(pr,0),(/1,2/)))*100
    ; copy_VarCoords(pr,pr_percent_mean)
    delete([/tmp1,tmp2,tmp3,tmp4,time1,time2,cld/])
    print_clock("  mrro"); ------------------
    var = "mrro"
    comp = "Lmon"
    ifile1 = pathin+"historical/"+comp+"/"+var+"/"+models_common(m)+"/r1i1p1/"+var+"_remapcon2_1x1.nc"
    a = addfile(ifile1(0),"r")
    tmp0    = a->$var$
    time0   = cd_calendar(tmp0&time,-1)
    ind0    = ind(time0 .eq. ende_hist)
    tmp1    = tmp0(0:ind0,:,:)
    delete([/tmp0,time0,ind0/])
    ifile1 = pathin+"rcp85/"+comp+"/"+var+"/"+models_common(m)+"/r1i1p1/"+var+"_remapcon2_1x1.nc"
    a = addfile(ifile1(0),"r")
    tmp2        = a->$var$
    tmp3        = array_append_record(tmp1,tmp2,0)
    time1       = cd_calendar(tmp1&time,-1)
    time2       = cd_calendar(tmp2&time,-1)
    cld         = array_append_record(time1,time2,0)
    ind_start   = ind(cld .eq. start)
    ind_ende    = ind(cld .eq. ende)
    tmp4        = tmp3(ind_start:ind_ende,:,:)
    mrro        = wateryear_avg(tmp4)
    copy_VarCoords(tmp4(18::12,:,:),mrro)
    ; mrro_percent_mean = (mrro/conform_dims(dimsizes(mrro),dim_avg_n_Wrap(mrro,0),(/1,2/)))*100
    ; copy_VarCoords(mrro,mrro_percent_mean)
    delete([/tmp1,tmp2,tmp3,tmp4,time1,time2,cld/])
    ; -- basins --
    print_clock("  loop over basins...")
    do r = 0,nbasins-1
      ; print_clock("  r = "+r)
      tmp1                = conform_dims(dimsizes(tas), basins(r,:,:), (/1,2/))
      tmp2                = tmp1 * tas
      tas_ts(m,r,:)       = tofloat(wgt_areaave(tmp2, clat, 1.0, 0))
      tmp2                = tmp1 * pr
      pr_ts(m,r,:)        = tofloat(wgt_areaave(tmp2, clat, 1.0, 0))
      tmp2                = tmp1 * mrro
      mrro_ts(m,r,:)      = tofloat(wgt_areaave(tmp2, clat, 1.0, 0))
      delete([/tmp1,tmp2/])
    end do ; end of basin loop
    delete([/pr,tas,mrro/])
  end do ; end of model loop
  print_clock("  save data...")
  fileout = pathout+"major_basins_tas_pr_mrro_ts.nc"
  system("rm -rf "+fileout)
  out = addfile(fileout,"c")
  out->tas_ts 	    = tas_ts
  out->pr_ts 	      = pr_ts
  out->mrro_ts 	    = mrro_ts
else
  print_clock("  load data...")
  fileout = pathout+"major_basins_tas_pr_mrro_ts.nc"
  a = addfile(fileout,"r")
  tas_ts 	      = a->tas_ts
  pr_ts 	      = a->pr_ts
  mrro_ts 	    = a->mrro_ts
end if



if recalc2 .eq. "yes"
  print_clock("  define sensitivities"); -----------------
  beta1 = new((/nmodels,nbasins,3/),float)
  beta2 = new((/nmodels,nbasins,4/),float)
  do m = 0,nmodels-1
    do b = 0,nbasins-1
      if ismissing(mrro_ts(m,b,0)) .or. mrro_ts(m,b,0).eq.0
        beta1(m,b,:) = mrro_ts@_FillValue
        beta2(m,b,:) = mrro_ts@_FillValue
      else
        mrro_percent_mean = (mrro_ts(m,b,refstart:refende)/dim_avg_n_Wrap(mrro_ts(m,b,refstart:refende),0))*100-100
        pr_percent_mean = (pr_ts(m,b,refstart:refende)/dim_avg_n_Wrap(pr_ts(m,b,refstart:refende),0))*100-100
        tas_anom = tas_ts(m,b,refstart:refende)-dim_avg_n_Wrap(tas_ts(m,b,refstart:refende),0)
        co2_anom = co2(refstart:refende)-dim_avg_n_Wrap(co2(refstart:refende),0)
        ; -- w/o CO2
        X1 = new((/3,(refende-refstart+1)/),float)
        X1(0,:) = pr_percent_mean
        X1(1,:) = tas_anom
        X1(2,:) = pr_percent_mean*tas_anom
        Y      = mrro_percent_mean
        beta1(m,b,:) = reg_multlin(Y,X1,False)
        ; -- w/ CO2
        X2 = new((/4,(refende-refstart+1)/),float)
        X2(0,:) = pr_percent_mean
        X2(1,:) = tas_anom
        X2(2,:) = tofloat(co2_anom)
        X2(3,:) = pr_percent_mean*tas_anom
        beta2(m,b,:) = reg_multlin(Y,X2,False)
      end if
    end do
  end do
  delete([/tas_anom,pr_percent_mean,mrro_percent_mean,co2_anom/])
  print_clock("  save data...")
  fileout = pathout+"beta.nc"
  system("rm -rf "+fileout)
  out = addfile(fileout,"c")
  out->beta1 	    = beta1
  out->beta2 	    = beta2
else
  print_clock("  load data...")
  fileout = pathout+"beta.nc"
  a = addfile(fileout,"r")
  beta1 	    = a->beta1
  beta2 	    = a->beta2
end  if


; -- find 2C warmer period and predict R change
futtas = 2
rmwl = 40
if recalc3 .eq. "yes"
  print_clock("  predict future runoff"); -----------------
  pr_delta         = new((/nmodels,nbasins/),float)
  ro_delta         = new((/nmodels,nbasins/),float)
  co2_delta        = new((/nmodels,nbasins/),float)
  ro_delta_predict1= new((/nmodels,nbasins/),float)
  ro_delta_predict2= new((/nmodels,nbasins/),float)
  corr1            = new((/nbasins/),float)
  corr2            = new((/nbasins/),float)
  do b = 0,nbasins-1
    do m = 0,nmodels-1
      if ismissing(mrro_ts(m,b,0)) .or. mrro_ts(m,b,0).eq.0
        ; print("basin with missing value")
      else
        tas_anom = tas_ts(m,b,:)-dim_avg_n_Wrap(tas_ts(m,b,refstart:refende),0)
        tas_anom_rm = runave_n_Wrap(tas_anom,rmwl,0,0)
        idx  = minind(abs(futtas-tas_anom_rm))
        pr_percent_mean = (pr_ts(m,b,:)/dim_avg_n_Wrap(pr_ts(m,b,refstart:refende),0))*100-100
        pr_percent_mean_rm = runave_n_Wrap(pr_percent_mean,rmwl,0,0)
        mrro_percent_mean = (mrro_ts(m,b,:)/dim_avg_n_Wrap(mrro_ts(m,b,refstart:refende),0))*100-100
        mrro_percent_mean_rm = runave_n_Wrap(mrro_percent_mean,rmwl,0,0)
        co2_anom = co2(refstart::)-dim_avg_n_Wrap(co2(refstart:refende),0)
        pr_delta(m,b) = pr_percent_mean_rm(idx)
        ro_delta(m,b) = mrro_percent_mean_rm(idx)
        co2_delta(m,b) = tofloat(co2_anom(idx))
        ro_delta_predict1(m,b) = beta1(m,b,0)*pr_delta(m,b) + beta1(m,b,1)*futtas
        ro_delta_predict2(m,b) = beta2(m,b,0)*pr_delta(m,b) + beta2(m,b,1)*futtas + beta2(m,b,2)*co2_delta(m,b)
      end if
    end do
    corr1(b) = abs(escorc(ro_delta(:,b),ro_delta_predict1(:,b)))
    corr2(b) = abs(escorc(ro_delta(:,b),ro_delta_predict2(:,b)))
  end do
  print_clock("  save data...")
  fileout = pathout+"prediction.nc"
  system("rm -rf "+fileout)
  out = addfile(fileout,"c")
  out->ro_delta 	       = ro_delta
  out->ro_delta_predict1 = ro_delta_predict1
  out->ro_delta_predict2 = ro_delta_predict2
  out->corr1              = corr1
  out->corr2              = corr2
else
  print_clock("  load data...")
  fileout = pathout+"prediction.nc"
  a = addfile(fileout,"r")
  ro_delta 	         = a->ro_delta
  ro_delta_predict1  = a->ro_delta_predict1
  ro_delta_predict2  = a->ro_delta_predict2
  corr1              = a->corr1
  corr2              = a->corr2
end  if

delete([/lat/])












; === PLOTTING ======================================================

ofile = "/home/flehner/publication/lehner18_runoff_efficiency_constraint/fig/"+get_script_prefix_name()
wks  = gsn_open_wks("eps",ofile) ; send graphics to PNG file

res                     = True

res@gsnDraw             = False       ; don't draw yet
res@gsnFrame            = False       ; don't advance frame yet
res@gsnMaximize         = True        ; maximize plot in frame

res@mpProjection        = "Robinson"    ; choose projection
res@mpLimitMode         = "LatLon"
res@mpMinLonF           = -180.
res@mpMaxLonF           = 180.
res@mpMinLatF           = -60.
res@mpMaxLatF           = 90.
; res@mpCenterLonF        = (res@mpMinLonF+res@mpMaxLonF)/2.0

; res@pmTickMarkDisplayMode = "Always"  ; turn on tickmarks

; res@cnLineLabelsOn      = True

res@tiMainString        = "Prediction of ~F33~D~F~R from runoff sensitivities"

plot = gsn_csm_map(wks,res)   ; Create map, but don't draw it yet.

; -- create colormap --
cmap = gsn_retrieve_colormap(wks)
hsv = rgbhsv(cmap)
hsv(10,1) = .47
hsv(8,1) = .45    ; 1c+1d
hsv(9,1) = .30    ; 1d
hsv(4,:) = (/ (hsv(3,:)+hsv(5,:))/2.0 /)
cmap2 = hsvrgb(hsv)
; gsn_define_colormap(wks,cmap2)
; -- merge colorbars/colormaps!
gsn_merge_colormaps(wks,cmap2,"precip_diff_12lev")


; -- load shapefile --
; load "/home/flehner/scripts/ncl/shapefile_utils.ncl"
shp_filename = "/project/cas/flehner/shapefiles/major_basins_of_the_world_0_0_0/Major_Basins_of_the_World.shp"
f = addfile(shp_filename, "r")   ; Open shapefile

;
; Read data off shapefile
;
segments = f->segments
geometry = f->geometry
segsDims = dimsizes(segments)
geomDims = dimsizes(geometry)
;
; Read global attributes
;
geom_segIndex = f@geom_segIndex
geom_numSegs  = f@geom_numSegs
segs_xyzIndex = f@segs_xyzIndex
segs_numPnts  = f@segs_numPnts

lines       = new(segsDims(0),graphic)   ; Array to hold polygons
numFeatures = geomDims(0)

data  = corr1^2 ; corr1 = w/o CO2 term; corr2 = w/ CO2 term
printVarSummary(data)
; print(data)
; exit
lon   = f->x
lat   = f->y

plbres             = True       ; resources for polylines
plbres@gsEdgesOn   = True       ; draw border around polygons
plbres@gsEdgeColor = "black"
plbres@cnLineLabelsOn      = True

; -- for correlation:
; levels = (/.38,.5,.7,.9/)
; labels = (/"0.38","0.5","0.7","0.9"/)
; -- for R^2:
levels = (/.14,.25,.5,.8/)
labels = (/"0.14","0.25","0.5","0.8"/)
nlevels = dimsizes(levels)
; colors = (/"blue","green","yellow","red"/)
colors = (/"(/1,1,1/)","(/1,.9,0/)","(/1,.7,0/)","(/1,.3,0/)","(/.85,0,0/)"/)

segNum = 0
do i=0, numFeatures-1

  if ismissing(data(i))
    plbres@gsFillColor = "black"
  else
    ; -- color assignment (probably a better way to do this?)
    ; if (data(i).ge.levels(0) .and. data(i).lt. levels(1)) then
    if (data(i).lt. levels(0)) then
       plbres@gsFillColor = colors(0)
    end if
    if (data(i).ge.levels(0) .and. data(i).lt. levels(1)) then
       plbres@gsFillColor = colors(1)
    end if
    if (data(i).ge.levels(1) .and. data(i).lt. levels(2)) then
       plbres@gsFillColor = colors(2)
    end if
    if (data(i).ge.levels(2) .and. data(i).lt. levels(3)) then
       plbres@gsFillColor = colors(3)
    end if
    if (data(i).ge.levels(3)) then
       plbres@gsFillColor = colors(4)
    end if
  end if

  startSegment = geometry(i, geom_segIndex)
  numSegments  = geometry(i, geom_numSegs)
  do seg=startSegment, startSegment+numSegments-1
    startPT = segments(seg, segs_xyzIndex)
    endPT = startPT + segments(seg, segs_numPnts) - 1
    lines(segNum) = gsn_add_polygon(wks, plot, lon(startPT:endPT),lat(startPT:endPT), plbres)
    segNum = segNum + 1
  end do
end do


lbres                    = True
lbres@vpWidthF           = 0.50             ; location
lbres@vpHeightF          = 0.07             ; " "
lbres@pmLabelBarOrthogonalPosF = -1.           ; move whole thing down
lbres@lbPerimOn          = False            ; Turn off perimeter.
lbres@lbOrientation      = "Horizontal"     ; Default is vertical.
lbres@lbLabelAlignment   = "InteriorEdges"     ; where to label
lbres@lbFillColors       = colors
lbres@lbMonoFillPattern  = True             ; Fill them all solid.
lbres@lbLabelFontHeightF = 0.016            ; label font height
lbres@lbTitleString      = "Correlation coefficient"        ; title
lbres@lbTitleString      = "R~S~2~N"        ; title
lbres@lbTitlePosition    = "Bottom"         ; location of title
lbres@lbTitleFontHeightF = 0.016            ; title font height


gsn_labelbar_ndc (wks,nlevels+1,labels,0.23,0.15,lbres)

; -- Maximize output in frame. This will draw everything and advance the frame
maximize_output(wks,False)


system("gv --orientation=seascape "+ofile+".eps &")
; system("gv "+ofile+".eps &")



end
