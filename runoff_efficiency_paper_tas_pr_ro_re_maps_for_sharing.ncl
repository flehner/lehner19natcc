; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ================================================;
begin

; ==============================================================================
;
; runoff_efficiency_paper_tas_pr_ro_re_maps_for_sharing.ncl
; author: flehner@ucar.edu
;
; Script to produce Fig. 1 in Lehner et al. 2019.
;
; Full citation:
; Lehner, F., A. W. Wood, J. A. Vano, D. M. Lawrence, M. P. Clark, J. S. Mankin (2019):
; The potential to reduce uncertainty in regional runoff projections from climate models
; Nature Climate Change, DOI: 10.1038/s41558-019-0639-x
; https://www.nature.com/articles/s41558-019-0639-x
;
;
; Note: this is a terribly long script with control parameters at the beginning.
; The "recalc[1..4]" parameters provide the option to redo certain calculation
; steps (which takes a while).
; If you have obtained the pre-processed data and adjusted the paths to your
; environment, all "recalc..." parameters can be set to "no" for it to
; run through and produce the figure.
;
; ==============================================================================

  print_clock("start script")

  recalc1 = "no" ; concatenate all models (historical and rcp85) and calculate wy mean values
  recalc2 = "no" ; calculate runoff efficiency on wy basis
  recalc3 = "no" ; calculate anomalies on reference period
  recalc4 = "no" ; calculate multi-model mean change, agreement on sign, and magnitude

  models = (/"bcc-csm1-1", "bcc-csm1-1-m", "BNU-ESM", "CanESM2", "CCSM4", "CESM1-BGC", "CESM1-CAM5", "CMCC-CESM", "CMCC-CM", "CMCC-CMS", "CNRM-CM5", "CSIRO-Mk3-6-0", "FGOALS-g2", "FIO-ESM", "GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M", "GISS-E2-H", "GISS-E2-H-CC", "GISS-E2-R", "GISS-E2-R-CC", "inmcm4", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "MIROC5", "MIROC-ESM", "MIROC-ESM-CHEM", "MPI-ESM-LR", "MPI-ESM-MR", "MRI-CGCM3", "MRI-ESM1", "NorESM1-M", "NorESM1-ME"/)

  pathin = "/project/cas/flehner/cmip5_ncar/"

  vars   = (/"mrro","tas","pr","re"/)

  ; -- fix parameters for file reading and concatenating --
  histstart   = 190101
  histende    = 200512
  rcpstart    = 200601
  rcpende     = 209912
  start       = round(histstart/100,3)
  ende        = round(rcpende/100,3)
  ; -- choose epochs for difference maps --
  refstart    = 1950; 1901 1929
  refende     = 2008; 1950 2008
  futtas_aim  = 2; 2 3 ; approx. warming in degC over the study basins
  if futtas_aim .eq. 2
    futstart    = 2000; 2001 2021 2021 2011
    futende     = 2058; 2050 2060
  else if futtas_aim .eq. 3
    futstart    = 2031; 2001 2021 2031 2050
    futende     = 2080; 2050 2050 2080 2099
  else if futtas_aim .eq. 4
    futstart    = 2050; 2001 2021 2031 2050
    futende     = 2099; 2050 2050 2080 2099
  end if
  end if
  end if

  region      = "north_america" ; north_america australia
  if region .eq. "north_america"
    minlat      = 22
    maxlat      = 74
    minlon      = -170 ; -145(default) -169
    maxlon      = -70
  end if


  if recalc1 .eq. "yes"
    print("recalc1")
    do v      = 0,dimsizes(vars)-2
      vari    = vars(v)
      print("variable = "+vari)

      if vari .eq. "mrro" .or. vari .eq. "re"
        comp      = "Lmon"
      else
        comp      = "Amon"
      end if

      ; -- create matrix
      ntime     = round(rcpende/100,3)-round(histstart/100,3)+1
      nlat      = 180
      nlon      = 360
      var_wy    = new((/ntime,dimsizes(models),nlat,nlon/),float)

      do m = 0,dimsizes(models)-1
        print("model = "+models(m))
        ifile = pathin+"historical/"+comp+"/"+vari+"/"+models(m)+"/r1i1p1/"+vari+"_remapcon2_1x1.nc"
        a     = addfile(ifile,"r")
        tmp0  = a->$vari$
        tim   = a->time
        cld   = cd_calendar(tim, -1)
        i_start = ind(cld .eq. histstart)
        i_ende = ind(cld .eq. histende)
        ; -- cut to common length
        tmp1 = tmp0(i_start:i_ende,:,:)
        delete([/tmp0,tim,cld,i_start,i_ende/])

        ifile = pathin+"rcp85/"+comp+"/"+vari+"/"+models(m)+"/r1i1p1/"+vari+"_remapcon2_1x1.nc"
        a     = addfile(ifile,"r")
        tmp0  = a->$vari$
        tim   = a->time
        cld   = cd_calendar(tim, -1)
        i_start = ind(cld .eq. rcpstart)
        i_ende = ind(cld .eq. rcpende)
        ; -- cut to common length
        tmp2 = tmp0(i_start:i_ende,:,:)
        tmp3 = array_append_record(tmp1,tmp2,0)
        ; -- make water year
        do i = 1,ntime-1
          print(""+i+"/"+ntime)
          var_wy(i,m,:,:) = tofloat( dim_avg_n_Wrap(tmp3(((i-1)*12)+9:(i*12)+8,:,:),0) )
        end do
        delete([/tmp0,tmp1,tmp2,tmp3,tim,cld,i_start,i_ende,ifile/])
      end do
      ; -- save data
      ifile = pathin+"historical/"+comp+"/"+vari+"/"+models(0)+"/r1i1p1/"+vari+"_remapcon2_1x1.nc"
      a     = addfile(ifile,"r")
      tmp0  = a->$vari$
      copy_VarCoords(tmp0(0,:,:),var_wy(0,0,:,:))
      var_wy!0 = "time"
      ofile   = pathin+"historical2rcp85/"+comp+"/"+vari+"/all_models_"+vari+"_remapcon2_1x1_wy.nc"
      system("rm -rf "+ofile)
      a         = addfile(ofile,"c")
      a->$vari$ = var_wy
    end do ; end of "vars" loop
    delete([/vari/])

  end if ; end of "recalc1" loop


  if recalc2 .eq. "yes"
    print("recalc2")
    ofile   = pathin+"historical2rcp85/Amon/pr/all_models_pr_remapcon2_1x1_wy.nc"
    a       = addfile(ofile,"r")
    pr      = a->pr
    ofile   = pathin+"historical2rcp85/Lmon/mrro/all_models_mrro_remapcon2_1x1_wy.nc"
    a       = addfile(ofile,"r")
    mrro    = a->mrro
    ; -- calculate runoff efficiency:
    pr0     = where(pr .eq. 0, pr@_FillValue, pr)
    re      = mrro/pr0
    copy_VarCoords(mrro,re)
    ofile   = pathin+"historical2rcp85/Lmon/re/all_models_re_remapcon2_1x1_wy.nc"
    system("rm -rf "+ofile)
    a         = addfile(ofile,"c")
    a->re = re
  end if



  if recalc3 .eq. "yes"
    print("recalc3")
    do v = 0,dimsizes(vars)-1
      vari          = vars(v)
      vari_anom     = vars(v)+"_anom"
      vari_anom_rel = vars(v)+"_anom_rel"
      print("variable = "+vari)
      if vari .eq. "mrro" .or. vari .eq. "re"
        comp      = "Lmon"
      else
        comp      = "Amon"
      end if
      ifile     = pathin+"historical2rcp85/"+comp+"/"+vari+"/all_models_"+vari+"_remapcon2_1x1_wy.nc"
      a         = addfile(ifile,"r")
      var_wy    = a->$vari$
      ; -- make anomalies to 1929-2008 mean --
      var_wy_anom = var_wy
      if vari .eq. "mrro" .or. vari .eq. "pr"
        var_wy_anom_rel = var_wy
      end if
      do m = 0,dimsizes(models)-1
        print("model = "+models(m))
        clim0 = dim_avg_n_Wrap(var_wy(refstart-start:refende-start,m,:,:),0)
        var_wy_anom(:,m,:,:) = var_wy(:,m,:,:)-conform_dims(dimsizes(var_wy(:,m,:,:)),clim0,(/1,2/))
        if vari .eq. "mrro" .or. vari .eq. "pr"
          clim1 = where(clim0 .eq. 0, clim0@_FillValue, clim0)
          var_wy_anom_rel(:,m,:,:) = var_wy_anom(:,m,:,:) / conform_dims(dimsizes(var_wy(:,m,:,:)),clim1,(/1,2/))
        end if
      end do
      copy_VarCoords(var_wy,var_wy_anom)
      ofile   = pathin+"historical2rcp85/"+comp+"/"+vari+"/all_models_"+vari+"_remapcon2_1x1_wy_anom"+refstart+"-"+refende+".nc"
      system("rm -rf "+ofile)
      a         = addfile(ofile,"c")
      a->$vari_anom$ = var_wy_anom
      if vari .eq. "mrro" .or. vari .eq. "pr"
        a->$vari_anom_rel$ = var_wy_anom_rel
      end if
    end do
  end if ; end of "recalc3" loop


  if recalc4 .eq. "yes"
    print("recalc4")
    do v = 0,dimsizes(vars)-1
      vari          = vars(v)
      vari_anom     = vars(v)+"_anom"
      vari_anom_rel = vars(v)+"_anom_rel"
      print("variable = "+vari)
      if vari .eq. "mrro" .or. vari .eq. "re"
        comp      = "Lmon"
      else
        comp      = "Amon"
      end if
      ofile   = pathin+"historical2rcp85/"+comp+"/"+vari+"/all_models_"+vari+"_remapcon2_1x1_wy.nc"
      a       = addfile(ofile,"r")
      tmp1    = a->$vari$

      ; -- multi-model mean (mmm) change
      mmm             = dim_avg_n_Wrap(tmp1,1)
      mean_change     = dim_avg_n_Wrap(tmp1(futstart-start:futende-start,:,:,:),0) - dim_avg_n_Wrap(tmp1(refstart-start:refende-start,:,:,:),0)
      mmm_change      = dim_avg_n_Wrap(mmm(futstart-start:futende-start,:,:),0) - dim_avg_n_Wrap(mmm(refstart-start:refende-start,:,:),0)
      if vari .eq. "mrro" .or. vari .eq. "pr"
        ; mean_change_rel = (mean_change/dim_avg_n_Wrap(tmp1(refstart-start:refende-start,:,:,:),0))*100
        mmm_change_rel  = (mmm_change/dim_avg_n_Wrap(mmm(refstart-start:refende-start,:,:),0))*100
      end if
      copy_VarCoords(tmp1(0,0,:,:),mean_change(0,:,:))
      copy_VarCoords(tmp1(0,0,:,:),mmm_change(:,:))
      if vari .eq. "mrro" .or. vari .eq. "pr"
        copy_VarCoords(tmp1(0,0,:,:),mmm_change_rel(:,:))
      end if
      ; -- agreement on sign (aos) of change - as fraction of models
      mean_change_binary = where(mean_change .gt. 0, 1, 0)
      mmm_change_binary = where(mmm_change .gt. 0, 1, 0)
      aos_tmp = where(mean_change_binary .eq. conform_dims(dimsizes(mean_change_binary), mmm_change_binary, (/1,2/)), 1, 0)
      aos_sum = dim_sum_n_Wrap(aos_tmp,0)
      tmp = new((/dimsizes(aos_sum(:,0)),dimsizes(aos_sum(0,:))/),float)
      tmp(:,:) = dimsizes(models)
      aos = aos_sum/tmp
      delete([/tmp,aos_sum,aos_tmp/])
      copy_VarCoords(tmp1(0,0,:,:),aos(:,:))
      ; -- magnitude (mag) of change
      mag = where(mmm_change .gt. 1*dim_stddev_n_Wrap(mean_change,0), 1, 0)
      copy_VarCoords(tmp1(0,0,:,:),mag(:,:))
      ofile   = pathin+"historical2rcp85/"+comp+"/"+vari+"/all_models_"+vari+"_remapcon2_1x1_wy_anom"+refstart+"-"+refende+"_mmm_change"+futstart+"-"+futende+"_aos_mag.nc"
      system("rm -rf "+ofile)
      a         = addfile(ofile,"c")
      a->mmm = mmm
      a->mmm_change = mmm_change
      if vari .eq. "mrro" .or. vari .eq. "pr"
        a->mmm_change_rel = mmm_change_rel
      end if
      a->mean_change = mean_change
      a->aos = aos
      a->mag = mag
    end do
  end if  ; end of "recalc4" loop






  ;**************************************************************************
  ; PLOTTING

  ofile           = pathin+"historical2rcp85/Amon/tas/all_models_tas_remapcon2_1x1_wy_anom"+refstart+"-"+refende+"_mmm_change"+futstart+"-"+futende+"_aos_mag.nc"
  a               = addfile(ofile,"r")
  tas_mmm_change  = a->mmm_change
  tas_aos         = a->aos
  tas_mag         = a->mag
  ofile           = pathin+"historical2rcp85/Amon/pr/all_models_pr_remapcon2_1x1_wy_anom"+refstart+"-"+refende+"_mmm_change"+futstart+"-"+futende+"_aos_mag.nc"
  a               = addfile(ofile,"r")
  pr_mmm_change   = a->mmm_change_rel
  pr_aos          = a->aos
  pr_mag          = a->mag
  ofile           = pathin+"historical2rcp85/Lmon/mrro/all_models_mrro_remapcon2_1x1_wy_anom"+refstart+"-"+refende+"_mmm_change"+futstart+"-"+futende+"_aos_mag.nc"
  a               = addfile(ofile,"r")
  mrro_mmm_change = a->mmm_change_rel
  mrro_aos        = a->aos
  mrro_mag        = a->mag
  ofile           = pathin+"historical2rcp85/Lmon/re/all_models_re_remapcon2_1x1_wy_anom"+refstart+"-"+refende+"_mmm_change"+futstart+"-"+futende+"_aos_mag.nc"
  a               = addfile(ofile,"r")
  tmp             = a->mmm_change
  re_mmm_change   = tmp*100
  copy_VarCoords(tmp,re_mmm_change)
  delete(tmp)
  re_aos          = a->aos
  re_mag          = a->mag

  varp 		    = "/home/flehner/publication/lehner18_runoff_efficiency_constraint/fig/cmip5_mmm_change_wy_"+futtas_aim+"C_warming_"+region
  wks_type 	  = "eps"
  wks 		    = gsn_open_wks(wks_type,varp)
  nbplot      = 4
  plot 		    = new(nbplot,graphic)
  plot_stip1  = new(nbplot,graphic)
  plot_stip2  = new(nbplot,graphic)

  ; -- merge Adam's colormap with precip_diff_12lev
  gsn_define_colormap(wks,"amwg_blueyellowred")
  cmap = gsn_retrieve_colormap(wks)
  hsv = rgbhsv(cmap)
  hsv(10,1) = .47
  hsv(8,1) = .45
  hsv(9,1) = .30
  hsv(4,:) = (/ (hsv(3,:)+hsv(5,:))/2.0 /)
  cmap2 = hsvrgb(hsv)
  gsn_define_colormap(wks,cmap2)
  newcolor = NhlNewColor(wks,0.42,0.42,0.42)
  nc0 = NhlNewColor(wks,0.95,0.95,0.95)
  nc1 = NhlNewColor(wks,0.7,0.7,0.7)
  nc3 = NhlNewColor(wks,0.5,0.5,0.5)
  nc4 = NhlNewColor(wks,0.30,0.30,0.30)
  nc2 = NhlNewColor(wks,139/256.,35/256.,35/256.)
  gsn_merge_colormaps(wks,cmap2,"precip_diff_12lev")

  res                     = True
  ; res@mpProjection        = "Robinson"
  res@mpFillOn = True                       ; default = True
  res@mpFillColors = (/-1,0,-1,0/)      ; map color fill ocean
  res@cnFillDrawOrder = "PreDraw"        ; color fill first
  res@mpOutlineDrawOrder = "Draw"      ; draw map outline after cnFill
  res@mpFillDrawOrder = "Draw"              ; draw map fill next
  res@mpPerimDrawOrder = "PostDraw"   ; draw map perimeter last
  res@gsnDraw              = False    ; do not draw picture
  res@gsnFrame             = False    ; do not advance frame
  res@gsnPaperOrientation  = "landscape"
  res@gsnAddCyclic = False;True
  res@mpLimitMode   = "LatLon"
  res@mpOutlineOn = True
  res@mpDataSetName = "Earth..4"
  res@mpDataBaseVersion = "MediumRes"
  res@mpGeophysicalLineColor = "gray42"
  res@mpCountyLineColor = "transparent"
  res@mpUSStateLineColor = "gray42"
  res@mpNationalLineColor = "gray22"
  res@mpGeophysicalLineThicknessF = 1.4
  res@mpOutlineBoundarySets = "National" ; "AllBoundaries" "GeophysicalAndUSStates";"USStates";"National"

  res@mpMinLatF            = minlat
  res@mpMaxLatF            = maxlat
  res@mpMinLonF            = minlon
  res@mpMaxLonF            = maxlon
;  res@mpCenterLonF         = 255.         ; defailt is 0 [GM]

; Define colors for filled contours
  res@cnFillOn             = True
  res@gsnSpreadColors      = True                  ; use full colormap
  res@cnLevelSelectionMode = "ExplicitLevels"

  res@cnLinesOn           = False         ; Turn lines off
  res@cnLineLabelsOn      = False         ; Turn labels off
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@mpGridLineColor      =  -1       ; Turn off grid lat/lon lines
  res@mpGridAndLimbOn      =  True
  res@lbLabelBarOn        = True; False
  res@lbTitlePosition         = "Bottom"
  res@lbLabelFontHeightF      = 0.02
  res@lbTitleFontHeightF      = 0.02
  res@pmLabelBarHeightF       = 0.08
  res@cnMissingValFillPattern = "SolidFill"
  res@cnMissingValFillColor =-1
  res@tmXBOn               = False
  res@tmXTOn               = False
  res@tmYLOn               = False
  res@tmYROn               = False

  ; -- tas resources
  rest = res
  if futtas_aim .eq. 2
    rest@cnFillColors     = (/10,11,12,13,14,15,16/) ; for lower warming period
    rest@cnLevels         = fspan(1.5,4,6) ; for lower warming period
  else if futtas_aim .eq. 3
    rest@cnFillColors     = (/10,11,12,13,14,15,16/)
    rest@cnLevels         = fspan(2.5,5,6)
  end if
  end if
  rest@lbTitleString    = "Temperature change (~F34~0~F~C)"
  ; -- pr resources
  resp = res
  delete(resp@cnFillColors)
  if futtas_aim .eq. 2
    resp@cnFillColors     = (/18,19,20,21,23,25,26,27,29,30/) ; for lower warming period
    resp@cnLevels         = fspan(-20,20,9) ; for lower warming period
  else if futtas_aim .eq. 3
    resp@cnFillColors     = (/18,19,20,21,22,23,25,26,27,28,29,30/)
    resp@cnLevels         = fspan(-25,25,11)
  end if
  end if
  resp@lbTitleString    = "Precipitation change (%)"
  ; -- mrro resources
  resm = resp
  resm@lbTitleString    = "Runoff change (%)"
  ; -- re resources
  resr = res
  if futtas_aim .eq. 2
    resr@cnFillColors     = (/18,19,20,21,23,25,26,27,29,30/) ; for lower warming period
    resr@cnLevels         = fspan(-2,2,9) ; for lower warming period
  else if futtas_aim .eq. 3
    resr@cnFillColors     = (/18,19,20,21,22,23,25,26,27,28,29,30/)
    resr@cnLevels         = fspan(-2.5,2.5,11)
  end if
  end if
  resr@lbTitleString    = "Runoff efficiency change (percentage points)"


  res2 = res
  res2@cnLinesOn           = True         ; Turn lines off
  res2@cnLineColor        = "blue";"gray10"
  res2@cnLineThicknessF      = 1.5
  res2@cnFillOn            = False
  res2@cnInfoLabelOn       = False
  res2@cnLevelSelectionMode = "ManualLevels"      ; manually set cn levels
  delete(res2@cnLevels)
  res2@cnMinLevelValF       = -10.25                 ; min level
  res2@cnMaxLevelValF       =  10.25                 ; max level
  res2@cnLevelSpacingF      =  .5                  ; contour level spacing
  res2@gsnContourNegLineDashPattern = 2

  res3=True
  res3@gsnDraw              = False        ; don't draw yet
  res3@gsnFrame             = False
  res3@cnFillOn            = True
  res3@cnLinesOn           = False
  res3@cnLineLabelsOn      = False
  res3@gsnLeftString = ""
  res3@cnLineThicknessF      = 3
  res3@cnMonoFillColor       = True
  res3@cnMonoFillPattern     = False
  res3@cnMonoFillScale       = True
  res3@cnFillColor = "black"
  res3cnFillDotSizeF = 0.003
  res3@cnLevelSelectionMode  = "ExplicitLevels"
  res3@cnFillPatterns        = (/8,-1/)
  res3@cnFillScaleF = 0.5
  res3@cnLevels              = fspan(-3,3,7)
  res3@cnInfoLabelOn         = False
  res3@lbLabelBarOn        = False

  res@tiMainFontHeightF     = 0.015
  res@tiYAxisFontHeightF    = 0.015

; ---------------------------------------------------------------------

  dum = new(nbplot,graphic)

  plot(0) 	= gsn_csm_contour_map_ce(wks,tas_mmm_change,rest)
  plot(1) 	= gsn_csm_contour_map_ce(wks,pr_mmm_change,resp)
  plot(2) 	= gsn_csm_contour_map_ce(wks,mrro_mmm_change,resm)
  plot(3) 	= gsn_csm_contour_map_ce(wks,re_mmm_change,resr)

  ; -- aos
  ; -- overlay hatching
  resov1       = True
  ; resov1@mpProjection        = "Robinson"
  resov1@mpFillColors = (/-1,0,-1,0/)      ; map color fill ocean
  resov1@cnFillDrawOrder = "PreDraw"        ; color fill first
  resov1@mpOutlineDrawOrder = "Draw"      ; draw map outline after cnFill
  resov1@mpFillDrawOrder = "Draw"              ; draw map fill next
  resov1@mpPerimDrawOrder = "PostDraw"   ; draw map perimeter last
  resov1@gsnDraw              = False        ; don't draw yet
  resov1@gsnFrame             = False
  resov1@cnFillOn            = True
  resov1@cnLinesOn           = False
  resov1@cnLineLabelsOn      = False
  resov1@gsnLeftString = ""
  resov1@cnLineThicknessF      = 3
  resov1@cnMonoFillColor       = True
  resov1@cnMonoFillPattern     = False
  resov1@cnMonoFillScale       = True
  resov1@cnFillColor = "gray40"
  resov1cnFillDotSizeF = 0.005
  resov1@cnLevelSelectionMode  = "ExplicitLevels"
  resov1@cnFillPatterns        = (/8,-1/)
  resov1@cnFillScaleF = 1.2
  resov1@cnLevels              = (/.666/)
  resov1@cnInfoLabelOn         = False
  resov1@lbLabelBarOn        = False
  plot_stip1(1) =  gsn_csm_contour(wks,pr_aos,resov1)
  overlay(plot(1),plot_stip1(1))
  plot_stip1(2) =  gsn_csm_contour(wks,mrro_aos,resov1)
  overlay(plot(2),plot_stip1(2))
  plot_stip1(3) =  gsn_csm_contour(wks,re_aos,resov1)
  overlay(plot(3),plot_stip1(3))

  ; -- mag
  ; -- overlay stippling
  resov2 = resov1
  resov2@cnFillPatterns        = (/17,-1/)
  plot_stip2(1) =  gsn_csm_contour(wks,pr_mag,resov2)
  overlay(plot(1),plot_stip2(1))
  plot_stip2(2) =  gsn_csm_contour(wks,mrro_mag,resov2)
  overlay(plot(2),plot_stip2(2))

  re_mag({lat|minlat+1},{lon|minlon+360+1}) = 1
  plot_stip2(3) =  gsn_csm_contour(wks,re_mag,resov2)
  overlay(plot(3),plot_stip2(3))


  ; -- region mask contours (shapefiles of watersheds) -------------------------
  load "/home/flehner/scripts/ncl/shapefile_utils.ncl"
  shp_filename = "/project/cas/flehner/shapefiles/NHDPlusv2_HUC4.DALLE.shp"
  f = addfile(shp_filename, "r")   ; Open shapefile
  segments = f->segments
  geometry = f->geometry
  lon      = f->x
  lat      = f->y
  segsDims = dimsizes(segments)
  geomDims = dimsizes(geometry)
  ; -- Read global attributes
  geom_segIndex = f@geom_segIndex
  geom_numSegs  = f@geom_numSegs
  segs_xyzIndex = f@segs_xyzIndex
  segs_numPnts  = f@segs_numPnts
  lines1_0       = new(segsDims(0),graphic)   ; Array to hold polygons
  lines1_1       = new(segsDims(0),graphic)
  lines1_2       = new(segsDims(0),graphic)
  lines1_3       = new(segsDims(0),graphic)
  numFeatures = geomDims(0)
  segNum = 0
  plres             = True       ; resources for polylines
  plres@gsEdgesOn   = True       ; draw border around polygons
  plres@gsEdgeColor = "blue"
  plres@gsFillColor = "Transparent"
  plres@gsFillBackgroundColor = "Transparent"
  plres@gsEdgeThicknessF = 2
  plres@gsLineThicknessF = 1
  plres@tfPolyDrawOrder = "PostDraw"
  do i = 0,numFeatures-1
    startSegment = geometry(i, geom_segIndex)
    numSegments  = geometry(i, geom_numSegs)
    do seg=startSegment, startSegment+numSegments-1
      startPT = segments(seg, segs_xyzIndex)
      endPT = startPT + segments(seg, segs_numPnts) - 1
      lines1_0(segNum) = gsn_add_polygon(wks, plot(0), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines1_1(segNum) = gsn_add_polygon(wks, plot(1), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines1_2(segNum) = gsn_add_polygon(wks, plot(2), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines1_3(segNum) = gsn_add_polygon(wks, plot(3), lon(startPT:endPT), lat(startPT:endPT), plres)
      segNum = segNum + 1
    end do
  end do
  delete([/segments,geometry,lat,lon,numFeatures/])

  shp_filename = "/project/cas/flehner/shapefiles/HUC02_conus.shp"
  f = addfile(shp_filename, "r")   ; Open shapefile
  segments = f->segments
  geometry = f->geometry
  lon      = f->x
  lat      = f->y
  segsDims = dimsizes(segments)
  geomDims = dimsizes(geometry)
  ; -- Read global attributes
  geom_segIndex = f@geom_segIndex
  geom_numSegs  = f@geom_numSegs
  segs_xyzIndex = f@segs_xyzIndex
  segs_numPnts  = f@segs_numPnts
  lines2_0       = new(segsDims(0),graphic)   ; Array to hold polygons
  lines2_1       = new(segsDims(0),graphic)
  lines2_2       = new(segsDims(0),graphic)
  lines2_3       = new(segsDims(0),graphic)
  numFeatures = geomDims(0)
  segNum = 0
  plres             = True       ; resources for polylines
  plres@gsEdgesOn   = True       ; draw border around polygons
  plres@gsEdgeColor = "blue"
  plres@gsFillColor = "Transparent"
  plres@gsFillBackgroundColor = "Transparent"
  plres@gsEdgeThicknessF = 2
  plres@gsLineThicknessF = 1
  plres@tfPolyDrawOrder = "PostDraw"
  do i = 8,8;0,numFeatures-1
    startSegment = geometry(i, geom_segIndex)
    numSegments  = geometry(i, geom_numSegs)
    do seg=startSegment, startSegment+numSegments-1
      startPT = segments(seg, segs_xyzIndex)
      endPT = startPT + segments(seg, segs_numPnts) - 1
      lines2_0(segNum) = gsn_add_polygon(wks, plot(0), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines2_1(segNum) = gsn_add_polygon(wks, plot(1), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines2_2(segNum) = gsn_add_polygon(wks, plot(2), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines2_3(segNum) = gsn_add_polygon(wks, plot(3), lon(startPT:endPT), lat(startPT:endPT), plres)
      segNum = segNum + 1
    end do
  end do
  delete([/segments,geometry,lat,lon,numFeatures/])

  shp_filename = "/project/cas/flehner/shapefiles/northern_sierras.shp"
  f = addfile(shp_filename, "r")   ; Open shapefile
  segments = f->segments
  geometry = f->geometry
  lon      = f->x
  lat      = f->y
  segsDims = dimsizes(segments)
  geomDims = dimsizes(geometry)
  ; -- Read global attributes
  geom_segIndex = f@geom_segIndex
  geom_numSegs  = f@geom_numSegs
  segs_xyzIndex = f@segs_xyzIndex
  segs_numPnts  = f@segs_numPnts
  lines3_0       = new(segsDims(0),graphic)   ; Array to hold polygons
  lines3_1       = new(segsDims(0),graphic)
  lines3_2       = new(segsDims(0),graphic)
  lines3_3       = new(segsDims(0),graphic)
  numFeatures = geomDims(0)
  segNum = 0
  plres             = True       ; resources for polylines
  plres@gsEdgesOn   = True       ; draw border around polygons
  plres@gsEdgeColor = "blue"
  plres@gsFillColor = "Transparent"
  plres@gsFillBackgroundColor = "Transparent"
  plres@gsEdgeThicknessF = 2
  plres@gsLineThicknessF = 1
  plres@tfPolyDrawOrder = "PostDraw"
  do i = 0,numFeatures-1
    startSegment = geometry(i, geom_segIndex)
    numSegments  = geometry(i, geom_numSegs)
    do seg=startSegment, startSegment+numSegments-1
      startPT = segments(seg, segs_xyzIndex)
      endPT = startPT + segments(seg, segs_numPnts) - 1
      lines3_0(segNum) = gsn_add_polygon(wks, plot(0), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines3_1(segNum) = gsn_add_polygon(wks, plot(1), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines3_2(segNum) = gsn_add_polygon(wks, plot(2), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines3_3(segNum) = gsn_add_polygon(wks, plot(3), lon(startPT:endPT), lat(startPT:endPT), plres)
      segNum = segNum + 1
    end do
  end do
  delete([/segments,geometry,lat,lon,numFeatures/])
  ; -- end of region mask contours (shapefiles of watersheds) ------------------

  ; -- add a box --
  ypts = (/ 31.0,  31.0,  42.0,   42.0, 39.0, 31.0/)
  xpts = (/-118.0, -101.0,-101.0, -125.0, -125.0, -118.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 2                      ; thickness of lines
  dum1 = new(5,graphic)

  ; ---------------------------------------------------------------------

  resP1                  = True                    ; modify the panel plot
  resP1@gsnMaximize     = True
  resP1@gsnPanelRowSpec = True
  resP1@gsnPanelBottom    = 0.1
  resP1@gsnPanelTop    = 1.0

  labels =  (/"(a)","(b)","(c)","(d)"/)
  resP1@gsnPanelFigureStrings = labels
  resP1@gsnPanelFigureStringsFontHeightF = 0.015; 0.011

  resP1@txString        = "CMIP5 multi-model mean change ("+futstart+"-"+futende+")-("+refstart+"-"+refende+")"
  resP1@txFontHeightF = 0.02

  resP1@gsnPanelYWhiteSpacePercent = 1
  resP1@gsnPanelXWhiteSpacePercent = 1
  resP1@pmLabelBarOrthogonalPosF = -.04
  resP1@lbOrientation     = "horizontal"



  resP1@gsnPanelFigureStringsPerimOn = False
  resP1@amJust          = "BottomLeft"

  resP1@gsnPaperOrientation = "portrait"
  resP1@lbLabelStride = 1
  resP1@lbBoxLineColor = "gray20"
  resP1@gsnPanelLabelBar    = False;True       ; common label bar
  resP1@gsnFrame            = False      ; don't advance frame yet
  resP1@lbLabelAutoStride   = True       ; nice label bar labels

  gsn_panel(wks,plot,(/2,2/),resP1)
  delete(wks)
  ; system("gv "+varp+".eps &")

  print_clock("script end")

end
