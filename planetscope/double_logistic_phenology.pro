PRO DOUBLE_LOGISTIC_PHENOLOGY
;**********************************************************************************************;
; This script is for calculate fractional cover of different PFTs/species from UAS
; classification image based on the resolution of Airborne/Satellite images

; Key processes:
; 1. adding/deleting columns and rows to match the size of UAS data with Airborne/Satellite data
;    for calculating Fcover
; 2. calculating Fcover of different PFTs or species within each Airborne/Satellite pixel

; Critial notes:
; 1. the numbers of adding/deleting columns & rows to add are needed, please manually compare the
;    UAS and Airborne/Satellite image to detemine how many collumns/rows are need to be added or
;    deleted. This is critial to minimize the spatial mismatched between UAS and Airborne/
;    Satellite
; 2. the resolution of Airborne/Satellite should be an integer fold of UAS

;           Latest updated on 2020/03/23 by Daryl Yang <<dediyang@bnl.gov>>
;**********************************************************************************************;

;*************************** THIS IS THE START OF THE PROGRAM**********************************;
Start_Time=systime(1)



;********************************* SET OUPUT DIRECTORY ****************************************;
out_dir = 'D:\Temp_Figure_forSS\Step6_Phenology'
FILE_MKDIR, out_dir
if out_dir then begin
  print, '--- Creating --- "', out_dir, '" --- Successful ---'
endif else begin
  print, '--- Creating --- "', out_dir, '" --- Failed ---'
endelse
;**********************************************************************************************;

;********************************* SET USER PARAMETER *****************************************;
; Set dates of data acquisition
date_stamp = [93,125,129,136,141,152,155,159,163,166,171,178,185,190,192,215,217,228,231,236,$
  256,267,273,284,286]

;**********************************************************************************************;

;************************************** LOAD DATA *********************************************;
;Load UAS classification as raster file
ndvi_series_dir = 'D:\Temp_Figure_forSS\Step5_NDVI\NDVI_TIME_SERIES.img'
ndvi_raster = envi.OpenRaster(ndvi_series_dir)
ndvi_img_data = ndvi_raster.GetData()

;**********************************************************************************************;

;********************************* SG filtering NDVI ******************************************;
dims_ndvi = size(ndvi_img_data, /dimensions)
ndvi_sg_filtered = fltarr(dims_ndvi)
for lon_col=0L, dims_ndvi[0]-1 do begin
  for lon_row = 0L, dims_ndvi[1]-1 do begin
    pixel_ndvi = ndvi_img_data[lon_col, lon_row, *]
    pixel_ndvi[where(pixel_ndvi lt 0 or pixel_ndvi gt 1)] = 0
    pixel_ndvi[where(~finite(pixel_ndvi))] = 0
    subscribe_zero = where(pixel_ndvi gt 0 and pixel_ndvi lt 1, value_count)
    if (dims_ndvi[2] - value_count) lt dims_ndvi[2] then begin
      pixel_ndvi_filled = fill_up(temporary(pixel_ndvi), 1)
      ndvi_sg_filtered[lon_col, lon_row, *] = sgfilter(pixel_ndvi_filled)
    endif
  endfor
endfor

print, '...... saving SG smoothed NDVI time series'
envi_open_file, ndvi_series_dir, r_fid=fid
map_info = envi_get_map_info(fid = fid)
file_name = getfilename(ndvi_series_dir)
out_name = out_dir + '\' + 'SG_Smoothed_' + file_name
envi_write_envi_file, ndvi_sg_filtered, out_name = out_name, map_info = map_info, $
  bnames = ndvi_raster.metadata['BAND NAMES']
;**********************************************************************************************;

;******************************* Double Logistic fitting **************************************;
parms = [0.1, 0.8, 1, 0.01, 2, 0.01]
yerr = randomn(seed, n_elements(date_stamp))/100
doublog_parms = fltarr(dims_ndvi[0], dims_ndvi[1], 6)
fitted_ndvi = fltarr(dims_ndvi)
for lon_col=0L, dims_ndvi[0]-1 do begin
  print, lon_col
  for lon_row = 0L, dims_ndvi[1]-1 do begin
    pixel_ndvi = reform(ndvi_sg_filtered[lon_col, lon_row, *])
    result = MPFITFUN('double_logistic', date_stamp, pixel_ndvi, yerr, parms, yfit = yfit, $
      STATUS=ST, ERRMSG=ERR, /QUIET) ;
    doublog_parms[lon_col, lon_row, *] = result
    fitted_ndvi[lon_col, lon_row, *] = yfit
  endfor
endfor

print, '...... saving Double Logistic results'
envi_open_file, ndvi_series_dir, r_fid=fid
map_info = envi_get_map_info(fid = fid)
file_name = getfilename(ndvi_series_dir)
out_name = out_dir + '\' + 'Double_Logistic_Parms_' + file_name
envi_write_envi_file, doublog_parms, out_name = out_name, map_info = map_info, $
  bnames = ['NDVImin', 'NDVImax', 'm1', 'n1', 'm2', 'n2']
  
out_name = out_dir + '\' + 'Double_Logistic_Fitted_NDVI_' + file_name
envi_write_envi_file, fitted_ndvi, out_name = out_name, map_info = map_info, $
  bnames = ndvi_raster.metadata['BAND NAMES']

END

Function double_logistic, t, parms
pred_y = parms[0] + parms[1]*(1/(1+exp(parms[2]-parms[3]*t))-1/(1+exp(parms[4]-parms[5]*t)))
return, pred_y
End







Function fill_up, vector_in, maxNDVI               ;remove the cloudy values

  num_elements=n_elements(vector_in)

  ;remove  continuous 0
  in_pixel = 0
  if vector_in[in_pixel] EQ 0 then begin
    pixel_start = in_pixel
    while vector_in[in_pixel] EQ 0 && in_pixel LT num_elements -2 do begin
      in_pixel = in_pixel +1
    endwhile
    pixel_end = in_pixel
    vector_in[pixel_start:pixel_end-1] = vector_in[pixel_end]
  endif

  in_pixel = num_elements-1
  if vector_in[in_pixel] EQ 0 then begin
    pixel_end = in_pixel
    while vector_in[in_pixel] EQ 0 && in_pixel GT 2 do begin
      in_pixel = in_pixel -1
    endwhile
    pixel_start = in_pixel
    vector_in[pixel_start+1:pixel_end] = vector_in[pixel_start]
  endif
  ;-----------
  in_pixel = 1
  while in_pixel LT num_elements - 2 do begin
    if vector_in[in_pixel] EQ 0 then begin
      pixel_start = in_pixel
      num_pixel = 1
      while vector_in[in_pixel] EQ 0 && in_pixel LT num_elements -2 do begin
        in_pixel = in_pixel +1
        num_pixel = num_pixel +1
      endwhile
      pixel_end = in_pixel
      temp = (vector_in[pixel_end] - vector_in[pixel_start-1]) / num_pixel
      for pixel = 0, num_pixel-2 do begin
        vector_in[pixel_start + pixel] = vector_in[pixel_start -1] + (pixel+1)*temp
      endfor
    endif
    in_pixel = in_pixel + 1
  endwhile

  in_pixel = 1
  while  in_pixel LT num_elements -2 do begin
    if (vector_in[in_pixel] - vector_in[in_pixel -1]) GE 0.2*maxNDVI $
      && (vector_in[in_pixel] - vector_in[in_pixel +1]) GE 0.2*maxNDVI then begin
      vector_in[in_pixel] = (vector_in[in_pixel -1] + vector_in[in_pixel +1]) /2.0
      in_pixel = in_pixel +2
    endif else begin
      in_pixel = in_pixel +1
    endelse
  endwhile

  return, vector_in

End

Function sgfilter,vector_in                 ;S-G filter

  num_elements = n_elements(vector_in)
  ; The first Savitzky-Golay fitting
  vector_in=reform(vector_in,num_elements)                         ; num_elements is the number of values of time-series
  savgolFilter = SAVGOL(4,4,0,2)                          ;set the window width(4,4) and degree (2) for computing trend curve
  rst = CONVOL(vector_in, savgolFilter, /EDGE_TRUNCATE)

  ; Calculate the threshold for loop control, so that the fit is maximize
  gdis = 0.0
  fl = IntARR(num_elements)

  for i =0,(num_elements-1) do begin
    fl[i] = (vector_in[i] ge rst[i])
    gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])*abs(vector_in[i]-rst[i])
  endfor

  ra4 = fltARR(num_elements)
  pre = fltARR(num_elements)

  ormax = gdis
  num   = 0

  loop_times = 0l
  while (gdis le ormax) && loop_times LT 10 do begin
    loop_times = loop_times +1
    for i =0,(num_elements-1) do begin
      ra4[i] = (vector_in[i] ge rst[i]) ? vector_in[i] : rst[i]
      pre[i] = rst[i]
    endfor

    ; The Savitzky-Golay fitting
    savgolFilter = SAVGOL(3, 3, 0, 3)        ;set the window width(4,4) and degree (6) for repetition
    rst = CONVOL(ra4, savgolFilter, /EDGE_TRUNCATE)
    ormax = gdis
    ; Calculate the fitting-effect index
    gdis = 0.0
    for i =0,(num_elements-1) do begin
      gdis = gdis + fl[i]*abs(vector_in[i]-rst[i])*abs(vector_in[i]-rst[i])
    endfor
  endwhile

  if loop_times GE 1000 then begin
    print, 'loop times is: ', loop_times
  endif


  return, pre

End ; of function sgfilter


Function getfilename, pathname
  ;**********************************************************************************************;
  ; This function is for extracting the file name from a file directory
  ;**********************************************************************************************;
  ; get the filename in a directory
  IF (N_PARAMS() NE 1) THEN RETURN, ''
  idx = STRPOS(pathname,'\', /REVERSE_SEARCH)
  filename = STRMID(pathname, idx + 1)
  RETURN, filename
End
