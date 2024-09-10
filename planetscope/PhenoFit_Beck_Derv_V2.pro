PRO PhenoFit_Beck_Derv_V2
  ;**********************************************************************************************;

  ;           Latest updated on 2023/02/11 by Daryl Yang <<dediyang@bnl.gov>>
  ;**********************************************************************************************;

  ;*************************** THIS IS THE START OF THE PROGRAM**********************************;
  Start_Time=systime(1)


  ;********************************* SET OUPUT DIRECTORY ****************************************;
  out_dir = '\\modex.bnl.gov\data2\dyang\projects\ngee_arctic\seward\analysis\phenology\kougarok\planetscope\2022\phenoR\pheno_idl_20230212'
  FILE_MKDIR, out_dir
  if out_dir then begin
    print, '--- Creating --- "', out_dir, '" --- Successful ---'
  endif else begin
    print, '--- Creating --- "', out_dir, '" --- Failed ---'
  endelse
  ;**********************************************************************************************;

  ;********************************* SET USER PARAMETER *****************************************;
  ; Set dates of data acquisition
  ;date_stamp = [93,125,129,136,141,152,155,159,163,166,171,178,185,190,192,215,217,228,231,236,$
  ;  256,267,273,284,286]

  date_stamp_orig = read_csv('\\modex.bnl.gov\data2\dyang\projects\ngee_arctic\seward\analysis\phenology\kougarok\planetscope\2022\VIs\doy.csv', $
    header = header)
  doy_stamp = date_stamp_orig.field2
  ; define experiment start and end dates (doy of year)
  expBEG = 130
  expEND = 275
  ;**********************************************************************************************;

  ;************************************** LOAD DATA *********************************************;
  ;Load UAS classification as raster file
  ndvi_series_dir = '\\modex.bnl.gov\data2\dyang\projects\ngee_arctic\seward\analysis\phenology\kougarok\planetscope\2022\VIs\ndgi_time_series_test_watershed.tif'
  ndvi_raster = envi.OpenRaster(ndvi_series_dir)
  ndvi_img_data = ndvi_raster.GetData()
  ndvi_img_data = transpose(ndvi_img_data, [2, 1, 0])
  ;**********************************************************************************************;

  ;date_stamp1 = doy_stamp[where(doy_stamp ge expBEG and doy_stamp le expEND)]
  ;mask = max(ndvi_img_data, dimension = 3)
  ;mask[where(mask lt 0.4)] = !VALUES.F_NAN
  ;mask[where(mask ge 0.4)] = 1
  ;********************************* SG filtering NDVI ******************************************;
  dims_ndvi = size(ndvi_img_data, /dimensions)

  ;date_stamp = date_stamp1
  ; define the first date extend to
  doy_start = 1
  ; define the last date extend to
  doy_end = 410
  ; define the interval for extending
  intv = 3

  pheno_pars = fltarr(dims_ndvi[0], dims_ndvi[1], 10)
  pheno_ndvi = fltarr(dims_ndvi[0], dims_ndvi[1], 365)
  for lon_col=0L, dims_ndvi[0]-1 do begin ; dims_ndvi[0]-1
    print, lon_col
    for lon_row = 0L, dims_ndvi[1]-1 do begin ;dims_ndvi[1]-1
      pixel_ndvi_in = reform(ndvi_img_data[lon_col, lon_row, *])
      ; remove observations before and after snow
      date_stamp = doy_stamp
      pixl_pheno = phenofit(pixel_ndvi_in, date_stamp, expBEG, expEND, doy_start, doy_end, intv, pheno_par=pheno_par, fitted_ndvi=fitted_ndvi)
      pheno_pars[lon_col,lon_row,*] = pheno_par
      pheno_ndvi[lon_col,lon_row,*] = fitted_ndvi
      delvar, pixel_ndvi_in
    endfor
  endfor

  ; clean the pheno results
  for i=0, 9 do begin
    band = pheno_pars[*,*,i]
    band_smoothed = simple_smooth(band, 5)
    band_smoothed = band_smoothed
    pheno_pars[*,*,i] = band_smoothed
  endfor

  pheno_pars_out = transpose(pheno_pars, [1, 0, 2])
  pheno_ndvi_out = transpose(pheno_ndvi, [1, 0, 2])

  print, '...... saving Double Logistic results'
  envi_open_file, ndvi_series_dir, r_fid=fid
  map_info = envi_get_map_info(fid = fid)
  file_name = file_basename(ndvi_series_dir, '.tif')
  out_name = out_dir + '\' + 'phenmetrics'
  envi_write_envi_file, pheno_pars_out, out_name = out_name, map_info = map_info, $
    bnames = ['ud', 'sd', 'dd', 'rd', 'ndvi_ud', 'ndvi_sd', 'ndvi_dd', 'ndvi_rd', 'ndvimin', 'ndvimax']

  out_name = out_dir + '\' + 'fitted_index'
  envi_write_envi_file, pheno_ndvi_out, out_name = out_name, map_info = map_info, $
    bnames = [1:365:1]

END

Function phenofit, pixel_ndvi, date_stamp, expBEG, expEND, doy_start, doy_end, intv, pheno_par=pheno_par, fitted_ndvi=fitted_ndvi
  ; replace data outside of desire doy range with baseline value
  subsnow = where(date_stamp le expBEG or date_stamp gt expEND)
  pixel_ndvi[subsnow] = 0.05 ; 0.05
  ; replace extreme low values using baseline values
  pixel_ndvi[where(pixel_ndvi lt 0.05)] = 0.05 ;0.05
  ; detect outliers leveraging sg filter
  ndvi_smooth = ts_smooth(pixel_ndvi, 7)
  ndvi_diff = pixel_ndvi - ndvi_smooth
  subout = where(abs(ndvi_diff) gt 0.1)
  pixel_ndvi[subout] = 0
  ; fill up outlies
  subscribe_zero = where(pixel_ndvi gt 0 and pixel_ndvi lt 1, value_count)
  if (n_elements(pixel_ndvi) - value_count) lt n_elements(pixel_ndvi) then begin
    pixel_ndvi_filled = fill_up(temporary(pixel_ndvi), 1)
  endif

  ;;; fit ndvi to get phenology
  date_stamp = date_stamp
  pixel_ndvi = pixel_ndvi_filled
  ; filter NDVI with sg
  num = n_elements(pixel_ndvi)
  if num gt 15 then begin
    ;pixel_ndvi_sgfilterred =  sgfilter(pixel_ndvi)
    ; initialize phenofit parameters
    front_extend_dates = [doy_start:(min(date_stamp)-1):intv]
    end_extend_dates = [doy_end:(max(date_stamp)+1):-intv]
    end_extend_dates = end_extend_dates[sort(end_extend_dates)]
    ext_doy = [front_extend_dates, date_stamp, end_extend_dates]
    ;doublog_parms = fltarr(dims_ndvi[0], dims_ndvi[1], 6)
    ;fitted_ndvi = fltarr(dims_ndvi[0], dims_ndvi[1], n_elements(ext_doy))
    pixel_t_ndvi = extend(date_stamp, pixel_ndvi, ext_doy, expBEG, expEND, ext_ndvi = ext_ndvi)
    pixel_ndvi_sgfilterred =  sgfilter(ext_ndvi)
    
    parms = [0.01, 1.0, 0.05, 0.02, 140, 260]
    yerr = randomn(seed, n_elements(ext_doy))/100
    result = MPFITFUN('DoubleLogBeck', ext_doy, pixel_ndvi_sgfilterred, yerr, parms, yfit = yfit, $
      STATUS=ST, ERRMSG=ERR, /QUIET) ;
    ;doublog_parms[lon_col, lon_row, *] = result

    ; generate smooth fitted ndvi time series
    t = [1:365:1]
    if n_elements(result) eq 6 then begin
      mn = result[0]
      mx = result[1]
      rsp = result[2]
      rau = result[3]
      sos = result[4]
      eos = result[5]
      fitted_ndvi = mn + (mx-mn)*(1/(1+exp(-rsp*(t-sos))) + 1/(1+exp(rau*(t-eos))))

      ; determine pheno dates based third-order deritive
      ndvimax = max(fitted_ndvi[expBEG:expEND])
      ndvimin = min(fitted_ndvi[expBEG:expEND])

      UD = 2*alog(sqrt(3)-sqrt(2))/rsp + sos
      SD = 2*alog(sqrt(3)+sqrt(2))/rsp + sos
      DD = 2*alog(sqrt(3)-sqrt(2))/rau + eos
      RD = 2*alog(sqrt(3)+sqrt(2))/rau + eos

      if UD gt expBEG and UD lt expEND then begin
        ndviUD = fitted_ndvi[UD-1];
      endif else begin
        ndviUD = !VALUES.F_NAN
      endelse
      if SD gt expBEG and SD lt expEND then begin
        ndviSD = fitted_ndvi[SD-1];
      endif else begin
        ndviSD = !VALUES.F_NAN
      endelse
      if DD gt expBEG and DD lt expEND then begin
        ndviDD = fitted_ndvi[DD-1];
      endif else begin
        ndviDD = !VALUES.F_NAN
      endelse
      if RD gt expBEG and RD lt expEND then begin
        ndviRD = fitted_ndvi[RD-1];
      endif else begin
        ndviRD = !VALUES.F_NAN
      endelse

      ; store pheno results
      pheno_par = [UD, SD, DD, RD, ndviUD, ndviSD, ndviDD, ndviRD, ndvimin, ndvimax]
      fitted_ndvi = fitted_ndvi
    endif
    if n_elements(result) ne 6 then begin
      pheno_par = replicate(0, 10)
      fitted_ndvi = replicate(0, 365)
    endif
  endif
  if num le 15 then begin
    ; if the number of valid ndvi is less than tha half of input then skip
    pheno_par = replicate(0, 10)
    fitted_ndvi = replicate(0, 365)
  endif

  delvar, pixel_ndvi, date_stamp
End

Function double_logistic, t, parms
  pred_y = parms[0] + parms[1]*(1/(1+exp(parms[2]-parms[3]*t))-1/(1+exp(parms[4]-parms[5]*t)))
  return, pred_y
End

Function DoubleLogBeck, t, parms
  mn = parms[0]
  mx = parms[1]
  rsp = parms[2]
  rau = parms[3]
  sos = parms[4]
  eos = parms[5]
  pred_y = mn + (mx-mn)*(1/(1+exp(-rsp*(t-sos))) + 1/(1+exp(rau*(t-eos))))
  return, pred_y
End

Function simple_smooth, pheno_pars, windsize
  dims = size(pheno_pars, /dimensions)
  ns = dims[0]
  nl = dims[1]
  num_point_sample= ns mod windsize
  num_point_line= ns mod windsize
  if num_point_sample ne 0 then begin
    num_count_sample=fix(ns/windsize)+1
  endif else begin
    num_count_sample=fix(ns/windsize)
  endelse

  if num_point_line ne 0 then begin
    num_count_line=fix(nl/windsize)+1
  endif else begin
    num_count_line=fix(nl/windsize)
  endelse

  pheno_par_corrected=fltarr(dims)

  for i=0L,num_count_sample-1 do begin
    for j=0L,num_count_line-1 do begin
      end_point_sample=(i+1)*windsize-1
      end_point_line=(j+1)*windsize-1

      if end_point_sample le ns-1 and end_point_line le nl-1 then begin
        cut_img=pheno_pars[[i*windsize]:[(i+1)*windsize-1],[j*windsize]:[(j+1)*windsize-1]]
        resistant_mean, cut_img, 1, meanv, meansig
        cut_img[where(abs(cut_img-meanv) gt 2*stddev(cut_img))] = meanv
        pheno_par_corrected[[i*windsize]:[(i+1)*windsize-1],[j*windsize]:[(j+1)*windsize-1]] = cut_img
      endif
      if end_point_sample le ns-1 and end_point_line gt nl-1 then begin
        cut_img=pheno_pars[[i*windsize]:[(i+1)*windsize-1],[j*windsize]:[nl-1]]
        resistant_mean, cut_img, 1, meanv, meansig
        cut_img[where(abs(cut_img-meanv) gt 2*stddev(cut_img))] = meanv
        pheno_par_corrected[[i*windsize]:[(i+1)*windsize-1],[j*windsize]:[nl-1]] = cut_img
      endif
      if end_point_sample gt ns-1 and end_point_line le nl-1 then begin
        cut_img=pheno_pars[[i*windsize]:[ns-1],[j*windsize]:[(j+1)*windsize-1]]
        resistant_mean, cut_img, 1, meanv, meansig
        cut_img[where(abs(cut_img-meanv) gt 2*stddev(cut_img))] = meanv
        pheno_par_corrected[[i*windsize]:[ns-1],[j*windsize]:[(j+1)*windsize-1]] = cut_img
      endif
      if end_point_sample gt ns-1 and end_point_line gt nl-1 then begin
        cut_img=pheno_pars[[i*windsize]:[ns-1],[j*windsize]:[nl-1]]
        resistant_mean, cut_img, 1, meanv, meansig
        cut_img[where(abs(cut_img-meanv) gt 2*stddev(cut_img))] = meanv
        pheno_par_corrected[[i*windsize]:[ns-1],[j*windsize]:[nl-1]] = cut_img
      endif
    endfor
  endfor
  return, pheno_par_corrected
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
  a
End

Function sgfilter,vector_in                 ;S-G filter

  num_elements = n_elements(vector_in)
  ; The first Savitzky-Golay fitting
  vector_in=reform(vector_in,num_elements)                         ; num_elements is the number of values of time-series
  savgolFilter = SAVGOL(5,5,0,2)                          ;set the window width(4,4) and degree (2) for computing trend curve
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
    savgolFilter = SAVGOL(5, 5, 0, 6)        ;set the window width(4,4) and degree (6) for repetition
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

Function extend, time_stamp, ndvi, ext_doy, expBEG, expEND, ext_ndvi = ext_ndvi
  ; define the first date extend to
  ;ndvi = ndvi[where(time_stamp ge expBEG and time_stamp le expEND)]
  ;time_stamp = time_stamp[where(time_stamp ge expBEG and time_stamp le expEND)]

  ndvi_ref = percentiles(ndvi, conflimit = 0.95)

  n_front_dates = n_elements(where(ext_doy lt time_stamp[0]))
  front_extend_ndvi = replicate(0.05, n_front_dates) ;+ randomu(seed, n_front_dates)*0.01 ;
  ;front_extend_ndvi = replicate(ndvi_ref[0], n_front_dates) + randomu(seed, n_front_dates)*0.02 ;
  n_end_dates = n_elements(where(ext_doy gt time_stamp[n_elements(time_stamp)-1]))
  end_extend_ndvi = replicate(0.05, n_end_dates) ;+ randomu(seed, n_end_dates)*0.01 
  ;end_extend_ndvi = replicate(ndvi_ref[0], n_end_dates) + randomu(seed, n_end_dates)*0.02 

  ext_ndvi = [front_extend_ndvi, ndvi, end_extend_ndvi]
End

function percentiles, values, conflimit=conflimit

  n = n_elements(values)
  if n_elements(conflimit) eq 0 then conflimit=0.68

  lowindex = long(((1.0-conflimit)/2)*n)
  highindex = n-lowindex-1
  sortvalues = values[sort(values)]

  return, [sortvalues[lowindex],sortvalues[highindex]]

end

