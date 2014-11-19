;+
; :NAME:
;   minmax
;
; :PURPOSE:
; Find maximum (positive and negative) change in time series.
;
; CATEGORY:
; Vegetation monitoring program
;
; :Returns:
;  Array: maximum positive (negative) change and year.
;
; :KEYWORD PARAMETERS:
;   DATA:     time series array (T, M, N).
;   TS_FILE:  time series file name.
;   OUT_FILE: set if you want to write the result to a file.
;
; :CALLING SEQUENCE:
;
;
; MODIFICATION HISTORY:
;    Written by:  Ron Drori 2014
;-
;



FUNCTION minmax,  ts_data=ts_data, $
                  ts_file=ts_file, $
                  out_file=out_file, $
                  base_year=base_year

IF KEYWORD_SET(ts_file) then begin
  ; Read MODIS yearly minimum NDVI
  ts_data = READ_TIFF(ts_file, GEOTIFF=gtag)
ENDIF

max_positive  = 0
max_year      = 1
max_negative  = 2
min_year      = 3
total_dif     = 4

dims = ts_data.DIM
ts = dims[0] & xs = dims[1] & ys=dims[2]

; array to store results: 2xchange,year
dif = fltarr(5, xs, ys)

; for each pixel find maximum change.
; last elemets is the change between start and end
for y=0, ys-1 do begin
  for x=0, xs-1 do begin
    ts_pixel = ts_data[* , x, y]
    s_ts_pixel = shift(ts_pixel,-1)
    dif_pixel = s_ts_pixel - ts_pixel
    resmax = MAX(dif_pixel[0:ts-2], smax, MIN=pmin, SUBSCRIPT_MIN=smin)
    
    if resmax gt 0. then begin
      dif[max_positive, x, y] = resmax
      dif[max_year, x, y] = smax 
    endif
    if pmin lt 0. then begin
      dif[max_negative, x, y] = pmin
      dif[min_year, x, y] = smin
    endif
    
    dif[total_dif, x, y] = ts_pixel[ts-1] - ts_pixel[0]
  endfor
endfor



; adjust to true year
if KEYWORD_SET(base_year) then begin
  base_year = base_year + 1
  dif[max_year,*,*] = dif[max_year,*,*] + base_year
  dif[min_year,*,*] = dif[min_year,*,*] + base_year
endif

IF KEYWORD_SET(out_file) then begin
  WRITE_TIFF, out_file, dif, GEOTIFF=gtag, /FLOAT
ENDIF

END


 