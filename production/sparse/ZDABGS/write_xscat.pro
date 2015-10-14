PRO write_xscat, filename, energy, xpos, rpos, data, $
    COMMENTS=comments, NANVALUE=nanvalue
;
;  Writes a XSCAT FITS file using the data in the passed structure.
;
;  The code first checks to see if filename exists.  If it doesn't, it
;  creates filename and then writes the primary header.  If it does,
;  the code checks to make sure this particular ion stage hasn't
;  already been written to the file.  Assuming all is OK, we are
;  positioned at the end of the file.  Then the header is created, and
;  the data is written to the file.
;
;
;  Check to see if file exists
;
result = findfile(filename,count=count)
if (count eq 0) then begin      ; must create file and primary HDU
    fxhmake, primary_header, /extend, /initialize, /date
    fxaddpar, primary_header, 'XSCAT','Dust scattering cross section'
    fxaddpar, primary_header, 'FILENAME','IDL_routine','Source of data'
    fxaddpar, primary_header, 'ORIGIN','AtomDB','http://www.atomdb.org'
    fxaddpar, primary_header, 'HDUCLASS','DSCAT','Dust Scattering data'
    fxaddpar, primary_header, 'HDUCLAS1','SIGMA',$
      'X-ray Dust scattering cross section data'
    version = '1.0.0'
    fxaddpar, primary_header, 'HDUVERS1', version, 'Version of datafile'
    fxwrite, filename, primary_header
endif

if (count gt 1) then begin
    print,'More than 1 file exists with this name: ',result
    stop
endif
;
;  Now we're all set, ready to create the new header.
;
integer = 1
float = 1.0
array8 = fltarr(8)
array4 = fltarr(4)
arrayN = fltarr(4, 8)
strings = '01234567890123456789'
nrows = n_elements(energy)*n_elements(array8)*n_elements(array4)
fxbhmake, header, nrows, 'XSCAT','Dust scattering cross section'
fxaddpar, header, 'HDUCLASS','DSCAT', 'Dust Scattering data'
fxaddpar, header, 'HDUCLAS1','SIGMA','X-ray Dust scattering cross section data'
version = '1.0.0'
fxaddpar, header, 'HDUVERS1', version, 'Version of datafile'
fxaddpar, header, 'NENERGY', n_elements(energy), 'Number of Energies'
fxaddpar, header, 'NXPOS'  , n_elements(array8), 'Number of Dust positions'
fxaddpar, header, 'NREXT'  , n_elements(array4), 'Number of Extraction Regions'
fxbaddcol, 1, header, float, 'Energy'
fxbaddcol, 2, header, float, 'DustPos'
fxbaddcol, 3, header, float, 'RadExtr'
fxbaddcol, 4, header, float, 'Sigma'
;
; Now add any additional comments
;
for i=0,n_elements(comments)-1 do begin
    if (strpos(comments(i),'COMMENT') eq 0) then begin
        commstr = strmid(comments(i),8)
        add_reference, header, commstr, /header
    endif else begin
        add_reference, header, comments[i], /header
    endelse
endfor
;
;  Now open the file and write the new header
;
fxbcreate, unit, filename, header
;
;  Write the data
;
;structnames = tag_names(data)
version = 1.0

row = 0
for i=0L,n_elements(energy)-1 do begin
  for j=0L, n_elements(array8)-1 do begin
     for k=0L, n_elements(array4)-1 do begin
     	row = row+1
        fxbwrite, unit, energy[i],   1, row
        fxbwrite, unit, xpos[j],     2, row
	fxbwrite, unit, rpos[k],     3, row
	fxbwrite, unit, data[k,j,i], 4, row
     endfor
  endfor
endfor
;
; Close the file
;
fxbfinish, unit

return
end
