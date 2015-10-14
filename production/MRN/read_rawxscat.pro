FUNCTION read_rawxscat, code

xpos = [0,0.2,0.5,0.75,0.9,0.95,0.990,0.999] ; 8 values
pos_str =['0.000','0.200','0.500','0.750','0.900','0.950','0.990','0.999']

rad  = [10., 30., 60., 120.]  ; 4 values
rad_str = ['010','030','060','120']

result = fltarr(4, 8, 1450)


for iX=0,n_elements(xpos)-1 do begin
   for iR=0,n_elements(rad)-1 do begin
      filename = 'xs_'+code+'_'+pos_str[iX]+'_'+rad_str[iR]+'.dat'
      print,'Reading ',filename
      read_array, filename, data
      energy = reform(data[0,*])
      Nelem = n_elements(reform(data[*,0]))
      values = reform(data[Nelem-1,*])
      for iVal=0,n_elements(energy)-1 do result[iR, iX, iVal]=values[iVal]
   endfor
endfor

return, result
end


FUNCTION read_energy, code

xpos = [0,0.2,0.5,0.75,0.9,0.95,0.990,0.999] ; 8 values
pos_str =['0.000','0.200','0.500','0.750','0.900','0.950','0.990','0.999']

rad  = [10., 30., 60., 120.]  ; 4 values
rad_str = ['010','030','060','120']

result = fltarr(4, 8, 49)


filename = 'xs_'+code+'_'+pos_str[0]+'_'+rad_str[0]+'.dat'
print,'Reading ',filename
read_array, filename, data
energy = reform(data[0,*])

return, energy
end


