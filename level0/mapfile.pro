pro mapfile,file,gl
;map the downloaded file to see the coverage in shell
if n_params() lt 1 then begin
	print,'Syntax - MPFILE, file, [gl]'
	return
endif
if n_params() lt 2 then gl = 80 else gl = fix(gl)
f=file_search(file)
xt=fix(strmid(f,0,4))
yt=fix(strmid(f,4,4))
sb=strpos(f,'U') ge 0

x = (xt - gl*10)/5
y = (yt + 55)/5
in = where(x ge 0 and x le 19, count)
if count eq 0 then begin
	print, 'Error - no file found'
	return
endif
x=x[in]
y=y[in]
sb=sb[in]

ucount = intarr(20,23)
lcount = intarr(20,23)
for i=0,n_elements(x)-1 do if sb[i] then ucount(x[i],y[i])++ else lcount(x[i],y[i])++
diff = strarr(20,23)
tmp = where(ucount eq lcount) & if tmp[0] ne -1 then diff[tmp]=' '
tmp = where(ucount ne lcount) & if tmp[0] ne -1 then diff[tmp]='#'

t = string(ucount,format='(i0)')+diff+string(lcount,format='(i0)')
empty = where(t eq '0 0')
if empty[0] ne -1 then t[empty] = '   '
t = reform(t,20,23)

i = indgen(20)
gl = string(i/2 + gl,format='(i3)')
gl[where(i mod 2 eq 1)] = '   '
print,'   |'+strjoin(reverse(gl),'|')+'|   '
seg = strjoin(replicate('---+',21))+'---'
print,seg;strjoin(replicate('-',89))
for i=22,0,-1 do begin
gb = (i mod 2 eq 1)?string((i-1)/2-5,format='(i+03)'):'   '
print,gb+'|'+strjoin(reverse(t[*,i]),'|')+'|'+gb
print,seg;strjoin(replicate('-',89))
endfor
print,'   |'+strjoin(reverse(gl),'|')+'|   '
end
