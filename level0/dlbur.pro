pro dlbur, range=range, append=append
;Generate download.sh to download bur data from website
;Steps:
;1. Download webpage "downloadclasslist.php.html" from http://www.radioast.nsdc.cn/
;2. Modify gl, gb range in this file
;3. run dlbur and get download.sh
;4. run download.sh in shell
if ~keyword_set(range) then range=[0,360,-6,6]
openr,lunr,'downloadclasslist.php.html',/get_lun
openw,lunw,'download.sh',/get_lun, append=append
line=''
WHILE ~ eof(lunr) DO BEGIN
  readf, lunr, line

  if strcmp(line, ' <td height="25"><div align="left">',35) then begin;get file zie
    v=strpos(line,'left">')
    t=strmid(line,v+6)
    v=strpos(t,'<')
    t=strmid(t,0,v)
    if strlen(t) le 7 then filesize=float(t)
    continue
  endif

  if ~strcmp(line, " <td><a href=",13) then continue
  v=strpos(line,'href="')
  t=strmid(line,v+6)
  v=strpos(t,'"')
  path=strmid(t,0,v);get filename and suffix
  
  name=file_basename(path)
  onlyname=gettok(name,'.')
  suffix=name
  if suffix ne 'bur' or strlen(onlyname) ne 18 then continue;throw none bur file

  name=file_basename(path)
  gl=float(strmid(name,0,4))/10
  gb=float(strmid(name,4,4))/10
  if gl lt range[0] or gl gt range[1] or gb lt range[2] or gb gt range[3] then continue;throw bur outside box
  print,name,filesize,(file_info(name)).size/1024000d
  if filesize - (file_info(name)).size/1024000d gt 0.1 or (file_info(name)).exists eq 0  then begin
    print,'Will download: '+ name
    if file_test(name) then print, 'file size is not right'
    printf, lunw, 'axel -n 8 http://www.radioast.nsdc.cn/'+path
  endif
ENDWHILE
close,lunr,lunw
free_lun,lunr,lunw
end
