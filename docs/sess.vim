let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/bgwas3/docs
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +38 proj.wf
badd +6 todo.wf
badd +0 inbox.wf
badd +0 report.Rmd
argglobal
silent! argdel *
$argadd proj.wf
edit report.Rmd
set splitbelow splitright
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
2wincmd h
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd w
wincmd w
wincmd t
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 119) / 238)
exe '2resize ' . ((&lines * 45 + 31) / 62)
exe 'vert 2resize ' . ((&columns * 67 + 119) / 238)
exe '3resize ' . ((&lines * 14 + 31) / 62)
exe 'vert 3resize ' . ((&columns * 67 + 119) / 238)
exe 'vert 4resize ' . ((&columns * 69 + 119) / 238)
argglobal
setlocal fdm=expr
setlocal fde=FoldLevelRmd(v:lnum)
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 7 - ((6 * winheight(0) + 30) / 60)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
7
normal! 0
wincmd w
argglobal
if bufexists('proj.wf') | buffer proj.wf | else | edit proj.wf | endif
setlocal fdm=expr
setlocal fde=FoldLevel(v:lnum)
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=4
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 1 - ((0 * winheight(0) + 22) / 45)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
argglobal
if bufexists('todo.wf') | buffer todo.wf | else | edit todo.wf | endif
setlocal fdm=expr
setlocal fde=FoldLevel(v:lnum)
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 6 - ((1 * winheight(0) + 7) / 14)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
6
normal! 075|
wincmd w
argglobal
if bufexists('inbox.wf') | buffer inbox.wf | else | edit inbox.wf | endif
setlocal fdm=expr
setlocal fde=FoldLevel(v:lnum)
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 1 - ((0 * winheight(0) + 30) / 60)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 119) / 238)
exe '2resize ' . ((&lines * 45 + 31) / 62)
exe 'vert 2resize ' . ((&columns * 67 + 119) / 238)
exe '3resize ' . ((&lines * 14 + 31) / 62)
exe 'vert 3resize ' . ((&columns * 67 + 119) / 238)
exe 'vert 4resize ' . ((&columns * 69 + 119) / 238)
tabnext 1
if exists('s:wipebuf') && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 winminheight=1 winminwidth=1 shortmess=filnxtToOF
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
