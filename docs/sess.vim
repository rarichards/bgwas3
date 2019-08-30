let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/bgwas3/code/docs
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +0 report.Rmd
badd +0 todo.wf
badd +0 bgwas3.bib
badd +0 inbox.wf
argglobal
silent! argdel *
$argadd report.Rmd
edit report.Rmd
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd w
wincmd t
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 84 + 116) / 233)
exe '2resize ' . ((&lines * 26 + 28) / 56)
exe 'vert 2resize ' . ((&columns * 74 + 116) / 233)
exe '3resize ' . ((&lines * 26 + 28) / 56)
exe 'vert 3resize ' . ((&columns * 73 + 116) / 233)
exe '4resize ' . ((&lines * 27 + 28) / 56)
exe 'vert 4resize ' . ((&columns * 148 + 116) / 233)
argglobal
setlocal fdm=expr
setlocal fde=FoldLevelRmd(v:lnum)
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
90
normal! zo
112
normal! zo
119
normal! zo
129
normal! zo
133
normal! zo
let s:l = 171 - ((74 * winheight(0) + 27) / 54)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
171
normal! 01|
wincmd w
argglobal
if bufexists('todo.wf') | buffer todo.wf | else | edit todo.wf | endif
setlocal fdm=expr
setlocal fde=FoldLevel(v:lnum)
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=3
setlocal fml=1
setlocal fdn=20
setlocal fen
1
normal! zo
let s:l = 9 - ((8 * winheight(0) + 13) / 26)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
9
normal! 05|
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
19
normal! zo
22
normal! zo
let s:l = 29 - ((17 * winheight(0) + 13) / 26)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
29
normal! 068|
wincmd w
argglobal
if bufexists('bgwas3.bib') | buffer bgwas3.bib | else | edit bgwas3.bib | endif
setlocal fdm=marker
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 1 - ((0 * winheight(0) + 13) / 27)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
3wincmd w
exe 'vert 1resize ' . ((&columns * 84 + 116) / 233)
exe '2resize ' . ((&lines * 26 + 28) / 56)
exe 'vert 2resize ' . ((&columns * 74 + 116) / 233)
exe '3resize ' . ((&lines * 26 + 28) / 56)
exe 'vert 3resize ' . ((&columns * 73 + 116) / 233)
exe '4resize ' . ((&lines * 27 + 28) / 56)
exe 'vert 4resize ' . ((&columns * 148 + 116) / 233)
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
