let SessionLoad = 1
if &cp | set nocp | endif
nnoremap  :echo 'word' expand('<cword>') 'has length' strlen(expand('<cword>'))
vnoremap ,u :norm ^xx
vnoremap ,k :norm I//
nnoremap ,d yy`v}kpi		ww2df:A;
nnoremap ,l ^2f:lyw/{%acC//end_0^%
nnoremap ,c /class =expand("<cword>")
nnoremap ,f /::=expand("<cword>")
let s:cpo_save=&cpo
set cpo&vim
nmap gx <Plug>NetrwBrowseX
nnoremap <silent> <Plug>NetrwBrowseX :call netrw#NetrwBrowseX(expand("<cWORD>"),0)
nnoremap <F4> :bn
nnoremap <F3> :b#
noremap <F2> :bp
noremap!  klyiwj$pa
let &cpo=s:cpo_save
unlet s:cpo_save
set backspace=indent,eol,start
set cindent
set cscopeprg=/usr/bin/cscope
set cscopetag
set cscopeverbose
set expandtab
set fileencodings=ucs-bom,utf-8,latin1
set guicursor=n-v-c:block,o:hor50,i-ci:hor15,r-cr:hor30,sm:block,a:blinkon0
set helplang=en
set hidden
set hlsearch
set iskeyword=48-57,a-z,A-Z,192-255,_
set ruler
set shiftwidth=2
set softtabstop=2
set viminfo='20,\"50
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/dealii-8.3.0/examples/darcy_p
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +489 sdd_simple_parallel.cc
badd +32 CMakeLists.txt
badd +188 ../step-22/step-22.cc
badd +310 stokes_stationary.w
badd +111 darcy_initialization.w
badd +2254 ../step-43/step-43.cc
badd +11 print_basic_stats.w
badd +247 stokes_initialization.w
badd +6 run_stationary.w
badd +34 lambda_initialization.w
badd +83 compute_flux.w
badd +494 ../step-40/step-40.cc
badd +13 flux_initialization.w
badd +95 map_linker.cc
badd +35 cg_lsq.w
badd +21 assemble_stokes_preconditioner.w
badd +6 flux_function_2.cc
badd +10 flux_function_2.h
badd +223 darcy_stationary.w
badd +23 build_stokes_preconditioner.w
badd +28 map_linker.h
badd +16 convergence_rates.cc
badd +2 convergence_rates.h
badd +55 N_operators.w
badd +25 darcy_K_inverse.cc
badd +2 darcy_K_inverse.h
badd +443 ../darcy_x/CNLSD/sdd_simple.cc
badd +1 _workers
badd +79 darcy_simple_linear_solvers.w
badd +732 ../step-20/step-20.cc
badd +77 stokes_solver.w
badd +34 stokes_linear_solvers.w
badd +92 darcy_simple_solver.w
badd +50 darcy_solver.w
badd +126 darcy_operator.w
badd +77 assemble_darcy_preconditioner.w
badd +61 darcy_linear_solvers.w
badd +3 darcy_minres_solver.w
badd +251 ../step-32/step-32.cc
badd +174 ../step-55/step-55.cc
badd +23 build_darcy_preconditioner.w
badd +6 assemble_darcy_system.w
badd +80 darcy_fgmres_solver.w
badd +8 darcy_block_schur_preconditioner.w
argglobal
silent! argdel *
argadd _workers
edit darcy_stationary.w
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
setlocal keymap=
setlocal noarabic
setlocal noautoindent
setlocal backupcopy=
setlocal nobinary
setlocal nobreakindent
setlocal breakindentopt=
setlocal bufhidden=
setlocal buflisted
setlocal buftype=
setlocal cindent
setlocal cinkeys=0{,0},0),:,0#,!^F,o,O,e
setlocal cinoptions=
setlocal cinwords=if,else,while,do,for,switch
setlocal colorcolumn=
setlocal comments=sO:*\ -,mO:*\ \ ,exO:*/,s1:/*,mb:*,ex:*/,://
setlocal commentstring=/*%s*/
setlocal complete=.,w,b,u,t,i
setlocal concealcursor=
setlocal conceallevel=0
setlocal completefunc=
setlocal nocopyindent
setlocal cryptmethod=
setlocal nocursorbind
setlocal nocursorcolumn
setlocal nocursorline
setlocal define=
setlocal dictionary=
setlocal nodiff
setlocal equalprg=
setlocal errorformat=
setlocal expandtab
if &filetype != 'cpp'
setlocal filetype=cpp
endif
setlocal foldcolumn=0
setlocal foldenable
setlocal foldexpr=0
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldmarker={{{,}}}
setlocal foldmethod=manual
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldtext=foldtext()
setlocal formatexpr=
setlocal formatoptions=croql
setlocal formatlistpat=^\\s*\\d\\+[\\]:.)}\\t\ ]\\s*
setlocal grepprg=
setlocal iminsert=0
setlocal imsearch=0
setlocal include=
setlocal includeexpr=
setlocal indentexpr=
setlocal indentkeys=0{,0},:,0#,!^F,o,O,e
setlocal noinfercase
setlocal iskeyword=48-57,a-z,A-Z,192-255,_
setlocal keywordprg=
setlocal nolinebreak
setlocal nolisp
setlocal lispwords=
setlocal nolist
setlocal makeprg=
setlocal matchpairs=(:),{:},[:]
setlocal modeline
setlocal modifiable
setlocal nrformats=octal,hex
setlocal nonumber
setlocal numberwidth=4
setlocal omnifunc=ccomplete#Complete
setlocal path=
setlocal nopreserveindent
setlocal nopreviewwindow
setlocal quoteescape=\\
setlocal noreadonly
setlocal norelativenumber
setlocal norightleft
setlocal rightleftcmd=search
setlocal noscrollbind
setlocal shiftwidth=2
setlocal noshortname
setlocal nosmartindent
setlocal softtabstop=2
setlocal nospell
setlocal spellcapcheck=[.?!]\\_[\\])'\"\	\ ]\\+
setlocal spellfile=
setlocal spelllang=en
setlocal statusline=
setlocal suffixesadd=
setlocal swapfile
setlocal synmaxcol=3000
if &syntax != 'cpp'
setlocal syntax=cpp
endif
setlocal tabstop=8
setlocal tags=
setlocal textwidth=0
setlocal thesaurus=
setlocal noundofile
setlocal undolevels=-123456
setlocal nowinfixheight
setlocal nowinfixwidth
setlocal wrap
setlocal wrapmargin=0
silent! normal! zE
let s:l = 218 - ((12 * winheight(0) + 17) / 35)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
218
normal! 0
tabnext 1
if exists('s:wipebuf')
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToO
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
