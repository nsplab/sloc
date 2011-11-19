" Vim command for searching the Deal.II source and include files,
" assuming you've installed it in ~/dev/deal.II/
"
" To use this script, source it directly from within vim
"
"   :source dealii.vim
"
" Now you can search the deal.II source code by passing an argument
" to any of the GD, GDH, GDS commands.
"
"   :GD Triangulation
"   :GDH DoFHandler
"   :GDS GridIn.*::read_ucd
"

let s:dealii_path = $HOME . '/dev/deal.II/'

function! GrepDeal2(regexp)
    let d = s:dealii_path
    execute ':vimgrep /' . a:regexp . '/gj ' . d . 'include/deal.II/*/*.h ' . d . 'source/*/*.{cc,in} '
    copen
endfunction
command! -nargs=1  GD  call GrepDeal2(<q-args>)

function! GrepDeal2Headers(regexp)
    let d = s:dealii_path
    execute ':vimgrep /' . a:regexp . '/gj ' . d . 'include/deal.II/*/*.h'
    copen
endfunction
command! -nargs=1  GDH  call GrepDeal2Headers(<q-args>)

function! GrepDeal2Sources(regexp)
    let d = s:dealii_path
    execute ':vimgrep /' . a:regexp . '/gj ' . d . 'source/*/*.{cc,in}'
    copen
endfunction
command! -nargs=1  GDS  call GrepDeal2Sources(<q-args>)

