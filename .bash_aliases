#!/bin/bash

alias ckanext='cd ~/ckanext-cfpb-extrafields/ckanext/cfpb_extrafields/'

#################### alias functions
printcode () { enscript --pretty-print --color --landscape --columns=2 --fancy-header -p code.ps $1 ; }
wiki () { dig +short txt $1.wp.dg.cx ; }
vic () { vi `/bin/ls -ltr ${*}* | tail -1|awk '{print($NF)}'`; }
calc () { echo "$1" | bc -l ; }
findh () { find ./ -name '*'"$1"'*' ; }
function pngtomov { for x in `/bin/ls *`; do new=`echo $x | sed 's/\.\([0-9][0-9]\)\./\1./'`; echo mv $x $new ; mv $x $new; done; }
function columns { awk '{print NF}' $1 | sort -nu | tail -n 1; }
#function latex { latex $1.tex; convert $1.dvi $1.ps; convert $1.ps $1.pdf; }
function dir_ps_to_eps { for fps in *.ps; do feps="`echo $fps|sed 's|\.ps|\.eps|g'`"; echo "$fps->$feps"; convert $fps $feps ;done }
function e2 { echo $1; }
function eg { egrep -in --exclude=*~ $@ ; }
function feg { egrep -in --exclude=*~ $1 `find "$2" -name '*'` ; }
function svnmerg { svn merge --dry-run -r $1:$2 $svndrudd/trunk/ .; }
function egf { egrep -in --exclude=*~ $@ ../*/*; }
# Find a file with a pattern in name:        
function ff()  { find . -type f -iname '*'$*'*' -ls ; }
function ffs() { find . -type f -iname "$*"'*' -ls ; }
function ffe() { find . -type f -iname '*'"$@" -ls ; }

# Find a file with pattern $1 in name and Execute $2 on it:
function fe()
{ find . -type f -iname '*'${1:-}'*' -exec ${2:-file} {} \;  ; }

function grepfind() { find . -type f -name "$2" -print0 | xargs -0 grep "$1" ; }

cat_pdfs () { python '/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py' "$@" ; }
function pdfcombine { gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=finished.pdf $1 $2 ;}

function line { sed -n "${1}p" ${2}; }
#might want to find { , but avoiding that in case of if statements that begin with {: grep "void|^{"
#find { when int double ... precedes on the same line or the previous line 
#finds first void before the line; $1=line# $2=filepath
function voidline { i=$1; func="" ; while [ -z "$func" ]; do i=$((i-1));  if [ $(($1-i)) -eq 200 ]; then func="fail";fi; func=`sed -n "${i}p" "$2"|grep void` ; done ; istart=$i; if [ "$func" != "fail" ]; then for i in {-10..10} ; do line $((i+istart)) $2; done; fi; }

#################### alias apps
alias ipynb="ipython notebook" 
alias egrep='egrep --exclude=*~ -s'
alias rm='rm -i'
alias rmt='rm -f *~ *.pyc'
alias mv='mv -i'
alias cp='cp -i'
alias cpl='cp -a'
alias ren='set noglob; ren ; unset noglob'
alias lss="/bin/ls"
alias ls='ls -G'
alias lst='ls -G -ltr -h'
alias gv='evince'

alias back='cd $OLDPWD'
alias diff='diff -a'

alias emacs='emacs'
alias vi='emacs -nw'
alias sm='sm -m ~/.sm2/start.sm'
alias make='make -j8'
alias ssh='ssh -o TCPKeepAlive=no -o ServerAliveInterval=15'
alias gnuplot='/usr/local/Cellar/gnuplot/4.6.3/bin/gnuplot'
alias kinit='kinit -r 1d'
alias gthumb='open -a preview'

export EDITOR='vi'
export SVN_EDITOR='vi'

export HISTCONTROL=ignoredups
export HISTSIZE=10000
export HISTIGNORE="&:clear:exit"
export PROMPT_COMMAND='history -a'
#PS1='[\u@\h$:\w]$  \[\033]0;\u@\h:\w\007\]'
#PS1='[@\h \w]$  \[\033]0;\u@\h:\w\007\]'

export PATH="~/.exec/:$PATH" 
export CPPFLAGS="$CPPFLAGS" #-I/opt/local/include 
export LDFLAGS="-L/Users/sleitner/local/lib/ $LDFLAGS" #-L/opt/local/lib 

########
alias act='source ../virtualenv/bin/activate'
alias killmysql='mysqladmin -u root -p shutdown'
alias startmysql='mysql.server start'

#
alias gitdontchange='git update-index --assume-unchanged' #file
alias gitdochange='git update-index --no-assume-unchanged' #file

