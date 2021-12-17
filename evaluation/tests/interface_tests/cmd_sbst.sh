echo $(Testvar=set
       unset Testvar
       echo $Testvar${Testvar-sh_352.10}${Testvar+set}
      )
x=$(set one two three; echo sh_352.11 $1 $2 $3 $# $* "$@"); echo "$x"
x=$(set one "twoA   twoB"; echo sh_352.12 $1 "$2" $3 $# $* "$@"); echo "$x"