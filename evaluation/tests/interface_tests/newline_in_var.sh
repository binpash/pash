x=$(for i in sh_352.09 one two
       do
         if [ "$i" != "one" ]
         then echo "$i"
         fi
       done
      ); echo "$x"

x=$(echo sh_352.18 line 1 > sh_352.18tmp && echo sh_352.18 line 2 >> sh_352.18tmp && cat sh_352.18tmp ); echo "$x"