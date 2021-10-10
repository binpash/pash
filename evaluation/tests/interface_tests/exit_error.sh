( ( true ) 3>/dev/null/abc; echo $?; false); echo $?
({ true; } 3>/dev/null/abc; echo $?; false); echo $?
(for i in 1; do true; done 3>/dev/null/abc; echo $?; false); echo $?
(case x in (x) true ;; esac 3>/dev/null/abc; echo $?; false); echo $?
(if true; then true; fi 3>/dev/null/abc; echo $?; false); echo $?
(while false; do true; done 3>/dev/null/abc; echo $?; false); echo $?
(until true; do true; done 3>/dev/null/abc; echo $?; false); echo $?
(func() { true; } 3>/dev/null/abc && func; echo $?; false); echo $?
func() { true; }; (func 3>/dev/null/abc; echo $?; false); echo $?
(name_of_a_command_that_will_not_be_found; echo $?; false); echo $?