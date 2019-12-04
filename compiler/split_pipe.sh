## This is a utility command that will be injected by the planner to split pipes
size=100
head -n $size $in_pipe > $out_pipe1 & cat $in_pipe > $out_pipe2
