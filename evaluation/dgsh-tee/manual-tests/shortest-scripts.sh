rm -f "#file2"
rm -f "#file4"
rm -f "#file6"
rm -f "#file8"
rm -f "#file10"
rm -f "#file12"
rm -f "#file14"
rm -f "#file17"
rm -f "#file18"
rm -f "#file19"
rm -f "#file20"
rm -f "#file21"
rm -f "#file22"
rm -f "#file23"
rm -f "#file24"
rm -f "#file25"
rm -f "#file26"
rm -f "#file27"
rm -f "#file28"
rm -f "#file29"
rm -f "#file30"
rm -f "#file31"
rm -f "#file32"
rm -f "#file33"
rm -f "#file34"
rm -f "#file35"
rm -f "#file36"
rm -f "#file37"
rm -f "#file38"
rm -f "#file39"
rm -f "#file40"
rm -f "#file41"
rm -f "#file42"
rm -f "#file46"
rm -f "#file43"
rm -f "#file47"
rm -f "#file44"
rm -f "#file48"
rm -f "#file45"
rm -f "#file49"
rm -f "#file50"
rm -f "#file51"
rm -f "#file52"
rm -f "#file53"
rm -f "#file54"
rm -f "#file55"
rm -f "#file56"
rm -f "#file57"
rm -f "#file58"
rm -f "#file59"
rm -f "#file60"
rm -f "#file61"
rm -f "#file62"
mkfifo "#file2"
mkfifo "#file4"
mkfifo "#file6"
mkfifo "#file8"
mkfifo "#file10"
mkfifo "#file12"
mkfifo "#file14"
mkfifo "#file17"
mkfifo "#file18"
mkfifo "#file19"
mkfifo "#file20"
mkfifo "#file21"
mkfifo "#file22"
mkfifo "#file23"
mkfifo "#file24"
mkfifo "#file25"
mkfifo "#file26"
mkfifo "#file27"
mkfifo "#file28"
mkfifo "#file29"
mkfifo "#file30"
mkfifo "#file31"
mkfifo "#file32"
mkfifo "#file33"
mkfifo "#file34"
mkfifo "#file35"
mkfifo "#file36"
mkfifo "#file37"
mkfifo "#file38"
mkfifo "#file39"
mkfifo "#file40"
mkfifo "#file41"
mkfifo "#file42"
mkfifo "#file46"
mkfifo "#file43"
mkfifo "#file47"
mkfifo "#file44"
mkfifo "#file48"
mkfifo "#file45"
mkfifo "#file49"
mkfifo "#file50"
mkfifo "#file51"
mkfifo "#file52"
mkfifo "#file53"
mkfifo "#file54"
mkfifo "#file55"
mkfifo "#file56"
mkfifo "#file57"
mkfifo "#file58"
mkfifo "#file59"
mkfifo "#file60"
mkfifo "#file61"
mkfifo "#file62"

mkfifo "#file63"
mkfifo "#file64"
mkfifo "#file65"
mkfifo "#file66"

{ cat $PASH_TOP/evaluation/scripts/input/1G.txt >"#file2" & }
{ $PASH_TOP/runtime/r_split "#file2" 10000000 "#file63" "#file64" "#file65" "#file66" & }

{ $PASH_TOP/runtime/dgsh_tee.sh "#file63" "#file17" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file64" "#file18" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file65" "#file19" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file66" "#file20" -I -f & }

{ $PASH_TOP/runtime/r_wrap xargs file <"#file17" >"#file22" & }
{ $PASH_TOP/runtime/r_wrap xargs file <"#file18" >"#file23" & }
{ $PASH_TOP/runtime/r_wrap xargs file <"#file19" >"#file24" & }
{ $PASH_TOP/runtime/r_wrap xargs file <"#file20" >"#file25" & }
{ $PASH_TOP/runtime/r_wrap grep "shell script" <"#file22" >"#file26" & }
{ $PASH_TOP/runtime/r_wrap grep "shell script" <"#file23" >"#file27" & }
{ $PASH_TOP/runtime/r_wrap grep "shell script" <"#file24" >"#file28" & }
{ $PASH_TOP/runtime/r_wrap grep "shell script" <"#file25" >"#file29" & }
{ $PASH_TOP/runtime/r_wrap cut -d: -f1 <"#file26" >"#file30" & }
{ $PASH_TOP/runtime/r_wrap cut -d: -f1 <"#file27" >"#file31" & }
{ $PASH_TOP/runtime/r_wrap cut -d: -f1 <"#file28" >"#file32" & }
{ $PASH_TOP/runtime/r_wrap cut -d: -f1 <"#file29" >"#file33" & }
{ $PASH_TOP/runtime/r_wrap xargs -L 1 wc -l <"#file30" >"#file34" & }
{ $PASH_TOP/runtime/r_wrap xargs -L 1 wc -l <"#file31" >"#file35" & }
{ $PASH_TOP/runtime/r_wrap xargs -L 1 wc -l <"#file32" >"#file36" & }
{ $PASH_TOP/runtime/r_wrap xargs -L 1 wc -l <"#file33" >"#file37" & }
{ $PASH_TOP/runtime/r_wrap grep -v "^0$" <"#file34" >"#file38" & }
{ $PASH_TOP/runtime/r_wrap grep -v "^0$" <"#file35" >"#file39" & }
{ $PASH_TOP/runtime/r_wrap grep -v "^0$" <"#file36" >"#file40" & }
{ $PASH_TOP/runtime/r_wrap grep -v "^0$" <"#file37" >"#file41" & }
{ $PASH_TOP/runtime/r_unwrap <"#file38" >"#file46" & }
{ sort -n <"#file53" >"#file42" & }
{ $PASH_TOP/runtime/r_unwrap <"#file39" >"#file47" & }
{ sort -n <"#file54" >"#file43" & }
{ $PASH_TOP/runtime/r_unwrap <"#file40" >"#file48" & }
{ sort -n <"#file55" >"#file44" & }
{ $PASH_TOP/runtime/r_unwrap <"#file41" >"#file49" & }
{ sort -n <"#file56" >"#file45" & }
{ sort -n -m "#file57" "#file58" >"#file50" & }
{ sort -n -m "#file59" "#file60" >"#file51" & }
{ sort -n -m "#file61" "#file62" >"#file14" & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file46" "#file53" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file47" "#file54" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file48" "#file55" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file49" "#file56" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file42" "#file57" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file43" "#file58" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file44" "#file59" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file45" "#file60" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file50" "#file61" -I -f & }
{ $PASH_TOP/runtime/dgsh_tee.sh "#file51" "#file62" -I -f & }
{ head -15 <"#file14" & }
source $PASH_TOP/runtime/wait_for_output_and_sigpipe_rest.sh ${!}
rm -f "#file2"
rm -f "#file4"
rm -f "#file6"
rm -f "#file8"
rm -f "#file10"
rm -f "#file12"
rm -f "#file14"
rm -f "#file17"
rm -f "#file18"
rm -f "#file19"
rm -f "#file20"
rm -f "#file21"
rm -f "#file22"
rm -f "#file23"
rm -f "#file24"
rm -f "#file25"
rm -f "#file26"
rm -f "#file27"
rm -f "#file28"
rm -f "#file29"
rm -f "#file30"
rm -f "#file31"
rm -f "#file32"
rm -f "#file33"
rm -f "#file34"
rm -f "#file35"
rm -f "#file36"
rm -f "#file37"
rm -f "#file38"
rm -f "#file39"
rm -f "#file40"
rm -f "#file41"
rm -f "#file42"
rm -f "#file46"
rm -f "#file43"
rm -f "#file47"
rm -f "#file44"
rm -f "#file48"
rm -f "#file45"
rm -f "#file49"
rm -f "#file50"
rm -f "#file51"
rm -f "#file52"
rm -f "#file53"
rm -f "#file54"
rm -f "#file55"
rm -f "#file56"
rm -f "#file57"
rm -f "#file58"
rm -f "#file59"
rm -f "#file60"
rm -f "#file61"
rm -f "#file62"

rm -f "#file63"
rm -f "#file64"
rm -f "#file65"
rm -f "#file66"