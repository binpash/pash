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
rm -f "#file32"
rm -f "#file31"
rm -f "#file33"
rm -f "#file34"
rm -f "#file35"
rm -f "#file36"
rm -f "#file37"
rm -f "#file38"
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
mkfifo "#file32"
mkfifo "#file31"
mkfifo "#file33"
mkfifo "#file34"
mkfifo "#file35"
mkfifo "#file36"
mkfifo "#file37"
mkfifo "#file38"
{ cat /home/tamlu/pash/evaluation/scripts/input/all_cmds_x100.txt >"#file2" & }
{ /home/tamlu/pash/runtime/r_split "#file2" 100000 "#file17" "#file18" & }
{ /home/tamlu/pash/runtime/r_wrap xargs file <"#file17" >"#file20" & }
{ /home/tamlu/pash/runtime/r_wrap xargs file <"#file18" >"#file21" & }
{ /home/tamlu/pash/runtime/r_wrap grep "shell script" <"#file20" >"#file22" & }
{ /home/tamlu/pash/runtime/r_wrap grep "shell script" <"#file21" >"#file23" & }
{ /home/tamlu/pash/runtime/r_wrap cut -d: -f1 <"#file22" >"#file24" & }
{ /home/tamlu/pash/runtime/r_wrap cut -d: -f1 <"#file23" >"#file25" & }
{ /home/tamlu/pash/runtime/r_wrap xargs -L 1 wc -l <"#file24" >"#file26" & }
{ /home/tamlu/pash/runtime/r_wrap xargs -L 1 wc -l <"#file25" >"#file27" & }
{ /home/tamlu/pash/runtime/r_wrap grep -v "^0$" <"#file26" >"#file28" & }
{ /home/tamlu/pash/runtime/r_wrap grep -v "^0$" <"#file27" >"#file29" & }
{ /home/tamlu/pash/runtime/r_unwrap <"#file28" >"#file32" & }
{ sort -n <"#file35" >"#file30" & }
{ /home/tamlu/pash/runtime/r_unwrap <"#file29" >"#file33" & }
{ sort -n <"#file36" >"#file31" & }
{ sort -n -m "#file37" "#file38" >"#file14" & }
{ /home/tamlu/pash/runtime/eager.sh "#file32" "#file35" "/tmp/pash_eager_intermediate_#file1" & }
{ /home/tamlu/pash/runtime/eager.sh "#file33" "#file36" "/tmp/pash_eager_intermediate_#file2" & }
{ /home/tamlu/pash/runtime/eager.sh "#file30" "#file37" "/tmp/pash_eager_intermediate_#file3" & }
{ /home/tamlu/pash/runtime/eager.sh "#file31" "#file38" "/tmp/pash_eager_intermediate_#file4" & }
{ head -15 <"#file14" & }
source /home/tamlu/pash/runtime/wait_for_output_and_sigpipe_rest.sh ${!}
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
rm -f "#file32"
rm -f "#file31"
rm -f "#file33"
rm -f "#file34"
rm -f "#file35"
rm -f "#file36"
rm -f "#file37"
rm -f "#file38"