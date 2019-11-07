(* This is straight-up copied from the libdash tests *)

let verbose = ref false
let input_src : string option ref = ref None

let parse_args () =
  Arg.parse
    [("-v",Arg.Set verbose,"verbose mode")]
    (function | "-" -> input_src := None | f -> input_src := Some f)
    "Final argument should be either a filename or empty (for STDIN); only the last such argument is used"

let read_channel chan =
let lines = ref [] in
try
  while true; do
    lines := input_line chan :: !lines
  done; !lines
with End_of_file ->
  close_in chan;
  List.rev !lines

let read_lines () =
  match !input_src with
  | None -> read_channel stdin
  | Some filename -> read_channel (open_in filename)

let parse_lines () : Ast.t list =
  let lines = read_lines () in
  List.map (fun line -> Ast_json.t_of_string line) lines


let main () = 
  parse_args ();
  let cs = parse_lines () in
  List.map (fun c -> print_endline (Ast.to_string c)) cs
;;

main ()
