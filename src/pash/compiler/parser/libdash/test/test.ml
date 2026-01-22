let verbose = ref false
let input_src : string option ref = ref None

let set_input_src () =
  match !input_src with
  | None -> Dash.setinputtostdin ()
  | Some f -> Dash.setinputfile f

let parse_args () =
  Arg.parse
    ["-v",Arg.Set verbose,"verbose mode"]
    (function | "-" -> input_src := None | f -> input_src := Some f)
    "Final argument should be either a filename or - (for STDIN); only the last such argument is used"

exception Parse_error

let rec parse_all () : Ast.t list =
  let stackmark = Dash.init_stack () in
  match Dash.parse_next ~interactive:false () with
  | Dash.Done -> Dash.pop_stack stackmark; []
  | Dash.Error -> Dash.pop_stack stackmark; raise Parse_error
  | Dash.Null -> Dash.pop_stack stackmark; parse_all ()
  | Dash.Parsed n -> 
     (* translate to our AST *)
     let c = Ast.of_node n in
     (* deallocate *)
     Dash.pop_stack stackmark;
     (* keep calm and carry on *)
     c::parse_all ()

let main () = 
  Dash.initialize ();
  parse_args ();
  set_input_src ();
  let cs = parse_all () in
  List.map (fun c -> print_endline (Ast.to_string c)) cs
;;

main ()
