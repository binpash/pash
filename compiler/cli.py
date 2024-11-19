import argparse
import os


class BaseParser(argparse.ArgumentParser):
    """
    Base class for all Argument Parsers used by PaSh. It has two configurable flags
    by default: debug and log_file.

    Other flags are available by classes which inherit BaseParser
    """

    @staticmethod
    def _get_width():
        cpus = os.cpu_count()
        assert cpus is not None
        return cpus // 8 if cpus >= 16 else 2

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_argument(
            "-t",
            "--output_time",  # FIXME: --time
            help="(obsolete, time is always logged now) output the time it took for every step",
            action="store_true",
        )
        self.add_argument(
            "-d",
            "--debug",
            type=int,
            help="configure debug level; defaults to 0",
            default=0,
        )
        self.add_argument(
            "--log_file",
            help="configure where to write the log; defaults to stderr.",
            default="",
        )

    def add_pash_args(self):
        self.add_argument(
            "-w",
            "--width",
            type=int,
            default=self._get_width(),
            help="set data-parallelism factor",
        )
        self.add_argument(
            "--no_optimize",
            help="not apply transformations over the DFG",
            action="store_true",
        )
        self.add_argument(
            "--dry_run_compiler",
            help="not execute the compiled script, even if the compiler succeeded",
            action="store_true",
        )
        self.add_argument(
            "--assert_compiler_success",
            help="assert that the compiler succeeded with no general error occuring",
            action="store_true",
        )
        self.add_argument(
            "--assert_all_regions_parallelizable",
            help="assert that the compiler succeeded with all regions being parallelizable and no general error occuring (used to make tests more robust); more strict than --assert_compiler_success flag",
            action="store_true",
        )
        self.add_argument(
            "--avoid_pash_runtime_completion",
            help="avoid the pash_runtime execution completion (only relevant when --debug > 0)",
            action="store_true",
        )
        self.add_argument(
            "-p",
            "--output_optimized",  # FIXME: --print
            help="output the parallel shell script for inspection",
            action="store_true",
        )
        self.add_argument(
            "--graphviz",
            help="generates graphical representations of the dataflow graphs. The option argument corresponds to the format. PaSh stores them in a timestamped directory in the argument of --graphviz_dir",
            choices=["no", "dot", "svg", "pdf", "png"],
            default="no",
        )
        ## TODO: To discuss: Do we maybe want to have graphviz to always be included
        ##       in the temp directory (under a graphviz subdirectory) instead of in its own?
        ##   kk: I think that ideally we want a log-directory where we can put logs, graphviz,
        ##       and other observability and monitoring info (instead of putting them in the temp).
        self.add_argument(
            "--graphviz_dir",
            help="the directory in which to store graphical representations",
            default="/tmp",
        )
        self.add_argument(
            "--no_parallel_pipelines",
            help="Disable parallel running of independent pipelines",
            action="store_true",
            default=False,
        )
        self.add_argument(
            "--parallel_pipelines_limit",
            help="Maximum number of parallel independent pipelines",
            type=int,
            default=2,
        )
        self.add_argument(
            "--r_split_batch_size",
            type=int,
            help="configure the batch size of r_split (default: 1MB)",
            default=1000000,
        )
        self.add_argument(
            "--config_path",
            help="determines the config file path. By default it is 'PASH_TOP/compiler/config.yaml'.",
            default="",
        )
        self.add_argument(
            "--version",
            action="version",
            version="%(prog)s {version}".format(
                version="0.12.2"
            ),  # What does this version mean?
        )

        self.add_experimental_args()
        self.add_obsolete_args()

    def add_obsolete_args(self):
        self.add_argument(
            "--no_daemon",
            help="(obsolete) does nothing -- Run the compiler everytime we need a compilation instead of using the daemon",
            action="store_true",
            default=False,
        )
        self.add_argument(
            "--parallel_pipelines",
            help="(obsolete) Run multiple pipelines in parallel if they are safe to run. Now true by default. See --no_parallel_pipelines.",
            action="store_true",
            default=True,
        )
        self.add_argument(
            "--r_split",
            help="(obsolete) does nothing -- only here for old interfaces (not used anywhere in the code)",
            action="store_true",
        )
        self.add_argument(
            "--dgsh_tee",
            help="(obsolete) does nothing -- only here for old interfaces (not used anywhere in the code)",
            action="store_true",
        )
        self.add_argument(
            "--speculation",
            help="(obsolete) does nothing -- run the original script during compilation; if compilation succeeds, abort the original and run only the parallel (quick_abort) (Default: no_spec)",
            choices=["no_spec", "quick_abort"],
            default="no_spec",
        )

    def add_experimental_args(self):
        self.add_argument(
            "--no_eager",
            help="(experimental) disable eager nodes before merging nodes",
            action="store_true",
        )
        self.add_argument(
            "--profile_driven",
            help="(experimental) use profiling information when optimizing",
            action="store_true",
        )
        self.add_argument(
            "--speculative",
            help="(experimental) use the speculative execution preprocessing and runtime (NOTE: this has nothing to do with --speculation, which is actually misnamed, and should be named concurrent compilation/execution and is now obsolete)",
            action="store_true",
            default=False,
        )
        self.add_argument(
            "--termination",
            help="(experimental) determine the termination behavior of the DFG. Defaults to cleanup after the last process dies, but can drain all streams until depletion",
            choices=["clean_up_graph", "drain_stream"],
            default="clean_up_graph",
        )
        self.add_argument(
            "--daemon_communicates_through_unix_pipes",
            help="(experimental) the daemon communicates through unix pipes instead of sockets",
            action="store_true",
        )
        self.add_argument(
            "--distributed_exec",
            help="(experimental) execute the script in a distributed environment. Remote machines should be configured and ready",
            action="store_true",
            default=False,
        )


class RunnerParser(BaseParser):
    """
    Parser for the PaSh Runner in compiler/pash.py
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_pash_args()

        self.add_argument(
            "input",
            nargs="*",
            help="the script to be compiled and executed (followed by any command-line arguments",
        )
        self.add_argument(
            "--preprocess_only",
            help="only preprocess the input script and not execute it",
            action="store_true",
        )
        self.add_argument(
            "--output_preprocessed",
            help=" output the preprocessed script",
            action="store_true",
        )
        self.add_argument(
            "--interactive",
            help="Executes the script using an interactive internal shell session (experimental)",
            action="store_true",
        )
        self.add_argument(
            "-c",
            "--command",
            help="Evaluate the following as a script, rather than a file",
            default=None,
        )
        ## This is not the correct way to parse these, because more than one option can be given together, e.g., -ae
        self.add_argument(
            "-a",
            help="Enabling the `allexport` shell option",
            action="store_true",
            default=False,
        )
        self.add_argument(
            "+a",
            help="Disabling the `allexport` shell option",
            action="store_false",
            default=False,
        )
        ## These two are here for compatibility with respect to bash
        self.add_argument(
            "-v",
            help="(experimental) prints shell input lines as they are read",
            action="store_true",
        )
        self.add_argument(
            "-x",
            help="(experimental) prints commands and their arguments as they execute",
            action="store_true",
        )
        self.set_defaults(preprocess_mode="pash")


class CompilerParser(BaseParser):
    """
    Parser for the PaSh compiler in compiler/pash_compiler.py
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_pash_args()

        self.add_argument(
            "compiled_script_file",
            help="the file in which to output the compiled script",
        )
        self.add_argument(
            "input_ir",
            help="the file containing the dataflow graph to be optimized and executed",
        )
        self.add_argument(
            "--var_file",
            help="determines the path of a file containing all shell variables.",
            default=None,
        )


class PreprocessorParser(BaseParser):
    """
    Parser for the preprocessor in compiler/preprocessor/preprocessor.py
    Generates two subparsers
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        subparser = self.add_subparsers(help="sub-command help")
        self.add_pash_subparser(subparser)
        self.add_spec_subparser(subparser)

    @staticmethod
    def add_pash_subparser(subparser):
        parser_pash = subparser.add_parser(
            "pash", help="Preprocess the script so that it can be run with PaSh"
        )
        parser_pash.add_pash_args()
        parser_pash.add_argument("input", help="the script to be preprocessed")
        parser_pash.set_defaults(preprocess_mode="pash")

    @staticmethod
    def add_spec_subparser(subparser):
        # create the parser for the "b" command
        parser_spec = subparser.add_parser(
            "spec", help="Preprocess the script so that it can be run with speculation"
        )
        parser_spec.add_argument("input", help="the script to be preprocessed")

        ## TODO: When we better integrate, this should be automatically set.
        parser_spec.add_argument(
            "partial_order_file",
            help="the file to store the partial order (currently just a sequence)",
        )
        parser_spec.set_defaults(preprocess_mode="spec")
