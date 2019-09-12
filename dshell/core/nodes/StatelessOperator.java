package dshell.core.nodes;

import dshell.core.Operator;
import dshell.core.misc.Utilities;
import dshell.core.worker.RemoteExecutionData;

import java.io.OutputStream;
import java.util.Arrays;

public class StatelessOperator extends Operator<Object, Object> {
    private final int parallelizationHint;

    public StatelessOperator(int inputArity, int outputArity, String program) {
        this(inputArity, outputArity, program, null, 1);
    }

    public StatelessOperator(int inputArity, int outputArity, String program, String[] commandLineArguments) {
        this(inputArity, outputArity, program, commandLineArguments, 1);
    }

    public StatelessOperator(int inputArity, int outputArity, String program, String[] commandLineArguments, int parallelizationHint) {
        super(inputArity, outputArity, program, commandLineArguments);

        this.parallelizationHint = parallelizationHint;
    }

    @Override
    public void next(int inputChannel, Object data) {
        try {
            // creating new process
            ProcessBuilder processBuilder = new ProcessBuilder(program, getArgumentsAsString());
            Process process = processBuilder.start();

            // if the following condition hold then either the program does not exists or the command line arguments
            // were specified wrong
            // NOTE: arguments passing syntax can hereby be different than it would be as if it was typed in raw linux terminal
            if (!process.isAlive() && process.exitValue() != 0)
                throw new RuntimeException("Execution of program '" + program + "' returned with exit value " + process.exitValue() + ".");

            // passing the data to standard input of the process
            if (data != null) {
                OutputStream os = process.getOutputStream();

                // write to standard input of the process
                os.write((byte[])data);
                os.flush();
                os.close();
            }

            // getting the standard output of the process
            byte[] systemOut = process.getInputStream().readAllBytes();

            // waiting for process to finish and terminating it
            process.waitFor();
            process.destroy();

            // sending computed output to the next subscribed operator
            consumers[0].next(0, systemOut);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int getParallelizationHint() {
        return parallelizationHint;
    }
}