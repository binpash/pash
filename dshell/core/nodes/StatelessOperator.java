package dshell.core.nodes;

import dshell.core.Operator;
import dshell.core.worker.RemoteExecutionData;

import java.io.OutputStream;

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

            if (!process.isAlive() && process.exitValue() != 0)
                throw new RuntimeException("Execution of program '" + program + "' returned with exit value " + process.exitValue() + ".");

            // getting standard input
            if (data != null) {
                OutputStream os = process.getOutputStream();
                RemoteExecutionData d = (RemoteExecutionData) data;

                // write to standard input of the process
                os.write(d.getInputData());
                os.flush();
                os.close();
            }

            byte[] systemOut = process.getInputStream().readAllBytes();

            // waiting for process to finish and terminating it
            process.waitFor();
            process.destroy();

            // send output of current operator to nextOperator that will output its output to localhost:4000
            RemoteExecutionData red = new RemoteExecutionData(nextOperator, systemOut, "localhost", 4000, "localhost", 3529);
            consumers[0].next(0, red);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int getParallelizationHint() {
        return parallelizationHint;
    }
}