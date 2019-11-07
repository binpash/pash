package dshell.core.nodes;

import dshell.core.Operator;
import dshell.core.OperatorType;
import dshell.core.misc.SystemMessage;

import java.io.InputStream;
import java.io.OutputStream;

public class StatelessOperator extends Operator<Object, Object> {
    private Object processBuilder;
    private Object process;
    private boolean initialized = false;
    private byte[] buffer = new byte[BUFFER_SIZE];

    // TODO: find out ideal packet length for TCP socket
    private static final int BUFFER_SIZE = 1024;

    public StatelessOperator(int inputArity, int outputArity, String program) {
        this(inputArity, outputArity, program, null, 1);
    }

    public StatelessOperator(int inputArity, int outputArity, String program, String[] commandLineArguments) {
        this(inputArity, outputArity, program, commandLineArguments, 1);
    }

    public StatelessOperator(int inputArity, int outputArity, String program, String[] commandLineArguments, int parallelizationHint) {
        super(OperatorType.STATELESS, inputArity, outputArity, program, commandLineArguments, parallelizationHint);
    }

    @Override
    public synchronized void next(int inputChannel, Object data) {
        if (!initialized) {
            try {
                processBuilder = new ProcessBuilder(program, getArgumentsAsString());
                process = ((ProcessBuilder) processBuilder).start();

                // if the following condition hold then either the program does not exists or the command line arguments
                // were specified wrong
                // NOTE: arguments passing syntax can here be different than it would be as if it was typed in raw linux terminal
                if (!((Process) process).isAlive() && ((Process) process).exitValue() != 0)
                    throw new RuntimeException("Execution of program '" + program + "' returned with exit value " + ((Process) process).exitValue() + ".");
            } catch (Exception ex) {
                throw new RuntimeException("Error creating process '" + program + "'.");
            }

            if (process == null)
                throw new RuntimeException("Child process was not created. Operator execution has been aborted.");

            initialized = true;
        }

        try {
            // passing the data to standard input of the process
            boolean endReceived = false;
            if (data != null) {
                OutputStream standardInput = ((Process) process).getOutputStream();

                if (data instanceof SystemMessage.EndOfData) {
                    endReceived = true;
                    //standardInput.close();
                } else {
                    // write to standard input of the process
                    standardInput.write((byte[]) data);
                    standardInput.close();
                    ;   //-> NOT NEEDED WHEN THE PROCESS IS CONSIDERED AS PIPELINE
                }
            }
            else if (data instanceof SystemMessage.EndOfData)
                endReceived = true;

            if (endReceived) {
                consumers[0].next(0, new SystemMessage.EndOfData());
                ((Process) process).destroy();
                return;
            }

            // getting the standard output of the process
            InputStream standardOutput = ((Process) process).getInputStream();

            byte[] raw = standardOutput.readAllBytes();

            // sending computed output to the next subscribed operator
            consumers[0].next(0, raw);
            /*while (standardOutput.available() != 0 && standardOutput.read(buffer, 0, BUFFER_SIZE) != -1)
                consumers[0].next(0, buffer);*/

            // send end of data signal
            if (data == null) {
                consumers[0].next(0, new SystemMessage.EndOfData());
                ((Process) process).destroy();
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e.getMessage());
        }
    }

    @Override
    protected Object clone() throws CloneNotSupportedException {
        StatelessOperator o = (StatelessOperator) super.clone();
        o.parallelizationHint = 1;

        return o;
    }
}