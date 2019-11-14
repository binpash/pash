package dshell.core.nodes;

import dshell.core.Operator;
import dshell.core.OperatorType;
import dshell.core.misc.SystemMessage;
import org.apache.commons.lang3.ArrayUtils;

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
        try {
            if (!initialized) {
                processBuilder = new ProcessBuilder(ArrayUtils.addAll(new String[]{program}, commandLineArguments));
                process = ((ProcessBuilder) processBuilder).start();

                ((ProcessBuilder) processBuilder).redirectError(ProcessBuilder.Redirect.INHERIT);

                if (process == null)
                    throw new Exception();

                initialized = true;
            }

            // passing the data to standard input of the process
            /*boolean endReceived = false;
            if (data != null) {
                OutputStream standardInput = ((Process) process).getOutputStream();

                if (data instanceof SystemMessage.EndOfData) {
                    endReceived = true;
                    standardInput.close();
                } else {
                    // write to standard input of the process
                    standardInput.write((byte[]) data);
                    standardInput.flush();
                    //return;   -> NOT NEEDED WHEN THE PROCESS IS CONSIDERED AS PIPELINE
                }
            } else {
                endReceived = true;
                ((Process) process).getOutputStream().close();
            }

            // getting the standard output of the process
            InputStream standardOutput = ((Process) process).getInputStream();

            // sending computed output to the next subscribed operator
            while (true) {
                if (standardOutput.available() >= BUFFER_SIZE) {
                    standardOutput.read(buffer, 0, BUFFER_SIZE);
                    consumers[0].next(0, buffer);
                } else if (standardOutput.available() > 0 && standardOutput.available() < BUFFER_SIZE) {
                    byte[] temp = new byte[standardOutput.available()];
                    standardOutput.read(temp, 0, temp.length);
                    consumers[0].next(0, temp);
                } else
                    break;
            }

            // send end of data signal
            if (endReceived) {
                consumers[0].next(0, new SystemMessage.EndOfData());
                ((Process) process).destroy();
            }*/



            boolean endReceived = false;
            if (data != null) {
                OutputStream standardInput = ((Process) process).getOutputStream();

                if (data instanceof SystemMessage.EndOfData) {
                    standardInput.close();
                    endReceived = true;
                } else
                    standardInput.write((byte[]) data);
            } else
                endReceived = true;

            if (endReceived) {
                byte[] standardOutput = ((Process) process).getInputStream().readAllBytes();
                consumers[0].next(0, standardOutput);
                consumers[0].next(0, new SystemMessage.EndOfData());
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            throw new RuntimeException(ex.getMessage());
        }
    }

    @Override
    protected Object clone() throws CloneNotSupportedException {
        StatelessOperator o = (StatelessOperator) super.clone();
        o.parallelizationHint = 1;

        return o;
    }
}