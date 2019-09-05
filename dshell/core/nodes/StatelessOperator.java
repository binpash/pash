package dshell.core.nodes;

import dshell.core.DFileSystem;
import dshell.core.Operator;

import java.io.OutputStream;

public class StatelessOperator extends Operator<byte[], byte[]> {
    private final int parallelizationHint;

    public StatelessOperator(String program) {
        super(program);
        this.parallelizationHint = 1;
    }

    public StatelessOperator(String program, int parallelizationHint) {
        super(program);
        this.parallelizationHint = parallelizationHint;
    }

    public StatelessOperator(String program, String[] commandLineArguments) {
        super(program, commandLineArguments);
        this.parallelizationHint = 1;
    }

    public StatelessOperator(String program, String[] commandLineArguments, int parallelizationHint) {
        super(program, commandLineArguments);
        this.parallelizationHint = parallelizationHint;
    }

    @Override
    public void next(byte[] data) {
        try {
            // creating new process
            ProcessBuilder processBuilder = new ProcessBuilder(program, getArgumentsAsString());
            Process process = processBuilder.start();

            // getting standard input
            if (data != null) {
                OutputStream os = process.getOutputStream();

                // write to standard input of the process
                os.write(data);
                os.flush();
                os.close();
            }

            // writing standard output to HDFS file
            String outputFilename = DFileSystem.generateFilename();
            byte[] systemOut = process.getInputStream().readAllBytes();

            // waiting for process to finish and terminating it
            process.waitFor();
            process.destroy();

            consumer.next(systemOut);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int getParallelizationHint() {
        return parallelizationHint;
    }
}