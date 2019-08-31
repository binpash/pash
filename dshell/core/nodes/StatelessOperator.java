package dshell.core.nodes;

import dshell.core.DFileSystem;
import dshell.core.Operator;

import java.io.OutputStream;

public class StatelessOperator extends Operator<String, String> {
    public StatelessOperator(String program) {
        super(program);
    }

    public StatelessOperator(String program, String[] commandLineArguments) {
        super(program, commandLineArguments);
    }

    @Override
    public void next(String inputFilename) {
        try {
            // creating new process
            ProcessBuilder processBuilder = new ProcessBuilder(program, getArgumentsAsString());
            Process process = processBuilder.start();

            // getting standard input
            if (inputFilename != null) {
                OutputStream os = process.getOutputStream();
                byte[] inputData = DFileSystem.downloadFile(inputFilename);

                // write to standard input of the process
                os.write(inputData);
                os.flush();
                os.close();
            }

            // writing standard output to HDFS file
            String outputFilename = DFileSystem.generateFilename();
            byte[] systemOut = process.getInputStream().readAllBytes();
            DFileSystem.uploadFile(outputFilename, systemOut);

            // waiting for process to finish and terminating it
            process.waitFor();
            process.destroy();

            consumer.next(outputFilename);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}