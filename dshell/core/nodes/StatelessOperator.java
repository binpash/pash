package dshell.core.nodes;

import dshell.core.Operator;

import java.io.OutputStream;

public class StatelessOperator extends Operator<byte[], byte[]> {
    public StatelessOperator(String program) {
        super(program);
    }

    public StatelessOperator(String program, String[] commandLineArguments) {
        super(program, commandLineArguments);
    }

    @Override
    public void next(byte[] data) {
        try {
            ProcessBuilder processBuilder = new ProcessBuilder(program, getArgumentsAsString());
            Process process = processBuilder.start();

            if (data != null) {
                OutputStream os = process.getOutputStream();
                os.write(data);
                os.flush();
                os.close();
            }

            byte[] systemOut = process.getInputStream().readAllBytes();

            process.waitFor();
            process.destroy();

            consumer.next(systemOut);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}