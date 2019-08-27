package dshell.core.nodes;

import dshell.core.Operator;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.OutputStream;

public class StatelessOperator extends Operator<String, String> {
    public StatelessOperator(String program) {
        super(program);
    }

    public StatelessOperator(String program, String[] commandLineArguments) {
        super(program, commandLineArguments);
    }

    @Override
    public void next(String data) {
        String execParameter = program + " " + getArgumentsAsString();
        StringBuilder processOutput = new StringBuilder();
        String s;

        try {
            ProcessBuilder processBuilder = new ProcessBuilder(execParameter);
            Process process = processBuilder.start();

            if (data != null) {
                OutputStream os = process.getOutputStream();
                os.write(data.getBytes());
            }

            BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            while ((s = bufferedReader.readLine()) != null)
                processOutput.append(s);

            process.waitFor();
            process.destroy();

            if (consumer != null)
                consumer.next(processOutput.toString());
            else
                System.out.println(processOutput.toString());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}