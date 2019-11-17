package dshell.paper;

import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class Executor {
    private static final String FILE_IDENTIFIERS = "fids";
    private static final String PIPE_NODES = "nodes";
    private static final String NODE_INPUT_FILES = "in";
    private static final String NODE_OUTPUT_FILES = "out";
    private static final String COMMAND = "command";
    private static final String INPUT_FILES = "in";
    private static final String OUTPUT_FILES = "out";

    private static List<ProcessBuilder> processBuilders = new ArrayList<>();
    private static List<Process> processes = new ArrayList<>();

    public static void main(String[] args) throws Exception {
        File jsonFile = new File(args[0]);
        int iterations = 1;
        long executionTime = 0;
        if (args.length > 1)
            iterations = Integer.parseInt(args[1]);

        StringBuilder jsonRaw = new StringBuilder();

        try (BufferedReader bufferedReader = new BufferedReader(new FileReader(jsonFile))) {
            String data;
            while ((data = bufferedReader.readLine()) != null)
                jsonRaw.append(data);
        }

        JSONObject rootElement = new JSONObject(jsonRaw.toString());
        JSONArray files = rootElement.getJSONArray(FILE_IDENTIFIERS);
        JSONObject nodes = rootElement.getJSONObject(PIPE_NODES);
        JSONArray inFiles = rootElement.getJSONArray(INPUT_FILES);
        JSONArray outFiles = rootElement.getJSONArray(OUTPUT_FILES);

        Iterator<String> keys = nodes.keys();
        while (keys.hasNext()) {
            String key = keys.next();
            JSONObject object = nodes.getJSONObject(key);

            JSONArray rawInput = object.getJSONArray(NODE_INPUT_FILES);
            JSONArray rawOutput = object.getJSONArray(NODE_OUTPUT_FILES);
            String command = object.getString(COMMAND);

            List<String> input = new ArrayList<>();
            for (int i = 0; i < rawInput.length(); i++)
                input.add(rawInput.getString(i).substring(1));
            List<String> output = new ArrayList<>();
            for (int i = 0; i < rawOutput.length(); i++)
                output.add(rawOutput.getString(i).substring(1));

            createNode(input, output, command);
        }

        for (int x = 0; x < iterations; x++) {
            for (int i = 0; i < files.length(); i++) {
                String file = files.getString(i).substring(1);

                remoteFIFOIfExists(file);
                createFIFO(file);
            }

            long startTime = System.currentTimeMillis();
            for (ProcessBuilder p : processBuilders)
                processes.add(p.start());


            for (int i = 0; i < outFiles.length(); i++) {
                String fileToWait = outFiles.getString(i).substring(1);
                executeBlocking(new String[]{"cat", fileToWait});
            }

            for (Process p : processes)
                p.waitFor();
            long endTime = System.currentTimeMillis();

            executionTime += endTime - startTime;
        }

        System.out.println("Script execution time is " + executionTime + " ms");
        if (iterations > 1)
            System.out.println("Script average execution time is " + executionTime / (double) iterations + " ms");
    }

    private static void createNode(List<String> inputFiles, List<String> outputFiles, String command) throws Exception {
        String[] args = compileArgumentList(inputFiles, outputFiles, command);
        createProcessBuilder(args);
    }

    private static String[] compileArgumentList(List<String> inputFiles, List<String> outputFiles, String command) {
        String[] argumentsToPass = new String[3];
        StringBuilder script = new StringBuilder();

        argumentsToPass[0] = "/bin/sh";
        argumentsToPass[1] = "-c";

        if (inputFiles.size() != 0) {
            script.append("cat ");
            for (String s : inputFiles)
                script.append(s + " ");
            script.append("| ");
        }

        script.append(command.trim());

        if (outputFiles.size() == 1) {
            script.append(" > ");
            script.append(outputFiles.get(0));
        } else if (outputFiles.size() > 1) {
            script.append(" | tee ");
            for (String s : inputFiles)
                script.append(s + " ");
        }

        argumentsToPass[2] = script.toString().trim();
        return argumentsToPass;
    }

    private static void remoteFIFOIfExists(String file) throws Exception {
        executeNonBlocking("rm", "-f", file);
    }

    private static void createFIFO(String file) throws Exception {
        executeNonBlocking("mkfifo", file);
    }

    private static void executeNonBlocking(String... command) throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder(command);
        Process process = processBuilder.start();
    }

    private static void createProcessBuilder(String... command) throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder(command);
        processBuilders.add(processBuilder);
    }

    private static void executeBlocking(String... command) throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder(command);
        Process process = processBuilder.start();

        byte[] data = process.getInputStream().readAllBytes();
        System.out.println(new String(data));

        process.waitFor();
    }
}