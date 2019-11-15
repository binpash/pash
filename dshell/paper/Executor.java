package dshell.paper;

import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class Executor {
    private static final String FILE_NODES = "fids";
    private static final String PIPE_NODES = "nodes";
    private static final String INPUT_FILES = "in";
    private static final String OUTPUT_FILES = "out";
    private static final String COMMAND = "command";

    /*private static class Node {
        String[] inputFiles;
        String[] outputFiles;
        String command;
    }*/

    public static void main(String[] args) throws Exception {
        File jsonFile = new File(args[0]);
        StringBuilder jsonRaw = new StringBuilder();

        try (BufferedReader bufferedReader = new BufferedReader(new FileReader(jsonFile))) {
            String data;
            while ((data = bufferedReader.readLine()) != null)
                jsonRaw.append(data);
        }

        JSONObject rootElement = new JSONObject(jsonRaw.toString());
        JSONArray files = rootElement.getJSONArray(FILE_NODES);
        JSONObject nodes = rootElement.getJSONObject(PIPE_NODES);

        for (int i = 0; i < files.length(); i++) {
            String file = files.getString(i);

            remoteFIFOIfExists(file);
            createFIFO(file);
        }

        //Map<Integer, Node> topology = new HashMap<>();
        Iterator<String> keys = nodes.keys();
        while (keys.hasNext()) {
            String key = keys.next();
            JSONObject object = nodes.getJSONObject(key);

            JSONArray rawInput = object.getJSONArray(INPUT_FILES);
            JSONArray rawOutput = object.getJSONArray(OUTPUT_FILES);
            String command = object.getString(COMMAND);

            List<String> input = new ArrayList<>();
            for (int i = 0; i < rawInput.length(); i++)
                input.add(rawInput.getString(i).substring(1));
            List<String> output = new ArrayList<>();
            for (int i = 0; i < rawOutput.length(); i++)
                output.add(rawOutput.getString(i).substring(1));


            createNode(input, output, command);
        }
    }

    private static void createNode(List<String> inputFiles, List<String> outputFiles, String command) throws Exception {
        List<String> args = compileArgumentList(inputFiles, outputFiles, command);
        executeNonBlocking(args);
    }

    private static List<String> compileArgumentList(List inputFiles, List outputFiles, String command) {
        List<String> argumentsToPass = new ArrayList<>();

        if (inputFiles.size() != 0) {
            argumentsToPass.add("cat");
            argumentsToPass.addAll(inputFiles);
            argumentsToPass.add("|");
        }

        argumentsToPass.addAll(splitCommand(command));

        if (outputFiles.size() != 0) {
            argumentsToPass.add("|");
            argumentsToPass.add("tee");
            argumentsToPass.addAll(outputFiles);
        }

        return argumentsToPass;
    }

    private static void remoteFIFOIfExists(String file) throws Exception {
        executeNonBlocking("rm", "-f", file);
    }

    private static void createFIFO(String file) throws Exception {
        executeNonBlocking("mkfifo", file.substring(1));
    }

    private static void executeBlocking(String... command) throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder(command);
        Process process = processBuilder.start();

        process.waitFor();
    }

    private static void executeBlocking(List<String> command) throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder(command);
        Process process = processBuilder.start();

        process.waitFor();
    }

    private static void executeNonBlocking(String... command) throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder(command);
        Process process = processBuilder.start();
    }

    private static void executeNonBlocking(List<String>command) throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder(command);
        Process process = processBuilder.start();
    }

    private static List<String> splitCommand(String command) {
        String[] t = command.split(" ");
        return Arrays.asList(t);
    }
}