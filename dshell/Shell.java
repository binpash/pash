package dshell;

import dshell.core.Graph;
import dshell.core.Operator;
import dshell.core.OperatorFactory;
import dshell.core.nodes.AtomicGraph;
import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.Sink;
import dshell.core.nodes.StatelessOperator;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

public class Shell {
    public static Graph createGraph(String command, String outputFile) {
        String[] splitedCommand = command.split("\\|");
        List<Operator> pipeline = new LinkedList<>();

        for (String s : splitedCommand) {
            s = s.trim();

            int firstSpace = s.indexOf(' ');
            String program;
            String[] arguments;

            if (firstSpace != -1) {
                program = s.substring(0, firstSpace);
                arguments = s.substring(firstSpace + 1).split(" ");
            } else {
                program = s;
                arguments = null;
            }

            pipeline.add(new StatelessOperator(program, arguments));
        }

        AtomicGraph[] graph = new AtomicGraph[pipeline.size() + 1];
        for (int i = 0; i < pipeline.size(); i++)
            graph[i] = new AtomicGraph(pipeline.get(i));
        graph[pipeline.size()] = new AtomicGraph(OperatorFactory.createHDFSFilePrinter(outputFile));

        return new SerialGraph(graph);
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 1)
            throw new RuntimeException("Incorrect invocation argument.");

        //setupSlaves();

        // TODO: fix main method
        /*String bashCommand = args[0];
        Graph graph = createGraph(bashCommand);

        graph.executeLocallySingleThreaded();*/
    }

    /*private static void setupSlaves() throws IOException {
        ProcessBuilder processBuilder = new ProcessBuilder("");
        Process p1 = processBuilder.start();
    }*/
}