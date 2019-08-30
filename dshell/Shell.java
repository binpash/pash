package dshell;

import dshell.core.Graph;
import dshell.core.Operator;
import dshell.core.nodes.AtomicGraph;
import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.Sink;
import dshell.core.nodes.StatelessOperator;

import java.util.LinkedList;
import java.util.List;

public class Shell {
    public static Graph createGraph(String command) {
        String[] splitedCommand = command.split("\\|");
        List<Operator> pipeline = new LinkedList<>();

        for (String s : splitedCommand) {
            s = s.trim();

            int firstSpace = s.indexOf(' ');
            String program;
            String[] arguments;

            if (firstSpace != -1) {
                program = s.substring(0, firstSpace);
                arguments = s.substring(firstSpace + 1, s.length()).split(" ");
            } else {
                program = s;
                arguments = null;
            }

            pipeline.add(new StatelessOperator(program, arguments));
        }

        AtomicGraph[] graph = new AtomicGraph[pipeline.size() + 1];
        for (int i = 0; i < pipeline.size(); i++) {
            graph[i] = new AtomicGraph(pipeline.get(i));
        }
        graph[pipeline.size()] = new AtomicGraph(Sink.createPrinter());

        return new SerialGraph(graph);
    }

    public static void main(String[] args) {
        if (args.length < 1)
            throw new RuntimeException("Incorrect invocation argument.");

        String bashCommand = args[0];
        Graph graph = createGraph(bashCommand);

        graph.executeLocally();
    }
}