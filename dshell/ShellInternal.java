package dshell;

import dshell.core.OperatorFactory;
import dshell.core.nodes.AtomicGraph;
import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.StatelessOperator;

public class ShellInternal {

    private static String INPUT_FILE = "/home/cvetkovic/Desktop/i1M.txt";
    private static String HDFS_OUTPUT_FILE = "output.txt";

    public static void main(String[] args) {
        AtomicGraph cat = new AtomicGraph(new StatelessOperator(0, 1, "cat", new String[]{INPUT_FILE}));
        AtomicGraph wc = new AtomicGraph(new StatelessOperator(1, 1, "wc", new String[]{"-m"}, 2));
        AtomicGraph hdfsPrinter = new AtomicGraph(OperatorFactory.createHDFSFilePrinter("output.txt"));
        SerialGraph graph = new SerialGraph(cat, wc, hdfsPrinter);

        long time = 0;
        int cases = 5;

        for (int i = 0; i < cases; i++) {
            long start = System.currentTimeMillis();
            graph.executeRemote(3529);
            if (i > cases / 5)
                time += System.currentTimeMillis() - start;
        }

        System.out.println(time / (cases * 0.8));
    }
}