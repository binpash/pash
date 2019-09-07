package dshell.test;

import dshell.Shell;
import dshell.core.DFileSystem;
import dshell.core.Graph;
import dshell.core.Operator;
import dshell.core.nodes.AtomicGraph;
import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.Sink;
import dshell.core.nodes.StatelessOperator;
import org.junit.Assert;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

public class GraphTest {
    private static String INPUT_FILE = "/home/cvetkovic/sdsh/scripts/input.txt";
    private static String HDFS_OUTPUT_FILE = "output.txt";

    /*
    @Test
    public void wordCount() {
        String command = "cat " + INPUT_FILE + " | wc -w";
        Graph graph = Shell.createGraph(command);
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        System.setOut(new PrintStream(baos));

        graph.executeLocallySingleThreaded();
        String output = baos.toString();
        Assert.assertTrue(output.equals("2923\n"));
    }*/

    @Test
    public void wordCountDistributed() throws Exception {
        String command = "cat " + INPUT_FILE + " | wc -w";
        Graph graph = Shell.createGraph(command, HDFS_OUTPUT_FILE);
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        System.setOut(new PrintStream(baos));

        graph.executeDistributed("output.txt");

        Thread.sleep(15000);

        byte[] file = DFileSystem.downloadFile("output.txt");
        Assert.assertTrue((new String(file)).equals("566311\n"));
    }

    /*
    @Test
    public void grep() {
        String command = "cat " + INPUT_FILE + " | grep [a-zA-Z0-9]\\+@[a-zA-Z0-9]\\+\\.[a-z]\\{2,\\}";
        Graph graph = Shell.createGraph(command);

        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        System.setOut(new PrintStream(baos));

        graph.executeLocallySingleThreaded();

        String output = baos.toString();
        Assert.assertTrue(output.equals("tutabugarin@gmail.com\r\n"));
    }*/

    /*
    @Test
    public void wordFrequencies() {
        Operator cat = new StatelessOperator("cat", new String[]{INPUT_FILE});
        Operator tr1 = new StatelessOperator("tr", new String[]{"-cs A-Za-z'\n'"});
        Operator tr2 = new StatelessOperator("tr", new String[]{"A-Z a-z"});
        Operator sort1 = new StatelessOperator("sort");
        Operator uniq = new StatelessOperator("uniq", new String[]{"-c"});
        Operator sort2 = new StatelessOperator("sort", new String[]{"-rn"});
        Operator sed = new StatelessOperator("sed", new String[]{"${1}q"});
        Sink sink = Sink.createPrinter();

        SerialGraph graph = new SerialGraph(new AtomicGraph(cat),
                new AtomicGraph(tr1),
                new AtomicGraph(tr2),
                new AtomicGraph(sort1),
                new AtomicGraph(uniq),
                new AtomicGraph(sort2),
                new AtomicGraph(sed),
                new AtomicGraph(sink));

        System.out.println(graph.toString());
        Assert.assertTrue(graph.toString().equals("cat $INPUT | tr -cs A-Za-z'\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed ${1}q"));
    }

    @Test
    public void topNTerms() {
        Operator cat = new StatelessOperator("cat", new String[]{INPUT_FILE});
        Operator tr1 = new StatelessOperator("tr", new String[]{"-cs A-Za-z '\n'"});
        Operator tr2 = new StatelessOperator("tr", new String[]{"A-Z a-z"});
        Operator sort1 = new StatelessOperator("sort");
        Operator uniq = new StatelessOperator("uniq", new String[]{"-c"});
        Operator sort2 = new StatelessOperator("sort", new String[]{"-rn"});
        Operator sed = new StatelessOperator("sed", new String[]{"1000q"});
        Sink sink = Sink.createPrinter();

        SerialGraph graph = new SerialGraph(new AtomicGraph(cat),
                new AtomicGraph(tr1),
                new AtomicGraph(tr2),
                new AtomicGraph(sort1),
                new AtomicGraph(uniq),
                new AtomicGraph(sort2),
                new AtomicGraph(sed),
                new AtomicGraph(sink));

        System.out.println(graph.toString());
        Assert.assertTrue(graph.toString().equals("cat $INPUT | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed 1000q"));
    }*/
}