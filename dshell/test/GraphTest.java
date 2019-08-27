package dshell.test;

import dshell.core.Operator;
import dshell.core.nodes.AtomicGraph;
import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.StatelessOperator;
import org.junit.Assert;
import org.junit.Test;

public class GraphTest {
    @Test
    public void wordCount() {
        Operator cat = new StatelessOperator("cat", new String[]{"~/Desktop/tekst.txt"});
        Operator wc = new StatelessOperator("wc", new String[]{"-w"});

        SerialGraph graph = new SerialGraph(new AtomicGraph(cat), new AtomicGraph(wc));

        //System.out.println(graph.toString());
        //Assert.assertTrue(graph.toString().equals("cat /home/cvetkovic/Desktop/example_text.txt | wc -w"));

        graph.executeLocally();
    }

    @Test
    public void grep() {
        Operator cat = new StatelessOperator("cat", new String[]{"$INPUT"});
        Operator wc = new StatelessOperator("grep", new String[]{"'[a-zA-Z0-9]\\+@[a-zA-Z0-9]\\+\\.[a-z]\\{2,\\}'"});

        SerialGraph graph = new SerialGraph(new AtomicGraph(cat), new AtomicGraph(wc));

        System.out.println(graph.toString());
        Assert.assertTrue(graph.toString().equals("cat $INPUT | grep '[a-zA-Z0-9]\\+@[a-zA-Z0-9]\\+\\.[a-z]\\{2,\\}'"));
    }

    @Test
    public void wordFrequencies() {
        Operator cat = new StatelessOperator("cat", new String[]{"$INPUT"});
        Operator tr1 = new StatelessOperator("tr", new String[]{"-cs A-Za-z'\n'"});
        Operator tr2 = new StatelessOperator("tr", new String[]{"A-Z a-z"});
        Operator sort1 = new StatelessOperator("sort");
        Operator uniq = new StatelessOperator("uniq", new String[]{"-c"});
        Operator sort2 = new StatelessOperator("sort", new String[]{"-rn"});
        Operator sed = new StatelessOperator("sed", new String[]{"${1}q"});

        SerialGraph graph = new SerialGraph(new AtomicGraph(cat),
                new AtomicGraph(tr1),
                new AtomicGraph(tr2),
                new AtomicGraph(sort1),
                new AtomicGraph(uniq),
                new AtomicGraph(sort2),
                new AtomicGraph(sed));

        System.out.println(graph.toString());
        Assert.assertTrue(graph.toString().equals("cat $INPUT | tr -cs A-Za-z'\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed ${1}q"));
    }

    @Test
    public void topNTerms() {
        Operator cat = new StatelessOperator("cat", new String[]{"$INPUT"});
        Operator tr1 = new StatelessOperator("tr", new String[]{"-cs A-Za-z '\n'"});
        Operator tr2 = new StatelessOperator("tr", new String[]{"A-Z a-z"});
        Operator sort1 = new StatelessOperator("sort");
        Operator uniq = new StatelessOperator("uniq", new String[]{"-c"});
        Operator sort2 = new StatelessOperator("sort", new String[]{"-rn"});
        Operator sed = new StatelessOperator("sed", new String[]{"1000q"});

        SerialGraph graph = new SerialGraph(new AtomicGraph(cat),
                new AtomicGraph(tr1),
                new AtomicGraph(tr2),
                new AtomicGraph(sort1),
                new AtomicGraph(uniq),
                new AtomicGraph(sort2),
                new AtomicGraph(sed));

        System.out.println(graph.toString());
        Assert.assertTrue(graph.toString().equals("cat $INPUT | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn | sed 1000q"));
    }
}