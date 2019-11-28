package dshell.core.nodes;

import dshell.core.Graph;
import dshell.core.Operator;

public class SerialGraph extends Graph {
    private AtomicGraph[] atomicGraphs;

    public SerialGraph(AtomicGraph... atomicGraphs) {
        this.atomicGraphs = atomicGraphs;
    }

    @Override
    public Operator getOperator() {
        if (atomicGraphs == null)
            throw new RuntimeException("The serial graph instance hasn't been initialized yet.");

        return atomicGraphs[0].getOperator();
    }

    public AtomicGraph[] getAtomicGraphs() {
        return atomicGraphs;
    }
}