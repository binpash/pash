package dshell.core.nodes;

import dshell.core.Graph;
import dshell.core.Operator;

public class ParallelGraph extends Graph {
    private AtomicGraph[] atomicGraphs;

    public ParallelGraph(AtomicGraph... atomicGraphs) {
        this.atomicGraphs = atomicGraphs;
    }

    @Override
    public Operator getOperator() {
        throw new RuntimeException("This method was never supposed to be called.");
    }

    public AtomicGraph[] getAtomicGraphs() {
        return atomicGraphs;
    }
}