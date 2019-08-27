package dshell.core.nodes;

import dshell.core.Graph;
import dshell.core.Operator;

public class SerialGraph extends Graph {
    private AtomicGraph[] atomicGraphs;

    public SerialGraph(AtomicGraph... atomicGraphs) {
        this.atomicGraphs = atomicGraphs;

        createSubscriptions();
    }

    @Override
    protected Operator getOperator() {
        if (atomicGraphs == null)
            throw new RuntimeException("The serial graph instance hasn't been initialized yet.");

        return atomicGraphs[0].getOperator();
    }

    private void createSubscriptions() {
        for (int i = 0; i < atomicGraphs.length - 1; i++) {
            atomicGraphs[i].getOperator().subscribe(atomicGraphs[i + 1].getOperator());
        }
    }
}