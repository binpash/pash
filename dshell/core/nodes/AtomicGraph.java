package dshell.core.nodes;

import dshell.core.Graph;
import dshell.core.Operator;

public class AtomicGraph extends Graph {
    private Operator operator;

    public AtomicGraph(Operator operator) {
        this.operator = operator;
    }

    @Override
    public Operator getOperator() {
        return operator;
    }
}