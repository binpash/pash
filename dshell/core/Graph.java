package dshell.core;

import dshell.core.nodes.Sink;

public abstract class Graph {
    protected abstract Operator getOperator();

    public void executeLocally() {
        getOperator().next(null);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        Operator op = getOperator();

        while (op != null && !(op instanceof Sink)) {
            sb.append(op.getProgram());
            if (op.getCommandLineArguments() != null) {
                for (String s : op.getCommandLineArguments()) {
                    sb.append(" " + s);
                }
            }

            if (op.getConsumer() != null && !(op.getConsumer() instanceof Sink))
                sb.append(" | ");

            op = (Operator) op.getConsumer();
        }

        return sb.toString();
    }
}