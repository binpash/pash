package dshell.core.nodes;

import dshell.core.Operator;

public abstract class Sink extends StatelessOperator {
    public Sink() {
        super(1, 0, null);
    }

    @Override
    public abstract void next(int inputChannel, Object data);
}
