package dshell.core;

import dshell.core.interfaces.Consumer;
import dshell.core.interfaces.Producer;

import java.io.Serializable;

public abstract class Operator<A, B> implements Consumer<A>, Producer<B>, Serializable, Cloneable {
    protected int inputArity;
    protected int outputArity;

    protected OperatorType operatorType;

    protected Consumer<B>[] consumers;

    protected String program;
    protected String[] commandLineArguments;

    protected int parallelizationHint;

    public Operator(OperatorType operatorType, int inputArity, int outputArity, String program, int parallelizationHint) {
        this(operatorType, inputArity, outputArity, program, null, parallelizationHint);
    }

    public Operator(OperatorType operatorType, int inputArity, int outputArity, String program, String[] commandLineArguments, int parallelizationHint) {
        this.operatorType = operatorType;
        this.program = program;
        this.commandLineArguments = commandLineArguments;
        this.inputArity = inputArity;
        this.outputArity = outputArity;
        this.parallelizationHint = parallelizationHint;
    }

    @Override
    public abstract void next(int inputChannel, A data);

    @Override
    public void subscribe(Consumer<B>... consumers) {
        this.consumers = consumers;
    }

    public String getProgram() {
        return program;
    }

    public String[] getCommandLineArguments() {
        return commandLineArguments;
    }

    public Consumer<B>[] getConsumers() {
        return consumers;
    }

    public String getArgumentsAsString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < commandLineArguments.length; i++) {
            sb.append(commandLineArguments[i]);

            if (i != commandLineArguments.length - 1)
                sb.append(" ");
        }

        return sb.toString();
    }

    public int getInputArity() {
        return inputArity;
    }

    public int getOutputArity() {
        return outputArity;
    }

    @Override
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    public OperatorType getOperatorType() {
        return operatorType;
    }

    public int getParallelizationHint() {
        return parallelizationHint;
    }
}