package dshell.core;

import dshell.core.interfaces.Consumer;
import dshell.core.interfaces.Producer;

import java.io.Closeable;
import java.io.Serializable;

public abstract class Operator<A, B> implements Consumer<A>, Producer<B>, Serializable, Cloneable {
    protected int inputArity;
    protected int outputArity;

    protected Consumer<B>[] consumers;

    protected String program;
    protected String[] commandLineArguments;

    public Operator(int inputArity, int outputArity, String program) {
        this(inputArity, outputArity, program, null);
    }

    public Operator(int inputArity, int outputArity, String program, String[] commandLineArguments) {
        this.program = program;
        this.commandLineArguments = commandLineArguments;
        this.inputArity = inputArity;
        this.outputArity = outputArity;
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
        for (String s : commandLineArguments)
            sb.append(s);

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
}