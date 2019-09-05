package dshell.core;

import dshell.core.interfaces.Consumer;
import dshell.core.interfaces.Producer;

public abstract class Operator<A, B> implements Consumer<A>, Producer<B> {
    protected Consumer<B> consumer;

    protected String program;
    protected String[] commandLineArguments;

    public Operator(String program) {
        this.program = program;
    }

    public Operator(String program, String[] commandLineArguments) {
        this.program = program;
        this.commandLineArguments = commandLineArguments;
    }

    @Override
    public abstract void next(A data);

    @Override
    public void subscribe(Consumer<B> consumer) {
        if (this.consumer != null)
            throw new RuntimeException("Operator is immutable object and hence cannot be assigned with a consumer again.");

        this.consumer = consumer;
    }

    public String getProgram() {
        return program;
    }

    public String[] getCommandLineArguments() {
        return commandLineArguments;
    }

    public Consumer<B> getConsumer() {
        return consumer;
    }

    public String getArgumentsAsString()
    {
        StringBuilder sb = new StringBuilder();
        for (String s : commandLineArguments)
            sb.append(s);

        return sb.toString();
    }
}