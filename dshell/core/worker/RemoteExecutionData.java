package dshell.core.worker;

import dshell.core.Operator;

import java.io.Serializable;

/**
 * Class that all the data needed for the computation on a remote machine
 */
public class RemoteExecutionData implements Serializable {
    // representation of computation
    private Operator operator;
    private boolean initialOperator;

    // URL for receiving the data from other operators
    private int inputPort;

    // URL for sending the data to other operators
    private String[] outputHost;
    private int[] outputPort;

    // URL for notifying the client that the computation was finished or that it threw an exception
    private String callbackHost;
    private int callBackPort;

    public RemoteExecutionData() {
    }

    public RemoteExecutionData(Operator operator, boolean initialOperator, int inputPort, String[] outputHost, int[] outputPort, String callbackHost, int callBackPort) {
        this.operator = operator;
        this.initialOperator = initialOperator;
        this.inputPort = inputPort;
        this.outputHost = outputHost;
        this.outputPort = outputPort;
        this.callbackHost = callbackHost;
        this.callBackPort = callBackPort;
    }

    public Operator getOperator() {
        return operator;
    }

    public void setOperator(Operator operator) {
        this.operator = operator;
    }

    public boolean isInitialOperator() {
        return initialOperator;
    }

    public void setInitialOperator(boolean initialOperator) {
        this.initialOperator = initialOperator;
    }

    public int getInputPort() {
        return inputPort;
    }

    public void setInputPort(int inputPort) {
        this.inputPort = inputPort;
    }

    public String[] getOutputHost() {
        return outputHost;
    }

    public void setOutputHost(String[] outputHost) {
        this.outputHost = outputHost;
    }

    public int[] getOutputPort() {
        return outputPort;
    }

    public void setOutputPort(int[] outputPort) {
        this.outputPort = outputPort;
    }

    public String getCallbackHost() {
        return callbackHost;
    }

    public void setCallbackHost(String callbackHost) {
        this.callbackHost = callbackHost;
    }

    public int getCallBackPort() {
        return callBackPort;
    }

    public void setCallBackPort(int callBackPort) {
        this.callBackPort = callBackPort;
    }
}