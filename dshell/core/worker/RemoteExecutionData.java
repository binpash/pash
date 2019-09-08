package dshell.core.worker;

import dshell.core.Operator;
import dshell.core.interfaces.Consumer;

import java.io.Serializable;

public class RemoteExecutionData implements Serializable {
    private Operator operator;
    private byte[] inputData;

    private String outputHost;
    private int outputPort;

    private String callbackHost;
    private int callBackPort;

    public RemoteExecutionData(Operator operator, byte[] inputData, String outputHost, int outputPort, String callbackHost, int callBackPort) {
        this.operator = operator;
        this.inputData = inputData;
        this.outputHost = outputHost;
        this.outputPort = outputPort;
        this.callbackHost = callbackHost;
        this.callBackPort = callBackPort;
    }

    public Operator getOperator() {
        return operator;
    }

    public byte[] getInputData() {
        return inputData;
    }

    public String getOutputHost() {
        return outputHost;
    }

    public int getOutputPort() {
        return outputPort;
    }

    public String getCallbackHost() {
        return callbackHost;
    }

    public int getCallBackPort() {
        return callBackPort;
    }
}