package dshell.core.worker;

import dshell.core.Operator;
import dshell.core.OperatorFactory;
import dshell.core.OperatorType;
import dshell.core.misc.SystemMessage;
import dshell.core.nodes.Sink;
import dshell.core.nodes.StatelessOperator;
import org.apache.hadoop.hdfs.server.datanode.FileIoProvider;
import org.apache.hadoop.yarn.webapp.RemoteExceptionData;

import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.rmi.Remote;

public class WorkerProcess {
    private static RemoteExecutionData red;

    public static void main(String[] args) {
        red = extractREDFromCommandLineArgs(args);

        try {
            // TODO: fix this once non-stateless operator are introduced
            Operator operator = red.getOperator();
            byte[][] data = null;

            // do not wait for data in case that the operator is the first one to execute
            if (!red.isInitialOperator()) {
                data = new byte[operator.getInputArity()][];

                // socket called 'inputDataSocket' is used to get the input data from another operator
                try (ServerSocket inputDataServerSocket = new ServerSocket(red.getInputPort())) {
                    for (int i = 0; i < operator.getInputArity(); i++) {
                        try (Socket inputDataSocket = inputDataServerSocket.accept();
                             ObjectInputStream ois = new ObjectInputStream(inputDataSocket.getInputStream())) {

                            // here order is not dependent because the operator is stateless
                            data[i] = (byte[]) ois.readObject();
                        }
                    }
                }
            }

            if (operator.getOperatorType() != OperatorType.HDFS_OUTPUT &&
                    operator.getOperatorType() != OperatorType.SOCKETED_OUTPUT) {
                // connecting output socket to current operator
                Operator[] socketedOutput = new Operator[operator.getOutputArity()];
                for (int i = 0; i < operator.getOutputArity(); i++)
                    socketedOutput[i] = OperatorFactory.createSocketedOutput(red.getOutputHost()[i], red.getOutputPort()[i]);
                operator.subscribe(socketedOutput);

                // invoking the operator's computation; after the computation, the data is sent via socket to next node
                // if this is split operator the data splitting will be done inside an operator and the data will be
                // outputted to the sockets that were created few lines before this
                boolean done = false;
                for (int i = 0; i < operator.getInputArity() || (!done && operator.getInputArity() == 0); i++) {
                    operator.next(i, data);
                    done = true;
                }
            } else // instance of sink
            {
                operator.next(0, data);

                // note: this operator is the last operator in the pipeline and therefore it sends signal back to client
                // that the computation has been completed
                try (Socket callbackSocket = new Socket(red.getCallbackHost(), red.getCallBackPort());
                     ObjectOutputStream callbackOOS = new ObjectOutputStream(callbackSocket.getOutputStream())) {

                    callbackOOS.writeObject(new SystemMessage.ComputationFinished());
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private static RemoteExecutionData extractREDFromCommandLineArgs(String[] args) {
        RemoteExecutionData red = new RemoteExecutionData();
        int readFrom = 0;

        int inputArity = Integer.parseInt(args[readFrom++]);
        int outputArity = Integer.parseInt(args[readFrom++]);
        OperatorType operatorType = OperatorType.parseInteger(Integer.parseInt(args[readFrom++]));
        if (operatorType == OperatorType.STATELESS) {
            String program = args[readFrom++];
            int numberOfArgs = Integer.parseInt(args[readFrom++]);
            String[] commandLineArgs = null;

            if (numberOfArgs >= 1) {
                commandLineArgs = new String[numberOfArgs];
                for (int i = 0; i < numberOfArgs; i++)
                    commandLineArgs[i] = args[readFrom++];
            }
            // PARALLELIZATION HINT IS INVALID PARAMETER HERE
            Operator operator = new StatelessOperator(inputArity,
                    outputArity,
                    program,
                    commandLineArgs);
            red.setOperator(operator);
        } else if (operatorType == OperatorType.MERGE)
            red.setOperator(OperatorFactory.createMerger(inputArity));
        else if (operatorType == OperatorType.SPLIT)
            red.setOperator(OperatorFactory.createSplitter(outputArity));
        else if (operatorType == OperatorType.HDFS_OUTPUT)
            red.setOperator(OperatorFactory.createHDFSFilePrinter("output.txt"));
        else
            throw new RuntimeException("Not supported type of operator");

        boolean initialOperator = Boolean.parseBoolean(args[readFrom++]);
        red.setInitialOperator(initialOperator);

        int inputPort = Integer.parseInt(args[readFrom++]);
        red.setInputPort(inputPort);

        int numberOfOutputHosts = Integer.parseInt(args[readFrom++]);
        String[] outputHosts = new String[numberOfOutputHosts];
        int[] outputPorts = new int[numberOfOutputHosts];
        for (int i = 0; i < numberOfOutputHosts; i++) {
            outputHosts[i] = args[readFrom++];
            outputPorts[i] = Integer.parseInt(args[readFrom++]);
        }
        red.setOutputHost(outputHosts);
        red.setOutputPort(outputPorts);

        String callbackHost = args[readFrom++];
        int callbackPort = Integer.parseInt(args[readFrom++]);
        red.setCallbackHost(callbackHost);
        red.setCallBackPort(callbackPort);

        return red;
    }
}