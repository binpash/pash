package dshell.core.worker;

import dshell.core.Operator;
import dshell.core.OperatorFactory;
import dshell.core.OperatorType;
import dshell.core.misc.SystemMessage;
import dshell.core.nodes.StatelessOperator;

import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class WorkerProcess implements Runnable {
    // if this is changed to threaded implementation remove keyword 'static'
    private RemoteExecutionData red;

    public WorkerProcess(RemoteExecutionData red) {
        this.red = red;
    }

    // This method will only be called in process mode
    public static void main(String[] args) {
        RemoteExecutionData rem = deserializeArgs(args);
        WorkerProcess workerProcess = new WorkerProcess(rem);

        // Will not be executed as a thread
        workerProcess.run();
    }

    // This method is used both for the execution in threaded and in process mode
    @Override
    public void run() {
        try {
            Operator operator = red.getOperator();

            InternalBuffer[] internalBuffers = new InternalBuffer[operator.getInputArity()];
            Thread[] internalThreads = new Thread[operator.getInputArity()];

            if (operator.getOperatorType() != OperatorType.HDFS_OUTPUT &&
                    operator.getOperatorType() != OperatorType.SOCKETED_OUTPUT) {

                // connecting output socket to current operator
                Operator[] socketedOutput = new Operator[operator.getOutputArity()];
                for (int i = 0; i < operator.getOutputArity(); i++)
                    socketedOutput[i] = OperatorFactory.createSocketedOutput(red.getOutputHost().get(i), red.getOutputPort().get(i));
                operator.subscribe(socketedOutput);
            }

            // do not wait for data in case that the operator is the first one to execute
            if (!red.isInitialOperator()) {
                // all the operators that outputting to this one share the same port
                try (ServerSocket inputDataServerSocket = new ServerSocket(red.getInputPort())) {
                    // receive data from all of them
                    ExecutorService inputThreadPool = Executors.newFixedThreadPool(operator.getInputArity());
                    List<Callable<Object>> inputGateThreads = new ArrayList<>(operator.getInputArity());

                    for (int i = 0; i < operator.getInputArity(); i++) {
                        final int inputChannelParameter = i;

                        internalBuffers[i] = new InternalBuffer();
                        internalThreads[i] = new Thread(new SocketToProcessThread(internalBuffers[i],
                                inputChannelParameter,
                                operator));
                        internalThreads[i].start();

                        inputGateThreads.add(Executors.callable(() -> {// connect by socket with all of them
                            try (Socket inputDataSocket = inputDataServerSocket.accept();
                                 ObjectInputStream ois = new ObjectInputStream(inputDataSocket.getInputStream())) {
                                byte[] data = null;

                                while (true) {
                                    Object received = ois.readObject();

                                    if (received instanceof SystemMessage.EndOfData) {
                                        break;
                                    } else {
                                        data = (byte[]) received;
                                        internalBuffers[inputChannelParameter].write(data);
                                    }
                                }
                            } catch (Exception ex) {
                                ex.printStackTrace();
                            }
                        }));
                    }

                    // wait for the listening threads to complete before proceeding any further
                    inputThreadPool.invokeAll(inputGateThreads);

                    // all the threads have now finished outputting
                    for (int j = 0; j < operator.getInputArity(); j++)
                        operator.next(j, new SystemMessage.EndOfData());
                }
            } else // this will only be called if operator is initial
                operator.next(0, null);

            if (operator.getOperatorType() == OperatorType.HDFS_OUTPUT || operator.getOperatorType() == OperatorType.SOCKETED_OUTPUT) {

                // note: this operator is the last operator in the pipeline and therefore it sends signal back to client
                // that the computation has been completed
                try (Socket callbackSocket = new Socket(red.getCallbackHost(), red.getCallBackPort());
                     ObjectOutputStream callbackOOS = new ObjectOutputStream(callbackSocket.getOutputStream())) {

                    callbackOOS.writeObject(new SystemMessage.ComputationFinished());
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            SystemMessage.RemoteException exception = new SystemMessage.RemoteException(red.getOperator().getProgram(), ex.getMessage());

            try (Socket callbackSocket = new Socket(red.getCallbackHost(), red.getCallBackPort());
                 ObjectOutputStream callbackOOS = new ObjectOutputStream(callbackSocket.getOutputStream())) {

                callbackOOS.writeObject(exception);
            } catch (Exception e) {
                System.err.println("Error communicating with client whose submitted job was aborted.");
            }
        }
    }

    private static RemoteExecutionData deserializeArgs(String[] args) {
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
        List<String> outputHosts = new ArrayList<>((numberOfOutputHosts));
        List<Integer> outputPorts = new ArrayList<>((numberOfOutputHosts));
        for (int i = 0; i < numberOfOutputHosts; i++) {
            outputHosts.add(args[readFrom++]);
            outputPorts.add(Integer.parseInt(args[readFrom++]));
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