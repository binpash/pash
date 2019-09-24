package dshell.core.worker;

import dshell.core.misc.SystemMessage;

import java.io.InputStream;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.rmi.Remote;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;

public class DistributedTask {
    private final static int NUMBER_OF_THREADS = 8;
    private final static String WORKER_NODE = "";

    public static void main(String[] args) throws Exception {
        if (args.length < 1)
            throw new RuntimeException("Invalid program call arguments.");

        int managerPort = Integer.parseInt(args[0]);

        try (ServerSocket serverSocket = new ServerSocket(managerPort)) {
            ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(NUMBER_OF_THREADS);
            System.out.println("Running distributed task executor on port " + managerPort);

            while (true) {
                try (Socket client = serverSocket.accept();
                     ObjectInputStream ois = new ObjectInputStream(client.getInputStream())) {

                    while (true) {
                        Object received = ois.readObject();

                        if (received instanceof RemoteExecutionData) {
                            RemoteExecutionData red = (RemoteExecutionData) received;
                            List<String> argsToSend = new ArrayList<>();

                            argsToSend.add("java");                             // JRE
                            argsToSend.add("-cp");                              // classpath
                            argsToSend.add("./dshell.jar");                     // executable location
                            argsToSend.add("dshell.core.worker.WorkerProcess"); // main method class container

                            argsToSend.add(Integer.toString(red.getOperator().getInputArity()));
                            argsToSend.add(Integer.toString(red.getOperator().getOutputArity()));
                            argsToSend.add(red.getOperator().getProgram());
                            argsToSend.add(Integer.toString(red.getOperator().getCommandLineArguments().length));
                            for (String commandLineArgument : red.getOperator().getCommandLineArguments())
                                argsToSend.add(commandLineArgument);

                            argsToSend.add(Boolean.toString(red.isInitialOperator()));

                            argsToSend.add(Integer.toString(red.getInputPort()));

                            argsToSend.add(Integer.toString(red.getOutputHost().length));
                            for (int i = 0; i < red.getOutputHost().length; i++) {
                                argsToSend.add(red.getOutputHost()[i]);
                                argsToSend.add(Integer.toString(red.getOutputPort()[i]));
                            }

                            argsToSend.add(red.getCallbackHost());
                            argsToSend.add(Integer.toString(red.getCallBackPort()));

                            ProcessBuilder builder = new ProcessBuilder(argsToSend);
                            builder.start();
                            // process will terminate itself upon completion of main method
                        } else if (received instanceof SystemMessage.EndOfREM)
                            break;
                        else
                            throw new RuntimeException("Not supported type of input data received by worker coordinator.");
                    }
                }
            }
        } catch (Exception e) {
            System.out.println(e.getStackTrace());
        }
    }
}