package dshell.core.worker;

import dshell.core.misc.SystemMessage;
import dshell.core.nodes.StatelessOperator;
import org.apache.hadoop.hdfs.security.token.delegation.DelegationTokenSecretManager;
import org.jboss.netty.channel.socket.Worker;

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
    private static String JAVA_PATH; //= "/usr/lib/jvm/java-11-openjdk-amd64/bin/java";
    private final static String WORKER_NODE = "/home/cvetkovic/Desktop/dshell_executable/dshell.jar";

    static {
        try {
            JAVA_PATH = System.getenv("JAVA_HOME") + "/java";
        } catch (Exception ex) {
            throw new RuntimeException("Could not find environmental variable 'JAVA_HOME' defined on this machine.");
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 1)
            throw new RuntimeException("Invalid program call arguments.");

        int managerPort = Integer.parseInt(args[0]);

        try (ServerSocket serverSocket = new ServerSocket(managerPort)) {
            System.out.println("Running distributed task executor on port " + managerPort);

            while (true) {
                try (Socket client = serverSocket.accept();
                     ObjectInputStream ois = new ObjectInputStream(client.getInputStream())) {

                    // TODO: do the graph traversal from the end so that initial opeators is sent last
                    // if it is not last it can make problems during the execution because some ports won't be opened for listening
                    ProcessBuilder pipelineInitiator = null;
                    RemoteExecutionData initial = null;

                    while (true) {
                        Object received = ois.readObject();

                        if (received instanceof RemoteExecutionData) {
                            RemoteExecutionData red = (RemoteExecutionData) received;
                            /*if (!red.isInitialOperator())
                                new Thread(new WorkerProcess(red)).start();
                            else
                                initial = red;*/

                            List<String> argsToSend = new ArrayList<>();

                            argsToSend.add(JAVA_PATH);
                            argsToSend.add("-cp");
                            argsToSend.add(WORKER_NODE);
                            argsToSend.add("dshell.core.worker.WorkerProcess");

                            argsToSend.add(Integer.toString(red.getOperator().getInputArity()));
                            argsToSend.add(Integer.toString(red.getOperator().getOutputArity()));
                            argsToSend.add(Integer.toString(red.getOperator().getOperatorType().getNumericalType()));
                            if (red.getOperator().getProgram() != null) {
                                argsToSend.add(red.getOperator().getProgram());
                                if (red.getOperator().getCommandLineArguments() != null) {
                                    argsToSend.add(Integer.toString(red.getOperator().getCommandLineArguments().length));
                                    for (String commandLineArgument : red.getOperator().getCommandLineArguments())
                                        argsToSend.add(commandLineArgument);
                                } else
                                    argsToSend.add(Integer.toString(0));
                            }

                            argsToSend.add(Boolean.toString(red.isInitialOperator()));
                            argsToSend.add(Integer.toString(red.getInputPort()));

                            if (red.getOutputHost() != null) {
                                argsToSend.add(Integer.toString(red.getOutputHost().size()));
                                for (int i = 0; i < red.getOutputHost().size(); i++) {
                                    argsToSend.add(red.getOutputHost().get(i));
                                    argsToSend.add(Integer.toString(red.getOutputPort().get(i)));
                                }
                            } else
                                argsToSend.add(Integer.toString(0));

                            argsToSend.add(red.getCallbackHost());
                            argsToSend.add(Integer.toString(red.getCallBackPort()));

                            ProcessBuilder builder = new ProcessBuilder(argsToSend.toArray(new String[0]));
                            builder.redirectError(ProcessBuilder.Redirect.INHERIT);
                            builder.redirectOutput(ProcessBuilder.Redirect.INHERIT);
                            if (!red.isInitialOperator())
                                builder.start();
                            else
                                pipelineInitiator = builder;

                            // NOTE: process will terminate itself upon completion of main method
                        } else if (received instanceof SystemMessage.EndOfREM)
                            break;
                        else
                            throw new RuntimeException("Not supported type of input data received by worker coordinator.");
                    }

                    //new Thread(new WorkerProcess(initial)).start();
                    pipelineInitiator.start();
                }
            }
        } catch (Exception e) {
            System.out.println(e.getStackTrace());
        }
    }
}