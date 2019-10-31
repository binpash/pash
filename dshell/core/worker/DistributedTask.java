package dshell.core.worker;

import dshell.core.misc.SystemMessage;

import java.io.ObjectInputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.*;
import java.util.concurrent.CountDownLatch;

public class DistributedTask {
    private static String JAVA_PATH = "/usr/lib/jvm/java-11-openjdk-amd64/bin/java";
    private final static String WORKER_NODE = "/home/cvetkovic/Desktop/dshell_executable/dshell.jar";

    private static int workerID = 0;
    private static int addPortWith = 0;
    private static Map<Integer, List<Object>> runningWorkers = new HashMap<>();

    static {
        try {
            // TODO: fix JAVA_HOME
            //JAVA_PATH = "/usr/bin/java";                        // System.getenv("JAVA_HOME") + "/java";
            //WORKER_NODE = "/home/cvetkovic/dshell/dshell.jar";  //


            /*File jarDir = new File(ClassLoader.getSystemClassLoader().getResource(".").getPath());
            String absolutePath = jarDir.getAbsolutePath();
            WORKER_NODE = absolutePath.substring(0, absolutePath.indexOf("/target/classes")) + ".jar";*/
        } catch (Exception ex) {
            throw new RuntimeException("Could not find environmental variable 'JAVA_HOME' defined on this machine.");
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 2)
            throw new RuntimeException("Invalid program call arguments.");

        CountDownLatch workerServerSocketBarrier = null;
        int managerPort = Integer.parseInt(args[0]);
        boolean threadedMode;
        if (args[1].equals("-t"))
            threadedMode = true;
        else if (args[1].equals("-p"))
            threadedMode = false;
        else
            throw new RuntimeException("Unsupported running method specified.");

        try (ServerSocket serverSocket = new ServerSocket(managerPort)) {
            System.out.println("Running distributed task executor on port " + managerPort);

            List<RemoteExecutionData> threadsToRun = new LinkedList<>();

            while (true) {
                try (Socket client = serverSocket.accept();
                     ObjectInputStream ois = new ObjectInputStream(client.getInputStream())) {

                    // generate new topology ID and List for exception handling
                    List<Object> exceptionHandlingList = new LinkedList<>();
                    workerID++;

                    // TODO: do the graph traversal from the end so that initial opeators is sent last
                    // if it is not last it can make problems during the execution because some ports won't be opened for listening
                    ProcessBuilder pipelineInitiator = null;
                    RemoteExecutionData initial = null;

                    while (true) {
                        Object received = ois.readObject();

                        if (received instanceof RemoteExecutionData) {
                            RemoteExecutionData red = (RemoteExecutionData) received;

                            if (!threadedMode) {
                                List<String> argsToSend = compileArgsForSending(red);

                                ProcessBuilder builder = new ProcessBuilder(argsToSend.toArray(new String[0]));
                                builder.redirectError(ProcessBuilder.Redirect.INHERIT);
                                builder.redirectOutput(ProcessBuilder.Redirect.INHERIT);
                                if (!red.isInitialOperator()) {
                                    exceptionHandlingList.add(builder);
                                    builder.start();
                                } else
                                    pipelineInitiator = builder;
                            } else {        // threaded mode
                                red.setInputPort(addPortWith + red.getInputPort());
                                if (red.getOutputPort() != null)
                                    for (int i = 0; i < red.getOutputPort().size(); i++)
                                        red.getOutputPort().set(i, addPortWith + red.getOutputPort().get(i));

                                if (!red.isInitialOperator())
                                    threadsToRun.add(red);
                                else
                                    initial = red;
                            }

                            // NOTE: process will terminate itself upon completion of main method
                        } else if (received instanceof SystemMessage.EndOfREM)
                            break;
                        else
                            throw new RuntimeException("Not supported type of input data received by worker coordinator.");
                    }

                    if (threadedMode) {
                        Thread newThread = new Thread(new WorkerProcess(initial, workerServerSocketBarrier, workerID));
                        exceptionHandlingList.add(newThread);
                        newThread.setUncaughtExceptionHandler(threadedModeExceptionHandler);

                        workerServerSocketBarrier = new CountDownLatch(threadsToRun.size());
                        for (RemoteExecutionData red : threadsToRun) {
                            Thread t = new Thread(new WorkerProcess(red, workerServerSocketBarrier, workerID));
                            exceptionHandlingList.add(t);
                            //t.setUncaughtExceptionHandler(threadedModeExceptionHandler);

                            t.start();
                        }

                        workerServerSocketBarrier.await();

                        //Thread.sleep(100);
                        newThread.start();
                    } else {
                        exceptionHandlingList.add(pipelineInitiator);
                        pipelineInitiator.start();
                    }

                    runningWorkers.put(workerID, exceptionHandlingList);
                    addPortWith += exceptionHandlingList.size();

                    threadsToRun.clear();
                }
            }
        } catch (Exception e) {
            System.out.println(e.getStackTrace());
        }
    }

    // TODO: do exception handling in process mode

    private static Thread.UncaughtExceptionHandler threadedModeExceptionHandler = new Thread.UncaughtExceptionHandler() {
        @Override
        public void uncaughtException(Thread thread, Throwable throwable) {
            if (throwable instanceof WorkerException) {
                int topologyToShutdown = ((WorkerException) throwable).getTopologyID();

                List<Object> workers = runningWorkers.get(topologyToShutdown);
                for (Object o : workers) {
                    Thread t = (Thread) o;

                    // deprecated -> safe in this context -> ???
                    t.stop();
                }
            } else {
                System.out.println("Error encountered in one of the workers. The whole node must be shutdown.");
                System.exit(1);
            }
        }
    };

    private static List<String> compileArgsForSending(RemoteExecutionData red) {
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
        argsToSend.add(Integer.toString(red.getInputPort() + addPortWith));

        if (red.getOutputHost() != null) {
            argsToSend.add(Integer.toString(red.getOutputHost().size()));
            for (int i = 0; i < red.getOutputHost().size(); i++) {
                argsToSend.add(red.getOutputHost().get(i));
                argsToSend.add(Integer.toString(red.getOutputPort().get(i) + addPortWith));
            }
        } else
            argsToSend.add(Integer.toString(0));

        argsToSend.add(red.getCallbackHost());
        argsToSend.add(Integer.toString(red.getCallBackPort()));
        argsToSend.add(Integer.toString(workerID));

        return argsToSend;
    }
}