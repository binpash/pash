package dshell.core;

import dshell.core.interfaces.Consumer;
import dshell.core.misc.SystemMessage;
import dshell.core.misc.Tuple;
import dshell.core.misc.Utilities;
import dshell.core.nodes.AtomicGraph;
import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.Sink;
import dshell.core.nodes.StatelessOperator;
import dshell.core.worker.DistributedTask;
import dshell.core.worker.RemoteExecutionData;
import dshell.core.worker.WorkloadDistributer;
import org.xbill.DNS.Serial;

import java.io.*;
import java.net.ServerSocket;
import java.net.Socket;
import java.nio.file.Files;
import java.rmi.Remote;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public abstract class Graph {
    public abstract Operator getOperator();

    public void executeLocallySingleThreaded() {
        compileTopology();

        getOperator().next(0, null);
    }

    public void executeRemote(int clientSocket) {
        try {
            // uploading has to be done before the graph expansion
            // TODO: uncomment HDFS
            //uploadInputFiles();
            compileTopology();
            delegateWorkers();

            // block here until the computation has finished distributively
            waitWhileNotFinished(clientSocket);
        } catch (Exception ex) {
            throw new RuntimeException(ex.getMessage());
        }
    }

    private void compileTopology() {
        if (this instanceof SerialGraph) {
            SerialGraph graph = (SerialGraph) this;
            AtomicGraph[] graphs = graph.getAtomicGraphs();

            Operator[] previous = null;

            for (int i = 0; i < graphs.length - 1; i++) {
                Operator current = graphs[i].getOperator();
                Operator next = graphs[i + 1].getOperator();

                if (current.getParallelizationHint() == 1 && next.getParallelizationHint() > 1) {
                    // add splitter here
                    // -------------------------------------------
                    // create replicas of the same operator
                    Operator[] replicas = replicateOperator(next, next.getParallelizationHint());
                    // create splitter
                    Operator splitter = OperatorFactory.createSplitter(next.getParallelizationHint());
                    // current -> splitter -> replicas
                    current.subscribe(splitter);
                    splitter.subscribe(replicas);
                    // keep previous replicas for linkage to next operator
                    previous = replicas;
                } else if (current.getParallelizationHint() == 1 && next.getParallelizationHint() == 1) {
                    // just link the operators
                    // -------------------------------------------
                    current.subscribe(next);
                    previous = null;
                } else if (current.getParallelizationHint() > 1 && next.getParallelizationHint() > 1 && current.getParallelizationHint() == next.getParallelizationHint()) {
                    // link each instance of current to next operator
                    // -------------------------------------------
                    Operator[] replicas = replicateOperator(current, current.getParallelizationHint());
                    for (int j = 0; j < current.getParallelizationHint(); j++)
                        previous[j].subscribe(replicas[j]);
                    previous = replicas;
                } else if (current.getParallelizationHint() > 1 && next.getParallelizationHint() > 1 && current.getParallelizationHint() != next.getParallelizationHint()) {
                    // TODO: solve this by adding merger and the splitter between 'current' and 'next'
                    // maybe another way (??? SEMANTICS ???) would be to link everything with everything
                    // -------------------------------------------
                    throw new RuntimeException("Could not resolve this type of graph.");
                } else if (current.getParallelizationHint() > 1 && next.getParallelizationHint() == 1) {
                    // add merger here
                    // -------------------------------------------
                    Operator merger = OperatorFactory.createMerger(previous.length);
                    for (int j = 0; j < previous.length; j++)
                        previous[j].subscribe(merger);
                    merger.subscribe(next);
                    previous = null;
                } else
                    throw new RuntimeException("This case of graph has not been covered with compiler analysis. The compilation has been aborted.");
            }
        }
    }

    private Operator[] replicateOperator(Operator operator, int numberOfReplicas) {
        Operator[] operators = new Operator[numberOfReplicas];

        for (int i = 0; i < numberOfReplicas; i++) {
            try {
                operators[i] = (Operator) operator.clone();
            } catch (CloneNotSupportedException ex) {
                // this would never happen
            }
        }

        return operators;
    }

    private void uploadInputFiles() throws IOException {
        Matcher matcher = Pattern.compile(Utilities.linuxFileRegex).matcher(this.toString());
        while (matcher.find()) {
            String uri = matcher.group();
            byte[] data = Files.readAllBytes(new File(uri).toPath());

            DFileSystem.uploadFile(uri, data);
        }
    }

    private void delegateWorkers() {
        // TODO: support multiple servers
        //List listOfWorkers = WorkloadDistributer.getAvailableNodes();
        Operator initial = null;

        // figuring out who is the initial operator in the graph
        if (this instanceof AtomicGraph)
            initial = this.getOperator();
        else if (this instanceof SerialGraph)
            initial = ((SerialGraph) this).getAtomicGraphs()[0].getOperator();
        else
            throw new RuntimeException("Type of graph is not supported by the compiler.");

        // creating data that will be sent to coordinator process
        List<RemoteExecutionData> remoteExecutionDataList = createWorkerMetadata(initial);

        // TODO: fix hardcoded values of server
        try (Socket coordinator = new Socket("localhost", 4000);
             ObjectOutputStream oos = new ObjectOutputStream(coordinator.getOutputStream())) {

            for (RemoteExecutionData rem : remoteExecutionDataList)
                oos.writeObject(rem);

            oos.writeObject(new SystemMessage.EndOfREM());
        } catch (Exception ex) {
            throw new RuntimeException(ex.getStackTrace().toString());
        }
    }

    private List<RemoteExecutionData> createWorkerMetadata(Operator root) {
        List<RemoteExecutionData> metadata = new ArrayList<>();
        depthTraversal(root, null, metadata);

        return metadata;
    }

    private static class TraversalPair {
        private Operator operator;
        private Operator previous;

        public TraversalPair(Operator operator, Operator previous) {
            this.operator = operator;
            this.previous = previous;
        }

        public Operator getOperator() {
            return operator;
        }

        public Operator getPrevious() {
            return previous;
        }

    }

    private void depthTraversal(Operator current, Operator previous, List<RemoteExecutionData> metadata) {
        int port = 4001;
        Set<Operator> visited = new HashSet<>();
        Map<Operator, Tuple<List<Integer>, Integer>> previousOperatorPorts = new HashMap<>();
        Map<Operator, RemoteExecutionData> newlyCreatedREMs = new HashMap<>();

        Stack<TraversalPair> stack = new Stack<>();
        stack.push(new TraversalPair(current, previous));

        while (!stack.empty()) {
            TraversalPair pair = stack.pop();
            current = pair.operator;
            previous = pair.previous;

            if (!visited.contains(current)) {
                // graph traversal data
                visited.add(current);
                for (int i = 0; current.getConsumers() != null && i < current.getConsumers().length; i++)
                    stack.push(new TraversalPair((Operator) current.getConsumers()[i], current));

                RemoteExecutionData rem = new RemoteExecutionData();
                rem.setOperator(current);

                // input
                if (previous == null)
                    rem.setInitialOperator(true);
                else {
                    Tuple<List<Integer>, Integer> portPair = previousOperatorPorts.get(previous);
                    int inputPort = portPair.getU().get(portPair.getV());
                    portPair.setV(portPair.getV() + 1);
                    rem.setInputPort(inputPort);
                }
                // ------------------------------------------------------------------------------

                if (current.getConsumers() != null) {
                    // output
                    List<String> outputHosts = new ArrayList<>(current.getOutputArity());
                    List<Integer> outputPorts = new ArrayList<>(current.getOutputArity());
                    for (int i = 0; i < current.getOutputArity(); i++) {
                        // TODO: hardcoded server location
                        if (current.getConsumers() == null || ((Operator) current.getConsumers()[i]).getInputArity() == 1 ||
                                newlyCreatedREMs.containsKey(current.getConsumers()[i]) == false) {
                            outputHosts.add("localhost");
                            outputPorts.add(port++);
                        } else {
                            // TODO: how to get address of the host?
                            outputHosts.add("localhost");
                            outputPorts.add(newlyCreatedREMs.get(current.getConsumers()[i]).getInputPort());
                        }
                    }
                    rem.setOutputHost(outputHosts);
                    rem.setOutputPort(outputPorts);
                    previousOperatorPorts.put(current, new Tuple<List<Integer>, Integer>(outputPorts, 0));
                }

                // callback to client that the graph execution has been completed
                rem.setCallbackHost("localhost");
                rem.setCallBackPort(3529);
                // ------------------------------------------------------------------------------

                metadata.add(rem);
                newlyCreatedREMs.put(current, rem);
            }
        }
    }

    private void waitWhileNotFinished(int clientSocket) throws Exception {
        ServerSocket serverSocket = new ServerSocket(clientSocket);
        try (Socket s = serverSocket.accept();
             ObjectInputStream ois = new ObjectInputStream(s.getInputStream())) {

            Object rx = ois.readObject();
            if (rx instanceof SystemMessage.RemoteException) {
                SystemMessage.RemoteException ex = (SystemMessage.RemoteException) rx;
                throw new RuntimeException("An error happened during the execution of operator '" + ex.getOperatorName() + "'.\n" + ex.getMessage());
            } else if (rx instanceof SystemMessage.ComputationFinished)
                System.out.println("Graph has finished computing successfully. The output was written into the provided file.");
        } finally {
            serverSocket.close();
        }
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        if (this instanceof AtomicGraph) {
            Operator op = getOperator();

            sb.append(op.getProgram());
            if (op.getCommandLineArguments() != null) {
                for (String s : op.getCommandLineArguments()) {
                    sb.append(" " + s);
                }
            }
        } else if (this instanceof SerialGraph) {
            AtomicGraph[] atomicGraphs = ((SerialGraph) this).getAtomicGraphs();

            for (int i = 0; i < atomicGraphs.length; i++) {
                Operator op = atomicGraphs[i].getOperator();

                if (atomicGraphs[i].getOperator() instanceof Sink)
                    continue;

                sb.append(op.getProgram());
                if (op.getCommandLineArguments() != null) {
                    for (String s : op.getCommandLineArguments()) {
                        sb.append(" " + s);
                    }
                }

                if (i < atomicGraphs.length - 1 &&
                        atomicGraphs[i + 1].getOperator() != null &&
                        !(atomicGraphs[i + 1].getOperator() instanceof Sink))
                    sb.append(" | ");
            }
        }

        return sb.toString();
    }
}