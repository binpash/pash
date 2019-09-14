package dshell.core;

import dshell.core.interfaces.Consumer;
import dshell.core.misc.SystemMessage;
import dshell.core.nodes.AtomicGraph;
import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.Sink;
import dshell.core.nodes.StatelessOperator;
import dshell.core.worker.DistributedTask;
import dshell.core.worker.RemoteExecutionData;
import dshell.core.worker.WorkloadDistributer;
import org.xbill.DNS.Serial;

import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.rmi.Remote;
import java.util.*;

public abstract class Graph {
    public abstract Operator getOperator();

    public void executeLocallySingleThreaded() {
        compileTopology();

        getOperator().next(0, null);
    }

    public void executeRemote(int clientSocket) {
        try {
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
                StatelessOperator current = (StatelessOperator) graphs[i].getOperator();
                StatelessOperator next = (StatelessOperator) graphs[i + 1].getOperator();

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
             ObjectOutputStream oos = new ObjectOutputStream(coordinator.getOutputStream());
             ObjectInputStream ois = new ObjectInputStream(coordinator.getInputStream())) {

            for (RemoteExecutionData rem : remoteExecutionDataList)
                oos.writeObject(rem);
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
        Map<Operator, Integer> previousOperatorPorts = new HashMap<>();

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
                else
                    rem.setInputPort(previousOperatorPorts.get(previous));
                // ------------------------------------------------------------------------------

                if (current.getConsumers() != null) {
                    // output
                    rem.setOutputHost("localhost");
                    rem.setOutputPort(port);
                    previousOperatorPorts.put(current, port++);
                } else {
                    // callback to client that the graph execution has been completed
                    rem.setCallbackHost("localhost");
                    rem.setCallBackPort(3529);
                }
                // ------------------------------------------------------------------------------

                metadata.add(rem);
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

    /*public String toString() {
        StringBuilder sb = new StringBuilder();
        Operator op = getOperator();

        while (op != null && !(op instanceof Sink)) {
            sb.append(op.getProgram());
            if (op.getCommandLineArguments() != null) {
                for (String s : op.getCommandLineArguments()) {
                    sb.append(" " + s);
                }
            }

            if (op.getConsumers() != null && !(op.getConsumer() instanceof Sink))
                sb.append(" | ");

            op = (Operator) op.getConsumer();
        }

        return sb.toString();
    }*/
}