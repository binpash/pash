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
import java.util.ArrayList;
import java.util.List;

public abstract class Graph {
    public abstract Operator getOperator();

    private void compileTopology() {
        if (this instanceof SerialGraph) {
            SerialGraph graph = (SerialGraph) this;
            AtomicGraph[] graphs = graph.getAtomicGraphs();

            Operator[] previous = null;

            for (int i = 0; i < graphs.length - 1; i++) {
                StatelessOperator current = (StatelessOperator) graphs[i].getOperator();
                StatelessOperator next = (StatelessOperator) graphs[i].getOperator();

                if (current.getParallelizationHint() == 1 && next.getParallelizationHint() > 1) {
                    // add splitter here
                    // -------------------------------------------
                    // create replicas of the same operator
                    Operator[] replicas = replicateOperator(current, next.getParallelizationHint());
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

    public void executeLocallySingleThreaded() {
        compileTopology();
        getOperator().next(0, null);
    }

    public void executeDistributed(int clientSocket) {
        try {
            /*if (this instanceof SerialGraph) {
                SerialGraph g = (SerialGraph) this;
                for (int i = 0; i < g.getAtomicGraphs().length; i++) {
                    Sink sink = null;
                    if (i != g.getAtomicGraphs().length - 1) {
                        sink = OperatorFactory.createSocketedOutput("localhost", 4000);

                        g.getAtomicGraphs()[i].getOperator().subscribe(sink);
                        g.getAtomicGraphs()[i].getOperator().setNextOperator(g.getAtomicGraphs()[i + 1].getOperator());
                    }
                }
            }

            System.out.println("The graph has started computing.");
            getOperator().next(0, null);
            System.out.println("Waiting for the graph to finish computing...");*/
            /*
            // TODO: support other types of graphs rather then just serial ones
            if (this instanceof SerialGraph) {
                SerialGraph serialGraph = (SerialGraph) this;
                AtomicGraph[] atomicGraphs = serialGraph.getAtomicGraphs();
                for (int i = 0; i < atomicGraphs.length - 1; i++) {
                    Operator current = atomicGraphs[i].getOperator();
                    Operator next = atomicGraphs[i + 1].getOperator();

                    if (current instanceof StatelessOperator && next instanceof StatelessOperator) {
                        if (((StatelessOperator) current).getParallelizationHint() != ((StatelessOperator) next).getParallelizationHint() &&
                                ((StatelessOperator) next).getParallelizationHint() != 1) {
                            // e.g. 3 replicas should output to two replicas of the following operator
                            // invalid case
                            throw new RuntimeException("Error defining the topology. Not implemented yet.");
                        }
                    }
                }
            }*/
            compileTopology();
            delegateWorkers();

            // block here until the computation has finished distributively
            waitWhileNotFinished(clientSocket);
        } catch (Exception ex) {
            throw new RuntimeException(ex.getMessage());
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
        List listOfWorkers = WorkloadDistributer.getAvailableNodes();
        Operator initial = null;

        // figuring out who is the initial operator in the graph
        if (this instanceof AtomicGraph)
            initial = this.getOperator();
        else if (this instanceof SerialGraph)
            initial = ((SerialGraph) this).getAtomicGraphs()[0].getOperator();
        else
            throw new RuntimeException("Type of graph is not supported by the compiler.");

        traverseGraph(initial, null);
    }

    private void traverseGraph(Operator current, Operator previous) {
        for (int i = 0; i < current.getConsumers().length; i++) {
            Operator next = (Operator) current.getConsumers()[i];

            RemoteExecutionData rem = new RemoteExecutionData();
            if (previous == null)
                rem.setInitialOperator(true);

            traverseGraph(next, current);
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