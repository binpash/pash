package dshell.core;

import dshell.core.misc.SystemMessage;
import dshell.core.nodes.AtomicGraph;
import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.Sink;
import dshell.core.nodes.StatelessOperator;

import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;

public abstract class Graph {
    public abstract Operator getOperator();

    public void executeLocallySingleThreaded() {
        if (this instanceof SerialGraph) {
            SerialGraph g = (SerialGraph) this;
            for (int i = 0; i < g.getAtomicGraphs().length - 1; i++) {
                g.getAtomicGraphs()[i].getOperator().subscribe(g.getAtomicGraphs()[i + 1].getOperator());
            }
        }

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
            // TODO: support other types of graphs rather then just serial ones
            if (this instanceof SerialGraph) {
                SerialGraph serialGraph = (SerialGraph) this;
                AtomicGraph[] atomicGraphs = serialGraph.getAtomicGraphs();
                for (int i = 0; i < atomicGraphs.length - 1; i++) {
                    Operator current = atomicGraphs[i].getOperator();
                    Operator next = atomicGraphs[i + 1].getOperator();

                    if (current instanceof StatelessOperator && next instanceof StatelessOperator)
                    {
                        if (((StatelessOperator)current).getParallelizationHint() != ((StatelessOperator)next).getParallelizationHint() &&
                                ((StatelessOperator)next).getParallelizationHint() != 1) {
                            // e.g. 3 replicas should output to two replicas of the following operator
                            // invalid case
                            throw new RuntimeException("Error defining the topology. Not implemented yet.");
                        }
                    }
                }
            }


            // TODO: ---------------------------------------------------------- EVERYTHING BELOW IS OK
            // wait here until the computation has finished
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
        } catch (Exception ex) {
            throw new RuntimeException(ex.getMessage());
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