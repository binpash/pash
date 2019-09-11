package dshell.core;

import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.Sink;

import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;

public abstract class Graph {
    public abstract Operator getOperator();

    public void executeLocallySingleThreaded() {
        // link the operators

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
            if (this instanceof SerialGraph) {
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
            System.out.println("Waiting for the graph to finish computing...");

            // wait here until the computation has finished
            ServerSocket serverSocket = new ServerSocket(clientSocket);
            try (Socket s = serverSocket.accept();
                 ObjectInputStream ois = new ObjectInputStream(s.getInputStream())) {

                boolean val = ois.readBoolean();
                if (!val)
                    throw new RuntimeException("Error receiving the end of computation marker");
            } finally {
                serverSocket.close();
            }
            System.out.println("The graph has finished computing.");
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