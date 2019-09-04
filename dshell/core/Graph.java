package dshell.core;

import dshell.core.nodes.AtomicGraph;
import dshell.core.nodes.SerialGraph;
import dshell.core.nodes.Sink;

import java.io.IOException;
import java.net.ServerSocket;
import java.net.Socket;

public abstract class Graph {
    private final static int SERVER_PORT = 56789;

    public abstract Operator getOperator();

    public void executeLocallySingleThreaded() {
        // link the operators

        if (this instanceof SerialGraph) {
            SerialGraph g = (SerialGraph) this;
            for (int i = 0; i < g.getAtomicGraphs().length - 1; i++) {
                g.getAtomicGraphs()[i].getOperator().subscribe(g.getAtomicGraphs()[i + 1].getOperator());
            }
        }

        getOperator().next(null);
    }

    private int i;

    public void executeLocallyDistributed() {
        try {
            Thread[] threads = null;

            if (this instanceof AtomicGraph) {
                threads = new Thread[1];
            } else if (this instanceof SerialGraph) {
                SerialGraph graph = (SerialGraph) this;
                threads = new Thread[graph.getAtomicGraphs().length];

                for (i = 0; i < threads.length; i++) {
                    graph.getAtomicGraphs()[i].getOperator().subscribe(Sink.createNetworkPrinter("localhost", SERVER_PORT + i + 1));

                    // this simulates the distributed environment
                    threads[i] = new Thread(new Runnable() {
                        private final int port = SERVER_PORT + i;

                        @Override
                        public void run() {
                            try {
                                ServerSocket ss = new ServerSocket(port);

                                if (port - SERVER_PORT == 0)
                                    getOperator().next(null);
                                else {
                                    Socket s = ss.accept();

                                    Object data = new String(s.getInputStream().readAllBytes());
                                    getOperator().next(data);

                                    s.close();
                                }
                            } catch (IOException e) {
                                e.printStackTrace();
                            }
                        }
                    });
                }

                for (i = 0; i < threads.length; i++)
                    threads[i].start();

                // TODO: put barrier here for synchronization and then close all the sockets
                while(true);
            } else
                throw new RuntimeException("Other types of graphs have not been implemented yet.");
        } catch (Exception ex) {
            throw new RuntimeException(ex.getMessage());
        }
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        Operator op = getOperator();

        while (op != null && !(op instanceof Sink)) {
            sb.append(op.getProgram());
            if (op.getCommandLineArguments() != null) {
                for (String s : op.getCommandLineArguments()) {
                    sb.append(" " + s);
                }
            }

            if (op.getConsumer() != null && !(op.getConsumer() instanceof Sink))
                sb.append(" | ");

            op = (Operator) op.getConsumer();
        }

        return sb.toString();
    }
}