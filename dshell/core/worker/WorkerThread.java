package dshell.core.worker;

import dshell.core.Operator;
import dshell.core.nodes.Sink;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.Socket;

public class WorkerThread extends Thread {
    private Socket client;
    private String outputHost;
    private int outputPort;

    public WorkerThread(Socket client, String outputHost, int outputPort) {
        this.client = client;
        this.outputHost = outputHost;
        this.outputPort = outputPort;
    }

    @Override
    public void run() {
        try (Socket socket = this.client;
             ObjectOutputStream outputStream = new ObjectOutputStream(socket.getOutputStream());
             ObjectInputStream inputStream = new ObjectInputStream(socket.getInputStream())) {

            Operator operator = (Operator) inputStream.readObject();
            byte[] data = (byte[]) inputStream.readObject();

            operator.subscribe(Sink.createNetworkPrinter(outputHost, outputPort));
            operator.next(data);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}