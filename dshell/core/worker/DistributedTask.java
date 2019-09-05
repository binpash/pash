package dshell.core.worker;

import java.io.InputStream;
import java.net.ServerSocket;
import java.net.Socket;

public class DistributedTask {
    public static void main(String[] args) throws Exception {
        if (args.length < 4)
            throw new RuntimeException("Invalid program call arguments.");

        // input
        String inputHost = args[0];
        int inputPort = Integer.parseInt(args[1]);
        // output
        String outputHost = args[2];
        int outputPort = Integer.parseInt(args[3]);

        try (ServerSocket serverSocket = new ServerSocket(inputPort)) {
            System.out.println("Node has become ready.");

            while (true) {
                Socket client = serverSocket.accept();
                new WorkerThread(client, outputHost, outputPort).start();
                // TODO: add mechanism that shutdowns this node
            }
        } catch (Exception e) {
            System.out.println(e.getStackTrace());
        }

        System.out.println("Node has been shutdown successfully.");
    }
}