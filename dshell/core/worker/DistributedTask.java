package dshell.core.worker;

import dshell.core.misc.SystemMessage;

import java.io.InputStream;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.rmi.Remote;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;

public class DistributedTask {
    private final static int NUMBER_OF_THREADS = 8;

    public static void main(String[] args) throws Exception {
        if (args.length < 1)
            throw new RuntimeException("Invalid program call arguments.");

        int managerPort = Integer.parseInt(args[0]);

        try (ServerSocket serverSocket = new ServerSocket(managerPort)) {
            ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(NUMBER_OF_THREADS);
            System.out.println("Running distributed task executor on port " + managerPort);

            while (true) {
                try (Socket client = serverSocket.accept();
                     ObjectInputStream ois = new ObjectInputStream(client.getInputStream())) {

                    while (true) {
                        Object received = ois.readObject();

                        if (received instanceof RemoteExecutionData)
                            executor.execute(new WorkerThread((RemoteExecutionData) received));
                        else if (received instanceof SystemMessage.EndOfREM)
                            break;
                        else
                            throw new RuntimeException("Not supported type of input data received by worker coordinator.");
                    }
                }
            }
        } catch (Exception e) {
            System.out.println(e.getStackTrace());
        }
    }
}