package dshell.core.worker;

import java.io.InputStream;
import java.net.ServerSocket;
import java.net.Socket;
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
                Socket client = serverSocket.accept();
                executor.execute(new WorkerThread(client));
            }
        } catch (Exception e) {
            System.out.println(e.getStackTrace());
        }

        System.out.println("Node has been shutdown successfully.");
    }
}