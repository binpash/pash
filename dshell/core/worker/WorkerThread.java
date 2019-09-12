package dshell.core.worker;

import dshell.core.Operator;
import dshell.core.OperatorFactory;
import dshell.core.misc.SystemMessage;
import dshell.core.nodes.Sink;

import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;

public class WorkerThread extends Thread {
    private Socket client;

    public WorkerThread(Socket client) {
        this.client = client;
    }

    @Override
    public void run() {
        // socket called 'socket' is used for communication with client from whom the worker will get all the information
        // about the operation it should execute
        try (Socket socket = this.client;
             ObjectOutputStream outputStream = new ObjectOutputStream(socket.getOutputStream());
             ObjectInputStream inputStream = new ObjectInputStream(socket.getInputStream())) {

            // operation to execute with all the additional info
            RemoteExecutionData red = (RemoteExecutionData) inputStream.readObject();
            byte[] data;

            // socket called 'inputDataSocket' is used to get the input data from another operator
            try (ServerSocket inputDataServerSocket = new ServerSocket(red.getInputPort())) {
                try (Socket inputDataSocket = inputDataServerSocket.accept();
                     ObjectOutputStream oos = new ObjectOutputStream(inputDataSocket.getOutputStream());
                     ObjectInputStream ois = new ObjectInputStream(inputDataSocket.getInputStream())) {
                    data = (byte[]) ois.readObject();
                }
            }

            // operation to execute
            Operator operator = red.getOperator();

            if (!(operator instanceof Sink)) {
                // connecting output socket to current operator
                operator.subscribe(OperatorFactory.createSocketedOutput(red.getOutputHost(), red.getOutputPort()));
                // invoking the operator's computation; after the computation, the data is sent via socket to next node
                operator.next(0, data);
            } else // instance of sink
            {
                operator.next(0, data);

                // note: this operator is the last operator in the pipeline and therefore it sends signal back to client
                // that the computation has been completed
                try (Socket callbackSocket = new Socket(red.getCallbackHost(), red.getCallBackPort());
                     ObjectOutputStream callbackOOS = new ObjectOutputStream(callbackSocket.getOutputStream())) {

                    callbackOOS.writeObject(new SystemMessage.ComputationFinished());
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}