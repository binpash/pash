package dshell.core.worker;

import dshell.core.Operator;
import dshell.core.OperatorFactory;
import dshell.core.misc.SystemMessage;
import dshell.core.nodes.Sink;
import dshell.core.nodes.StatelessOperator;

import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.net.ServerSocket;
import java.net.Socket;

public class WorkerThread extends Thread {
    private RemoteExecutionData red;

    public WorkerThread(RemoteExecutionData red) {
        this.red = red;
    }

    @Override
    public void run() {
        try {
            // TODO: fix this once non-stateless operator are introduced
            Operator operator = red.getOperator();

            byte[][] data = null;

            // do not wait for data in case that the operator is the first one to execute
            if (!red.isInitialOperator()) {
                data = new byte[operator.getInputArity()][];

                // socket called 'inputDataSocket' is used to get the input data from another operator
                try (ServerSocket inputDataServerSocket = new ServerSocket(red.getInputPort())) {
                    for (int i = 0; i < operator.getInputArity(); i++) {
                        try (Socket inputDataSocket = inputDataServerSocket.accept();
                             ObjectOutputStream oos = new ObjectOutputStream(inputDataSocket.getOutputStream());
                             ObjectInputStream ois = new ObjectInputStream(inputDataSocket.getInputStream())) {

                            // here order is not dependent because the operator is stateless
                            data[i] = (byte[]) ois.readObject();
                        }
                    }
                }
            }

            if (!(operator instanceof Sink)) {
                // connecting output socket to current operator
                Operator[] socketedOutput = new Operator[operator.getConsumers().length];
                for (int i = 0; i < operator.getConsumers().length; i++)
                    socketedOutput[i] = OperatorFactory.createSocketedOutput(red.getOutputHost()[i], red.getOutputPort()[i]);
                operator.subscribe(socketedOutput);

                // invoking the operator's computation; after the computation, the data is sent via socket to next node
                // if this is split operator the data splitting will be done inside an operator and the data will be
                // outputted to the sockets that were created few lines before this
                boolean done = false;
                for (int i = 0; i < operator.getInputArity() || (!done && operator.getInputArity() == 0); i++) {
                    operator.next(i, data);
                    done = true;
                }
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
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}