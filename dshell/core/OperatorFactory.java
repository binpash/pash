package dshell.core;

import dshell.core.misc.Utilities;
import dshell.core.nodes.Sink;
import org.apache.hadoop.yarn.webapp.hamlet2.Hamlet;

import java.io.*;
import java.net.Socket;

public class OperatorFactory {

    public static Sink createSystemOutPrinter() {
        return new Sink() {
            @Override
            public void next(int inputChannel, Object data) {
                System.out.print(data.toString());
            }
        };
    }

    public static Operator<Object, Object> createHDFSFilePrinter(String filename) {
        return new Operator(OperatorType.HDFS_OUTPUT, 1, 0, null, 1) {
            @Override
            public void next(int inputChannel, Object data) {
                //DFileSystem.uploadFile(filename, ((byte[][]) data)[0]);
            }
        };
    }

    public static Operator<Object, Object> createSocketedOutput(String address, int port) {
        return new Operator(OperatorType.SOCKETED_OUTPUT, 1, 0, null, 1) {
            @Override
            public void next(int inputChannel, Object data) {
                try (Socket socket = new Socket(address, port);
                     ObjectOutputStream outputStream = new ObjectOutputStream(socket.getOutputStream())) {

                    outputStream.writeObject(data);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        };
    }

    public static Operator<Object, Object> createSplitter(int outputArity) {
        return new Operator<>(OperatorType.SPLIT, 1, outputArity, null, 1) {
            @Override
            public void next(int inputChannel, Object data) {
                byte[][] split = Utilities.splitData(((byte[][]) data)[0], outputArity);

                for (int i = 0; i < this.outputArity; i++)
                    consumers[i].next(0, split[i]);
            }
        };
    }

    public static Operator<Object, Object> createMerger(int inputArity) {
        return new Operator<>(OperatorType.MERGE, inputArity, 1, null, 1) {
            private int received = 0;
            private byte[][] buffer = new byte[inputArity][];

            @Override
            public void next(int inputChannel, Object data) {
                if (buffer[inputChannel] == null) {
                    buffer[inputChannel] = ((byte[][]) data)[0];
                    received++;
                } else
                    throw new RuntimeException("The data for the specified channel has been received already.");

                if (received == inputArity) {
                    byte[] input = Utilities.mergeData(buffer);
                    consumers[0].next(0, input);
                }
            }
        };
    }
}