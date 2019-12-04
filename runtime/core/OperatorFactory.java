package dshell.core;

import dshell.core.misc.SystemMessage;
import dshell.core.misc.Utilities;
import dshell.core.nodes.Sink;

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

    public static Operator createHDFSFilePrinter(String filename) {
        return new Operator(OperatorType.HDFS_OUTPUT, 1, 0, null, 1) {
            @Override
            public void next(int inputChannel, Object data) {
                //DFileSystem.uploadFile(filename, ((byte[][]) data)[0]);

                if (data instanceof byte[]) {
                    byte[] d = (byte[]) data;
                    for (int i = 0; i < d.length; i++)
                        if ((char) d[i] != 0)
                            System.out.print((char) d[i]);
                } else    // discard all the other data (SystemMessages)
                    return;
            }
        };
    }

    public static Operator createSocketedOutput(String address, int port) {
        return new Operator(OperatorType.SOCKETED_OUTPUT, 1, 0, null, 1) {
            private Socket socket;
            private ObjectOutputStream outputStream;

            @Override
            public void next(int inputChannel, Object data) {
                try {
                    if (socket == null) {
                        socket = new Socket(address, port);
                        outputStream = new ObjectOutputStream(socket.getOutputStream());
                    }

                    // here writeUnshared() is used because of the usage of the same buffer reference
                    // in other cause it wouldn't work on other side of the socket
                    outputStream.writeUnshared(data);
                    outputStream.flush();

                    // close socket on signal
                    if (data instanceof SystemMessage.EndOfData)
                        closeOutput();
                } catch (Exception e) {
                    e.printStackTrace();

                    try {
                        closeOutput();
                    } catch (Exception ex) {
                        e.printStackTrace();
                    }
                }
            }

            public void closeOutput() throws Exception {
                /*if (socket != null) {
                    outputStream.close();
                    socket.close();
                }*/
            }
        };
    }

    public static Operator<Object, Object> createSplitter(int outputArity) {
        return new Operator<>(OperatorType.SPLIT, 1, outputArity, null, 1) {
            @Override
            public void next(int inputChannel, Object data) {
                if (data instanceof SystemMessage.Payload)
                    for (int i = 0; i < this.outputArity; i++)
                        consumers[i].next(0, data);
                else {
                    byte[][] split = Utilities.splitData((byte[]) data, outputArity);

                    for (int i = 0; i < this.outputArity; i++)
                        consumers[i].next(0, split[i]);
                }
            }
        };
    }

    public static Operator<Object, Object> createMerger(int inputArity) {
        return new Operator<>(OperatorType.MERGE, inputArity, 1, null, 1) {
            private int receivedEOD = 0;

            // NOTE: this method has to be declared as synchronized
            @Override
            public synchronized void next(int inputChannel, Object data) {
                if (data instanceof SystemMessage.EndOfData) {
                    receivedEOD++;
                    if (receivedEOD == inputArity)
                        consumers[0].next(0, new SystemMessage.EndOfData());
                } else
                    consumers[0].next(0, data);
            }
        };
    }
}