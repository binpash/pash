package dshell.core;

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

    public static Sink createHDFSFilePrinter(String filename) {
        return new Sink() {
            @Override
            public void next(int inputChannel, Object data) {
                DFileSystem.uploadFile(filename, (byte[]) data);
            }
        };
    }

    public static Sink createSocketedOutput(String address, int port) {
        return new Sink() {
            @Override
            public void next(int inputChannel, Object data) {
                try (Socket socket = new Socket(address, port);
                     ObjectOutputStream outputStream = new ObjectOutputStream(socket.getOutputStream());
                     ObjectInputStream inputStream = new ObjectInputStream(socket.getInputStream())) {

                    outputStream.writeObject(data);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        };
    }

    public static Operator<byte[], byte[]> createSplitter(int outputArity) {
        return new Operator<>(1, outputArity, null) {
            @Override
            public void next(int inputChannel, byte[] data) {
                byte[][] split = Utilities.splitData(data, outputArity);

                for (int i = 0; i < this.outputArity; i++)
                    consumers[i].next(0, split[i]);
            }
        };
    }

    public static Operator<byte[], byte[]> createMerger(int inputArity) {
        return new Operator<>(inputArity, 1, null) {
            private int received = 0;

            @Override
            public void next(int inputChannel, byte[] data) {
                byte[][] buffer = new byte[inputArity][];

                if (buffer[inputChannel] == null) {
                    buffer[inputChannel] = data;
                    received++;
                } else
                    throw new RuntimeException("The data for the specified channel has been received already.");

                if (received == inputArity)
                    consumers[0].next(0, Utilities.mergeData(buffer));
            }
        };
    }
}