package dshell.core;

import dshell.core.nodes.Sink;
import dshell.core.worker.RemoteExecutionData;

import java.io.*;
import java.net.Socket;

public class OperatorFactory {

    public static Sink createSystemOutPrinter() {
        return new Sink() {
            @Override
            public void next(Object data) {
                System.out.print(data.toString());
            }
        };
    }

    public static Sink createHDFSFilePrinter(String filename) {
        return new Sink() {
            @Override
            public void next(Object data) {
                DFileSystem.uploadFile(filename, (byte[]) data);
            }
        };
    }

    public static Sink createSocketedOutput(String address, int port) {
        return new Sink() {
            @Override
            public void next(Object data) {
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
}