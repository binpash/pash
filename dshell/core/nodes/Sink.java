package dshell.core.nodes;

import dshell.core.DFileSystem;
import dshell.core.Operator;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.Socket;

public abstract class Sink extends Operator<String, String> {
    public static final int NETWORK_PRINTER_PORT = 38137;

    private Sink(String program) {
        super(program);
    }

    @Override
    public abstract void next(String data);

    public static Sink createPrinter() {
        return new Sink(null) {
            @Override
            public void next(String data) {
                System.out.print(new String(DFileSystem.downloadFile(data)));
            }
        };
    }

    public static Sink createNetworkPrinter(String address, int port) {
        return new Sink(null) {
            private Socket socket;
            private OutputStream os;
            private DataOutputStream dos;

            @Override
            public void next(String data) {
                try {
                    if (socket == null) {
                        socket = new Socket(address, port);
                        os = socket.getOutputStream();
                        dos = new DataOutputStream(os);
                    }

                    dos.write(data.getBytes());
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        };
    }
}
