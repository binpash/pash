package dshell.core.nodes;

import dshell.core.Operator;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.Socket;

public abstract class Sink extends Operator<byte[], byte[]> {
    public static final int NETWORK_PRINTER_PORT = 38137;

    private Sink(String program) {
        super(program);
    }

    @Override
    public abstract void next(byte[] data);

    public static Sink createPrinter() {
        return new Sink(null) {
            @Override
            public void next(byte[] data) {
                System.out.print(new String(data));
            }
        };
    }

    public static Sink createNetworkPrinter() {
        return new Sink(null) {
            private Socket socket;
            private OutputStream os;
            private DataOutputStream dos;

            @Override
            public void next(byte[] data) {
                try {
                    if (socket == null) {
                        socket = new Socket("localhost", NETWORK_PRINTER_PORT);
                        os = socket.getOutputStream();
                        dos = new DataOutputStream(os);
                    }

                    dos.write(data);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        };
    }
}
