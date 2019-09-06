package dshell.core.nodes;

import dshell.core.DFileSystem;
import dshell.core.Operator;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.net.Socket;

public abstract class Sink extends Operator<Object, Object> {
    public Sink() {
        super(null);
    }

    @Override
    public abstract void next(Object data);
}
