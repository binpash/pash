package dshell.core.worker;

import java.util.ArrayDeque;
import java.util.LinkedList;
import java.util.Queue;
import java.util.concurrent.Semaphore;

public class InternalBuffer {
    private final static int BUFFER_CAPACITY = 100;

    private ArrayDeque<Object> buffer = new ArrayDeque<>(BUFFER_CAPACITY);
    private Semaphore mutex, emptySlots, fullSlots;

    public InternalBuffer() {
        mutex = new Semaphore(1);
        emptySlots = new Semaphore(BUFFER_CAPACITY);
        fullSlots = new Semaphore(0);
    }

    public void write(Object data) {
        emptySlots.acquireUninterruptibly();
        mutex.acquireUninterruptibly();
        buffer.add(data);
        mutex.release();
        fullSlots.release();
    }

    public Object read() {
        Object d;

        fullSlots.acquireUninterruptibly();
        mutex.acquireUninterruptibly();
        d = buffer.poll();
        mutex.release();
        emptySlots.release();

        return d;
    }
}