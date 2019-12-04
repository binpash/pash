package dshell.core.worker;

import java.util.concurrent.ArrayBlockingQueue;

public class InternalBuffer {
    private final static int BUFFER_CAPACITY = 1024;

    private volatile ArrayBlockingQueue<Object> buffer = new ArrayBlockingQueue<>(BUFFER_CAPACITY, true);
    //private Semaphore mutex, emptySlots, fullSlots;

    /*public InternalBuffer() {
        mutex = new Semaphore(1);
        emptySlots = new Semaphore(BUFFER_CAPACITY, true);
        fullSlots = new Semaphore(0, true);
    }*/

    public void write(Object data) {
        try {
            buffer.put(data);
        } catch (Exception ex) {
            // won't happen
        }
        /*emptySlots.acquireUninterruptibly();
        mutex.acquireUninterruptibly();
        buffer.add(data);
        mutex.release();
        fullSlots.release();*/
    }

    public Object read() {
        Object d = null;
        try {
            d = buffer.take();
        } catch (Exception ex) {
            // won't happen
        }

        /*fullSlots.acquireUninterruptibly();
        mutex.acquireUninterruptibly();
        d = buffer.poll();
        mutex.release();
        emptySlots.release();*/

        return d;
    }
}