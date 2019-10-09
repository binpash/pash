package dshell.core.misc;

import java.util.concurrent.Semaphore;

public class ReaderWriter {
    private Semaphore mutex, common;
    private volatile int readerCount = 0;

    public ReaderWriter() {
        mutex = new Semaphore(1);
        common = new Semaphore(1);
    }

    public void startRead() {
        mutex.acquireUninterruptibly();
        readerCount++;
        if (readerCount == 1)
            common.acquireUninterruptibly();
        mutex.release();
    }

    public void endRead() {
        mutex.acquireUninterruptibly();
        readerCount--;
        if (readerCount == 0)
            common.release();
        mutex.release();
    }

    public void startWrite() {
        common.acquireUninterruptibly();
    }

    public void endWrite() {
        common.release();
    }
}
