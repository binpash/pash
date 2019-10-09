package dshell.core.worker;

import dshell.core.interfaces.Consumer;
import dshell.core.misc.ReaderWriter;
import dshell.core.misc.SystemMessage;

import java.util.Queue;

public class SocketToProcessThread implements Runnable {
    private ReaderWriter readerWriter;
    private int inputChannel;
    private Consumer nextOperator;
    private Queue<Object> commonQueue;

    public SocketToProcessThread(ReaderWriter readerWriter, int inputChannel, Consumer nextOperator, Queue<Object> commonQueue) {
        this.readerWriter = readerWriter;
        this.inputChannel = inputChannel;
        this.nextOperator = nextOperator;
        this.commonQueue = commonQueue;
    }

    @Override
    public void run() {
        Object data;
        boolean endNotReceived = false;

        while (endNotReceived) {
            readerWriter.startRead();
            data = commonQueue.poll();
            readerWriter.endRead();

            if (data instanceof SystemMessage.EndOfData) {
                endNotReceived = false;
                return;
            }
            else {
                // invoking the operator's computation; after the computation, the data is sent via socket to next node
                // if this is split operator the data splitting will be done inside an operator and the data will be
                // outputted to the sockets that were created few lines before this
                nextOperator.next(inputChannel, data);
            }
        }
    }
}
