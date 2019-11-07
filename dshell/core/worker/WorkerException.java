package dshell.core.worker;

public class WorkerException extends RuntimeException {
    private int topologyID;

    public WorkerException(String message, int topologyID) {
        super(message);
        this.topologyID = topologyID;
    }

    public int getTopologyID() {
        return topologyID;
    }
}
