package dshell.core.misc;

import java.io.Serializable;

public class SystemMessage implements Serializable {
    public static class Payload implements Serializable {
        // an empty class
    }

    public static class RemoteException extends Payload {
        String operatorName;
        String message;

        public RemoteException(String operatorName, String message) {
            this.operatorName = operatorName;
            this.message = message;
        }

        public String getOperatorName() {
            return operatorName;
        }

        public String getMessage() {
            return message;
        }
    }

    public static class ComputationFinished extends Payload {

    }

    public static class EndOfREM extends Payload {

    }

    public static class EndOfData extends Payload {

    }
}