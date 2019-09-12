package dshell.core.misc;

public class SystemMessage {
    public static class Payload {
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
}