package dshell.core;

public enum OperatorType {
    MERGE(0),
    SPLIT(1),
    STATELESS(2),
    STATEFULL(3),
    HDFS_OUTPUT(4),
    SOCKETED_OUTPUT(5);

    private int type;

    OperatorType(int type) {
        this.type = type;
    }

    public int getNumericalType() {
        return type;
    }

    public static OperatorType parseInteger(int x) {
        switch (x) {
            case 0:
                return MERGE;
            case 1:
                return SPLIT;
            case 2:
                return STATELESS;
            case 3:
                return STATEFULL;
            case 4:
                return HDFS_OUTPUT;
            case 5:
                return SOCKETED_OUTPUT;
            default:
                throw new RuntimeException("Not supported type.");
        }
    }
}