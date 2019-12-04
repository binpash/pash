package dshell.core.misc;

public class Tuple<U,V> {
    private U u;
    private V v;

    public Tuple(U u, V v) {
        this.u = u;
        this.v = v;
    }

    public U getU() {
        return u;
    }

    public V getV() {
        return v;
    }

    public void setU(U u) {
        this.u = u;
    }

    public void setV(V v) {
        this.v = v;
    }
}