package dshell.core.interfaces;

public interface Producer<T> {
    void subscribe(Consumer<T>... consumers);
}