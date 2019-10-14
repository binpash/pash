package dshell.core.interfaces;

public interface Consumer<T> {
    void next(int inputChannel, T data);
}