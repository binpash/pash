package dshell.core.interfaces;

public interface Consumer<T> {
    void next(T data);
}