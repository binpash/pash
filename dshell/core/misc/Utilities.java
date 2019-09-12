package dshell.core.misc;

import jdk.jshell.spi.ExecutionControl;

import java.util.Arrays;

public class Utilities {
    public static byte[][] splitData(byte[] input, int numberOfChunks) {
        // TODO: test data splitting
        int chunkSize = input.length / numberOfChunks;
        byte[][] result = new byte[numberOfChunks][];

        for (int i = 0; i < numberOfChunks; i++)
            result[i] = Arrays.copyOfRange(input, chunkSize * i, chunkSize * (i + 1));

        return result;

    }

    public static byte[] mergeData(byte[][] buffer) {
        throw new RuntimeException("Not implemented yet.");
    }
}