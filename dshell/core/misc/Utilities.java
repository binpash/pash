package dshell.core.misc;

import jdk.jshell.spi.ExecutionControl;

import java.util.Arrays;

public class Utilities {

    public static String linuxFileRegex = "(\\/[^\\/ ]*)+\\/?";

    public static byte[][] splitData(byte[] input, int numberOfChunks) {
        int chunkSize = (int) Math.ceil(input.length / (double) numberOfChunks);
        byte[][] result = new byte[numberOfChunks][];

        for (int i = 0; i < numberOfChunks; i++)
            result[i] = Arrays.copyOfRange(input, chunkSize * i, Math.min(chunkSize * (i + 1), input.length));

        return result;
    }

    public static byte[] mergeData(byte[][] buffer) {
        int integralSize = 0;
        for (byte[] b : buffer)
            integralSize += b.length;

        byte[] result = new byte[integralSize];
        int writeAt = 0;
        for (int i = 0; i < buffer.length; i++) {
            for (byte b : buffer[i])
                result[writeAt++] = b;
        }

        return result;
    }
}