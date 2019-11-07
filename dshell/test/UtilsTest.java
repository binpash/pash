package dshell.test;

import dshell.core.misc.Utilities;
import dshell.core.worker.InternalBuffer;
import org.junit.Assert;
import org.junit.Test;

public class UtilsTest {
    @Test
    public void splitTest() {
        byte[] full = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05};
        byte[] partial = {0x00, 0x01, 0x02, 0x03, 0x04};

        byte[][] s1 = Utilities.splitData(full, 3);
        Assert.assertTrue(s1.length == 3);
        Assert.assertTrue(s1[0].length == 2);
        Assert.assertTrue(s1[1].length == 2);
        Assert.assertTrue(s1[2].length == 2);
        Assert.assertTrue(s1[0][0] == 0x00);
        Assert.assertTrue(s1[0][1] == 0x01);
        Assert.assertTrue(s1[1][0] == 0x02);
        Assert.assertTrue(s1[1][1] == 0x03);
        Assert.assertTrue(s1[2][0] == 0x04);
        Assert.assertTrue(s1[2][1] == 0x05);

        byte[][] s2 = Utilities.splitData(partial, 3);
        Assert.assertTrue(s2.length == 3);
        Assert.assertTrue(s2[0].length == 2);
        Assert.assertTrue(s2[1].length == 2);
        Assert.assertTrue(s2[2].length == 1);
        Assert.assertTrue(s2[0][0] == 0x00);
        Assert.assertTrue(s2[0][1] == 0x01);
        Assert.assertTrue(s2[1][0] == 0x02);
        Assert.assertTrue(s2[1][1] == 0x03);
        Assert.assertTrue(s2[2][0] == 0x04);


        byte[] twoPartsA = {0x00, 0x01, 0x02, 0x03, 0x04};
        byte[] twoPartsB = {0x00, 0x01, 0x02, 0x03};

        byte[][] s3 = Utilities.splitData(twoPartsA, 2);
        byte[][] s4 = Utilities.splitData(twoPartsB, 2);
        Assert.assertTrue(s3.length == 2);
        Assert.assertTrue(s3[0].length == 3);
        Assert.assertTrue(s3[1].length == 2);

        Assert.assertTrue(s4.length == 2);
        Assert.assertTrue(s4[0].length == 2);
        Assert.assertTrue(s4[1].length == 2);

        Assert.assertTrue(s3[0][0] == 0x00);
        Assert.assertTrue(s3[0][1] == 0x01);
        Assert.assertTrue(s3[0][2] == 0x02);
        Assert.assertTrue(s3[1][0] == 0x03);
        Assert.assertTrue(s3[1][1] == 0x04);

        Assert.assertTrue(s4[0][0] == 0x00);
        Assert.assertTrue(s4[0][1] == 0x01);
        Assert.assertTrue(s4[1][0] == 0x02);
        Assert.assertTrue(s4[1][1] == 0x03);
    }

    @Test
    public void mergeTest() {
        byte[][] data = {{0x00, 0x01}, {0x02, 0x03}, {0x04, 0x05}};
        byte[] merger = Utilities.mergeData(data);

        Assert.assertTrue(merger[0] == 0x00);
        Assert.assertTrue(merger[1] == 0x01);
        Assert.assertTrue(merger[2] == 0x02);
        Assert.assertTrue(merger[3] == 0x03);
        Assert.assertTrue(merger[4] == 0x04);
        Assert.assertTrue(merger[5] == 0x05);

        data = new byte[][]{{0x00, 0x01}, {0x02, 0x03}, {0x04}};
        merger = Utilities.mergeData(data);

        Assert.assertTrue(merger[0] == 0x00);
        Assert.assertTrue(merger[1] == 0x01);
        Assert.assertTrue(merger[2] == 0x02);
        Assert.assertTrue(merger[3] == 0x03);
        Assert.assertTrue(merger[4] == 0x04);
    }

    @Test
    public void readerWriterTest() throws Exception {
        InternalBuffer rw = new InternalBuffer();

        Thread reader = new Thread(new Runnable() {
            @Override
            public void run() {
                for (int i = 0; i < 1000; i++) {
                    System.out.println(rw.read());
                }
            }
        });
        Thread writer = new Thread(new Runnable() {
            @Override
            public void run() {
                for (int i = 0; i < 1000; i++) {
                    rw.write(i + " WRITER");
                }
            }
        });

        reader.start();
        writer.start();

        reader.join();
        writer.join();
    }
}