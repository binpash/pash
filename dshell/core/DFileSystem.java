package dshell.core;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IOUtils;

import java.util.UUID;

public class DFileSystem {
    private static final String HDFS_SERVER = "hdfs://";

    public static String generateFilename() {
        UUID uuid = UUID.randomUUID();
        return uuid.toString();
    }

    public static void uploadFile(String filename, byte[] rawData) {
        Configuration configuration = new Configuration();
        configuration.set("fs.default.name", HDFS_SERVER);
        FSDataOutputStream dos = null;

        try {
            FileSystem fileSystem = FileSystem.get(configuration);
            Path outputFile = new Path(filename);

            if (fileSystem.exists(outputFile))
                fileSystem.delete(outputFile);
            fileSystem.createFile(outputFile);

            dos = fileSystem.create(outputFile);
            dos.write(rawData);
        } catch (Exception ex) {
            ex.printStackTrace();
        } finally {
            IOUtils.closeStream(dos);
        }
    }

    public static byte[] downloadFile(String filename) {
        Configuration configuration = new Configuration();
        configuration.set("fs.default.name", HDFS_SERVER);
        FSDataInputStream dis = null;
        byte[] result = null;

        try {
            FileSystem fileSystem = FileSystem.get(configuration);
            Path inputFile = new Path(filename);

            if (!fileSystem.exists(inputFile))
                throw new RuntimeException("File with the given filename does not exists within current HDFS context.");

            dis = fileSystem.open(inputFile);
            result = dis.readAllBytes();
        } catch (Exception ex) {
            ex.printStackTrace();
        } finally {
            IOUtils.closeStream(dis);
        }

        return result;
    }
}