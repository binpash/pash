package dshell.core;

import dshell.core.misc.Tuple;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.*;
import org.apache.hadoop.io.IOUtils;

import java.util.LinkedList;
import java.util.List;
import java.util.UUID;

public class DFileSystem {
    private static final String HDFS_SERVER = "hdfs://localhost:9000";
    private static final Configuration configuration;
    private static FileSystem fileSystem;

    static {
        configuration = new Configuration();
        configuration.set("fs.defaultFS", HDFS_SERVER);
        configuration.set("fs.hdfs.impl", org.apache.hadoop.hdfs.DistributedFileSystem.class.getName());
        configuration.set("fs.file.impl", org.apache.hadoop.fs.LocalFileSystem.class.getName());

        try {
            fileSystem = FileSystem.get(configuration);
        } catch (Exception ex) {
            throw new RuntimeException("File system failed to initialize.");
        }
    }

    public static String generateFilename() {
        UUID uuid = UUID.randomUUID();
        return uuid.toString();
    }

    public static void uploadFile(String filename, byte[] rawData) {
        FSDataOutputStream dos = null;

        try {
            Path outputFile = new Path(filename);

            if (fileSystem.exists(outputFile))
                fileSystem.delete(outputFile);

            dos = fileSystem.create(outputFile);
            dos.write(rawData);
        } catch (Exception ex) {
            ex.printStackTrace();
        } finally {
            IOUtils.closeStream(dos);
        }
    }

    public static byte[] downloadFile(String filename) {
        FSDataInputStream dis = null;
        byte[] result = null;

        try {
            Path inputFile = new Path(filename);

            if (!fileSystem.exists(inputFile))
                throw new RuntimeException("File with the given filename does not exist within current HDFS context.");

            dis = fileSystem.open(inputFile);
            result = dis.readAllBytes();
        } catch (Exception ex) {
            ex.printStackTrace();
        } finally {
            IOUtils.closeStream(dis);
        }

        return result;
    }

    public static List<Tuple<String, Integer>> getFileLocations(String filename) {
        List<Tuple<String, Integer>> result = new LinkedList<>();
        Configuration configuration = new Configuration();

        try {
            Path file = new Path(filename);

            if (!fileSystem.exists(file))
                throw new RuntimeException("File with the given filename does not exist within current HDFS context.");

            BlockLocation[] locations = fileSystem.getFileBlockLocations(file, 0, fileSystem.getLength(file));
            for (BlockLocation b : locations) {
                String[] names = b.getNames();
                for (String n : names) {
                    String host = n.substring(0, n.indexOf(':') - 1);
                    Integer port = Integer.parseInt(n.substring(n.indexOf(':') + 1));

                    result.add(new Tuple(host, port));
                }
            }
        } catch (
                Exception ex) {
            ex.printStackTrace();
        }

        return result;
    }
}