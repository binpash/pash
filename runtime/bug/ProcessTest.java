//package dshell;
// To run this
// 1. create a FIFO: `mkfifo s1`
// 2. in a different terminal window run `cat s1`; it will block
// 3. in this terminal window run `javac ProcessTest.java && java ProcessTest`
//
// The problem is that the other process did not flush the stream; readAll waits
// for the End-of-Stream / EOF, which never came; this reads characters
// incrementally

import java.io.InputStream;
import java.io.OutputStream;

public class ProcessTest {
    public static void main(String[] args) throws Exception {
        ProcessBuilder builder = new ProcessBuilder("tee", "s1");
        Process process = builder.start();

        for (int i = 0; i<10; i++) {
            String sentence = "Lorem ipsum dolor sit amen";

            OutputStream os = process.getOutputStream();
            InputStream is = process.getInputStream();

            System.out.println("one");
            os.write(sentence.getBytes());
            os.flush();

            System.out.println("two");
            //os.close();

            String o = new String("chars: " + is.read());
            System.out.println(o);
            System.out.println(is.available());
        }
    }
}
