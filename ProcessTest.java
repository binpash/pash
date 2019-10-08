package dshell;

import java.io.InputStream;
import java.io.OutputStream;

public class ProcessTest {
    public static void main(String[] args) throws Exception {
        ProcessBuilder builder = new ProcessBuilder("wc", "-m");
        Process process = builder.start();

        while (true) {
            String sentence = "Lorem ipsum dolor sit amen";

            OutputStream os = process.getOutputStream();
            InputStream is = process.getInputStream();

            os.write(sentence.getBytes());
            os.flush();

            os.close();

            String o = new String(is.readAllBytes());
            System.out.println(o);
            System.out.println(is.available());
        }
    }
}
