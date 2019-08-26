import java.io.*;
import java.util.*;

public class CommandAnalyzer
{
    private static Map<String, Integer> commands = new LinkedHashMap();
    private static int numberOfEncounteredCommands = 0;

    private static class KV
    {
        public String k;
        public Double v;

        public KV(String k, Double v)
        {
            this.k = k;
            this.v = v;
        }
    }

    public static void main(String[] args) throws Exception
    {
        if (args.length < 3)
        {
            System.out.println("Call the program with two parameters. First should be the file that " +
                    "contains linux command list. The second should be directory with shell scripts, " +
                    " and the third should be output.");

            return;
        }

        String linuxCommandList = args[0];
        String shellScriptsDirectory = args[1];
        String outputFile = args[2];

        populateCommandMap(linuxCommandList);

        List<String> listOfFiles = new LinkedList<>();
        getAllFiles(shellScriptsDirectory, listOfFiles);

        processFiles(listOfFiles);

        deliverOutput(outputFile);
    }

    private static void populateCommandMap(String fileLocation) throws Exception
    {
        try (BufferedReader reader = new BufferedReader(new FileReader(fileLocation)))
        {
            String line;

            while ((line = reader.readLine()) != null)
                commands.put(line.trim().toLowerCase(), 0);
        }
    }

    private static void getAllFiles(String root, List<String> result)
    {
        File directory = new File(root);

        File[] files = directory.listFiles();
        if (files != null)
        {
            for (File f : files)
            {
                if (f.isFile())
                    result.add(f.getAbsolutePath());
                else
                    getAllFiles(f.getAbsolutePath(), result);
            }
        }
    }

    private static void processFiles(List<String> files) throws Exception
    {
        for (int i = 0; i < files.size(); i++)
        {
            System.out.println("Processing - " + (i + 1) + "/" + files.size());

            try (BufferedReader reader = new BufferedReader(new FileReader(files.get(i))))
            {
                String line;

                while ((line = reader.readLine()) != null)
                {
                    if (line.startsWith("#"))
                        continue;

                    String[] words = line.split(" ");
                    for (String s : words)
                    {
                        Integer count = commands.get(s);
                        if (count != null)
                        {
                            commands.replace(s, count + 1);
                            numberOfEncounteredCommands++;
                        }
                    }
                }
            }
        }
    }

    private static void deliverOutput(String outputFile) throws Exception
    {
        List<KV> statistics = new LinkedList<>();

        for (Map.Entry<String, Integer> entry : commands.entrySet())
        {
            double counted = entry.getValue();
            double p = counted / numberOfEncounteredCommands;

            statistics.add(new KV(entry.getKey(), p));
        }

        statistics.sort(Comparator.comparing(kv -> kv.v, Comparator.reverseOrder()));

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile)))
        {
            for (KV kv : statistics)
                writer.write(String.format("%-15s - %s\n", kv.k, kv.v));
        }

        System.out.println("The program has successfully finished the computation.");
    }
}