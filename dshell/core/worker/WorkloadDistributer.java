package dshell.core.worker;

import dshell.core.misc.Tuple;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

public class WorkloadDistributer {
    private static List<Tuple<String, Integer>> getAvailableNodes() {
        List<Tuple<String, Integer>> result = new ArrayList<>();

        try (BufferedReader bufferedReader = new BufferedReader(new FileReader("/home/cvetkovic/IdeaProjects/nodes_config.txt"))) {
            String line;

            while ((line = bufferedReader.readLine()) != null) {
                String host = line.substring(0, line.indexOf(':') - 1);
                Integer port = Integer.parseInt(line.substring(line.indexOf(':') + 1));

                result.add(new Tuple(host, port));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return result;
    }
}