# GRPC File Reader Service

## Components
-  Server: Accepts read requests and forwards local files to caller

- Client: The `dfs_split_reader` client takes in a path for a config file (containing file blocks and their hosts) and a split number determining the logical split to read.

### Config File 
The config file looks as follows
```
Config {
    Blocks : Array[Block]
}

Block {
    Path : str
    Hosts : Array[str] (e.g 127.0.0.1)
}
```

Read more [here](https://tammammustafa.notion.site/HDFS-newline-ff2aabde3f9e45c0914760c24f164154)