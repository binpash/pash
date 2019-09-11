# Find markdown files, compile them, and expose them over the network

find . '*.md' | # Parallelizable, given a distributed FS
    xargs mdc | # xargs is higher-order, trivially parallelizable; mdc is third-party
    nc -l 80    # netcat could default-but-configurably parallelizable


