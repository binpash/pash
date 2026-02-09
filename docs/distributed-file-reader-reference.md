# File Reader Architecture

## Overview

This project implements a **distributed file reading system** inspired by HDFS (Hadoop Distributed File System) concepts, using gRPC for remote procedure calls. The system is designed to read files that are split across multiple machines in a distributed manner, with support for fault tolerance through replica management and intelligent boundary handling for record-based data processing.

The architecture follows a hybrid design pattern combining:
- **Client-Server**: For remote file reading operations
- **Service Discovery**: For dynamic peer discovery and coordination
- **Distributed Systems**: HDFS-inspired block-based file distribution with newline-delimited split handling

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                         Distributed File System                      │
├─────────────────────────────────────────────────────────────────────┤
│                                                                       │
│  ┌──────────────┐      ┌──────────────┐      ┌──────────────┐      │
│  │   Machine 1  │      │   Machine 2  │      │   Machine 3  │      │
│  │              │      │              │      │              │      │
│  │ FileReader   │      │ FileReader   │      │ FileReader   │      │
│  │ Server       │      │ Server       │      │ Server       │      │
│  │ :50051       │      │ :50051       │      │ :50051       │      │
│  │              │      │              │      │              │      │
│  │ [block1.txt] │      │ [block2.txt] │      │ [block3.txt] │      │
│  └──────────────┘      └──────────────┘      └──────────────┘      │
│         ▲                      ▲                      ▲             │
│         │                      │                      │             │
│         │     gRPC Streaming   │                      │             │
│         └──────────────────────┼──────────────────────┘             │
│                                │                                     │
│                      ┌─────────▼────────┐                           │
│                      │  DFS Split       │                           │
│                      │  Reader Client   │                           │
│                      │                  │                           │
│                      │  - Reads split N │                           │
│                      │  - Newline       │                           │
│                      │    boundary      │                           │
│                      │    handling      │                           │
│                      │  - Replica       │                           │
│                      │    failover      │                           │
│                      └──────────────────┘                           │
│                                                                       │
├─────────────────────────────────────────────────────────────────────┤
│                      Service Discovery Layer                         │
├─────────────────────────────────────────────────────────────────────┤
│                                                                       │
│                      ┌──────────────────┐                           │
│                      │  Discovery       │                           │
│                      │  Server :50052   │                           │
│                      │                  │                           │
│                      │  - PutAddr       │                           │
│                      │  - GetAddr       │                           │
│                      │  - RemoveAddr    │                           │
│                      │  - Streams map   │                           │
│                      └────────┬─────────┘                           │
│                               │                                      │
│                    ┌──────────┴──────────┐                          │
│                    │                     │                          │
│              ┌─────▼──────┐       ┌─────▼──────┐                   │
│              │ Datastream │       │ Datastream │                   │
│              │ Writer     │──────▶│ Reader     │                   │
│              │            │  TCP  │            │                   │
│              │ - Register │       │ - Query    │                   │
│              │ - Listen   │       │ - Connect  │                   │
│              └────────────┘       └────────────┘                   │
│                                                                       │
└─────────────────────────────────────────────────────────────────────┘
```

**Key Communication Patterns:**
- **gRPC Streaming**: FileReader servers stream file chunks to clients
- **Direct TCP**: Peer-to-peer data streaming bypasses gRPC for efficiency
- **Service Discovery**: Central coordination for dynamic peer discovery

---

## File Structure

### Protocol Definitions
```
filereader/
├── file_reader.proto           # FileReader service definition
└── file_reader.pb.go          # Generated Go code
└── file_reader_grpc.pb.go     # Generated gRPC code

datastream/
├── data_stream.proto          # Discovery service definition
└── data_stream.pb.go          # Generated Go code
└── data_stream_grpc.pb.go     # Generated gRPC code
```

### Server Implementations
```
server/
├── server.go                  # FileReader server (port 50051)
└── discovery_server.go        # Discovery server (port 50052)
```

### Client Implementations
```
client/
├── client.go                  # Simple file reader client
├── dfs_split_reader.go        # Distributed split reader (primary client)
└── datastream.go              # Peer-to-peer streaming client
```

### Utilities
```
filereader/
└── utils.go                   # Path resolution utilities
```

### Configuration
```
go.mod                         # Go module definition
docs/README.md                 # Project documentation
```

### Future Extensions
```
aws/                           # (Empty - planned AWS integration)
```

---

## Protocol Definitions

### FileReader Service

The FileReader service provides remote file reading with streaming support.

**Protocol Definition** (`filereader/file_reader.proto`):
```protobuf
syntax = "proto3";

option go_package = "dspash/filereader";

service FileReader {
    rpc ReadFile(FileRequest) returns (stream ReadReply) {}
}

message FileRequest {
    string path = 1;
}

message ReadReply {
    bytes buffer = 1;
}
```

**Key Characteristics:**
- **Streaming RPC**: Returns a stream of `ReadReply` messages for efficient large file transfer
- **Simple Interface**: Single RPC method takes a file path and streams chunks
- **Chunk-based**: File is split into configurable chunks (default 4KB)

---

### Discovery Service

The Discovery service enables service discovery and peer-to-peer data streaming coordination.

**Protocol Definition** (`datastream/data_stream.proto`):
```protobuf
syntax = "proto3";

option go_package = "dspash/datastream";

service Discovery {
    rpc PutAddr(PutAddrMsg) returns (Status) {}
    rpc GetAddr(AddrReq) returns(GetAddrReply) {}
    rpc RemoveAddr(AddrReq) returns (Status) {}
    rpc readStream(AddrReq) returns (stream Data) {}
    rpc writeStream(stream Data) returns (Status) {}
}

message PutAddrMsg {
    string Id = 1;
    string Addr = 2;
}

message AddrReq {
    string Id = 1;
}

message GetAddrReply {
    bool Success = 1;
    string Addr = 2;
}

message Status {
    bool Success = 1;
}

message Data {
    bytes buffer = 1;
    string Id = 2; // only sent with first message
}
```

**Key Characteristics:**
- **Service Registry**: `PutAddr`, `GetAddr`, `RemoveAddr` manage node addresses
- **Stream Coordination**: `readStream` and `writeStream` enable data streaming (note: currently marked as buggy and unused in favor of direct TCP)
- **ID-based Discovery**: Nodes register with unique IDs and peers discover each other by ID

---

## Core Components

### Server Components

#### FileReader Server

**Location**: `server/server.go`
**Port**: 50051 (default, configurable)
**Purpose**: Reads files from the local filesystem and streams them to remote clients via gRPC

**Key Implementation** (server/server.go:26-56):
```go
func (s *fileReaderServer) ReadFile(req *pb.FileRequest, stream pb.FileReader_ReadFileServer) error {
    filename, err := pb.GetAbsPath(req.Path)
    if err != nil {
        log.Println(err)
        return err
    }

    file, err := os.Open(filename)
    if err != nil {
        log.Println(err)
        return err
    }
    defer file.Close()

    reader := bufio.NewReader(file)
    buffer := make([]byte, *chunkSize)
    for {
        _, err := reader.Read(buffer)
        stream.Send(&pb.ReadReply{Buffer: buffer})

        if err == io.EOF {
            break
        }

        if err != nil {
            return err
        }
    }

    return nil
}
```

**Features:**
- Path expansion using `GetAbsPath()` utility (supports `~/` expansion)
- Buffered reading with configurable chunk size (default 4KB, flag: `--chunk_size`)
- Streams file chunks until EOF
- Automatic error propagation to client

---

#### Discovery Server

**Location**: `server/discovery_server.go`
**Port**: 50052 (default, configurable)
**Purpose**: Central coordination service for peer discovery and data streaming

**Internal Data Structures** (server/discovery_server.go:25-30):
```go
type DiscoveryServer struct {
    pb.UnimplementedDiscoveryServer
    addrs   map[string]string           // Maps IDs to network addresses
    streams map[string]chan []byte      // Maps IDs to data channels
    mu      sync.Mutex                  // protects addrs and streams
}
```

**Key Implementation - Service Registration** (server/discovery_server.go:32-43):
```go
func (s *DiscoveryServer) PutAddr(ctx context.Context, msg *pb.PutAddrMsg) (*pb.Status, error) {
    s.mu.Lock()
    defer s.mu.Unlock()

    addr, id := msg.Addr, msg.Id
    if _, ok := s.addrs[id]; ok {
        return &pb.Status{Success: false}, errors.New("PutAddr: id already inserted\n")
    }

    s.addrs[id] = addr
    return &pb.Status{Success: true}, nil
}
```

**Key Implementation - Service Discovery** (server/discovery_server.go:45-55):
```go
func (s *DiscoveryServer) GetAddr(ctx context.Context, msg *pb.AddrReq) (*pb.GetAddrReply, error) {
    s.mu.Lock()
    defer s.mu.Unlock()

    addr, ok := s.addrs[msg.Id]
    if !ok {
        return &pb.GetAddrReply{Success: false}, errors.New("GetAddr: id not found, retry in a little bit\n")
    }

    return &pb.GetAddrReply{Success: true, Addr: addr}, nil
}
```

**Features:**
- Thread-safe operations with mutex protection
- ID-based service registry (`addrs` map)
- Stream coordination via channels (`streams` map)
- Timeout mechanism for stream operations (default 10 seconds, flag: `-t`)
- Support for registration, discovery, and cleanup operations

---

### Client Components

#### Simple Client

**Location**: `client/client.go`
**Purpose**: Basic file reader client for reading a single file from a FileReader server

**Functionality:**
- Connects to a specified FileReader server address
- Requests a file by path
- Streams file contents to stdout

**Use Case**: Direct file reading from a known server

---

#### DFS Split Reader (Primary Client)

**Location**: `client/dfs_split_reader.go`
**Purpose**: Sophisticated distributed file system client that reads logical splits from files distributed across multiple machines

**Configuration Format** (JSON):
```json
{
  "Blocks": [
    {"Path": "/data/part1.txt", "Hosts": ["192.168.1.10", "192.168.1.11"]},
    {"Path": "/data/part2.txt", "Hosts": ["192.168.1.12"]},
    {"Path": "/data/part3.txt", "Hosts": ["192.168.1.13", "192.168.1.14"]}
  ]
}
```

**Key Implementation - Boundary Handling** (client/dfs_split_reader.go:38-85):
```go
func readFirstLine(block DFSBlock, writer *bufio.Writer) (ok bool, e error) {
    var opts []grpc.DialOption
    opts = append(opts, grpc.WithTransportCredentials(insecure.NewCredentials()))
    ctx, cancel := context.WithCancel(context.Background())
    defer cancel()

    ok = false
    e = errors.New("Failed to read newline from all replicas")
    for _, host := range block.Hosts {
        addr := fmt.Sprintf("%s:%d", host, *serverPort)
        conn, err := grpc.Dial(addr, opts...)

        if err != nil {
            continue // try next addr
        }
        defer conn.Close()

        client := pb.NewFileReaderClient(conn)
        stream, err := client.ReadFile(ctx, &pb.FileRequest{Path: block.Path})
        if err != nil {
            continue
        }

        for {
            reply, err := stream.Recv()
            if err == io.EOF {
                return ok, err
            }
            if err != nil {
                return ok, err
            }
            for _, byt := range reply.Buffer {
                err := writer.WriteByte(byt)
                if err != nil {
                    return
                }
                if byt == '\n' {
                    return true, nil  // Found newline boundary
                }
            }
        }
    }
    return
}
```

**Key Implementation - Split Reading Logic** (client/dfs_split_reader.go:111-146):
```go
func readDFSLogicalSplit(conf DFSConfig, split int) error {

    skipFirstLine := true
    writer := bufio.NewWriter(os.Stdout)
    defer writer.Flush()

    if split == 0 {
        skipFirstLine = false
    }

    filepath, err := pb.GetAbsPath(conf.Blocks[split].Path)
    if err != nil {
        return err
    }

    err = readLocalFile(filepath, skipFirstLine, writer)
    if err != nil {
        return err
    }

    // Read until newline from subsequent blocks
    for _, block := range conf.Blocks[split+1:] {
        done, err := readFirstLine(block, writer)
        if !done {
            if err == io.EOF {
                continue // read next block if first one didn't contain newline
            } else {
                return err
            }
        } else {
            break
        }
    }
    return nil
}
```

**Features:**
- **Newline-based Split Handling**: Ensures each split contains complete records
  - Skips first line for splits > 0 (handled by previous split)
  - Reads until newline from subsequent blocks to complete the last record
- **Replica Failover**: Tries multiple hosts for each block if primary fails
- **Local File Reading**: Optimized for reading local blocks directly from filesystem
- **Remote Boundary Fetching**: Uses gRPC to fetch continuation lines from remote blocks

**Flags:**
- `--config`: JSON configuration file path
- `--split`: Logical split number to read (0-indexed)
- `--port`: Server port (default 50051)

---

#### Datastream Client

**Location**: `client/datastream.go`
**Purpose**: Enables peer-to-peer data streaming using the Discovery service for coordination

**Architecture:**
- **Writer Mode**: Starts TCP listener, registers address with Discovery server, streams stdin data
- **Reader Mode**: Queries Discovery server for writer's address, connects directly via TCP, streams to stdout

**Key Feature:**
- Direct TCP connection between peers after discovery, bypassing gRPC overhead
- The gRPC stream methods (`readStream`/`writeStream`) in the proto are marked as buggy and unused

**Use Case**: High-throughput peer-to-peer data transfer with dynamic service discovery

---

### Utilities

#### Path Resolution Utility

**Location**: `filereader/utils.go`
**Purpose**: Resolves file paths with shell expansion support

**Implementation** (filereader/utils.go:8-15):
```go
func GetAbsPath(s string) (string, error) {
    // Hacky but should work for now
    out, err := exec.Command("bash", "-c", fmt.Sprintf("echo -n %s", s)).Output()
    if err != nil {
        return "", err
    }
    return string(out), nil
}
```

**Features:**
- Shell path expansion (e.g., `~/file.txt` → `/home/user/file.txt`)
- Environment variable expansion (e.g., `$HOME/file.txt`)
- Uses bash for consistent expansion behavior

**Note**: Comment indicates this is a temporary implementation; future versions may use a more robust approach.

---

## Architecture Patterns & Component Interaction

### Pattern 1: Distributed File Reading

This is the primary use case, enabling MapReduce-style parallel processing of distributed datasets.

```
Configuration File (JSON)                  Multiple Machines Running FileReader Servers
┌──────────────────────┐
│ {                    │                   Machine 1         Machine 2         Machine 3
│   "Blocks": [        │                   ┌──────────┐      ┌──────────┐      ┌──────────┐
│     {                │                   │ Server   │      │ Server   │      │ Server   │
│       "Path": "p1",  │                   │ :50051   │      │ :50051   │      │ :50051   │
│       "Hosts": [...] │                   │          │      │          │      │          │
│     },               │                   │ part1.txt│      │ part2.txt│      │ part3.txt│
│     ...              │                   │ [AAAA\n] │      │ BBBB\n   │      │ CCCC\n]  │
│   ]                  │                   │ [DDDD]   │      │ [EEEE\n] │      │ [FFF]    │
│ }                    │                   └──────────┘      └──────────┘      └──────────┘
└──────────────────────┘                         │                 │                 │
                                                  │                 │                 │
         ┌────────────────────────────────────────┴─────────────────┴─────────────────┘
         │                                        │                 │
         │  DFS Split Reader                      │                 │
         │  --split 0                             │                 │
         │                                        │                 │
         │  Step 1: Read local block 0            │                 │
         │          Output: AAAA\nDDDD            │                 │
         │                                        │                 │
         │  Step 2: Fetch first line from block 1 │                 │
         │          (to complete record)     ────────▶              │
         │          Read until '\n'               │                 │
         │          Output: (nothing, starts with \n)               │
         │                                        │                 │
         │  Final Output: AAAA\nDDDD              │                 │
         └────────────────────────────────────────┘                 │
                                                                     │
                                                                     │
         ┌────────────────────────────────────────────────────────────┘
         │  DFS Split Reader
         │  --split 1
         │
         │  Step 1: Read local block 1, skip first line
         │          (first line "BBBB\n" was part of split 0)
         │          Output: EEEE\n
         │
         │  Step 2: Fetch continuation from block 2
         │          Read until '\n': FFF
         │          Output: FFF (then stop at \n)
         │
         │  Final Output: EEEE\nFFF
         └────────────────────────────────────────
```

**Process Flow:**

1. **Configuration**: JSON file defines distributed blocks with paths and replica hosts
2. **Split Assignment**: Each worker/client is assigned a split number (e.g., `--split 0`)
3. **Local Read**: Client reads the local block file at the split index
   - If split > 0: Skips first line (was completed by previous split)
4. **Boundary Handling**: Client reads from subsequent blocks via gRPC until finding newline
   - Ensures the last record in the split is complete
   - Tries replica hosts if primary fails
5. **Output**: Streams complete records to stdout (ready for processing)

**Benefits:**
- **Record Integrity**: Newline-based splits guarantee complete records
- **Fault Tolerance**: Replica hosts provide failover capability
- **Parallel Processing**: Multiple splits can be processed concurrently
- **Data Locality**: Reads local blocks directly from filesystem (no network overhead)

---

### Pattern 2: Peer-to-Peer Data Streaming

Enables dynamic peer-to-peer communication through central service discovery.

```
┌───────────────────────────────────────────────────────────────────┐
│                      Discovery Server :50052                       │
│                                                                     │
│  State:                                                            │
│    addrs   = {"job123": "192.168.1.50:34567"}                    │
│    streams = {}                                                    │
└────────────────────────┬──────────────────────┬───────────────────┘
                         │                      │
                         │ 1. PutAddr          │ 2. GetAddr
                         │    (register)        │    (discover)
                         │                      │
                ┌────────▼──────────┐  ┌────────▼──────────┐
                │  Datastream       │  │  Datastream       │
                │  Writer           │  │  Reader           │
                │                   │  │                   │
                │  1. Listen on     │  │  1. Query for     │
                │     random port   │  │     "job123"      │
                │     (34567)       │  │                   │
                │                   │  │  2. Get address   │
                │  2. Register ID   │  │     192.168.1.50  │
                │     "job123" with │  │     :34567        │
                │     addr          │  │                   │
                │                   │  │  3. Connect       │
                │  3. Accept conn ◀─┼──┼─────directly      │
                │                   │  │                   │
                │  4. Stream stdin──┼──┼────▶ stdout       │
                │     via TCP       │  │     via TCP       │
                └───────────────────┘  └───────────────────┘
```

**Process Flow:**

**Writer Side:**
1. Start TCP listener on random available port
2. Register with Discovery server: `PutAddr(id="job123", addr="192.168.1.50:34567")`
3. Wait for connection
4. Stream stdin data directly over TCP

**Reader Side:**
1. Query Discovery server: `GetAddr(id="job123")`
2. Receive writer's address: `"192.168.1.50:34567"`
3. Establish direct TCP connection to writer
4. Stream received data to stdout

**Benefits:**
- **Dynamic Discovery**: Peers find each other without hardcoded addresses
- **Direct Connection**: After discovery, data flows peer-to-peer (no intermediary)
- **Efficient Streaming**: Direct TCP bypasses gRPC overhead for maximum throughput
- **Decoupled Lifecycle**: Writer and reader can start in any order

**Note**: The gRPC streaming methods (`readStream`/`writeStream`) defined in the proto are currently unused in favor of this direct TCP approach due to performance and reliability considerations.

---

## Technology Stack

**Language**: Go 1.17

**RPC Framework**:
- gRPC (google.golang.org/grpc)
- Protocol Buffers (protobuf v3)

**Module**: `dspash` (Distributed Spark-like data processing framework)

**Key Dependencies**:
- `google.golang.org/grpc`: gRPC framework
- `google.golang.org/grpc/credentials/insecure`: Insecure credentials for testing/internal networks

**Build & Generation**:
- `protoc` compiler for generating Go code from `.proto` files
- Generated files: `*_pb.go` (message types), `*_grpc.pb.go` (service interfaces)

---

## Key Design Decisions

### 1. Newline-based Split Boundaries

**Rationale**: Ensures record integrity across distributed blocks for line-oriented data (logs, CSV, TSV, etc.)

**Implementation**:
- Splits skip their first line (handled by previous split)
- Splits read continuation from next blocks until finding newline
- Handles edge case where blocks don't end with newline

**Trade-off**:
- ✅ Guarantees complete records per split
- ❌ Requires network calls to subsequent blocks for boundary handling
- ❌ Not suitable for non-line-oriented formats (binary, JSON, etc.)

---

### 2. gRPC Streaming for File Transfer

**Rationale**: Efficient large file transfer with flow control and backpressure

**Implementation**:
- Default 4KB chunks (configurable via `--chunk_size`)
- Server streams chunks until EOF
- Client receives and processes incrementally

**Benefits**:
- Memory-efficient (no need to buffer entire file)
- Network-efficient (flow control prevents overwhelming receiver)
- Standard protocol (interoperable with other gRPC clients)

---

### 3. Replica Support with Failover

**Rationale**: Fault tolerance in distributed environments where machines may be unavailable

**Implementation**:
- Each block can have multiple replica hosts
- Client tries hosts in order until successful
- Graceful fallback for connection/read errors

**Benefits**:
- Increases availability (system works even if some nodes fail)
- No single point of failure for data blocks
- Transparent to caller (automatic retry logic)

---

### 4. Direct TCP for Datastream (Instead of gRPC)

**Rationale**: Bypass gRPC overhead for high-throughput peer-to-peer streaming

**Implementation**:
- Discovery server coordinates via gRPC (lightweight operations)
- Actual data transfer uses raw TCP sockets
- Writer starts TCP listener, reader connects directly

**Trade-off**:
- ✅ Maximum throughput (no protobuf encoding overhead)
- ✅ Lower latency (no gRPC framing)
- ❌ Less portable (no language-agnostic protocol)
- ❌ Manual connection management

**Note**: The proto defines gRPC streaming methods (`readStream`/`writeStream`) but they are marked as buggy and currently unused.

---

### 5. Shell Path Expansion via Bash

**Rationale**: Support user-friendly paths with `~/` and environment variables

**Implementation**: `GetAbsPath()` uses `bash -c "echo -n <path>"` for expansion

**Trade-off**:
- ✅ Simple implementation
- ✅ Consistent with shell behavior
- ❌ Security risk (shell injection if paths are user-controlled)
- ❌ Platform-dependent (requires bash)
- ❌ Performance overhead (spawns process per call)

**Future Improvement**: Code comment indicates this is temporary; production systems should use Go's filepath package with manual home directory expansion.

---

### 6. Configurable Chunk Size

**Rationale**: Tune performance based on network characteristics and file sizes

**Default**: 4KB (typical filesystem page size)

**Tuning Guidance**:
- **Smaller chunks** (1-4KB): Lower latency, higher overhead per chunk
- **Larger chunks** (64KB-1MB): Higher throughput, more memory usage, higher latency

**Use Case Examples**:
- Small files, low latency network: 1-4KB
- Large files, high bandwidth network: 256KB-1MB

---

### 7. Hybrid Architecture Pattern

**Design**: Combines multiple architectural patterns for different use cases

**Patterns Used**:
- **Client-Server**: FileReader service (traditional RPC)
- **Service Discovery**: Discovery server (publish-subscribe registry)
- **Peer-to-Peer**: Datastream direct connections

**Rationale**:
- Different use cases benefit from different patterns
- Centralized discovery with decentralized data transfer
- Flexibility for various deployment scenarios

---

## Summary

This distributed file reading system provides a robust foundation for parallel data processing across multiple machines. Its key strengths are:

1. **HDFS-inspired design** with block-based file distribution
2. **Intelligent boundary handling** ensuring record integrity for line-oriented data
3. **Fault tolerance** through replica management and automatic failover
4. **Efficient streaming** using gRPC and direct TCP based on use case
5. **Service discovery** enabling dynamic peer coordination

The system is designed as part of the larger `dspash` distributed data processing framework, providing the foundational I/O layer for MapReduce-style computation patterns.
