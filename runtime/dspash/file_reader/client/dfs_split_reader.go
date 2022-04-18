package main

import (
	"bufio"
	"context"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"log"
	"os"

	"google.golang.org/grpc"
	"google.golang.org/grpc/credentials/insecure"

	pb "dspash/filereader"
)

var (
	config     = flag.String("config", "", "File to read")
	splitNum   = flag.Int("split", 0, "The logical split number")
	serverPort = flag.Int("port", 50051, "The server port, all machines should use same port")
	prefix     = flag.String("prefix", "", "The top directory")
)

// Distrubted file system block
type DFSBlock struct {
	Path  string
	Hosts []string
}

// Distributed file system config
type DFSConfig struct {
	Blocks []DFSBlock
}

// TODO: improve so that we don't use network if block is replicated on the same machine
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

		stream, err := client.ReadFile(ctx, &pb.FileRequest{Path: *prefix + block.Path})
		if err != nil {
			continue
		}

		for {
			reply, err := stream.Recv()
			if err == io.EOF {
				return ok, err
			}
			if err != nil {
				// Can't recover because we already wrote some bytes.
				// TODO: recover by using intermediate buffer or adding rpcs to allow
				// 		discarding on server side
				return ok, err
			}
			for _, byt := range reply.Buffer {
				err := writer.WriteByte(byt)
				if err != nil {
					return
				}
				if byt == '\n' {
					return true, nil
				}
			}
		}
	}
	return
}

func getAbsPath(s string) string {
	return os.ExpandEnv(*prefix + s)
}

func readLocalFile(p string, skipFirstLine bool, writer *bufio.Writer) error {
	file, err := os.Open(p)
	if err != nil {
		return err
	}
	defer file.Close()

	reader := bufio.NewReader(file)

	if skipFirstLine {
		_, err = reader.ReadString('\n') //discarded
		if err != nil {
			return err
		}
	}

	io.Copy(writer, reader)

	return nil
}

func readUntilDelim(reader *bufio.Reader, writer *bufio.Writer) {
}

func readDFSLogicalSplit(conf DFSConfig, split int) error {

	skipFirstLine := true
	writer := bufio.NewWriter(os.Stdout)
	defer writer.Flush()

	if split == 0 {
		skipFirstLine = false
	}

	err := readLocalFile(getAbsPath(conf.Blocks[split].Path), skipFirstLine, writer)
	if err != nil {
		return err
	}

	// Read until newline
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

func serialize_conf(p string) DFSConfig {
	conf := DFSConfig{}
	byt, err := os.ReadFile(p)
	if err != nil {
		log.Fatalln(err)
	}
	if err := json.Unmarshal(byt, &conf); err != nil {
		log.Fatalln(err)
	}
	return conf
}

func main() {
	flag.Parse()
	if flag.NArg() < 1 && *config == "" {
		flag.Usage()
		os.Exit(0)
	} else if *config == "" {
		*config = flag.Arg(0)
	}

	conf := serialize_conf(*config)
	err := readDFSLogicalSplit(conf, *splitNum)
	if err != nil {
		log.Fatalln(err)
	}
}
