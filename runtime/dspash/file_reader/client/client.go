package main

import (
	"context"
	"flag"
	"io"
	"log"
	"os"

	"google.golang.org/grpc"
	"google.golang.org/grpc/credentials/insecure"

	pb "dspash/filereader"
)

var (
	filepath   = flag.String("file", "go.sum", "file to read")
	serverAddr = flag.String("addr", "localhost:50051", "The server address in the format of host:port")
)

func readFile(client pb.FileReaderClient, path string) {
	ctx := context.Background()
	stream, err := client.ReadFile(ctx, &pb.FileRequest{Path: path})
	if err != nil {
		log.Fatalf("%v.readFile(_) = _, %v", client, err)
	}

	for {
		reply, err := stream.Recv()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("%v.readFile(_) = _, %v", client, err)
		}

		os.Stdout.Write(reply.Buffer)
	}

}

func main() {
	flag.Parse()
	var opts []grpc.DialOption
	opts = append(opts, grpc.WithTransportCredentials(insecure.NewCredentials()))

	conn, err := grpc.Dial(*serverAddr, opts...)
	if err != nil {
		log.Fatalf("fail to dial: %v", err)
	}
	defer conn.Close()
	client := pb.NewFileReaderClient(conn)
	readFile(client, *filepath)

}
