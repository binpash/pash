package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"net"
	"os"

	"google.golang.org/grpc"

	pb "dspash/filereader"
)

var (
	port      = flag.Int("port", 50051, "The server port")
	chunkSize = flag.Int("chunk_size", 4*1024, "The chunk size for the rpc stream")
)

type fileReaderServer struct {
	pb.UnimplementedFileReaderServer
}

func (s *fileReaderServer) ReadFile(req *pb.FileRequest, stream pb.FileReader_ReadFileServer) error {
	file, err := os.Open(os.ExpandEnv(req.Path))
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

func newServer() *fileReaderServer {
	s := &fileReaderServer{}
	return s
}

func main() {
	flag.Parse()
	lis, err := net.Listen("tcp", fmt.Sprintf("0.0.0.0:%d", *port))
	if err != nil {
		log.Fatalf("failed to listen: %v", err)
	}
	var opts []grpc.ServerOption
	grpcServer := grpc.NewServer(opts...)
	pb.RegisterFileReaderServer(grpcServer, newServer())
	fmt.Printf("File server running on %v\n", lis.Addr())
	grpcServer.Serve(lis)
}
