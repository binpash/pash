package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"net"
	"os"
	"os/exec"

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

func getAbsPath(s string) (string, error) {
	// Hacky but should work for now
	out, err := exec.Command("bash", "-c", fmt.Sprintf("echo %s", s)).Output()
	if err != nil {
		return "", err
	}
	return string(out), nil
}

func (s *fileReaderServer) ReadFile(req *pb.FileRequest, stream pb.FileReader_ReadFileServer) error {
	filename, err := getAbsPath(req.path)
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
