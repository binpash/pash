package main

import (
	"context"
	"errors"
	"flag"
	"fmt"
	"io"
	"log"
	"net"
	"sync"
	"time"

	"google.golang.org/grpc"

	pb "dspash/datastream"
)

var (
	port          = flag.Int("port", 50052, "The server port")
	chunkSize     = flag.Int("chunk_size", 4*1024, "The chunk size for the rpc stream")
	streamTimeout = flag.Int("t", 10, "Wait period in seconds before we give up")
)

type DiscoveryServer struct {
	pb.UnimplementedDiscoveryServer
	addrs   map[string]string
	streams map[string]chan []byte
	mu      sync.Mutex // protects addrs
}

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

func (s *DiscoveryServer) GetAddr(ctx context.Context, msg *pb.AddrReq) (*pb.GetAddrReply, error) {
	s.mu.Lock()
	defer s.mu.Unlock()

	addr, ok := s.addrs[msg.Id]
	if !ok {
		return &pb.GetAddrReply{Success: false}, errors.New("GetAddr: id not found, retry in a little bit\n")
	}

	return &pb.GetAddrReply{Success: true, Addr: addr}, nil
}

func (s *DiscoveryServer) RemoveAddr(ctx context.Context, msg *pb.AddrReq) (*pb.Status, error) {
	s.mu.Lock()
	defer s.mu.Unlock()

	_, ok := s.addrs[msg.Id]
	if !ok {
		return &pb.Status{Success: false}, errors.New("RemoveAddr: id not found\n")
	}

	delete(s.addrs, msg.Id)
	return &pb.Status{Success: true}, nil
}

func (s *DiscoveryServer) ReadStream(req *pb.AddrReq, stream pb.Discovery_ReadStreamServer) error {
	totimer := time.NewTimer(time.Duration(*streamTimeout) * time.Second)
	defer totimer.Stop()
	var ch chan []byte
	for {
		s.mu.Lock()
		val, ok := s.streams[req.Id]
		s.mu.Unlock()
		if ok {
			ch = val
			break
		}
		select {
		case <-time.After(time.Millisecond * 100):
			continue
		case <-totimer.C:
			return errors.New("No writer subscribed in timeout period")
		}
	}

	for buf := range ch {
		stream.Send(&pb.Data{Buffer: buf})
	}

	return nil

}

func (s *DiscoveryServer) WriteStream(stream pb.Discovery_WriteStreamServer) error {
	data, err := stream.Recv() // first message contains id
	if err != nil {
		return err
	}

	ch := make(chan []byte)
	s.mu.Lock()
	s.streams[data.Id] = ch
	s.mu.Unlock()
	defer delete(s.streams, data.Id)

	for {
		data, err := stream.Recv()
		if err == io.EOF {
			close(ch)
			return stream.SendAndClose(&pb.Status{Success: true})
		}
		if err != nil {
			return err
		}
		ch <- data.Buffer
	}
}

func newServer() *DiscoveryServer {
	s := &DiscoveryServer{}
	s.addrs = map[string]string{}
	s.streams = map[string]chan []byte{}
	s.mu = sync.Mutex{}
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
	pb.RegisterDiscoveryServer(grpcServer, newServer())
	fmt.Printf("Discovery server running on %v\n", lis.Addr())
	grpcServer.Serve(lis)
}
