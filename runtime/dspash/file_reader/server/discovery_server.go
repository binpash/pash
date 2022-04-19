package main

import (
	"context"
	"errors"
	"flag"
	"fmt"
	"log"
	"net"
	"sync"

	"google.golang.org/grpc"

	pb "dspash/datastream"
)

var (
	port = flag.Int("port", 50052, "The server port")
)

type DiscoveryServer struct {
	pb.UnimplementedDiscoveryServer
	addrs map[string]string
	mu    sync.Mutex // protects addrs
}

func (s *DiscoveryServer) PutAddr(ctx context.Context, msg *pb.PutAddrMsg) (*pb.Status, error) {
	s.mu.Lock()
	defer s.mu.Unlock()

	addr, id := msg.Addr, msg.Id
	if _, ok := s.addrs[id]; ok {
		return &pb.Status{Success: false}, errors.New("PutAddr: id already inserted")
	}

	s.addrs[id] = addr
	return &pb.Status{Success: true}, nil
}

func (s *DiscoveryServer) GetAddr(ctx context.Context, msg *pb.AddrReq) (*pb.GetAddrReply, error) {
	s.mu.Lock()
	defer s.mu.Unlock()

	addr, ok := s.addrs[msg.Id]
	if !ok {
		return &pb.GetAddrReply{Success: false}, errors.New("GetAddr: id not found, retry in a little bit")
	}

	return &pb.GetAddrReply{Success: true, Addr: addr}, nil
}

func (s *DiscoveryServer) RemoveAddr(ctx context.Context, msg *pb.AddrReq) (*pb.Status, error) {
	s.mu.Lock()
	defer s.mu.Unlock()

	_, ok := s.addrs[msg.Id]
	if !ok {
		return &pb.Status{Success: false}, errors.New("RemoveAddr: id not found")
	}

	delete(s.addrs, msg.Id)
	return &pb.Status{Success: true}, nil
}

func newServer() *DiscoveryServer {
	s := &DiscoveryServer{}
	s.addrs = map[string]string{}
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
	fmt.Printf("File server running on %v\n", lis.Addr())
	grpcServer.Serve(lis)
}
