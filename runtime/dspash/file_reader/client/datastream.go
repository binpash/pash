package main

import (
	"context"
	"errors"
	"flag"
	"io"
	"log"
	"net"
	"os"
	"time"

	pb "dspash/datastream"

	"google.golang.org/grpc"
	"google.golang.org/grpc/credentials/insecure"
)

var (
	streamType = flag.String("type", "", "Either read/write")
	serverAddr = flag.String("addr", "localhost:50052", "The server address in the format of host:port")
	streamId   = flag.String("id", "", "The id of the stream")
	debug      = flag.Bool("d", false, "Turn on debugging messages")
)

func getAddr(client pb.DiscoveryClient, timeout time.Duration) (string, error) {
	ctx, cancel := context.WithTimeout(context.Background(), timeout)
	defer cancel()

	totimer := time.NewTimer(timeout)
	defer totimer.Stop()
	for {
		reply, err := client.GetAddr(ctx, &pb.AddrReq{Id: *streamId})
		if err == nil {
			return reply.Addr, nil
		}
		select {
		case <-time.After(time.Millisecond * 100):
			log.Printf("%s retrying to connect\n", err)
			continue
		case <-totimer.C:
			return "", err
		}
	}
}

func removeAddr(client pb.DiscoveryClient) {
	ctx, cancel := context.WithTimeout(context.Background(), time.Second)
	defer cancel()

	client.RemoveAddr(ctx, &pb.AddrReq{Id: *streamId})
}

func read(client pb.DiscoveryClient) (int64, error) {
	timeout := 10 * time.Second
	addr, err := getAddr(client, timeout)
	if err != nil {
		return 0, err
	}

	conn, err := net.Dial("tcp4", addr)
	if err != nil {
		return 0, err
	}
	defer conn.Close()

	n, err := io.Copy(os.Stdout, conn)
	return n, err
}

func write(client pb.DiscoveryClient) (int64, error) {
	ctx, cancel := context.WithTimeout(context.Background(), time.Second)
	defer cancel()

	host, err := getMyIP()
	if err != nil {
		return 0, err
	}

	ln, err := net.Listen("tcp4", host+":0") // connect to random port
	if err != nil {
		return 0, err
	}
	defer ln.Close()

	addr := ln.Addr().String()
	_, err = client.PutAddr(ctx, &pb.PutAddrMsg{Id: *streamId, Addr: addr})
	if err != nil {
		return 0, err
	}
	defer removeAddr(client)

	conn, err := ln.Accept()
	if err != nil {
		return 0, err
	}
	defer conn.Close()

	n, err := io.Copy(conn, os.Stdin)

	return n, err
}

func getMyIP() (string, error) {
	addrs, err := net.InterfaceAddrs()
	if err != nil {
		return "", err
	}

	for _, a := range addrs {
		if ipnet, ok := a.(*net.IPNet); ok && !ipnet.IP.IsLoopback() {
			if ipnet.IP.To4() != nil {
				return ipnet.IP.String(), nil
			}
		}
	}

	return "", errors.New("Couldn't figure local machine ip")
}

func main() {
	flag.Parse()
	arg_idx := 0
	if *streamType == "" {
		*streamType = flag.Arg(arg_idx)
		arg_idx += 1
	}

	if *streamId == "" {
		*streamId = flag.Arg(arg_idx)
		arg_idx += 1
	}

	if !*debug {
		log.SetOutput(io.Discard)
	}

	var opts []grpc.DialOption
	opts = append(opts, grpc.WithTransportCredentials(insecure.NewCredentials()))

	conn, err := grpc.Dial(*serverAddr, opts...)
	if err != nil {
		log.Fatalf("Failed to connect to grpc server: %v", err)
	}
	defer conn.Close()
	client := pb.NewDiscoveryClient(conn)

	var reqerr error
	var n int64
	if *streamType == "read" {
		n, reqerr = read(client)
	} else if *streamType == "write" {
		n, reqerr = write(client)
	} else {
		flag.Usage()
		os.Exit(1)
	}

	if reqerr != nil {
		log.Fatalln(reqerr)
	}

	log.Printf("Success %s %d bytes\n", *streamType, n)
}
