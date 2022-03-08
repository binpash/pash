package main

import (
	"io"
	"log"
	"net"
	"os"
	"strings"
	"time"

	"github.com/urfave/cli/v2"
)

var RETRY_INTERVAL time.Duration = 10 * time.Millisecond
var DEFAULT_TIMEOUT int = 0 // In seconds

func Port(ln net.Listener) string {
	x := strings.Split(ln.Addr().String(), ":")
	return x[1]
}

func failOnError(err error, msg string, args ...interface{}) {
	if err != nil {
		log.Printf(msg, args...)
		log.Fatalln(err)
	}
}

func listen(protocol string, address string) (net.Conn, error) {
	ln, err := net.Listen(protocol, address)
	if err != nil {
		return nil, err
	}
	defer ln.Close()

	log.Printf("Connected to %s waiting for a connection to accept\n", ln.Addr())

	conn, err := ln.Accept()
	if err != nil {
		return nil, err
	}

	return conn, nil
}

func retry(method func(string, string) (net.Conn, error), protocol string, address string, timeout time.Duration) net.Conn {
	totimer := time.NewTimer(timeout)
	defer totimer.Stop()

	for {
		conn, err := method(protocol, address)
		if err == nil {
			return conn
		}
		select {
		case <-time.After(RETRY_INTERVAL):
			log.Printf("%s retrying to connect\n", err)
			continue
		case <-totimer.C:
			failOnError(err, "Failed to establish connection after %s ms\n", timeout)
		}
	}
}

func getConn(ctx *cli.Context) net.Conn {
	protocol := "tcp4"
	host := "0.0.0.0"
	port := "0"
	listen_flag := ctx.Bool("listen")
	var timeout time.Duration = time.Duration(ctx.Int("timeout")) * time.Second

	if ctx.IsSet("host") {
		host = ctx.String("host")
	} else if ctx.Args().Get(0) != "" {
		host = ctx.Args().Get(0)
	}

	if ctx.IsSet("port") {
		port = ctx.String("port")
	} else if ctx.Args().Get(1) != "" {
		port = ctx.Args().Get(1)
	}

	address := host + ":" + port

	log.Printf("Will try connecting to %s with %s protocol\n", address, protocol)

	var conn net.Conn
	if listen_flag {
		conn = retry(listen, protocol, address, timeout)
	} else {
		conn = retry(net.Dial, protocol, address, timeout)
	}
	return conn
}

func Read(conn net.Conn) {
	//read
	n, err := io.Copy(os.Stdout, conn)
	failOnError(err, "Failed while reading socket data [%s]\n", conn.RemoteAddr())

	log.Printf("Read %d bytes", n)
}

func Write(conn net.Conn) {
	//write
	n, err := io.Copy(conn, os.Stdin)
	failOnError(err, "Failed while writing socket data [%s]", conn.RemoteAddr())

	log.Printf("Wrote %d bytes", n)
}

func main() {
	app := cli.NewApp()
	app.Name = "dspash_pipe"
	app.Usage = "Connect dspash remote reads and writes (simulate unix pipe on network)"
	app.Version = "0.1"
	app.EnableBashCompletion = true

	// globalFlags := []cli.Flag{}
	// app.Flags = globalFlags

	commandFlags := []cli.Flag{
		&cli.IntFlag{
			Name:    "timeout",
			Aliases: []string{"t"},
			Value:   DEFAULT_TIMEOUT,
			Usage:   "Timeout to wait for messages",
		},
		&cli.BoolFlag{
			Name:    "listen",
			Value:   false,
			Aliases: []string{"l"},
			Usage:   "listen to host:port",
		},
		&cli.StringFlag{
			Name:  "host",
			Value: "0.0.0.0",
			Usage: "Host name",
		},
		&cli.StringFlag{
			Name:    "port",
			Aliases: []string{"p"},
			Value:   "0",
			Usage:   "Port number",
		},
		&cli.BoolFlag{
			Name:    "debug",
			Aliases: []string{"d"},
			Value:   false,
			Usage:   "Turn on debugging, use --log_file to specify file",
		},
		&cli.StringFlag{
			Name:    "log_file",
			Aliases: []string{"lg"},
			Value:   "",
			EnvVars: []string{"LOG_FILE"},
			Usage:   "Turn on debugging, use --log_file to specify file",
		},
	}

	prepare := func(ctx *cli.Context) {
		if !ctx.Bool("debug") {
			log.SetOutput(io.Discard)
		} else {
			if ctx.IsSet("log_file") {
				file, err := os.OpenFile(ctx.String("log_file"), os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0666)
				if err != nil {
					log.Fatal(err)
				}
				log.SetOutput(file)
			}
		}
	}

	app.Commands = []*cli.Command{
		{
			Name:    "read",
			Aliases: []string{"r"},
			Usage:   "Read messages from socket",
			Action: func(c *cli.Context) error {
				prepare(c)
				conn := getConn(c)
				defer conn.Close()
				Read(conn)
				return nil
			},
			Flags: commandFlags,
		},
		{
			Name:    "write",
			Aliases: []string{"w"},
			Usage:   "Write to socket",
			Action: func(c *cli.Context) error {
				prepare(c)
				conn := getConn(c)
				defer conn.Close()
				Write(conn)
				return nil
			},
			Flags: commandFlags,
		},
	}

	err := app.Run(os.Args)
	if err != nil {
		log.Fatal(err)
	}
}
