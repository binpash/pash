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
		log.Fatalf(msg, args...)
	}
}

func getConn(ctx *cli.Context) net.Conn {
	protocol := "tcp4"
	host := "0.0.0.0"
	port := "0"
	listen := ctx.Bool("listen")
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

	log.Printf("connecting to %s:%s with protocol %s\n", host, port, protocol)

	if listen {
		ln, err := net.Listen(protocol, host+":"+port)
		failOnError(err, "Failed to open listening socket on %s:%s\n", host, port)
		defer ln.Close()

		log.Printf("Connected to %s waiting for a connection\n", ln.Addr())

		conn, err := ln.Accept()
		failOnError(err, "Failed at accepting connection [%s]\n", ln.Addr())

		return conn
	} else {
		totimer := time.NewTimer(timeout)
		defer totimer.Stop()

		for {
			conn, err := net.Dial(protocol, host+":"+port)
			if err == nil {
				return conn
			}
			select {
			case <-time.After(RETRY_INTERVAL):
				log.Printf("Failed connecting to %s:%s\n retrying to connect\n", host, port)
				continue
			case <-totimer.C:
				failOnError(err, "Failed connecting to %s:%s\n after %d ms", port, host, timeout)
			}

		}
	}
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
				Read(conn)
				conn.Close()
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
				Write(conn)
				conn.Close()
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
