package main

import (
	"io"
	"log"
	"net"
	"os"
	"strings"

	"github.com/urfave/cli/v2"
)

func Port(ln net.Listener) string {
	x := strings.Split(ln.Addr().String(), ":")
	return x[1]
}

func get_conn(host string, port string, listen bool) net.Conn {
	if listen {
		ln, err := net.Listen("tcp4", host+":"+port)
		defer ln.Close()
		if err != nil {
			log.Fatalln("Trying to connect to [%s]: ERROR: %s", ln.Addr(), err)
		}
		log.Println("Connected to %s waiting for a connection", ln.Addr())

		con, err := ln.Accept()
		if err != nil {
			log.Fatalln("Trying to connect to [%s]: ERROR: %s", con.RemoteAddr(), err)
		}
		return con
	} else {
		con, err := net.Dial("tcp4", host+":"+port)
		if err != nil {
			log.Fatalln("Trying to connect to [%s]: ERROR: %s", con.RemoteAddr(), err)
		}
		return con
	}
}

func Read(host string, port string, listen bool) {
	//connect
	con := get_conn(host, port, listen)
	defer con.Close()

	//read
	n, err := io.Copy(os.Stdout, con)
	if err != nil {
		log.Fatalln("Failed while reading socket data [%s]: ERROR: %s", con.RemoteAddr(), err)
	}

	log.Printf("Read %d bytes", n)
}

func Write(host string, port string, listen bool) {
	//connect
	con := get_conn(host, port, listen)
	defer con.Close()

	//write
	n, err := io.Copy(con, os.Stdin)
	if err != nil {
		log.Fatalln("Failed while writing socket data [%s]: ERROR: %s", con.RemoteAddr(), err)
	}

	log.Printf("Wrote %d bytes", n)
}

func main() {
	app := cli.NewApp()
	app.Name = "dspash_pipe"
	app.Usage = "Connect dspash remote reads and writes (simulate unix pipe on network)"
	app.Version = "0.1"
	app.EnableBashCompletion = true

	globalFlags := []cli.Flag{
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
	app.Flags = globalFlags

	commandFlags := []cli.Flag{
		&cli.IntFlag{
			Name:  "timeout",
			Value: 1,
			Usage: "Timeout to wait for messages",
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
	}

	app.Before = func(ctx *cli.Context) error {
		if !ctx.Bool("debug") {
			log.SetOutput(io.Discard)
		}

		if ctx.IsSet("log_file") {
			file, err := os.OpenFile(ctx.String("log_file"), os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0666)
			if err != nil {
				log.Fatal(err)
			}
			log.SetOutput(file)
		}
		return nil
	}

	app.Commands = []*cli.Command{
		{
			Name:    "read",
			Aliases: []string{"r"},
			Usage:   "Read messages from socket",
			Action: func(c *cli.Context) error {
				Read(c.String("host"), c.String("port"), c.Bool("listen"))
				return nil
			},
			Flags: commandFlags,
		},
		{
			Name:    "write",
			Aliases: []string{"w"},
			Usage:   "Write to socket",
			Action: func(c *cli.Context) error {
				Write(c.String("host"), c.String("port"), c.Bool("listen"))
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
