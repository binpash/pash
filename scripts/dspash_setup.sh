
# TODO: install any extra needed python debs

# TODO: install golang
rm -rf /usr/local/go && tar -C /usr/local -xzf go1.17.7.linux-amd64.tar.gz
echo -e '\nexport PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
export PATH=$PATH:/usr/local/go/bin

# TODO: install cli tool
GO111MODULE=on go get github.com/urfave/cli/v2

# TODO: compile socket_pipe
cd ../runtime/dspash
go build socket_pipe.go
