
# TODO: install any extra needed python debs

# Get PASH_TOP
if git rev-parse --git-dir > /dev/null 2>&1; then
    # set PASH_TOP
    PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
else
    # set PASH_TOP to the root folder of the project if it is not available
    PASH_TOP=${PASH_TOP:-$PWD/..}
fi

# Install Go
wget https://go.dev/dl/go1.17.7.linux-amd64.tar.gz
rm -rf /usr/local/go && tar -C /usr/local -xzf go1.17.7.linux-amd64.tar.gz
echo -e '\nexport PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
export PATH=$PATH:/usr/local/go/bin
rm go1.17.7.linux-amd64.tar.gz

# Install deps
GO111MODULE=on go get github.com/urfave/cli/v2

# Compile runtime
cd $PASH_TOP/runtime/dspash
go build socket_pipe.go
