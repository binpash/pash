install_docker() {
    # Add Docker's official GPG key:
    sudo apt-get update
    sudo apt-get install ca-certificates curl
    sudo install -m 0755 -d /etc/apt/keyrings
    sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
    sudo chmod a+r /etc/apt/keyrings/docker.asc

    # Add the repository to Apt sources:
    echo \
    "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
    $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
    sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
    sudo apt-get update

    # To install the latest version
    sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin -y
}

setup_stun_lib() {
    git clone https://yizhengx:tokennnnn@github.com/binpash/splash-stun-lib.git # TODO: make it public
    cd splash-stun-lib
    sudo docker build --progress=plain -f ./pash/Dockerfile -t pash -o out .
    sudo chmod 777 out
    sudo mv out/pashlib /opt/pashlib
    cd ..
}

setup_pash()  {
    git clone https://github.com/binpash/pash.git
    cd pash
    git switch sls-wip
    ./scripts/distro-deps.sh
    ./scripts/setup-pash.sh
    export PASH_TOP=$(pwd)
    bash ./runtime/serverless/binaries.sh
    cd ..
    cp -r $PASH_TOP/runtime/serverless/runtime ./runtime
    cp -r $PASH_TOP/runtime/serverless/aws ./aws
}

pip_install() {
    pip install boto3
}

install_docker
setup_stun_lib
setup_pash
pip_install
export AWS_ACCOUNT_ID="192165654483"
export AWS_QUEUE="queue"
export AWS_BUCKET="yizhengx"
export AWS_DEFAULT_REGION=us-east-1
export AWS_ACCESS_KEY_ID=""
export AWS_SECRET_ACCESS_KEY=""
export AWS_SESSION_TOKEN=""