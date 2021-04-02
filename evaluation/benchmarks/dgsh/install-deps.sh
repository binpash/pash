sudo dpkg --add-architecture i386
sudo apt-get update
echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | sudo debconf-set-selections
sudo apt-get install -y csh
sudo apt-get install -y default-jdk
sudo apt-get install -y default-jre
sudo apt-get install -y libc6:i386
sudo apt-get install -y libstdc++6:i386
sudo apt-get install -y libx11-6:i386
sudo apt-get install -y libxext6:i386
sudo apt-get install -y msttcorefonts
