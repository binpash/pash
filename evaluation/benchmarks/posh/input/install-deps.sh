pkgs='imagemagick'
if ! dpkg -s $pkgs >/dev/null 2>&1; then
  sudo apt-get install $pkgs
fi
