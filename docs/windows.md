# Windows Installation
To setup and install PaSh on a Windows machine, `Windows Subsystem for Linux (WSL)` needs to be installed.

## What is WSL
The Windows Subsystem for Linux lets developers run a GNU/Linux environment 
-- including most command-line tools, utilities, and applications -- 
directly on Windows, unmodified, without the overhead of a traditional virtual 
machine or dualboot setup.

## Setup WSL

- Open PowerShell as Administrator and run:

    `dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart`

- Enable Linux Filesystem for Windows

    `dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart`

- Install the distro of your choice:

  - [Ubuntu 20.04-lts](https://www.microsoft.com/el-gr/p/ubuntu-2004-lts/9n6svws3rx71?rtc=1&activetab=pivot:overviewtab)
  - [Debian GNU/Linux](https://www.microsoft.com/el-gr/p/debian/9msvkqc78pk6?rtc=1&activetab=pivot:overviewtab)
  - [Fedora Remix](https://www.microsoft.com/el-gr/p/fedora-remix-for-wsl/9n6gdm4k2hnc?rtc=1&activetab=pivot:overviewtab)
  - [Available Distros](https://docs.microsoft.com/en-us/windows/wsl/install-win10#step-6---install-your-linux-distribution-of-choice)

- Reboot your computer
- On Windows search menu, find wsl
- It will take some time to initialize and then it will ask a username/password
- After this point, you may run full Linux applications on your Windows machine!
- Continue the PaSh installation from [here](https://github.com/andromeda/pash/blob/main/docs/tutorial.md#installation)


## Docker TODO
