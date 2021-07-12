# PaSh installation on Windows
To setup and install PaSh on a Windows machine, `Windows Subsystem for Linux (WSL)` needs to be installed first.

## What is WSL
The Windows Subsystem for Linux lets developers run a GNU/Linux environment 
— including most command-line tools, utilities, and applications — 
directly on Windows, unmodified, without the need to setup a traditional virtual 
machine or use dualboot.

## Setup WSL
You can either install WSL or WSL2. WSL2 has more features (notably, it includes the full Linux kernel),
but confilcts with other virtualization software such as VirtualBox or VMware, as it requires Hyper-V to work.
Here are the short instructions that assume installing WSL2 on a compatible Windows 10 machine running on x64.

In case of errors, consult the full installation guide on [Microsoft Docs](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
and the [WSL Troubleshooting Guide](https://github.com/MicrosoftDocs/WSL/blob/master/WSL/troubleshooting.md).

### Installing WSL2 — short guide
All commands below should be run in an Administrator Poweshell prompt.

- If Running Windows 10 Preview build 20262 or higher:
    - List the available Linux Distributions
        ```powershell
        wsl --list --online
        ```
    - Install WSL2 with with a given distribution
        ```powershell
        wsl --install -d <distribution-name-from-above>
        ```
    - Restart your machine
        ```powershell
        Restart-Computer
        ```
- Else (most people):
    - Enable WSL
        ```powershell
        dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
        ```
    - Enable Virtual Machine Platform (required for Hyper-V and WSL2)
        ```powershell
        dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
        ```
    - Restart your machine
       ```powershell
       Restart-Computer
       ```
    - Update WSL to WSL2
        - Download the WSL to WSL2 updater and run it
            ```powershell
            (New-Object System.Net.WebClient).DownloadFile("https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi", "C:/Windows/Temp/wsl2.msi");
            msiexec /i "C:/Windows/Temp/wsl2.msi";
            Remove-Item "C:/Windows/Temp/wsl2.msi"
            ```
        - Set WSL2 as default WSL version
            ```powershell
            wsl --set-default-version 2
            ```
    - Install a Linux distribution from the Microsoft Store (cannot be done using command line, 
      unless using something like Chocolatey, but in that case the whole process is automated):
        - [Ubuntu 18.04 LTS](https://www.microsoft.com/store/apps/9N9TNGVNDL3Q)
        - [Ubuntu 20.04 LTS](https://www.microsoft.com/store/apps/9n6svws3rx71)
        - [openSUSE Leap 15.1](https://www.microsoft.com/store/apps/9NJFZK00FGKV)
        - [SUSE Linux Enterprise Server 12 SP5](https://www.microsoft.com/store/apps/9MZ3D1TRP8T1)
        - [SUSE Linux Enterprise Server 15 SP1](https://www.microsoft.com/store/apps/9PN498VPMF3Z)
        - [Kali Linux](https://www.microsoft.com/store/apps/9PKR34TNCV07)
        - [Debian GNU/Linux](https://www.microsoft.com/store/apps/9MSVKQC78PK6)
        - [Fedora Remix for WSL](https://www.microsoft.com/store/apps/9n6gdm4k2hnc)
        - [Pengwin](https://www.microsoft.com/store/apps/9NV1GV1PXZ6P)
        - [Pengwin Enterprise](https://www.microsoft.com/store/apps/9N8LP0X93VCP)
        - [Alpine WSL](https://www.microsoft.com/store/apps/9p804crf0395)
    - Restart your machine
       ```powershell
       Restart-Computer
       ```
- Run the `wsl` command or find the installed Linux distribution in Windows Start menu and run it.
- After a few minutes of installation, enter a username and password for the internal WSL account to be created.

Continue the PaSh installation process from [here](https://github.com/andromeda/pash/blob/main/docs/tutorial.md#installation)
inside the WSL installation.


## Docker TODO
