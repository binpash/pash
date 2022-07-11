For this version of auto-checker, I already did following OS's compatibility check and pass them.

``
ubuntu:latest ubuntu:20.04 ubuntu:18.04 debian:latest fedora:latest archlinux:latest
``

Script files are all inside folder:
``$PATH_OF_PASH/scripts/offline_testing``

For compability check, please run ``testing-script.sh`` inside the above folder.

I already wrote comments inside these files, please check.



Here are some minor changes for distro-deps.sh and setup-pash.sh files:

1. Already added loops for installing packages for different distros in ``testing-script.sh``

2. Added ``which`` and ``words`` packages for ``fedora:latest`` in ``testing-script.sh``
3. Added ``which`` package for ``archlinux:latest`` in  ``testing-script.sh``

4. Added ``python3 -m pip install --upgrade && pip python3 -m pip install --upgrade Pillow`` for  ``setup-pash.sh``. This is because for ``ubuntu:18.04``, ``Pillow`` for some reasons not being installed properly for the execution of Pash.

