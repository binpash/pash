import setuptools
from setuptools.command.install import install
import subprocess
import os
import time
import shlex
import sys

class CustomInstallCommand(install):
    """Custom install setup to help run shell commands (outside shell) before installation"""

    def run(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        os.environ["PASH_TOP"] = dir_path

        print("Make sure to have the following packages installed libtool m4 automake pkg-config libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils", file=sys.stderr)
        cmd = shlex.split("bash scripts/setup-pip.sh")
        ret = subprocess.run(cmd, stdout=f, universal_newlines=True)
        if ret.returncode != 0:
            print(f"Setup script failed", file=sys.stderr)
           

        install.run(self)


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    cmdclass={'install': CustomInstallCommand},
    name="pash",
    version='0.0.1',
    author="binpash",
    author_email="",
    description='Parallelize POSIX shell scripts.',
    url='https://github.com/binpash/pash',
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        "Bug Tracker": "https://github.com/binpash/pash/issues",
    },
    package_dir={"pash": "."},
    python_requires=">=3.8",
    packages=setuptools.find_packages('.'),
    install_requires=[
        "jsonpickle",
        "pexpect",
        "numpy",
        "matplotlib",
        "wheel"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        "Operating System :: Unix"
    ],
#     entry_points={
#         'console_scripts': [
#             "pash = 'pa.sh'",
#         ],
# },
)
