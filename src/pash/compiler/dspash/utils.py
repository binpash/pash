import subprocess
import os
import tempfile
import uuid


def read_file(file, mode="r"):
    with open(file, mode) as f:
        data = f.read()
    return data


def write_file(file, data, mode="w"):
    with open(file, mode) as f:
        n = f.write(data)
    return n


def create_filename(dir, prefix="", temp=False):
    if temp:
        return tempfile.mkstemp(dir=dir, prefix=prefix)
    else:
        return os.path.join(dir, f"{prefix}_{uuid.uuid4()}")
