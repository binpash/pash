import os
from subprocess import PIPE, Popen
import time
import shlex
import sys

class DebugNode():
    def __init__(self, name, cmd, server):
        self.server = server
        self.cmd = cmd
        self.name = name

    def exec(self):
        proc = Popen([self.cmd], stdout=sys.stdout, stdin=sys.stdin, stderr=sys.stderr, env=os.environ, shell=True)
        return proc

    def monitor(self):
        proc = self.exec()
        while proc.poll() is None:
            print(f"{self.name} still alive", file=sys.stderr)
            time.sleep(0.01)
        print(f"{self.name} is dead", file=sys.stderr)

def main():
    name = sys.argv[1]
    cmd = sys.argv[2]
    server = sys.argv[3]
    d = DebugNode(name, cmd, server)
    d.monitor()

if __name__ == "__main__":
    main()