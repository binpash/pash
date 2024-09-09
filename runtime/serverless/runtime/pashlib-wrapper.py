# write a wrapper to run a command read from stdin and write to stdout

import sys
import subprocess


def run_command(command):
    """
    Run a command reading from stdin and writing to stdout.

    :param command: The command to run as a list of arguments.
    """
    # Run the command with stdin, stdout, and stderr connected directly
    # try: 
    #     result = subprocess.run(command, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr)
    # except Exception as e:
    #     print(f"[pashlib-wrapper.py] Error running command {command}: {e}")
    counter = 0
    if "send" in command:
        try:
            process = subprocess.Popen(command, stdin=subprocess.PIPE, text=False)     
            while True:
                chunk = sys.stdin.buffer.read(128)
                if not chunk:
                    break
                process.stdin.write(chunk)
                process.stdin.flush()  # Ensure the data is sent to the subprocess
                print(f"[pashlib-wrapper.py] {command} - Sent chunk {counter}", flush=True)
                counter += 1
            return 0
        except Exception as e:
            print(f"[pashlib-wrapper.py] Error sending data to command {command}: {e}")
            return 1
    else:
        try:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, text=False)     
            while True:
                chunk = process.stdout.read(128)
                if not chunk:
                    break
                sys.stdout.buffer.write(chunk)
                sys.stdout.flush()  # Ensure the data is sent to the subprocess
                print(f"[pashlib-wrapper.py] {command} - Received chunk {counter}", flush=True)
                counter += 1
            return 0
        except Exception as e:
            print(f"[pashlib-wrapper.py] Error receiving data from command {command}: {e}")
            return 1

# Example usage
if __name__ == "__main__":
    # Read the command from arguments
    print(sys.argv)
    command = sys.argv[1:]
    print(f"[pashlib-wrapper.py] Start running command: {command}")
    exit_code = run_command(command)
    print(f"[pashlib-wrapper.py] {command} Command exited with code {exit_code}")