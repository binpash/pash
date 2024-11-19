#!/usr/bin/env python3
import sys
import random

# Function to generate a random IP address
def generate_ip():
    return ".".join(str(random.randint(0, 255)) for _ in range(4))

# Function to generate random data with IP addresses and datacenter numbers
def generate_data():
    ip = generate_ip()
    num = random.randint(1, 200)
    line = f"{ip} {num}"
    return line
    
def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <number of data points>", file=sys.stderr)
        sys.exit(1)

    random.seed(42)
    n = sys.argv[1]
    n = int(n)
    for _ in range(n):
        print(generate_data())

if __name__ == "__main__":
    random.seed(42)
    main()
