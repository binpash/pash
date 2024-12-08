import re
from enum import Enum

def parse_diff(file_path):
    with open(file_path, 'r') as f:
        diff_lines = f.readlines()

    results = []
    current_file = None
    current_diff = {"original": [], "new": []}
    command_regex = re.compile(r"Daemon:.*?class: ([\w-]+) found for: (.+)")

    for line in diff_lines:
        if line.startswith("Comparing:"):
            # Save the current diff block
            if current_file:
                results.append((current_file, current_diff))
            # Start a new diff block
            current_file = line.strip().split(" ")[1]
            current_diff = {"original": [], "new": []}
        elif line.startswith("<"):
            match = command_regex.search(line)
            if match:
                current_diff["original"].append((match.group(2), match.group(1)))
        elif line.startswith(">"):
            match = command_regex.search(line)
            if match:
                current_diff["new"].append((match.group(2), match.group(1)))

    # Add the last diff block
    if current_file:
        results.append((current_file, current_diff))

    return results

class Parallelizability(Enum):
    S = "stateless"
    P = "pure"
    N = "non-pure"
    E = "side-effectful"

    @classmethod
    def _missing_(cls, value):
        if value == "parallelizable_pure":
            return Parallelizability.P

    def __gt__(self, other):
        if other is self:
            return False
        else:
            return self >= other

    def __ge__(self, other):
        match self:
            case Parallelizability.S:
                return True
            case Parallelizability.P:
                return other is not Parallelizability.S
            case Parallelizability.N:
                return other is self or other is Parallelizability.E
            case Parallelizability.E:
                return other is self



def compare_classes(results):
    for file, diff in results:
        original_commands = {cmd: cls for cmd, cls in diff["original"]}
        new_commands = {cmd: cls for cmd, cls in diff["new"]}

        all_commands = set(original_commands.keys()).union(new_commands.keys())
        found_commands = set(original_commands.keys()).intersection(new_commands.keys())
        for cmd in all_commands:
            original_cls = original_commands.get(cmd)
            new_cls = new_commands.get(cmd)

            if original_cls == "pure": 
                original_cls = "non-pure"

            if original_cls and new_cls:
                print(f"Command: {cmd}")
                print(f"  Original Class: {original_cls}")
                print(f"  New Class: {new_cls}")
                print(f"  Is okay: {Parallelizability(new_cls) <= Parallelizability(original_cls)}")
                print(f"  Is same: {Parallelizability(new_cls) == Parallelizability(original_cls)}")
                print(f"  Is too low: {Parallelizability(new_cls) == Parallelizability.E and Parallelizability(original_cls) != Parallelizability.E}")
            # elif original_cls:
            #     print(f"Command: {cmd} only in Original")
            #     print(f"  Original Class: {original_cls}")
            # elif new_cls:
            #     print(f"Command: {cmd} only in New")
            #     print(f"  New Class: {new_cls}")
        print("-" * 40)

if __name__ == "__main__":
    diff_file = "diffs"  
    diff_results = parse_diff(diff_file)
    compare_classes(diff_results)

