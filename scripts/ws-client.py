import argparse
import json
import sys
import time

from websocket import create_connection

RESULT_POLLING_FREQUENCY = 60


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-b", "--target_branch", help="the target branch to fork and run the tests on"
    )
    parser.add_argument(
        "-c",
        "--target_commit",
        help="the target commit to checkout to run the tests on",
    )
    parser.add_argument(
        "-m",
        "--mode",
        help="the execution mode. `run` runs and waits until the results are there, `wait` just waits, and `check` just returns the current task",
        choices=["run", "wait", "check"],
        default="run",
    )
    args = parser.parse_args()
    return args


def issue_test_run(websocket, target_commit, target_branch):
    run_tests_req_data = {
        "cmd": {
            "job": "issue",
            "benchmark": "CORRECTNESS",
            "commit": target_commit,
            "branch": target_branch,
        }
    }
    msg = json.dumps(run_tests_req_data)
    websocket.send(msg)
    print(
        "POSIX Tests request made for branch:",
        target_branch,
        "and commit:",
        target_commit,
        file=sys.stderr,
    )


def fetch_runs(websocket):
    data = {"cmd": {"job": "/fetch_runs", "count": 50}}
    msg = json.dumps(data)
    # print("Sending:", msg, file=sys.stderr)
    websocket.send(msg)
    # print("Sent!", file=sys.stderr)
    res = websocket.recv()
    runs_data = json.loads(res)
    return runs_data


def current_task(websocket):
    data = {"cmd": {"job": "/current_task"}}
    msg = json.dumps(data)
    # print("Sending:", msg, file=sys.stderr)
    websocket.send(msg)
    # print("Sent!", file=sys.stderr)
    res = websocket.recv()
    res_data = json.loads(res)
    return res_data


def wait_for_result(websocket, target_commit):
    found = False
    sleep_duration = RESULT_POLLING_FREQUENCY

    while not found:
        ## Fetch all runs
        runs_data = fetch_runs(websocket)
        result_rows = runs_data["data"]["rows"]

        ## Check if target commit exists in latest runs
        for row in result_rows:
            # print(row["commit"], row["bench"], file=sys.stderr)
            if row["commit"] == target_commit and row["bench"] == "CORRECTNESS":
                print("FOUND COMMIT:", target_commit, file=sys.stderr)
                found = True
                result_row = json.dumps(row)

        if not found:
            print("Results not present for commit:", target_commit, file=sys.stderr)
            print("Sleeping for", sleep_duration, "seconds", file=sys.stderr)
            time.sleep(sleep_duration)
    return result_row


## Parse the arguments
args = parse_args()

## Create a websocket connection
ws = create_connection("ws://antikythera.csail.mit.edu:52277/")

## Just check and exit
if args.mode == "check":
    ## Debug only
    print(current_task(ws))

elif args.mode in ["run", "wait"]:
    # target_commit = "e5937bca7644c6cae0a1b8214fa980b12145bbed"
    # target_branch = "pip-libdash"
    target_commit = args.target_commit
    target_branch = args.target_branch

    if args.mode == "run":
        ## Issue the POSIX tests requests
        issue_test_run(ws, target_commit, target_branch)

    ##
    ## Wait until we have the POSIX test results
    ##
    result_row = wait_for_result(ws, target_commit)

    ## Output the relevant test
    print(result_row)


ws.close()
