#!/usr/bin/env python3

import argparse
import sys
from typing import Optional

import boto3
from botocore.exceptions import BotoCoreError, ClientError


DEFAULT_LOG_GROUP = "/aws/lambda/lambda"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Delete all CloudWatch log streams in a log group."
    )
    parser.add_argument(
        "--log-group",
        default=DEFAULT_LOG_GROUP,
        help=f"CloudWatch log group name (default: {DEFAULT_LOG_GROUP})",
    )
    parser.add_argument(
        "--stream-prefix",
        help="Only delete streams whose names start with this prefix.",
    )
    parser.add_argument(
        "--region",
        help="AWS region (uses default AWS config if omitted).",
    )
    parser.add_argument(
        "--profile",
        help="AWS profile name (uses default AWS config if omitted).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List matching streams without deleting them.",
    )
    return parser.parse_args()


def make_logs_client(profile: Optional[str], region: Optional[str]):
    session_kwargs = {}
    client_kwargs = {}

    if profile:
        session_kwargs["profile_name"] = profile
    if region:
        client_kwargs["region_name"] = region

    session = boto3.session.Session(**session_kwargs)
    return session.client("logs", **client_kwargs)


def describe_first_page(logs_client, log_group: str, stream_prefix: Optional[str]):
    request = {"logGroupName": log_group}
    if stream_prefix:
        request["logStreamNamePrefix"] = stream_prefix
    return logs_client.describe_log_streams(**request)


def main() -> int:
    args = parse_args()

    try:
        logs_client = make_logs_client(args.profile, args.region)
    except Exception as exc:  # boto3 session/client setup failures
        print(f"Failed to create CloudWatch Logs client: {exc}", file=sys.stderr)
        return 1

    deleted = 0
    rounds = 0

    print(f"Log group: {args.log_group}")
    if args.stream_prefix:
        print(f"Stream prefix filter: {args.stream_prefix}")
    if args.dry_run:
        print("Dry run enabled: no streams will be deleted.")

    while True:
        try:
            response = describe_first_page(logs_client, args.log_group, args.stream_prefix)
        except logs_client.exceptions.ResourceNotFoundException:
            print(f"Log group not found: {args.log_group}", file=sys.stderr)
            return 1
        except (ClientError, BotoCoreError) as exc:
            print(f"Failed to list log streams: {exc}", file=sys.stderr)
            return 1

        streams = response.get("logStreams", [])
        if not streams:
            break

        rounds += 1
        print(f"Round {rounds}: {len(streams)} stream(s)")

        for stream in streams:
            stream_name = stream["logStreamName"]

            if args.dry_run:
                print(f"[dry-run] {stream_name}")
                deleted += 1
                continue

            try:
                logs_client.delete_log_stream(
                    logGroupName=args.log_group,
                    logStreamName=stream_name,
                )
                deleted += 1
            except (ClientError, BotoCoreError) as exc:
                print(f"Failed to delete stream '{stream_name}': {exc}", file=sys.stderr)
                return 1

    action = "Matched" if args.dry_run else "Deleted"
    print(f"{action} {deleted} log stream(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
