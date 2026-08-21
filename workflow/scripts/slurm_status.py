#!/usr/bin/env python3
"""Translate a SLURM job state for Snakemake's cluster-generic executor."""

from __future__ import annotations

import argparse
import subprocess


SUCCESS_STATES = {"COMPLETED"}
FAILURE_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "SPECIAL_EXIT",
    "TIMEOUT",
}


def normalize_state(state: str) -> str:
    """Remove details such as ``CANCELLED by 123`` and trailing ``+``."""

    return state.strip().upper().split()[0].rstrip("+") if state.strip() else ""


def classify_state(state: str, exit_code: str = "") -> str:
    """Return one of the three values required by cluster-generic."""

    normalized = normalize_state(state)
    if normalized in SUCCESS_STATES:
        return "success" if not exit_code or exit_code == "0:0" else "failed"
    if normalized in FAILURE_STATES:
        return "failed"
    return "running"


def run_command(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, capture_output=True, text=True, check=False)


def current_queue_states(job_id: str) -> list[str]:
    result = run_command(
        ["squeue", "--noheader", "--jobs", job_id, "--format", "%T"]
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def accounting_record(job_id: str) -> tuple[str, str] | None:
    result = run_command(
        [
            "sacct",
            "--noheader",
            "--parsable2",
            "--jobs",
            job_id,
            "--format",
            "JobIDRaw,State,ExitCode",
        ]
    )
    fallback: tuple[str, str] | None = None
    for line in result.stdout.splitlines():
        columns = line.split("|")
        if len(columns) < 3:
            continue
        record_id, state, exit_code = (column.strip() for column in columns[:3])
        if not state:
            continue
        if fallback is None and "." not in record_id:
            fallback = (state, exit_code)
        if record_id == job_id:
            return state, exit_code
    return fallback


def job_status(job_id: str) -> str:
    # --parsable can append a cluster name as ``jobid;cluster``.
    job_id = job_id.split(";", 1)[0]

    queue_states = current_queue_states(job_id)
    if queue_states:
        statuses = [classify_state(state) for state in queue_states]
        if "failed" in statuses:
            return "failed"
        if all(status == "success" for status in statuses):
            return "success"
        return "running"

    record = accounting_record(job_id)
    if record is None:
        # Accounting data can lag briefly after a job leaves squeue. Reporting
        # running makes Snakemake poll again instead of declaring a false error.
        return "running"
    return classify_state(*record)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("job_id")
    args = parser.parse_args()
    print(job_status(args.job_id))


if __name__ == "__main__":
    main()
