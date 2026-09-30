import argparse
import os
import subprocess
import time

import config

# One shared, authenticated SSH connection for every transfer in this pipeline.
#
# Without this, each rsync opens its own SSH session and asks for the password again --
# 10 prompts for 5 samples (filt + full), spread over hours because the full files are
# fetched lazily between MSMuTect runs. Instead we open a single OpenSSH "control master"
# up front (one password prompt) and point every later rsync at its control socket, so
# they piggyback on the already-authenticated connection.
#
# The master lingers in the background for ControlPersist, so re-running a script -- or
# running the steps as separate commands -- reuses the same login. Tear it down early with
#     python ssh_connection.py --close

CONTROL_DIR = os.path.expanduser("~/.ssh/controlmasters")
CONTROL_PATH = os.path.join(CONTROL_DIR, f"{config.REMOTE_USER}@{config.REMOTE_HOST}:22")
CONTROL_PERSIST = "4h"  # long enough to outlive a full 5-sample run, short enough to self-clean


def target() -> str:
    return f"{config.REMOTE_USER}@{config.REMOTE_HOST}"


def master_is_alive() -> bool:
    result = subprocess.run(["ssh", "-S", CONTROL_PATH, "-O", "check", target()],
                            capture_output=True, text=True)
    return result.returncode == 0


def ensure_master():
    # Opens the shared connection if it isn't already up. This is the ONLY place that can
    # ask for a password; everything afterwards rides the resulting socket.
    if master_is_alive():
        return
    os.makedirs(CONTROL_DIR, mode=0o700, exist_ok=True)
    print(f"opening a shared SSH connection to {target()} -- you will be asked for your password once")
    command = ["ssh", "-M", "-S", CONTROL_PATH,
               "-o", f"ControlPersist={CONTROL_PERSIST}",
               "-N", "-f", target()]
    # deliberately NOT capture_output: the password prompt has to reach your terminal
    result = subprocess.run(command)
    # `-f` backgrounds ssh once the socket is up, but give it a moment rather than racing it
    if result.returncode == 0:
        for _ in range(10):
            if master_is_alive():
                break
            time.sleep(0.2)
    if result.returncode != 0 or not master_is_alive():
        raise RuntimeError(
            f"could not open a shared SSH connection to {target()} (exit code {result.returncode}).\n"
            f"Run this from an interactive terminal so you can type your password."
        )
    print(f"shared SSH connection established (persists ~{CONTROL_PERSIST} after last use)")


def rsync_over_master(remote_path: str, local_path: str) -> subprocess.CompletedProcess:
    ensure_master()
    remote = f"{target()}:{remote_path}"
    # ControlMaster=no => "use the existing socket, don't try to become a master yourself"
    ssh_command = f"ssh -S {CONTROL_PATH} -o ControlMaster=no"
    return subprocess.run(["rsync", "-avP", "--partial", "-e", ssh_command, remote, local_path])


def close_master():
    if not master_is_alive():
        print("no shared SSH connection is open")
        return
    subprocess.run(["ssh", "-S", CONTROL_PATH, "-O", "exit", target()], capture_output=True)
    print(f"closed the shared SSH connection to {target()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--close", action="store_true", help="close the shared SSH connection and exit")
    parser.add_argument("--check", action="store_true", help="report whether the shared connection is open")
    args = parser.parse_args()

    if args.close:
        close_master()
    elif args.check:
        print("open" if master_is_alive() else "closed")
    else:
        ensure_master()


if __name__ == "__main__":
    main()
