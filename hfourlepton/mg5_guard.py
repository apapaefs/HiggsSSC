"""Run an owned MG5 process group and fail promptly on explicit fatal output.

Quiet integration is allowed. Only a reported failure, a nonzero exit, an
explicit diagnostic time limit, or a user interrupt stops the child group.
"""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import signal
import socket
import subprocess
import sys
import time


FATAL = re.compile(
    r"\b0 points have been tried|Error detected in|interrupted with error:|"
    r"Poles do not cancel|Some tests failed, run cannot continue|"
    r"Stop all computation|ends with non zero status:\s*\S+|"
    r"STOP No (?:dp|qp) loop library|double free or corruption|"
    r"Program received signal SIG(?:ABRT|SEGV|FPE)|Segmentation fault",
    re.IGNORECASE)


def fatal_message(text):
    text = re.sub(r"\x1b\[[0-9;]*m", "", text)
    match = FATAL.search(text)
    return match[0] if match else None


def current_log_failure(path):
    """Inspect the most recent attempt in an appended managed-stage log."""
    path = Path(path)
    if not path.is_file():
        return None
    with path.open('rb') as stream:
        stream.seek(max(0, path.stat().st_size - 1024 * 1024))
        text = stream.read().decode(errors='replace')
    return fatal_message(text.rsplit('\nCOMMAND:', 1)[-1])


def save_state(path, state):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(state, indent=2) + '\n')
    temporary.replace(path)


def stop_owned_group(child):
    """Never signal a process group that this invocation did not create."""
    def send(signum):
        try:
            os.killpg(child.pid, signum)
        except ProcessLookupError:
            pass
        except PermissionError:
            # Some sandboxed hosts report EPERM for an already exited group.
            # A running child must still be treated as a cleanup failure.
            if child.poll() is None:
                raise

    send(signal.SIGTERM)
    try:
        child.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    # The parent may exit before its workers. Clean up that same owned group.
    send(signal.SIGKILL)
    child.wait()


def run_guarded(command, log, state_path, stage, *, cwd=None, timeout=None, program='MG5'):
    if timeout is not None and (not math.isfinite(timeout) or timeout <= 0):
        raise ValueError('diagnostic timeout must be positive and finite')
    log, state_path = Path(log), Path(state_path)
    log.parent.mkdir(parents=True, exist_ok=True)
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state = dict(stage=stage, status='starting', host=socket.gethostname(),
                 supervisor_pid=os.getpid(), command=[str(x) for x in command],
                 log=str(log.resolve()), started_at=datetime.now(timezone.utc).isoformat())
    started = time.monotonic()
    child = None
    pending = ''
    try:
        # Each export/compile/launch has its own log; never replace an old attempt.
        with log.open('x') as output, log.open(errors='replace') as reader:
            child = subprocess.Popen(state['command'], cwd=cwd, stdout=output,
                                     stderr=subprocess.STDOUT, start_new_session=True)
            state.update(status='running', child_pid=child.pid)
            save_state(state_path, state)
            while True:
                if timeout is not None and time.monotonic() - started >= timeout:
                    raise RuntimeError(f'diagnostic time limit ({timeout:g} seconds) reached')
                chunk = reader.read(65536)
                if chunk:
                    print(chunk, end='', flush=True)
                    pending = pending[-1024:] + chunk
                    failure = fatal_message(pending)
                    if failure:
                        raise RuntimeError(program + ' reported: ' + failure)
                elif child.poll() is not None:
                    break
                else:
                    time.sleep(.2)
            if child.returncode:
                raise RuntimeError(f'{program} exited with status {child.returncode}')
        state.update(status='complete', returncode=child.returncode)
    except BaseException as error:
        state.update(status='failed', error=str(error) or 'command interrupted',
                     returncode=child.returncode if child else None)
        if child is not None:
            try:
                stop_owned_group(child)
            except OSError as cleanup_error:
                state['cleanup_error'] = str(cleanup_error)
            state['returncode'] = child.returncode
        raise
    finally:
        state.update(finished_at=datetime.now(timezone.utc).isoformat(),
                     elapsed_seconds=time.monotonic() - started)
        save_state(state_path, state)
    return state


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', required=True, type=Path)
    parser.add_argument('--state', required=True, type=Path)
    parser.add_argument('--stage', required=True)
    parser.add_argument('--timeout', type=float)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('a command is required after --')

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f'interrupted by signal {signum}')

    signal.signal(signal.SIGTERM, interrupted)
    try:
        run_guarded(command, args.log, args.state, args.stage, timeout=args.timeout)
    except (Exception, KeyboardInterrupt) as error:
        print(f'ERROR: {args.stage} failed: {error}; log: {args.log}', file=sys.stderr)
        return 130 if isinstance(error, KeyboardInterrupt) else 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
