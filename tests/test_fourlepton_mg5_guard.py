"""Exercise real process cleanup with disposable Python children, never MG5."""
import contextlib
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

from hfourlepton import mg5_guard as guard


class MG5GuardTests(unittest.TestCase):
    def run_child(self, folder, code, **kwargs):
        with contextlib.redirect_stdout(io.StringIO()):
            return guard.run_guarded([sys.executable, '-u', '-c', code],
                folder / 'stage.log', folder / 'state.json', 'compile', **kwargs)

    def test_success_and_quiet_work_are_allowed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self.run_child(root, "import time; time.sleep(.3); print('Poles successfully cancel for 20 points')")
            self.assertEqual(result['status'], 'complete')
            self.assertEqual(result['returncode'], 0)
            self.assertIsNone(guard.fatal_message('20 points have been tried'))

    def test_failure_stops_grandchild_in_the_owned_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            heartbeat = root / 'heartbeat'
            worker = ("import pathlib,time\np=pathlib.Path(" + repr(str(heartbeat)) + ")\n"
                      "while True:\n p.write_text(str(time.monotonic())); time.sleep(.02)")
            parent = ("import subprocess,sys,time,pathlib\n"
                      "subprocess.Popen([sys.executable, '-c', " + repr(worker) + "])\n"
                      "p=pathlib.Path(" + repr(str(heartbeat)) + ")\n"
                      "while not p.exists(): time.sleep(.01)\n"
                      "print('Stop all computation', flush=True)\ntime.sleep(30)")
            with self.assertRaisesRegex(RuntimeError, 'MG5 reported'):
                self.run_child(root, parent)
            time.sleep(.1)
            stopped = heartbeat.read_text()
            time.sleep(.2)
            self.assertEqual(heartbeat.read_text(), stopped, 'owned worker still running')

    def test_explicit_failure_stops_child_even_when_mg5_would_hang(self):
        for message in ('WARNING: 0 points have been tried',
                        'ends with non zero status: 134. Stop all computation',
                        'Error detected in launch', 'double free or corruption (out)'):
            with self.subTest(message=message), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                other = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)'], start_new_session=True)
                try:
                    started = time.monotonic()
                    with self.assertRaisesRegex(RuntimeError, 'MG5 reported:'):
                        self.run_child(root, f"import time; print({message!r}); time.sleep(30)")
                    self.assertLess(time.monotonic() - started, 3)
                    result = json.loads((root / 'state.json').read_text())
                    self.assertEqual(result['status'], 'failed')
                    self.assertIsNotNone(result['returncode'])
                    self.assertIsNone(other.poll(), 'unrelated process was stopped')
                finally:
                    other.terminate()
                    other.wait()

    def test_error_with_success_exit_and_nonzero_exit_are_rejected(self):
        for code in ("print('Error detected in launch')", 'raise SystemExit(7)'):
            with self.subTest(code=code), tempfile.TemporaryDirectory() as tmp:
                with self.assertRaises(RuntimeError):
                    self.run_child(Path(tmp), code)

    def test_timeout_also_applies_to_continuous_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            started = time.monotonic()
            with self.assertRaisesRegex(RuntimeError, 'time limit'):
                self.run_child(Path(tmp), "import time\nwhile True:\n print('working'); time.sleep(.01)", timeout=.4)
            self.assertLess(time.monotonic() - started, 3)

    def test_existing_attempt_log_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'stage.log').write_text('previous evidence')
            with self.assertRaises(FileExistsError):
                self.run_child(root, "print('new work')")
            self.assertEqual((root / 'stage.log').read_text(), 'previous evidence')

    def test_user_interrupt_cleans_up_owned_child(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = [sys.executable, guard.__file__, '--stage', 'compile',
                       '--log', str(root / 'stage.log'), '--state', str(root / 'state.json'),
                       '--', sys.executable, '-c', 'import time; time.sleep(30)']
            with (root / 'supervisor.log').open('w') as output:
                supervisor = subprocess.Popen(command, stdout=output, stderr=output, start_new_session=True)
                try:
                    deadline = time.monotonic() + 5
                    while not (root / 'state.json').exists():
                        if supervisor.poll() is not None or time.monotonic() > deadline:
                            self.fail('supervisor did not start')
                        time.sleep(.02)
                    os.kill(supervisor.pid, signal.SIGINT)
                    self.assertEqual(supervisor.wait(timeout=7), 130)
                    state = json.loads((root / 'state.json').read_text())
                    self.assertEqual(state['status'], 'failed')
                    with self.assertRaises(ProcessLookupError):
                        os.kill(state['child_pid'], 0)
                finally:
                    if supervisor.poll() is None:
                        supervisor.kill()
                        supervisor.wait()

    def test_status_ignores_errors_from_an_earlier_attempt(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'managed.log'
            path.write_text('\nCOMMAND: old\nError detected in launch\nCOMMAND: new\ncompiling\n')
            self.assertIsNone(guard.current_log_failure(path))
            with path.open('a') as stream:
                stream.write('WARNING: 0 points have been tried\n')
            self.assertEqual(guard.current_log_failure(path), '0 points have been tried')


if __name__ == '__main__':
    unittest.main()
