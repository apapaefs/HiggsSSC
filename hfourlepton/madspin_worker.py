"""Run the installed MadSpin bridge with a bounded, private core allocation."""
import argparse
import logging.config
from pathlib import Path
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mg5-dir", required=True, type=Path)
    parser.add_argument("--cores", required=True, type=int)
    parser.add_argument("--check", action="store_true", help="check startup and private core settings without generating events")
    parser.add_argument("card", type=Path, nargs="?")
    args = parser.parse_args(argv)
    if args.cores < 1:
        parser.error("--cores must be positive")
    if (args.card is None) != args.check:
        parser.error("provide either a MadSpin card or --check")
    sys.path.insert(0, str(args.mg5_dir.resolve()))
    # The stock MadSpin frontend registers logging.ColorFormatter this way.
    import madgraph.interface.coloring_logging  # noqa: F401
    logging.config.fileConfig(args.mg5_dir / "madgraph/interface/.mg5_logging.conf")
    from MadSpin.interface_madspin import MadSpinInterface

    class BoundedMadSpin(MadSpinInterface):
        def configure_cores(self):
            # Apply after importing the production banner: it may contain its
            # own nb_core setting. Never modify the installation configuration.
            for command in (f"set nb_core {args.cores} --no_save",
                            "set run_mode 2 --no_save", "set automatic_html_opening False --no_save"):
                self.mg5cmd.exec_cmd(command, printcmd=False)

        def do_launch(self, line):
            self.configure_cores()
            return super().do_launch(line)

    cmd = BoundedMadSpin()
    cmd.use_rawinput = False
    cmd.haspiping = False
    if args.check:
        cmd.configure_cores()
        options = cmd.mg5cmd.options
        if (int(options["nb_core"]) != args.cores or int(options["run_mode"]) != 2
                or options["automatic_html_opening"] is not False):
            raise RuntimeError("MadSpin did not accept the requested private core settings")
        print(f"MadSpin startup passed; {args.cores} private cores; no events generated.")
        cmd.run_cmd("quit")
        return
    cmd.import_command_file(str(args.card.resolve()))
    cmd.run_cmd("quit")


if __name__ == "__main__":
    main()
