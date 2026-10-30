"""Run the installed MadSpin bridge with a bounded, private core allocation."""
import argparse
import logging.config
from pathlib import Path
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mg5-dir", required=True, type=Path)
    parser.add_argument("--cores", required=True, type=int)
    parser.add_argument("card", type=Path)
    args = parser.parse_args(argv)
    if args.cores < 1:
        parser.error("--cores must be positive")
    sys.path.insert(0, str(args.mg5_dir.resolve()))
    from MadSpin.interface_madspin import MadSpinInterface
    logging.config.fileConfig(args.mg5_dir / "madgraph/interface/.mg5_logging.conf")

    class BoundedMadSpin(MadSpinInterface):
        def do_launch(self, line):
            # Apply after importing the production banner: it may contain its
            # own nb_core setting. Never modify the installation configuration.
            for command in (f"set nb_core {args.cores} --no_save",
                            "set run_mode 2 --no_save", "set automatic_html_opening False --no_save"):
                self.mg5cmd.exec_cmd(command, printcmd=False)
            return super().do_launch(line)

    cmd = BoundedMadSpin()
    cmd.use_rawinput = False
    cmd.haspiping = False
    cmd.import_command_file(str(args.card.resolve()))
    cmd.run_cmd("quit")


if __name__ == "__main__":
    main()
