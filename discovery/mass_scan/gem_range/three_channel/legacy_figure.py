"""Read the 125 GeV ordinate of both printed curves in the June 1993 GEM Figure 1."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from PIL import Image


SOURCE = "https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=13"


def tick_center(line, start, stop):
    ink = np.flatnonzero(line[start:stop])
    if not 1 <= len(ink) <= 10:
        raise ValueError("Tick missing or too wide in this rendering")
    return start + (ink[0] + ink[-1]) / 2


def read_graph(image):
    ink = np.asarray(image) < 128
    if ink.shape != (2210, 1710):
        raise ValueError("Use the original TN-93-373 page at 200 dpi")
    # Locate both vertical decades at the top and bottom. The ordinates use
    # the 10 and 1 fb^-1 major ticks at both sides of the slightly skewed plot.
    ys = [382, 1500]
    xs = [350, 1370]
    x100 = [tick_center(ink[y], 466, 489) for y in ys]
    x1000 = [tick_center(ink[y], 1379, 1394) for y in ys]
    y10 = [tick_center(ink[:, x], 735, 760) for x in xs]
    y1 = [tick_center(ink[:, x], 1122, 1146) for x in xs]
    y, x = np.nonzero(ink[915:1120, 550:580])
    y, x = y + 915, x + 550
    left, right = np.interp(y, ys, x100), np.interp(y, ys, x1000)
    upper, lower = np.interp(x, xs, y10), np.interp(x, xs, y1)
    mass = 100. * 10.**((x-left)/(right-left))
    log_lumi = (lower-y)/(lower-upper)
    measurements = {}
    for name, y_center, slope in [
        ("solid_default_bottom_mass", 972., 3.2),
        ("dashed_alternative_running_bottom_mass", 1065., 2.7),
    ]:
        # These strips identify the two ink traces, separated by ~90 pixels.
        # The regression uses the actual pixels; the strips only label them.
        curve = np.abs(y - (y_center + slope * (x-564))) < 35
        fits = {}
        for width in (.75, 1., 1.5):
            within = curve & (np.abs(mass - 125.) <= width)
            centered = mass[within] - 125.
            if len(np.unique(x[within])) < 5:
                raise ValueError("Insufficient ink around the mass of interest")
            fit = np.polyfit(centered, log_lumi[within], 2)
            fits[str(width)] = float(10.**fit[-1])
        luminosity = fits["1.0"]
        measurements[name] = dict(luminosity_fb_inverse=luminosity,
                                  design_months=1.2*luminosity,
                                  fits_by_half_width_gev=fits)
    return dict(source=SOURCE, figure=1, mass_gev=125, design_fb_inverse_per_year=10,
                interpretation="Graph reading of the historical sensitivity at the specified mass. GEM gives no calibrated global five-sigma threshold across masses or three-channel breakdown at 125 GeV.",
                x_100_gev_at_two_heights=x100, x_1000_gev_at_two_heights=x1000,
                y_10_fb_inverse_at_two_abscissae=y10, y_1_fb_inverse_at_two_abscissae=y1,
                curves=measurements)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_pdf", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        prefix = Path(tmp) / "figure1"
        subprocess.run(["pdftoppm", "-f", "13", "-l", "13", "-r", "200", "-singlefile",
                        "-gray", str(args.source_pdf), str(prefix)], check=True)
        with Image.open(prefix.with_suffix(".pgm")) as image:
            output = read_graph(image)
    output["source_sha256"] = hashlib.sha256(args.source_pdf.read_bytes()).hexdigest()
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
