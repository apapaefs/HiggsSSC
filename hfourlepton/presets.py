"""Versioned event budgets for the managed student campaigns (no generators)."""
from decimal import Decimal, InvalidOperation, ROUND_CEILING
import hashlib

VERSION = 1
CHANNELS = ("4e", "4mu", "2e2mu")
MINIMAL_EVENTS = {
    "signal_gg_h": 20000,
    "signal_vbf_h": 10000,
    "signal_wp_h": 5000,
    "signal_wm_h": 5000,
    "signal_z_h": 5000,
    "signal_gg_z_h": 5000,
    "signal_tt_h": 5000,
    "bkg_qq": 20000,
    "bkg_gg": 2000,
}


def multiplier(value):
    try:
        result = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError("--event-multiplier must be positive and finite") from exc
    if not result.is_finite() or result <= 0:
        raise ValueError("--event-multiplier must be positive and finite")
    if result > Decimal(2147483647 // 3) / 20000:
        raise ValueError("--event-multiplier would exceed the generator event-count limit")
    return result


def counts(value=1):
    scale = multiplier(value)
    return {f"{group}_{channel}": int((scale * count).to_integral_value(rounding=ROUND_CEILING))
            for group, count in MINIMAL_EVENTS.items() for channel in CHANNELS}


def seed_base(run_tag):
    # New tags are independent campaigns, not overlapping prefixes to concatenate.
    return 1 + int(hashlib.sha256(run_tag.encode()).hexdigest()[:12], 16) % 800000000


def description(args):
    if getattr(args, "preset", "standard") != "minimal":
        return None
    jobs = getattr(args, "jobs", getattr(args, "nb_core", 1))
    workers = min(8, jobs)
    return dict(version=VERSION, preset="minimal", event_multiplier=str(multiplier(args.event_multiplier).normalize()),
                sample_nevents=counts(args.event_multiplier), total_events=sum(counts(args.event_multiplier).values()),
                req_acc=0.03, loop_survey_accuracy=0.05, survey_splitting=jobs, nhel=1,
                loop_reduction_gg4l="7|6|1", loop_reduction_other="6|7|1",
                continuum_m4l_gev=None if args.wide_continuum else [90., 180.],
                postprocess_workers=workers, cores_per_decay=max(1, jobs // workers),
                shared_shower=True, seed_base=seed_base(args.run_tag) if args.run_tag else None)


def arguments(args):
    if getattr(args, "preset", "standard") != "minimal":
        return []
    return ["--preset", "minimal", "--event-multiplier", str(args.event_multiplier),
            *(["--wide-continuum"] if args.wide_continuum else [])]


def max_events(args):
    return max(counts(args.event_multiplier).values()) if getattr(args, "preset", "standard") == "minimal" else args.nevents
