"""Version-checked MG5 3.5.15 loop-integration bookkeeping fix.

Only the private managed installation is patched. Previously exported processes
keep their original Python sources. No matrix elements or numerical thresholds
are changed; the upstream unstable-subjob exclusion is made internally consistent.
"""
import ast
import hashlib
import json
from pathlib import Path

PATCH_ID = "higgsssc-loop-recovery-v1"
UPSTREAM_SHA256 = "40c9a676418204332373a74a50cd0eb07330ed4bef140c4aa8fea89e40492616"
SOURCE = Path("madgraph/madevent/gen_ximprove.py")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def corrected_source(original):
    """Patch the pinned source, refusing an unreviewed MG5 version."""
    if digest(original) != UPSTREAM_SHA256:
        raise ValueError("unrecognized MG5 gen_ximprove.py; loop recovery fix requires the pinned 3.5.15 source")
    text = original.decode()

    def change(old, new, count=1):
        nonlocal text
        if text.count(old) != count:
            raise ValueError("MG5 loop recovery patch does not match its expected source context")
        text = text.replace(old, new)

    # A recursive retry needs the same results and grids. Delete neither until
    # all exclusions have been resolved and the combined result is available.
    change("            os.remove(pjoin(path, 'results.dat'))\n",
           "            # " + PATCH_ID + ": cleanup follows the final combination.\n")
    change("        # Clean up grid_information to avoid border effects in case of a crash\n",
           "        # Save original subjob indices for the event-file collectors.\n"
           "        grid_calculator.higgsssc_job_indices = [i for i in\n"
           "            range(self.splitted_for_dir(Pdir, G)) if i not in exclude_sub_jobs]\n"
           "        # Clean up grid_information to avoid border effects in case of a crash\n")
    change("                os.remove(pjoin(path, 'grid_information'))\n",
           "                for name in ('results.dat', 'grid_information'):\n"
           "                    if os.path.exists(pjoin(path, name)):\n"
           "                        os.remove(pjoin(path, name))\n")
    # After one exclusion the compact results list no longer has the original
    # job numbering. Keep the mapping for any further outlier in this iteration.
    change("            th_maxwgt = [(r.th_maxwgt,i) for i,r in enumerate(grid_calculator.results)]\n",
           "            job_indices = [i for i in range(self.splitted_for_dir(Pdir, G))\n"
           "                           if i not in exclude_sub_jobs]\n"
           "            th_maxwgt = [(r.th_maxwgt,i) for i,r in zip(job_indices, grid_calculator.results)]\n")
    change("                grid_calculator.results.run_statistics['skipped_subchannel'] += 1\n", "")
    change("        grid_calculator.results.compute_values(update_statistics=True)\n",
           "        grid_calculator.results.compute_values(update_statistics=True)\n"
           "        grid_calculator.results.run_statistics['skipped_subchannel'] += len(exclude_sub_jobs)\n")
    change("                            evtRecord.write(pjoin(gPath,'events.lhe').read())\n",
           "                            with open(pjoin(gPath,'events.lhe')) as stream:\n"
           "                                evtRecord.write(stream.read())\n")
    # Survey must still resubmit every allocated job, but final event merging
    # must use exactly the jobs that contributed to the integral.
    change("            for Gdir in Gdirs:\n                fsock.write(open(pjoin(Gdir, 'events.lhe')).read()) \n",
           "            for i in grid_calculator.higgsssc_job_indices:\n"
           "                with open(pjoin(Gdirs[i], 'events.lhe')) as stream:\n"
           "                    fsock.write(stream.read())\n"
           "            fsock.close()\n")
    change("        # collect all the generated_event\n"
           "        Gdirs = [] #build the the list of directory\n"
           "        for i in range(self.splitted_for_dir(Pdir, G)):\n"
           "            path = pjoin(Pdir, \"G%s_%s\" % (G, i+1))\n"
           "            Gdirs.append(path)\n"
           "        assert len(grid_calculator.results) == len(Gdirs) == self.splitted_for_dir(Pdir, G)\n",
           "        # Merge only events from the subjobs used in combine_grid.\n"
           "        Gdirs = [pjoin(Pdir, \"G%s_%s\" % (G, i+1))\n"
           "                 for i in grid_calculator.higgsssc_job_indices]\n"
           "        assert len(grid_calculator.results) == len(Gdirs)\n")
    ast.parse(text)
    return text.encode()


def install(mg5_dir):
    """Install once, with original source and a checksum receipt for auditing."""
    mg5_dir = Path(mg5_dir)
    source = mg5_dir / SOURCE
    backup = source.with_name(source.name + ".before-" + PATCH_ID)
    receipt = mg5_dir / (PATCH_ID + ".json")
    if source.is_symlink() or not source.resolve().is_relative_to(mg5_dir.resolve()):
        raise ValueError(f"refusing to patch a shared MG5 source: {source}")
    current = source.read_bytes()
    original = backup.read_bytes() if backup.exists() else current
    corrected = corrected_source(original)
    record = {"id": PATCH_ID, "original_sha256": digest(original),
              "patched_sha256": digest(corrected), "source": str(source.resolve()),
              "backup": str(backup.resolve())}
    if receipt.exists() and json.loads(receipt.read_text()) != record:
        raise ValueError(f"MG5 loop recovery receipt differs: {receipt}")
    if current not in (original, corrected):
        raise ValueError(f"MG5 source changed since loop recovery backup: {source}")
    if current != corrected:
        if not backup.exists():
            with backup.open("xb") as stream:
                stream.write(original)
        temporary = source.with_name(source.name + "." + PATCH_ID + ".tmp")
        with temporary.open("xb") as stream:
            stream.write(corrected)
        temporary.chmod(source.stat().st_mode)
        temporary.replace(source)
    if not receipt.exists():
        with receipt.open("x") as stream:
            stream.write(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def prepare_campaign(mg5_dir, directory):
    """Record the hotfix separately from unchanged completed-stage identities."""
    directory = Path(directory)
    record = install(mg5_dir)
    path = directory / "mg5_loop_recovery.json"
    if path.exists():
        if json.loads(path.read_text()).get("installation") != record:
            raise ValueError(f"campaign MG5 loop recovery receipt differs: {path}")
        return
    exports = {}
    for source in sorted((directory / "production").glob("*/mg5_process/bin/internal/gen_ximprove.py")):
        exports[source.parents[3].name] = {"path": str(source), "sha256": digest(source.read_bytes())}
    directory.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump({"installation": record, "exports_already_present": exports,
                   "scope": "New MG5 exports only; completed exports and event files are unchanged."}, stream, indent=2)
        stream.write("\n")
