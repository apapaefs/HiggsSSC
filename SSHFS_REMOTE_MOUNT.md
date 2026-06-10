# Viewing Remote Results With SSHFS

This note is for viewing campaign outputs that were produced on
`timur.kennesaw.edu` from a local Linux machine.  SSHFS mounts a remote
directory over SSH, so your local browser can open the generated HTML report
and its image assets without copying the whole campaign directory first.

## 1. Install SSHFS Locally

On your local Linux machine, install SSHFS if it is not already available.

For Ubuntu or Debian:

```bash
sudo apt update
sudo apt install sshfs
```

For Fedora, Red Hat, or similar systems:

```bash
sudo dnf install sshfs
```

If the package is not found on a managed machine, ask the system administrator
whether SSHFS/FUSE is available.

## 2. Create A Local Mount Point

Choose a local directory where the remote repository will appear:

```bash
mkdir -p ~/mnt/timur-higgsssc
```

## 3. Mount The Remote Repository

If the repository was cloned on Timur as `~/HiggsSSC`, run:

```bash
sshfs USERNAME@timur.kennesaw.edu:~/HiggsSSC ~/mnt/timur-higgsssc \
  -o reconnect,ServerAliveInterval=15,ServerAliveCountMax=3,follow_symlinks
```

Replace `USERNAME` with your Timur username.

If the repository is somewhere else on Timur, replace `~/HiggsSSC` with the
correct remote path.  You can check the path by logging in:

```bash
ssh USERNAME@timur.kennesaw.edu
pwd
ls
```

## 4. Open The Gamma-Gamma Report

After mounting, the remote files look local.  For the default `run_01` report,
open:

```bash
xdg-open ~/mnt/timur-higgsssc/hgammagamma/LOAnalysis/plots/gammagamma_run_01/index.html
```

You can also paste this path into a browser's file-open dialog:

```text
~/mnt/timur-higgsssc/hgammagamma/LOAnalysis/plots/gammagamma_run_01/index.html
```

The report is a static HTML page, so the PNG, SVG, CSV, and zip download links
should work through the mounted directory.

## 5. Copy Results Locally If Needed

If you want a local copy of one report directory:

```bash
cp -r ~/mnt/timur-higgsssc/hgammagamma/LOAnalysis/plots/gammagamma_run_01 \
  ~/gammagamma_run_01
```

Then open:

```bash
xdg-open ~/gammagamma_run_01/index.html
```

## 6. Unmount When Finished

When you are done viewing files:

```bash
fusermount3 -u ~/mnt/timur-higgsssc
```

On older Linux systems, the command may be:

```bash
fusermount -u ~/mnt/timur-higgsssc
```

## Practical Notes

- Run MG5, Herwig, and the analysis on Timur, not through the SSHFS mount.
- Use SSHFS mainly for viewing plots, downloading report assets, or checking
  small text files.
- If the mount becomes stale and you see `Transport endpoint is not connected`,
  unmount it with `fusermount3 -u ~/mnt/timur-higgsssc` and mount it again.
