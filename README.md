# HiggsSSC

What if? Higgs Boson production at the SuperConducting SuperCollider.

This repository contains the course/project files, with the paper tracked in
`paper/` as a Git submodule.

## Start Here: Clone The Repository

If you are a student, this is the starting point. Clone the repository with
submodules included:

```bash
git clone --recurse-submodules git@github.com:apapaefs/HiggsSSC.git
cd HiggsSSC
```

The `--recurse-submodules` part is important because the paper lives in a
separate repository inside `paper/`.

If you already cloned the repository without that option, run this from inside
`HiggsSSC`:

```bash
git submodule update --init --recursive
```

Students should not run `git submodule add`. The repository setup has already
been done once by the instructor/maintainer. For students, the job is only to
clone the existing repository, update it when needed, and work with the files.

## Layout

```text
HiggsSSC/
  paper/   # paper submodule -> HiggsAtSSC
```

## Running Herwig On Timur

To log in to Timur, open a terminal and run:

```bash
ssh -Y YourUsername@timur.kennesaw.edu
```

Replace `YourUsername` with your Timur/KSU username.

Load Herwig after each login:

```bash
module load herwig/stable
```

The shared testing directory is:

```text
/home/shared/HerwigTesting/
```

Before you start working, copy the shared testing area into your home
directory. Work from your own copy instead of editing files in the shared
location:

```bash
mkdir -p ~/HerwigWork
cp -r /home/shared/HerwigTesting ~/HerwigWork/
cd ~/HerwigWork/HerwigTesting
```

Useful Linux commands:

| Command | What it does | Example |
| --- | --- | --- |
| `ls` | Lists files in the current directory. | `ls` |
| `ls -l` | Lists files with more detail. | `ls -l` |
| `cd` | Changes into a directory. | `cd ~/HerwigWork/HerwigTesting` |
| `cp` | Copies a file or directory. | `cp LHC.in LHC_backup.in` |
| `mv` | Moves or renames a file. | `mv old_name.in new_name.in` |
| `less` | Opens a file for reading in the terminal. | `less LHC.log` |

For example, `less LHC.log` opens the file for reading. Press `q` to quit.

Before editing an input file, it is a good idea to make a backup copy:

```bash
cp LHC.in LHC_original.in
```

If you make a mistake and want to restore the original version from your
backup, copy it back:

```bash
cp LHC_original.in LHC.in
```

You can also restore a fresh original copy from the shared directory:

```bash
cp /home/shared/HerwigTesting/LHC.in ~/HerwigWork/HerwigTesting/LHC.in
```

## Editing Files With Vim On Timur

Vim is a terminal text editor. To edit a file, run:

```bash
vim LHC.in
```

Basic Vim workflow:

1. Press `i` to enter insert mode and start editing.
2. Use the arrow keys to move around.
3. Press `Esc` when you are done typing.
4. Type `:w` and press `Enter` to save.
5. Type `:q` and press `Enter` to quit.

Common Vim commands:

```text
i       start editing text
Esc     leave editing mode
:w      save the file
:q      quit Vim
:wq     save and quit
:q!     quit without saving
/word   search for "word"
n       go to the next search result
```

If Vim feels stuck, press `Esc` once or twice, then type one of the commands
above. For example, `:q!` quits without saving changes.

## Basic Herwig Workflow

Read an input file:

```bash
Herwig read LHC.in
```

Run Herwig using the generated run file:

```bash
Herwig run LHC.run -N1000
```

This example generates `1000` events.

After the run, inspect:

- `LHC.out` for the processes that were generated.
- `LHC.log` for run and event information.

## Changing The Process In `LHC.in`

In `LHC.in`, the active physics process is controlled in the section called
`Matrix Elements for hadron-hadron collisions`.

Lines that begin with `#` are comments, so Herwig ignores them. A line without
`#` is active. The default active process is:

```text
insert SubProcess:MatrixElements[0] MEqq2gZ2ff
```

That corresponds to Drell-Yan `Z/gamma` production. To change the process,
comment out the current active line by adding `#` at the front, then remove
the `#` from one new process line.

For example, to switch from Drell-Yan `Z/gamma` to Drell-Yan `W` production:

```text
# insert SubProcess:MatrixElements[0] MEqq2gZ2ff
insert SubProcess:MatrixElements[0] MEqq2W2ff
```

Change one process at a time at first. This makes it much easier to understand
which change caused the output you see in the `.out` and `.log` files.

Some good first processes to try are:

- `MEqq2W2ff` for Drell-Yan `W` production.
- `MEWJet` or `MEZJet` for a vector boson produced with a jet.
- `MEQCD2to2` for ordinary QCD two-to-two scattering.
- `MEHeavyQuark` for top-antitop production.
- `MEHiggs` for inclusive Higgs production.
- `MEHiggsJet` for Higgs plus jet production.
- `MEPP2HiggsVBF` for vector-boson-fusion Higgs production.
- `MEPP2ttbarH` for Higgs production with a top-antitop pair.

For the Higgs examples, read the comments already in `LHC.in` carefully. Some
Higgs-associated processes include extra suggested settings, such as setting
the jet `pT` cut to zero:

```text
set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV
```

After editing `LHC.in`, read the file again before running:

```bash
Herwig read LHC.in
Herwig run LHC.run -N1000
```

Then compare the new `LHC.out` and `LHC.log` files with the previous run.
Look for the reported cross section in the output, and keep a short record of
which process you ran and what cross section Herwig reported.

## Estimating The Expected Number Of Events

Once you have a cross section, estimate how many events an experiment would
expect with:

```text
N = L * sigma
```

Here `N` is the expected number of events, `sigma` is the cross section, and
`L` is the integrated luminosity. Integrated luminosity measures the total
amount of collision data collected: a larger `L` means more chances for a
process to occur.

For this exercise, assume the experiment collected:

```text
L = 100 fb^-1
```

If the cross section is in femtobarns (`fb`), then:

```text
N = 100 * sigma
```

If the cross section is in picobarns (`pb`), use:

```text
100 fb^-1 = 100000 pb^-1
N = 100000 * sigma
```

For example, if Herwig reports `sigma = 5 pb`, then:

```text
N = 100000 * 5 = 500000 events
```

## Student Exercise

1. Read and run `LEP.in`.
2. Check the corresponding `.out` file.
3. Read the corresponding `.log` file.
4. Modify `LHC.in` to run a different process.
5. For each process you try, document the process name, the number of events,
   the reported cross section, and the units.
6. Using `L = 100 fb^-1`, calculate the expected number of events for each
   process.

## Instructor/Maintainer Notes

The repository is intended to act as the superproject, with the paper tracked
as a Git submodule in `paper/` using the GitHub repository
`git@github.com:apapaefs/HiggsAtSSC.git`.

The submodule setup below is a one-time setup step for the
instructor/maintainer. It is not something students need to do.

Because the paper is already synced with Overleaf through GitHub, using the
GitHub repository as the submodule remote is cleaner than pointing the
submodule directly at Overleaf. The important caveat is that Overleaf and
GitHub do not sync automatically: changes made in Overleaf must be pushed from
Overleaf to GitHub before this repository can see them locally.

Also remember that a submodule is recorded in the superproject as a pointer to
a specific commit of the submodule repository. After updating the paper, commit
the paper changes inside `paper/` first, then commit the updated submodule
pointer in `HiggsSSC`.

### One-Time Repository Setup

For this repository specifically, the intended setup was:

```bash
git clone git@github.com:apapaefs/HiggsSSC.git
cd HiggsSSC
git submodule add -b main git@github.com:apapaefs/HiggsAtSSC.git paper
git commit -m "Add paper submodule"
git push origin main
```

That is the standard `git submodule add <repository> <path>` workflow. Again,
students do not need to run these commands.

If `HiggsAtSSC` does not use `main`, replace `main` with its actual default
branch.

### When You Edit The Paper Locally

```bash
cd HiggsSSC/paper
git pull origin main

# edit files

git add -A
git commit -m "Revise paper"
git push origin main

cd ..
git add paper
git commit -m "Update paper submodule pointer"
git push origin main
```

### When You Edit The Paper In Overleaf

1. In Overleaf, use the GitHub sync dialog to push Overleaf changes to GitHub.
   Overleaf does not sync automatically.
2. Then locally:

```bash
cd HiggsSSC/paper
git pull origin main
cd ..
git add paper
git commit -m "Advance paper submodule after Overleaf changes"
git push origin main
```

One important note: this works well because the paper repository is the
submodule. Overleaf allows Overleaf-backed projects to be used as submodules
inside another repository, while Overleaf projects themselves cannot contain
submodules.

### Quality Of Life

A small setting that helps with submodules is:

```bash
git config --global submodule.recurse true
```

That tells Git to recurse into submodules for many common commands by default.

## Operational Discipline

- Code changes live in `HiggsSSC`.
- Paper changes live in `paper/`.
- After changing `paper/`, always commit inside `paper/` first, then commit
  the updated submodule pointer in `HiggsSSC`.
