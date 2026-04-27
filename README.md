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

Two useful Linux commands are:

```bash
cd directory_name
```

Use `cd` to change into a directory. For example:

```bash
cd ~/HerwigWork/HerwigTesting
```

To view a file in the terminal, use:

```bash
less filename
```

For example, `less LHC.log` opens the file for reading. Press `q` to quit.

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

## Student Exercise

1. Read and run `LEP.in`.
2. Check the corresponding `.out` file.
3. Read the corresponding `.log` file.
4. Modify `LHC.in` to run a different process.

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
