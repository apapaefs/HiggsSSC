# HiggsSSC

Higgs Boson production at the SuperConducting SuperCollider.

This repository is intended to act as the superproject, with the paper tracked
as a Git submodule in `paper/` using the GitHub repository
`git@github.com:apapaefs/HiggsAtSSC.git`.

Because the paper is already synced with Overleaf through GitHub, using the
GitHub repository as the submodule remote is cleaner than pointing the
submodule directly at Overleaf. The important caveat is that Overleaf and
GitHub do not sync automatically: changes made in Overleaf must be pushed from
Overleaf to GitHub before this repository can see them locally.

Also remember that a submodule is recorded in the superproject as a pointer to
a specific commit of the submodule repository. After updating the paper, you
should commit the paper changes inside `paper/` first, then commit the updated
submodule pointer in `HiggsSSC`.

## Layout

```text
HiggsSSC/
  src/
  ...
  paper/   # submodule -> HiggsAtSSC
```

## Setup From Scratch

To recreate the submodule setup from a fresh checkout:

```bash
git clone git@github.com:apapaefs/HiggsSSC.git
cd HiggsSSC

git submodule add git@github.com:apapaefs/HiggsAtSSC.git paper
git commit -m "Add HiggsAtSSC as paper submodule"
git push origin main
```

That is the standard `git submodule add <repository> <path>` workflow.

If `HiggsAtSSC` uses `main` and you want the submodule metadata to remember
that branch, use:

```bash
git clone git@github.com:apapaefs/HiggsSSC.git
cd HiggsSSC

git submodule add -b main git@github.com:apapaefs/HiggsAtSSC.git paper
git commit -m "Add HiggsAtSSC as paper submodule"
git push origin main
```

For this repository specifically, the intended setup is:

```bash
git clone git@github.com:apapaefs/HiggsSSC.git
cd HiggsSSC
git submodule add -b main git@github.com:apapaefs/HiggsAtSSC.git paper
git commit -m "Add paper submodule"
git push origin main
```

If `HiggsAtSSC` does not use `main`, replace `main` with its actual default
branch.

## Cloning A Repository With Submodules

Anyone cloning the code repository should do:

```bash
git clone --recurse-submodules git@github.com:apapaefs/HiggsSSC.git
```

Or, if they already cloned it:

```bash
git submodule update --init --recursive
```

That is the standard Git behavior for repositories containing submodules.

## Day-To-Day Workflow

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

## Quality Of Life

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
