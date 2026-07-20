# Pushing this repository to GitHub

The repo is already git-initialized with an initial commit on branch `main`.
GitHub does not allow spaces in repository names, so "growth-at-risk energy"
becomes **growth-at-risk-energy** (the display name can contain anything).

## Option A — GitHub CLI (fastest)

```sh
cd growth-at-risk-energy
gh auth login                      # once, if not already logged in
gh repo create growth-at-risk-energy --private --source=. --push
```

## Option B — web UI + git

1. Go to https://github.com/new
2. Repository name: `growth-at-risk-energy`; choose Private (recommended
   while the paper is under review); do NOT initialize with a README
   (this repo already has one).
3. Then:

```sh
cd growth-at-risk-energy
git remote add origin git@github.com:<your-username>/growth-at-risk-energy.git
git push -u origin main
```

(Use `https://github.com/<your-username>/growth-at-risk-energy.git` instead
if you authenticate with a token rather than SSH keys.)

## Notes

- `output/`, `data/`, caches and PDFs are gitignored — only code and the
  input CSVs are tracked (33 files, ~1 MB).
- Keep the repo private until you have checked the redistribution terms of
  the Eurostat/TTF/Brent input data.
- Suggested description: "R and Python pipelines for a Growth-at-Risk
  analysis of energy price shocks on European industrial production".
