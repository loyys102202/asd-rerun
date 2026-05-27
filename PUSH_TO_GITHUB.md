# How to push this repository to GitHub

This file documents the exact commands to publish this code to
`https://github.com/loyys102202/asd-rerun`, the URL declared in the
manuscript Resource availability statement.

You only need to do this once. After the first push, normal `git`
workflow applies for any subsequent edits.

## Prerequisites

- A GitHub account at `loyys102202`. If the account already exists,
  skip ahead. Otherwise sign up at <https://github.com/signup>.
- `git` installed locally. Verify with `git --version`.
- Either an SSH key registered with GitHub (recommended) or a
  Personal Access Token (PAT) you can paste when prompted for a
  password.

## Step 1 — Create the empty repository on github.com

1. Sign in to GitHub.
2. Top right corner → **`+` → New repository**.
3. Owner: `loyys102202`. Repository name: **`asd-rerun`** (exactly,
   case-sensitive; the manuscript references this string).
4. Description: paste this single line:
   > End-to-end reproducibility code for the iScience manuscript
   > "A monocyte-tracking blood expression score reveals platform-
   > specific dissociation between cellular composition and autism
   > diagnosis" (ISCIENCE-D-26-03472, Jing Wen).
5. Visibility: **Public** (the manuscript states the URL is publicly
   accessible as of the date of publication).
6. **Do NOT** tick "Initialize this repository with a README",
   "Add .gitignore", or "Choose a license". The local folder
   already contains all three. Initialising on the server side will
   force you to do an extra merge step.
7. Click **Create repository**.

GitHub then shows a page beginning *"Quick setup — if you've done this
kind of thing before"*. Leave this page open; you'll paste the URL it
shows into Step 2.

## Step 2 — Push the local folder

Open a terminal in the folder that contains this `PUSH_TO_GITHUB.md`
file (i.e., the `asd-rerun/` root). Then run:

```bash
# initialise local repo
git init -b main

# stage everything except items listed in .gitignore
git add .

# first commit
git commit -m "Initial release: code accompanying ISCIENCE-D-26-03472 revision"

# connect to your empty GitHub repo (HTTPS form; substitute SSH if preferred)
git remote add origin https://github.com/loyys102202/asd-rerun.git

# push and set upstream
git push -u origin main
```

When prompted for credentials:
- **Username**: `loyys102202`
- **Password**: paste a GitHub Personal Access Token (PAT). Account
  passwords are no longer accepted for HTTPS pushes. Create a PAT at
  <https://github.com/settings/tokens?type=beta> with the `repo`
  scope.

## Step 3 — Tag the v1.0.0 release (recommended)

A tagged release lets Zenodo capture a clean archival snapshot.

```bash
git tag -a v1.0.0 -m "v1.0.0 — submission code release for ISCIENCE-D-26-03472"
git push origin v1.0.0
```

On github.com, navigate to **Releases → Draft a new release →
Choose tag → v1.0.0 → Publish release**.

## Step 4 — Connect Zenodo to GitHub for automatic DOI minting

This is the path that produces a real Zenodo DOI without uploading the
zip manually.

1. Sign in to Zenodo with your GitHub account at
   <https://zenodo.org/account/settings/github/>.
2. In the list of repositories, locate `loyys102202/asd-rerun` and
   flip the toggle to **ON**.
3. Back on GitHub, **publish a new release** (e.g. `v1.0.0` from
   Step 3). Zenodo will detect it within a minute and mint a fresh
   DOI of the form `10.5281/zenodo.NNNNNNN`.
4. Copy that DOI. Replace every occurrence of `XXXXXXXX` in:
   - `README.md` (badge line and `## Citation` section)
   - `CITATION.cff` (no DOI field by default — add a `doi:` line)
   - `Manuscript_revised.docx` (Resource availability)
   - `Title_Page_revised.docx` (Data and code availability)
   - `Cover_Letter_revised.docx`

The `.zenodo.json` file in this repo will be picked up by Zenodo
automatically to populate the deposit metadata (title, authors,
keywords, license, related identifiers).

## Alternative: manual Zenodo upload (if you skip the GitHub link)

If you prefer not to link Zenodo to GitHub at all, you can upload
the standalone archive `asd-rerun-zenodo.zip` directly at
<https://zenodo.org/uploads/new>. Use the metadata in
`.zenodo.json` (or in `ZENODO_METADATA.md` in that zip) to fill
in the deposit form by hand. You will still get a DOI; you just
won't get automatic re-minting on each future GitHub release.

## Sanity check after pushing

After Step 2 completes, open
<https://github.com/loyys102202/asd-rerun> in a browser and confirm:

- The repository is public.
- The README renders with badges at the top.
- The file list shows `01_download_geo.py` through `12_make_results_table.py`,
  the `figures/` and `tables/` subdirectories, `requirements.txt`,
  `run_all.sh`, `LICENSE`, `CITATION.cff`, `.zenodo.json`, `.gitignore`.
- The license is detected as MIT (top right of the repo page).
- The `CITATION.cff` is detected (a "Cite this repository" button
  appears near the top of the repo).

If any of these are missing, re-check `.gitignore` did not accidentally
exclude a file you wanted committed: `git status --ignored` will show
ignored items.
