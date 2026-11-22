# Publishing this project as a GitHub repository

This archive is a complete **local repository source tree** and contains `.github/workflows/tests.yml`, a license, README, docs, reproducible source code, and tests.

GitHub repository creation is a separate account operation. The archive itself is **not** a remotely published GitHub URL.

## With the GitHub CLI

Install [GitHub CLI](https://cli.github.com/) and authenticate. From the extracted project directory:

```bash
python -m unittest discover -s tests -v

git init -b main
git add .
git commit -m "Initialize conic billiard research lab"
gh auth login
gh repo create conic-billiard-crypto-lab --private --source=. --remote=origin --push
```

If `git` reports an unset identity, set `git config user.name` and `git config user.email` to **your own** GitHub commit identity. To make the repository public later, change visibility in GitHub settings after reviewing all content.

Alternatively, create a **private, empty** GitHub repository via the website and push this local source tree to its URL with `git remote add origin <your-repository-url>` and `git push -u origin main`.

## Security

This is an experimental, variable-time finite-field arithmetic implementation, not a production DH library. The bytes in `results/` and test vectors are public mathematical examples, never real private keys.
