#!/usr/bin/env bash
set -euo pipefail

echo "[devcontainer] Configuring Git for SSH commit signing"
GIT_NAME="${GIT_AUTHOR_NAME:-${GIT_COMMITTER_NAME:-VCode}}"
GIT_EMAIL="${GIT_AUTHOR_EMAIL:-${GIT_COMMITTER_EMAIL:-vscode@users.noreply.github.com}}"

git config --global user.name "$GIT_NAME"
git config --global user.email "$GIT_EMAIL"
git config --global gpg.format ssh
git config --global commit.gpgsign true

# Get the SSH key from the forwarded agent and set it as signing key.
SSH_KEY=$(ssh-add -L 2>/dev/null | head -1)
if [ -n "$SSH_KEY" ]; then
  git config --global user.signingkey "$SSH_KEY"
  echo "[devcontainer] SSH signing key configured: ${SSH_KEY:0:50}..."
else
  echo "[devcontainer] WARNING: No SSH key found in agent. Signing may not work."
fi
