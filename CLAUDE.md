# Notes for Claude Code

## Commits and pull requests

- Never add a `Claude-Session:` trailer, a claude.ai session link, or any
  other session/conversation URL to commit messages, PR descriptions, or
  code. This repository is public. (`.claude/settings.json` also disables
  it via `attribution.sessionUrl: false`.)
- A `Co-Authored-By: Claude ... <noreply@anthropic.com>` trailer is fine.
