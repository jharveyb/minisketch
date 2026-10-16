---
name: bitcoin-knowledge-base
description: Search the Bitcoin and Lightning development record (BIPs, BOLTs, bLIPs, GitHub issues/PRs/commits, mailing lists, IRC, Delving Bitcoin, Optech) through the bkb MCP tools. Use for spec or history questions - e.g. BIP-330 (Erlay), bitcoin-core/minisketch or bitcoin/bitcoin PRs and their discussion - instead of answering from memory.
---

<!-- Vendored from tnull/bitcoin-knowledge-base@3d917d6 SKILL.md (the rev
     flake.nix pins for bkb-mcp); only this frontmatter and the Setup section
     are local. -->

# Bitcoin Knowledge Base (BKB) -- MCP Skill

Search across the Bitcoin and Lightning development ecosystem: BIPs,
BOLTs, bLIPs, GitHub issues/PRs/commits, mailing lists, IRC logs,
Delving Bitcoin, and Optech newsletters.

## Setup

Nothing to install by hand: the repo's `.mcp.json` defines the `bkb` server
as `nix run .#bkb-mcp` (built from the pinned rev in `flake.nix`), talking to
https://bitcoinknowledge.dev. Approve the project MCP server when Claude Code
asks; a user-scope `bkb` entry in `~/.claude.json`, if any, is overridden by
the project one.

## Available Tools

- `bkb_search` -- Full-text search across all sources (filters: `source_type`, `source_repo`, `author`, `after`, `before`, `limit`). Supports wildcard queries: use `query: "*"` with at least one filter to retrieve all documents matching the filters (e.g., all commits in a repo within a date range).
- `bkb_get_document` -- Get full document by ID with content, cross-references, and concept tags
- `bkb_get_references` -- Find all documents referencing an entity (e.g. `BIP-340`, `bitcoin/bitcoin#12345`)
- `bkb_lookup_bip` -- BIP spec with all referencing discussions and PRs
- `bkb_lookup_bolt` -- BOLT spec with all referencing documents
- `bkb_lookup_blip` -- bLIP spec with all referencing documents
- `bkb_timeline` -- Chronological timeline of a concept across all sources
- `bkb_find_commit` -- Find commits/PRs matching a description
