#!/usr/bin/env bash
# Refuse to push credentials or bulk data.
#
# Run by hand before pushing:
#     ./sandbox/fastq_to_bam/check_no_secrets.sh
#
# Or install as a hook so it runs automatically (git does not track hooks,
# so each clone installs its own):
#     ln -sf ../../sandbox/fastq_to_bam/check_no_secrets.sh .git/hooks/pre-push
#     chmod +x .git/hooks/pre-push
#
# Exits non-zero on any finding, which aborts the push. It inspects tracked
# and staged files only -- gitignored files cannot reach the remote, and
# saying so is the whole point of the check.

set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || exit 1

FAIL=0
note() { printf '  %s\n' "$*"; }
fail() { FAIL=1; printf '\nFAIL: %s\n' "$1"; }

# Files git would actually publish: tracked at HEAD plus anything staged.
mapfile -t FILES < <( { git ls-files; git diff --cached --name-only; } \
                      | sort -u | while read -r f; do [ -f "$f" ] && echo "$f"; done )

echo "checking ${#FILES[@]} tracked/staged files"

# --- 1. credential files that must never be tracked ----------------------
for pat in '.env' '.env.*' '.netrc' '_netrc' 'id_rsa' 'id_ed25519' '*.pem' '*.p12' '*.pfx' '*.keystore'; do
    hits=$(printf '%s\n' "${FILES[@]}" | grep -E "(^|/)$(echo "$pat" | sed 's/\./\\./g; s/\*/[^\/]*/g')$" || true)
    # .env.example is a value-free template and is meant to be committed.
    hits=$(echo "$hits" | grep -v '\.env\.example$' || true)
    if [ -n "$hits" ]; then
        fail "credential file is tracked: $hits"
        note "untrack with: git rm --cached <file>   (the local copy is kept)"
    fi
done

# --- 2. secret-shaped content --------------------------------------------
# Assignments with a non-empty, non-placeholder value. Empty values and
# <angle-bracket> or ${VAR} placeholders are how the templates are written.
SECRET_RE='(pass|passwd|password|pwd|secret|credential|api[_-]?key|access[_-]?token|auth[_-]?token|client[_-]?secret)[[:space:]]*[:=][[:space:]]*[^[:space:]<$"'"'"']'
KEY_RE='BEGIN [A-Z ]*PRIVATE KEY|AKIA[0-9A-Z]{16}|ghp_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{20,}|xox[baprs]-[A-Za-z0-9-]{10,}|glpat-[A-Za-z0-9_-]{15,}'
for f in "${FILES[@]}"; do
    case "$f" in
        *.png|*.jpg|*.pdf|*.gz|*.bam|*.bai|*.mmi|*.lock|*/check_no_secrets.sh) continue ;;
    esac
    if hit=$(grep -nEiH "$SECRET_RE" "$f" 2>/dev/null | head -3); then
        [ -n "$hit" ] && { fail "possible credential assignment"; note "${hit%%:*}: see line $(echo "$hit" | cut -d: -f2)"; }
    fi
    if hit=$(grep -nEH "$KEY_RE" "$f" 2>/dev/null | head -3); then
        [ -n "$hit" ] && { fail "private key or provider token pattern in $f"; }
    fi
done

# --- 3. bulk data ---------------------------------------------------------
# A single 9 GiB FASTQ will be rejected by GitHub anyway, but only after the
# push has spent your upload bandwidth getting there.
for f in "${FILES[@]}"; do
    sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
    if [ "$sz" -gt 52428800 ]; then
        fail "$f is $((sz/1048576)) MiB -- GitHub warns over 50 MiB, rejects over 100 MiB"
        note "add it to .gitignore, then: git rm --cached '$f'"
    fi
done

# --- 4. did the ignore rules actually take effect? -----------------------
for probe in .env sandbox/big_data/probe.fastq; do
    git check-ignore -q "$probe" || { fail "$probe is NOT ignored"; note "check .gitignore"; }
done

if [ "$FAIL" -eq 0 ]; then
    echo "OK: no credential files, no secret-shaped content, no bulk data."
else
    echo
    echo "Push aborted. Note that removing a secret in a NEW commit does not"
    echo "help if it is already in an earlier one -- history has to be rewritten"
    echo "(git filter-repo) and the credential rotated regardless."
fi
exit "$FAIL"
