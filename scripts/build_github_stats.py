#!/usr/bin/env python3
"""Build src/assets/github.json: repository stats + 12-month star trend + 52-week commits.

Runs daily in CI (GITHUB_TOKEN) and locally (falls back to `gh auth token`).
"""
from __future__ import annotations
import argparse, datetime as dt, json, os, re, subprocess, sys, time, urllib.request, urllib.error

API = "https://api.github.com"


def token() -> str | None:
    t = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if t:
        return t
    try:
        return subprocess.run(["gh", "auth", "token"], capture_output=True, text=True, check=True).stdout.strip() or None
    except Exception:
        return None


def get(url: str, tok: str | None, accept: str = "application/vnd.github+json"):
    req = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "omicverse-org-home"})
    if tok:
        req.add_header("Authorization", f"Bearer {tok}")
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                body = r.read()
                return (json.loads(body) if body else None), dict(r.headers)
        except urllib.error.HTTPError as e:
            if e.code == 202 and attempt < 3:   # stats being computed
                time.sleep(2 + attempt * 2); continue
            raise
    return None, {}


def last_page(link: str | None) -> int | None:
    if not link:
        return None
    m = re.search(r'[?&]page=(\d+)>; rel="last"', link)
    return int(m.group(1)) if m else None


def month_key(d: dt.datetime) -> str:
    return f"{d.year:04d}-{d.month:02d}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="omicverse/omicverse")
    ap.add_argument("--out", default="src/assets/github.json")
    ap.add_argument("--months", type=int, default=12)
    args = ap.parse_args()
    tok = token()

    repo, _ = get(f"{API}/repos/{args.repo}", tok)
    out = {
        "generated_at": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "repo": args.repo,
        "stars": repo["stargazers_count"], "forks": repo["forks_count"], "open_issues": repo["open_issues_count"],
        "watchers": repo.get("subscribers_count"), "license": (repo.get("license") or {}).get("spdx_id"),
        "pushed_at": repo["pushed_at"], "created_at": repo["created_at"],
    }
    try:
        rel, _ = get(f"{API}/repos/{args.repo}/releases/latest", tok)
        out["latest_release"] = {"tag": rel["tag_name"], "published_at": rel["published_at"]}
    except Exception:
        out["latest_release"] = None
    try:
        _, h = get(f"{API}/repos/{args.repo}/contributors?per_page=1&anon=1", tok)
        out["contributors"] = last_page(h.get("Link")) or None
    except Exception:
        out["contributors"] = None
    try:
        part, _ = get(f"{API}/repos/{args.repo}/stats/participation", tok)
        out["commits_by_week"] = (part or {}).get("all") or []
    except Exception:
        out["commits_by_week"] = []

    # Star history: every stargazer with its timestamp (100 per page).
    counts: dict[str, int] = {}
    total = 0
    page = 1
    while True:
        rows, h = get(f"{API}/repos/{args.repo}/stargazers?per_page=100&page={page}", tok, "application/vnd.github.star+json")
        if not rows:
            break
        for r in rows:
            ts = r.get("starred_at")
            if ts:
                counts[ts[:7]] = counts.get(ts[:7], 0) + 1
                total += 1
        if "next" not in (h.get("Link") or ""):
            break
        page += 1
        if page > 400:
            break
    now = dt.datetime.now(dt.timezone.utc)
    months = []
    y, m = now.year, now.month
    for _ in range(args.months):
        months.append(f"{y:04d}-{m:02d}")
        m -= 1
        if m == 0:
            m, y = 12, y - 1
    months.reverse()
    before = sum(v for k, v in counts.items() if k < months[0])
    series, running = [], before
    for k in months:
        running += counts.get(k, 0)
        series.append({"month": k, "added": counts.get(k, 0), "total": running})
    out["stars_by_month"] = series
    out["stars_counted"] = total

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
        f.write("\n")
    print(f"wrote {args.out}: {out['stars']} stars, {len(series)} months, {len(out['commits_by_week'])} weeks, {out.get('contributors')} contributors")
    return 0


if __name__ == "__main__":
    sys.exit(main())
