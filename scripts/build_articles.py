#!/usr/bin/env python3
"""Build the Articles section of omicverse.org from articles/<slug>/.

Each article directory holds
    meta.yaml     title, authors, affiliations, abstract, keywords, links, dates
    article.md    the body in a small Markdown dialect (see docs/ARTICLES.md)
    figures/      figN.jpg / figSN.jpg referenced by :::figure blocks

    python3 scripts/build_articles.py          # HTML pages + web figures
    python3 scripts/build_articles.py --pdf    # also typeset the two-column PDF with the
                                               # OmicVerse LaTeX template (needs latexmk)

Outputs: src/articles.html, src/articles/<slug>.html, src/assets/articles/<slug>/*.
Run scripts/sync_chrome.py afterwards (the build does it for you unless --no-chrome).
"""
from __future__ import annotations

import argparse
import datetime as dt
import html
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse

import yaml
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sync_chrome import ASSET_VERSION, footer_html, header_html  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ARTICLES = ROOT / "articles"
SRC = ROOT / "src"
OUT_PAGES = SRC / "articles"
OUT_ASSETS = SRC / "assets" / "articles"
TEMPLATE_REPO = "https://github.com/omicverse/omicverse-latex-template.git"
VENUE = "OmicVerse Technical Reports"
SITE = "https://omicverse.org"
KATEX = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist"

MONTHS = ["January", "February", "March", "April", "May", "June", "July",
          "August", "September", "October", "November", "December"]


# ═══════════════════════════════════════════════════════════════ parsing
def parse_article(md: str) -> list[dict]:
    """Split article.md into blocks: heading, para, math, list, figure, table, refs."""
    lines = md.splitlines()
    blocks: list[dict] = []
    para: list[str] = []
    i = 0
    section_kind = "body"

    def flush():
        nonlocal para
        if para:
            blocks.append({"t": "p", "text": " ".join(x.strip() for x in para)})
            para = []

    while i < len(lines):
        line = lines[i].rstrip()
        s = line.strip()
        if not s:
            flush(); i += 1; continue
        if s.startswith(":::figure") or s.startswith(":::table"):
            flush()
            kind, ident, *rest = s[3:].split()
            body = []
            i += 1
            while lines[i].strip() != ":::":
                body.append(lines[i].rstrip()); i += 1
            i += 1
            if kind == "figure":
                blocks.append({"t": "figure", "id": ident, "file": rest[0], "title": body[0].strip(),
                               "legend": " ".join(b.strip() for b in body[1:] if b.strip()),
                               "supp": ident.startswith("S")})
            else:
                rows = [[c.strip() for c in r.strip().strip("|").split("|")] for r in body[1:] if r.strip().startswith("|")]
                rows = [r for r in rows if not all(set(c) <= set("-: ") for c in r)]
                blocks.append({"t": "table", "id": ident, "title": body[0].strip(), "head": rows[0], "rows": rows[1:]})
            continue
        m = re.match(r"^(#\*?|##)\s+(.*)$", s)
        if m:
            flush()
            mark, text = m.groups()
            if mark == "##":
                blocks.append({"t": "h2", "text": text})
            else:
                section_kind = ("back" if mark == "#*" else
                                "refs" if text == "References" else
                                "supp" if text.startswith("Supplementary") else "body")
                blocks.append({"t": "h1", "text": text, "kind": section_kind})
            i += 1
            continue
        if section_kind == "refs":
            m = re.match(r"^(\d+)\.\s+(.*)$", s)
            if m:
                flush()
                if not blocks or blocks[-1]["t"] != "refs":
                    blocks.append({"t": "refs", "items": []})
                blocks[-1]["items"].append((int(m.group(1)), m.group(2)))
                i += 1
                continue
        m = re.match(r"^\$\$\s*(.*?)\s*\\tag\{([^}]*)\}\s*\$\$$", s)
        if m:
            flush()
            blocks.append({"t": "math", "tex": m.group(1), "tag": m.group(2)})
            i += 1
            continue
        if re.match(r"^\d+\.\s+", s) and not line.startswith(" "):
            flush()
            items = []
            while i < len(lines) and (lines[i].strip() == "" or re.match(r"^(\d+\.\s+|\s{2,}-\s+)", lines[i])):
                cur = lines[i]
                if not cur.strip():
                    j = i + 1
                    while j < len(lines) and not lines[j].strip():
                        j += 1
                    if j < len(lines) and re.match(r"^(\d+\.\s+|\s{2,}-\s+)", lines[j]):
                        i = j; continue
                    break
                if re.match(r"^\d+\.\s+", cur):
                    items.append([re.sub(r"^\d+\.\s+", "", cur).strip(), []])
                else:
                    items[-1][1].append(re.sub(r"^\s+-\s+", "", cur).strip())
                i += 1
            blocks.append({"t": "list", "items": items})
            continue
        para.append(line)
        i += 1
    flush()
    return blocks


def expand_cites(spec: str) -> list[int]:
    out: list[int] = []
    for part in spec.split(","):
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


INLINE = re.compile(
    r"(?P<math>\$[^$\n]+?\$)"
    r"|(?P<code>`[^`]+`)"
    r"|(?P<cite>\[@[0-9,\-]+\])"
    r"|(?P<link>\[(?P<ltext>[^\]]+)\]\((?P<lurl>[^)\s]+)\))"
    r"|(?P<bold>\*\*(?P<btext>.+?)\*\*)"
    r"|(?P<ital>\*(?P<itext>[^*\s][^*]*?)\*)"
    r"|(?P<url>https?://[^\s<>()\[\]]+[^\s<>()\[\].,;:])"
    r"|(?P<fig>\bFig\. (?P<fnum>S?\d+)(?P<fsub>[a-z](?:-[a-z0-9]+)?)?)"
)


# ═══════════════════════════════════════════════════════════════ HTML
def h_inline(text: str) -> str:
    out, pos = [], 0
    for m in INLINE.finditer(text):
        out.append(html.escape(text[pos:m.start()], quote=False))
        pos = m.end()
        g = m.lastgroup
        if g == "math":
            out.append(html.escape(m.group(0), quote=False))
        elif g == "code":
            out.append(f"<code>{html.escape(m.group(0)[1:-1])}</code>")
        elif g == "cite":
            spec = m.group(0)[2:-1]
            parts = []
            for chunk in spec.split(","):
                if "-" in chunk:
                    a, b = chunk.split("-")
                    parts.append(f'<a href="#ref-{a}">{a}</a>–<a href="#ref-{b}">{b}</a>')
                else:
                    parts.append(f'<a href="#ref-{chunk}">{chunk}</a>')
            out.append(f'<sup class="cite">{",".join(parts)}</sup>')
        elif g == "link":
            out.append(f'<a href="{html.escape(m.group("lurl"))}" target="_blank" rel="noopener">{h_inline(m.group("ltext"))}</a>')
        elif g == "bold":
            out.append(f"<strong>{h_inline(m.group('btext'))}</strong>")
        elif g == "ital":
            out.append(f"<em>{h_inline(m.group('itext'))}</em>")
        elif g == "url":
            u = m.group(0)
            out.append(f'<a href="{html.escape(u)}" target="_blank" rel="noopener">{html.escape(u)}</a>')
        elif g == "fig":
            out.append(f'<a class="figref" href="#fig-{m.group("fnum")}">{html.escape(m.group(0))}</a>')
    out.append(html.escape(text[pos:], quote=False))
    return "".join(out)


def ref_html(text: str) -> str:
    body = h_inline(text)
    m = re.search(r"(10\.\d{4,9}/[^\s]+?)\.?$", text)
    if m:
        doi = m.group(1)
        body = body.replace(html.escape(doi, quote=False), f'<a href="https://doi.org/{html.escape(doi)}" target="_blank" rel="noopener">{html.escape(doi)}</a>')
    return body


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def fmt_date(d: dt.date, short=False) -> str:
    return f"{d.day:02d} {MONTHS[d.month - 1][:3] if short else MONTHS[d.month - 1]} {d.year}"


def author_names(meta, fmt="full") -> list[str]:
    return [a["name"] for a in meta["authors"]]


def cite_text(meta) -> str:
    def short(n):
        parts = n.split()
        return f"{parts[-1]}, {'. '.join(p[0] for p in parts[:-1])}."
    names = [short(a["name"]) for a in meta["authors"]]
    who = ", ".join(names[:-1]) + " & " + names[-1] if len(names) > 1 else names[0]
    where = f"https://doi.org/{meta['doi']}" if meta.get("doi") else f"{SITE}/articles/{meta['slug']}"
    return f"{who} {meta['title']}. {VENUE} ({meta['published'].year}). {where}"


def bibtex(meta) -> str:
    authors = " and ".join(f"{a['name'].split()[-1]}, {' '.join(a['name'].split()[:-1])}" for a in meta["authors"])
    key = f"{meta['authors'][0]['name'].split()[-1].lower()}{meta['published'].year}{slugify(meta['short_title']).replace('-', '')}"
    return (f"@techreport{{{key},\n  title       = {{{meta['title']}}},\n  author      = {{{authors}}},\n"
            f"  institution = {{OmicVerse}},\n  type        = {{{meta['type']}}},\n  year        = {{{meta['published'].year}}},\n"
            + (f"  doi         = {{{meta['doi']}}},\n" if meta.get("doi") else "")
            + f"  url         = {{{SITE}/articles/{meta['slug']}}}\n}}")


def render_blocks_html(blocks, meta, fig_url) -> tuple[str, list[tuple[str, str]], list[dict]]:
    """Return (html, toc entries, figures)."""
    out, toc, figs = [], [], []
    open_section = False
    for b in blocks:
        t = b["t"]
        if t == "h1":
            if open_section:
                out.append("</section>")
            sid = "sec-" + slugify(b["text"])
            if b["kind"] == "refs":
                sid = "references"
            elif b["kind"] == "supp":
                sid = "supplementary"
            cls = f"paper-section is-{b['kind']}"
            out.append(f'<section class="{cls}" id="{sid}" aria-labelledby="{sid}-h">')
            out.append(f'<h2 id="{sid}-h">{h_inline(b["text"])}</h2>')
            if b["kind"] != "back":
                toc.append((sid, b["text"]))
            open_section = True
        elif t == "h2":
            out.append(f'<h3 id="{slugify(b["text"])}">{h_inline(b["text"])}</h3>')
        elif t == "p":
            out.append(f"<p>{h_inline(b['text'])}</p>")
        elif t == "math":
            out.append(f'<div class="paper-eq" id="eq-{b["tag"]}">$${html.escape(b["tex"], quote=False)} \\tag{{{b["tag"]}}}$$</div>')
        elif t == "list":
            lis = []
            for text, kids in b["items"]:
                sub = "".join(f"<li>{h_inline(k)}</li>" for k in kids)
                lis.append(f"<li>{h_inline(text)}{f'<ul>{sub}</ul>' if sub else ''}</li>")
            out.append(f"<ol>{''.join(lis)}</ol>")
        elif t == "figure":
            full, thumb = fig_url(b["file"])
            label = f"Fig. {b['id']}"
            figs.append({"id": b["id"], "title": b["title"], "thumb": thumb, "full": full})
            out.append(
                f'<figure class="paper-fig" id="fig-{b["id"]}">'
                f'<figcaption class="fig-title"><b>{label}:</b> {h_inline(b["title"])}</figcaption>'
                f'<a class="fig-img" href="{full}" target="_blank" rel="noopener"><img src="{full}" alt="{html.escape(label + ": " + b["title"])}" loading="lazy" /></a>'
                f'<div class="fig-legend">{h_inline(b["legend"])}</div>'
                f'<a class="fig-full" href="{full}" target="_blank" rel="noopener">Full size image</a>'
                f"</figure>")
        elif t == "table":
            head = "".join(f"<th>{h_inline(c)}</th>" for c in b["head"])
            rows = []
            for r in b["rows"]:
                if all(not c for c in r[1:]):
                    rows.append(f'<tr class="group"><th colspan="{len(b["head"])}">{h_inline(r[0])}</th></tr>')
                else:
                    rows.append("<tr>" + "".join(f"<td>{h_inline(c)}</td>" for c in r) + "</tr>")
            out.append(f'<figure class="paper-table" id="table-{b["id"]}"><figcaption class="fig-title"><b>Table {b["id"]}:</b> {h_inline(b["title"])}</figcaption>'
                       f'<div class="table-scroll"><table><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div></figure>')
        elif t == "refs":
            lis = "".join(f'<li id="ref-{n}"><span class="ref-n">{n}.</span><span class="ref-text">{ref_html(txt)}</span></li>' for n, txt in b["items"])
            out.append(f'<ol class="paper-refs">{lis}</ol>')
    if open_section:
        out.append("</section>")
    return "\n".join(out), toc, figs


ICON_MAIL = '<svg viewBox="0 0 24 24" width="15" height="15" aria-hidden="true"><path fill="none" stroke="currentColor" stroke-width="1.8" d="M3 6h18v12H3z M3 6l9 7 9-7"/></svg>'
ICON_DL = '<svg viewBox="0 0 24 24" width="18" height="18" aria-hidden="true"><path fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" d="M12 4v11m0 0l-4.5-4.5M12 15l4.5-4.5M5 19h14"/></svg>'


def article_page(meta, body_html, toc, figs, pdf_info, n_refs) -> str:
    slug, title = meta["slug"], meta["title"]
    pub = meta["published"]
    authors = []
    for a in meta["authors"]:
        sup = ",".join(str(x) for x in a["affiliations"])
        mail = (f' <a class="au-mail" href="mailto:{html.escape(a["email"])}" aria-label="Email {html.escape(a["name"])}">{ICON_MAIL}</a>'
                if a.get("email") else "")
        authors.append(f'<li><span class="au-name">{html.escape(a["name"])}</span><sup>{sup}</sup>{mail}</li>')
    affs = "".join(f'<li><sup>{k}</sup> {html.escape(v)}</li>' for k, v in meta["affiliations"].items())
    corr = ", ".join(f'<a href="mailto:{html.escape(a["email"])}">{html.escape(a["name"])}</a>' for a in meta["authors"] if a.get("email"))
    kw = "".join(f"<li>{html.escape(k)}</li>" for k in meta.get("keywords", []))
    toc_items = [("abstract", "Abstract")] + toc + [("cite", "Cite this article")]
    toc_html = "".join(f'<li><a href="#{sid}">{html.escape(label)}</a></li>' for sid, label in toc_items)
    fig_html = "".join(f'<li><a href="#fig-{f["id"]}"><img src="{f["thumb"]}" alt="" loading="lazy" /><span>Fig. {f["id"]}</span></a></li>' for f in figs)
    pdf_href = f"../assets/articles/{slug}/{slug}.pdf"
    pdf_note = f"Two-column PDF · {pdf_info['pages']} pages · {pdf_info['size']}" if pdf_info else "Two-column PDF"
    links = meta.get("links", {})
    side_links = []
    if links.get("code"):
        side_links.append(f'<li><span>Code</span><a href="{html.escape(links["code"])}" target="_blank" rel="noopener">{html.escape(links["code"].replace("https://github.com/", "github.com/"))}</a></li>')
    if links.get("api"):
        side_links.append(f'<li><span>In OmicVerse</span><code>{html.escape(links["api"])}</code></li>')
    if links.get("tutorial"):
        host = urlparse(links["tutorial"]).netloc or links["tutorial"]      # show where the link goes
        side_links.append(f'<li><span>Tutorial</span><a href="{html.escape(links["tutorial"])}" target="_blank" rel="noopener">{html.escape(host)}</a></li>')
    doi = meta.get("doi")
    pre = meta.get("preprint") or {}
    doi_url = f"https://doi.org/{doi}" if doi else ""
    doi_line, cite_meta = "", [("Published", fmt_date(pub))]
    if pre:
        cite_meta.append(("Preprint", f"{html.escape(pre['server'])}, posted {fmt_date(pre['posted'])}"))
    if doi:
        posted = f'<span class="sep" aria-hidden="true">·</span>{html.escape(pre["server"])} preprint, posted {fmt_date(pre["posted"])}' if pre else ""
        doi_line = f'<p class="paper-doi">DOI <a href="{doi_url}" target="_blank" rel="noopener">{doi_url}</a>{posted}</p>'
        cite_meta.append(("DOI", f'<a href="{doi_url}" target="_blank" rel="noopener">{doi_url}</a>'))
        side_links.insert(0, f'<li><span>DOI{" · " + html.escape(pre["server"]) if pre else ""}</span><a href="{doi_url}" target="_blank" rel="noopener">{html.escape(doi)}</a></li>')
    cite_meta_html = "".join(f"<div><dt>{k}</dt><dd>{v}</dd></div>" for k, v in cite_meta)
    desc = html.escape(re.sub(r"\s+", " ", meta["abstract"])[:300], quote=True)
    citation_meta = "\n".join(
        [f'  <meta name="citation_title" content="{html.escape(title, quote=True)}" />']
        + [f'  <meta name="citation_author" content="{html.escape(a["name"], quote=True)}" />' for a in meta["authors"]]
        + [f'  <meta name="citation_publication_date" content="{pub:%Y/%m/%d}" />',
           f'  <meta name="citation_journal_title" content="{VENUE}" />',
           f'  <meta name="citation_pdf_url" content="{SITE}/assets/articles/{slug}/{slug}.pdf" />']
        + ([f'  <meta name="citation_doi" content="{html.escape(doi, quote=True)}" />'] if doi else []))
    bib = html.escape(bibtex(meta))
    cite = html.escape(cite_text(meta))
    abstract_html = h_inline(re.sub(r"\s+", " ", meta["abstract"]).strip())
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(title)} — OmicVerse Articles</title>
  <meta name="description" content="{desc}" />
  <meta property="og:title" content="{html.escape(title, quote=True)}" />
  <meta property="og:description" content="{desc}" />
{citation_meta}
  <link rel="icon" type="image/png" sizes="32x32" href="../assets/brand/favicon-32.png" />
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500&family=Source+Serif+4:opsz,wght@8..60,400;8..60,600;8..60,700&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="{KATEX}/katex.min.css" />
  <link rel="stylesheet" href="../styles.css?v={ASSET_VERSION}" />
  <script src="../script.js?v={ASSET_VERSION}" defer></script>
  <script src="{KATEX}/katex.min.js" defer></script>
  <script src="{KATEX}/contrib/auto-render.min.js" defer></script>
  <script src="../articles.js?v={ASSET_VERSION}" defer></script>
</head>
<body class="paper-page">
  <div class="bg-grid" aria-hidden="true"></div>
{header_html("articles", "../")}

  <main class="paper page-width wide">
    <nav class="paper-crumbs" aria-label="Breadcrumb"><a href="../index.html">OmicVerse</a><span aria-hidden="true">›</span><a href="../articles.html">Articles</a><span aria-hidden="true">›</span><span>{html.escape(meta["type"])}</span></nav>
    <div class="paper-grid">
      <article class="paper-main">
        <p class="paper-meta"><span>{html.escape(meta["type"])}</span><span>Published: <time datetime="{pub.isoformat()}">{fmt_date(pub)}</time></span><span>Version {meta.get("version", 1)}</span></p>
        <h1 class="paper-title">{html.escape(title)}</h1>
        <ul class="paper-authors">{"".join(authors)}</ul>
        <details class="paper-affil"><summary>Affiliations</summary><ol>{affs}</ol></details>
        <p class="paper-journal"><i>{VENUE}</i> ({pub.year})<span class="sep" aria-hidden="true">|</span><a href="#cite">Cite this article</a></p>
        {doi_line}
        <a class="paper-pdf-inline btn btn-primary" href="{pdf_href}" download>{ICON_DL}Download PDF</a>

        <section class="paper-section" id="abstract" aria-labelledby="abstract-h">
          <h2 id="abstract-h">Abstract</h2>
          <p class="paper-abstract">{abstract_html}</p>
          <ul class="paper-keywords" aria-label="Keywords">{kw}</ul>
        </section>

{body_html}

        <section class="paper-section" id="author-information" aria-labelledby="author-information-h">
          <h2 id="author-information-h">Author information</h2>
          <h3>Affiliations</h3>
          <ol class="paper-affil-list">{affs}</ol>
          <h3>Corresponding authors</h3>
          <p>Correspondence to {corr}.</p>
        </section>

        <section class="paper-section" id="cite" aria-labelledby="cite-h">
          <h2 id="cite-h">Cite this article</h2>
          <p class="cite-text" id="cite-text">{cite}</p>
          <dl class="cite-meta">{cite_meta_html}</dl>
          <div class="cite-actions"><button class="btn btn-ghost" type="button" data-copy-target="#cite-text">Copy citation</button><button class="btn btn-ghost" type="button" data-copy-target="#cite-bib">Copy BibTeX</button></div>
          <pre class="cite-bib" id="cite-bib">{bib}</pre>
        </section>
      </article>

      <aside class="paper-side" aria-label="Article tools">
        <div class="side-card side-pdf">
          <a class="btn-pdf" href="{pdf_href}" download>Download PDF {ICON_DL}</a>
          <p>{pdf_note} · <a href="{pdf_href}" target="_blank" rel="noopener">open in browser</a></p>
        </div>
        <div class="side-tabs" data-side-tabs>
          <div class="side-tablist" role="tablist">
            <button type="button" role="tab" aria-selected="true" data-tab="sections">Sections</button>
            <button type="button" role="tab" aria-selected="false" data-tab="figures">Figures</button>
            <button type="button" role="tab" aria-selected="false" data-tab="references">References</button>
          </div>
          <ol class="side-panel side-toc" data-panel="sections">{toc_html}</ol>
          <ul class="side-panel side-figs" data-panel="figures" hidden>{fig_html}</ul>
          <div class="side-panel" data-panel="references" hidden><p>{n_refs} references, numbered in order of citation.</p><a href="#references">Go to references</a></div>
        </div>
        <ul class="side-card side-links">{"".join(side_links)}</ul>
      </aside>
    </div>
  </main>

{footer_html("../")}
</body>
</html>
"""


def index_page(entries) -> str:
    types = sorted({e["meta"]["type"] for e in entries})
    years = sorted({e["meta"]["published"].year for e in entries}, reverse=True)
    rows = []
    for e in entries:
        m = e["meta"]
        names = [a["name"] for a in m["authors"]]
        who = ", ".join(names) if len(names) <= 3 else f"{names[0]}, {names[1]} … {names[-1]}"
        abstract = re.sub(r"\s+", " ", m["abstract"]).strip()
        abstract = abstract if len(abstract) < 260 else abstract[:abstract.rfind(" ", 0, 250)] + " …"
        thumb = f'<a class="ar-thumb" href="articles/{m["slug"]}.html" tabindex="-1" aria-hidden="true"><img src="{e["thumb"]}" alt="" loading="lazy" /></a>' if e.get("thumb") else ""
        rows.append(f"""      <li class="ar-row" data-type="{html.escape(m['type'])}" data-year="{m['published'].year}">
        <div class="ar-meta"><span class="ar-type">{html.escape(m['type'])}</span><time datetime="{m['published'].isoformat()}">{fmt_date(m['published'], short=True)}</time></div>
        <div class="ar-body">
          <h2 class="ar-title"><a href="articles/{m['slug']}.html">{html.escape(m['title'])}</a></h2>
          <p class="ar-authors">{html.escape(who)}</p>
          <p class="ar-abstract">{html.escape(abstract)}</p>
        </div>
        {thumb}
      </li>""")
    opt_types = "".join(f'<option value="{html.escape(t)}">{html.escape(t)}</option>' for t in types)
    opt_years = "".join(f'<option value="{y}">{y}</option>' for y in years)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Articles — OmicVerse</title>
  <meta name="description" content="Technical reports for OmicVerse methods and rebuilt R packages: full methods, benchmarks and figures, typeset as papers." />
  <link rel="icon" type="image/png" sizes="32x32" href="assets/brand/favicon-32.png" />
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500&family=Source+Serif+4:opsz,wght@8..60,400;8..60,600;8..60,700&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="styles.css?v={ASSET_VERSION}" />
  <script src="script.js?v={ASSET_VERSION}" defer></script>
  <script src="articles.js?v={ASSET_VERSION}" defer></script>
</head>
<body class="paper-page">
  <div class="bg-grid" aria-hidden="true"></div>
{header_html("articles")}

  <main class="articles-index page-width wide">
    <nav class="paper-crumbs" aria-label="Breadcrumb"><a href="index.html">OmicVerse</a><span aria-hidden="true">›</span><span>Articles</span></nav>
    <h1 class="articles-h1">Articles</h1>
    <p class="articles-lede">Technical reports for OmicVerse methods and rebuilt R packages, with full methods, benchmarks and figures, typeset as papers with the <a href="https://github.com/omicverse/omicverse-latex-template" target="_blank" rel="noopener">OmicVerse LaTeX template</a>. Reports are released by their authors and are not peer reviewed.</p>
    <form class="articles-filter" data-articles-filter>
      <label><span>Article Type</span><select name="type"><option value="">All</option>{opt_types}</select></label>
      <label><span>Year</span><select name="year"><option value="">All</option>{opt_years}</select></label>
    </form>
    <ol class="ar-list" data-articles-list>
{chr(10).join(rows)}
    </ol>
    <p class="ar-empty" data-articles-empty hidden>No articles match these filters.</p>
    <aside class="articles-submit">
      <h2>Publish a technical report</h2>
      <p>Write the report with the <a href="https://github.com/omicverse/omicverse-latex-template" target="_blank" rel="noopener">OmicVerse LaTeX template</a> in the two-column layout, then open a pull request that adds <code>articles/&lt;name&gt;/</code> to <a href="https://github.com/omicverse/omicverse-org-home" target="_blank" rel="noopener">omicverse-org-home</a>. The page and the PDF are built from that folder.</p>
    </aside>
  </main>

{footer_html()}
</body>
</html>
"""


# ═══════════════════════════════════════════════════════════════ LaTeX
TEX_UNICODE = {
    "τ": r"$\tau$", "²": r"\textsuperscript{2}", "±": r"$\pm$", "≥": r"$\geq$", "≤": r"$\leq$", "×": r"$\times$",
    "→": r"$\rightarrow$", "’": "'", "‘": "`", "“": "``", "”": "''", "–": "--", "—": "---", "Å": r"\AA{}",
    "ϵ": r"$\epsilon$", "α": r"$\alpha$", "β": r"$\beta$", "γ": r"$\gamma$", "μ": r"$\mu$", "…": r"\ldots{}",
    "\u00a0": "~", "✉": "", "⁻": r"\textsuperscript{$-$}",
    **{c: r"\textsuperscript{%d}" % i for i, c in enumerate("⁰¹²³⁴⁵⁶⁷⁸⁹")},
}
TEX_SPECIAL = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
               "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}", "<": r"\textless{}", ">": r"\textgreater{}"}


def t_escape(s: str) -> str:
    s = "".join(TEX_SPECIAL.get(c, c) for c in s)
    return "".join(TEX_UNICODE.get(c, c) for c in s)


def t_url(u: str) -> str:
    return u.replace("\\", "/").replace("%", r"\%").replace("#", r"\#")


def t_inline(text: str) -> str:
    out, pos = [], 0
    for m in INLINE.finditer(text):
        out.append(t_escape(text[pos:m.start()]))
        pos = m.end()
        g = m.lastgroup
        if g == "math":
            out.append(m.group(0))
        elif g == "code":
            out.append(r"\texttt{" + t_escape(m.group(0)[1:-1]) + "}")
        elif g == "cite":
            out.append(r"\cite{" + ",".join(f"r{n}" for n in expand_cites(m.group(0)[2:-1])) + "}")
        elif g == "link":
            out.append(r"\href{" + t_url(m.group("lurl")) + "}{" + t_inline(m.group("ltext")) + "}")
        elif g == "bold":
            out.append(r"\textbf{" + t_inline(m.group("btext")) + "}")
        elif g == "ital":
            out.append(r"\emph{" + t_inline(m.group("itext")) + "}")
        elif g == "url":
            out.append(r"\url{" + t_url(m.group(0)) + "}")
        elif g == "fig":
            out.append(r"\figlink{fig:" + m.group("fnum") + "}{" + t_escape(m.group(0)).replace(". ", ".~", 1) + "}")
    out.append(t_escape(text[pos:]))
    return "".join(out)


def render_blocks_tex(blocks, meta, fig_dir: str) -> str:
    out = []
    supp_started = False
    for b in blocks:
        t = b["t"]
        if t == "h1":
            if b["kind"] == "refs":
                continue
            if b["kind"] == "supp":
                supp_started = True
                out.append("\\onecolumn\n\\section*{" + t_escape(b["text"]) + "}\n\\setcounter{figure}{0}\n\\renewcommand{\\thefigure}{S\\arabic{figure}}")
                continue
            star = "*" if b["kind"] == "back" else ""
            out.append(f"\\section{star}{{{t_inline(b['text'])}}}")
        elif t == "h2":
            out.append(f"\\subsection{{{t_inline(b['text'])}}}")
        elif t == "p":
            out.append(t_inline(b["text"]) + "\n")
        elif t == "math":
            out.append("\\begin{equation*}\n" + b["tex"] + " \\tag{" + b["tag"] + "}\n\\end{equation*}")
        elif t == "list":
            items = []
            for text, kids in b["items"]:
                sub = ("\n  \\begin{itemize}\n" + "\n".join(f"    \\item {t_inline(k)}" for k in kids) + "\n  \\end{itemize}") if kids else ""
                items.append(f"  \\item {t_inline(text)}{sub}")
            out.append("\\begin{enumerate}\n" + "\n".join(items) + "\n\\end{enumerate}")
        elif t == "figure":
            env, place = ("figure", "!htbp") if b["supp"] else ("figure*", "tp")
            out.append(f"\\begin{{{env}}}[{place}]\n\\centering\n"
                       f"\\includegraphics[width=\\textwidth,height=0.66\\textheight,keepaspectratio]{{{fig_dir}/{b['file']}}}\n"
                       f"\\pdcaption{{{t_inline(b['title'])}}}{{{t_inline(b['legend'])}}}\n"
                       f"\\label{{fig:{b['id']}}}\n\\end{{{env}}}")
        elif t == "table":
            n = len(b["head"])
            rows = []
            for r in b["rows"]:
                if all(not c for c in r[1:]):
                    rows.append(f"\\multicolumn{{{n}}}{{@{{}}l}}{{\\textbf{{{t_inline(r[0])}}}}}\\\\")
                else:
                    rows.append(" & ".join(t_inline(c) for c in r) + "\\\\")
            head = "\\tablehead " + " & ".join(f"\\textbf{{{t_inline(c)}}}" for c in b["head"]) + "\\\\"
            out.append("\\begin{table*}[tp]\n\\pdcaption{" + t_inline(b["title"]) + "}{}\n\\label{table:" + b["id"] + "}\n"
                       "\\footnotesize\\renewcommand{\\arraystretch}{1.15}\n"
                       "\\begin{tabular}{@{}p{0.28\\textwidth}p{0.24\\textwidth}p{0.42\\textwidth}@{}}\n\\toprule\n"
                       + head + "\n\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}")
        elif t == "refs":
            items = "\n".join(f"\\bibitem{{r{n}}} {t_inline(txt)}" for n, txt in b["items"])
            out.append("\\begin{thebibliography}{" + str(len(b["items"])) + "}\n" + items + "\n\\end{thebibliography}")
    body = "\n\n".join(out)
    # the bibliography must precede the supplementary figures
    if supp_started:
        refs = re.search(r"\\begin\{thebibliography\}.*?\\end\{thebibliography\}", body, re.S).group(0)
        body = body.replace(refs, "")
        body = body.replace("\\onecolumn\n\\section*{Supplementary", refs + "\n\n\\onecolumn\n\\section*{Supplementary", 1)
    return body


def tex_document(meta, body: str, variant: str = "site") -> str:
    pub = meta["published"]
    authors = []
    for i, a in enumerate(meta["authors"]):
        marks = ",".join(str(x) for x in a["affiliations"]) + (",*" if a.get("email") else "")
        sep = "" if i == len(meta["authors"]) - 1 else ","
        prefix = "and " if i == len(meta["authors"]) - 1 and len(meta["authors"]) > 1 else ""
        authors.append(f"{prefix}{t_escape(a['name'])}{sep}\\affmark{{{marks}}}")
    affs = "\n".join(f"\\affiliation{{{k}}}{{{t_escape(v)}.}}" for k, v in meta["affiliations"].items())
    corr = ", ".join(f"\\href{{mailto:{a['email']}}}{{\\textcolor{{pdprimary}}{{{t_escape(a['email'])}}}}}" for a in meta["authors"] if a.get("email"))
    summary = t_inline(re.sub(r"\s+", " ", meta["abstract"]).strip())
    doi_tex = f"\\articledoi{{{meta['doi']}}}\n" if meta.get("doi") else ""
    if variant == "preprint":
        # Copy for preprint servers: single column, and no logo, corner graphic or venue name
        # (Preprints.org asks for journal branding to be removed). CC BY 4.0, as preprints are posted.
        head = f"""\\documentclass[brand=omicverse]{{pdarticle}}

\\makeatletter
\\def\\pd@logo{{}}\\def\\pd@logosmall{{}}
\\makeatother
\\brandgraphic{{}}
\\journalname{{Preprint}}
\\articletype{{{t_escape(meta.get('preprint_type', meta['type']))}}}
\\journalinfo{{Preprint, {MONTHS[pub.month - 1]} {pub.year}}}
\\copyrightholder{{{pub.year} The Authors}}
\\licenseline{{This work is licensed under a Creative Commons Attribution 4.0 International License (CC BY 4.0).}}
\\renewcommand{{\\topfraction}}{{0.9}}
\\renewcommand{{\\floatpagefraction}}{{0.8}}
"""
        doi_tex = ""
    else:
        head = f"""\\documentclass[brand=omicverse, twocolumn]{{pdarticle}}

\\journalname{{OmicVerse}}
\\articletype{{{t_escape(meta['type'])}}}
\\journalinfo{{{VENUE}, {MONTHS[pub.month - 1]} {pub.year}}}
\\copyrightholder{{{pub.year} The Authors}}
\\licenseline{{Published at \\href{{{SITE}/articles/{meta['slug']}}}{{omicverse.org/articles/{meta['slug']}}}, version {meta.get('version', 1)}.}}
"""
    return f"""%% Generated by scripts/build_articles.py from articles/{meta['slug']}/ — edit article.md / meta.yaml, not this file.
{head}
{doi_tex}\\title{{{t_escape(meta['title'])}}}
\\author{{{' '.join(authors)}}}
{affs}
% full-width figures and tables: allow tall top floats, and start float-only pages at the top
\\renewcommand{{\\dbltopfraction}}{{0.9}}
\\renewcommand{{\\dblfloatpagefraction}}{{0.75}}
\\renewcommand{{\\textfraction}}{{0.08}}
\\setcounter{{dbltopnumber}}{{2}}
\\makeatletter
\\setlength{{\\@fptop}}{{0pt}}
\\setlength{{\\@dblfptop}}{{0pt}}
\\makeatother
\\makeatletter
\\def\\pd@corresp{{*Correspondence: {corr}\\par}}
\\makeatother

\\summary{{{summary}}}
\\keywords{{{t_escape(', '.join(meta.get('keywords', [])))}}}

\\begin{{document}}
\\maketitle

{body}

\\end{{document}}
"""


# ═══════════════════════════════════════════════════════════════ build
def ensure_template(path: Path) -> Path:
    if not (path / "pdarticle.cls").exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            shutil.rmtree(path)
        subprocess.run(["git", "clone", "--depth", "1", "-q", TEMPLATE_REPO, str(path)], check=True)
    return path


def build_pdf(meta, blocks, art_dir: Path, template: Path) -> dict:
    slug = meta["slug"]
    bdir = art_dir / "build"
    bdir.mkdir(exist_ok=True)
    for f in ("pdarticle.cls", "numbered.bst"):
        shutil.copy2(template / f, bdir / f)
    if (bdir / "logos").exists():
        shutil.rmtree(bdir / "logos")
    shutil.copytree(template / "logos", bdir / "logos")
    tex = tex_document(meta, render_blocks_tex(blocks, meta, "../figures"))
    (bdir / f"{slug}.tex").write_text(tex)
    r = subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", f"{slug}.tex"],
                       cwd=bdir, capture_output=True, text=True)
    if r.returncode != 0:
        log = (bdir / f"{slug}.log").read_text(errors="replace")
        errs = [l for l in log.splitlines() if l.startswith("!") or "Error" in l][:20]
        raise SystemExit(f"[articles] LaTeX failed for {slug}:\n" + "\n".join(errs))
    dest = OUT_ASSETS / slug
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(bdir / f"{slug}.pdf", dest / f"{slug}.pdf")
    pages = int(re.search(r"Pages:\s+(\d+)", subprocess.run(["pdfinfo", str(dest / f"{slug}.pdf")], capture_output=True, text=True).stdout).group(1))
    size = (dest / f"{slug}.pdf").stat().st_size
    info = {"pages": pages, "size": f"{size / 1e6:.1f} MB"}
    (art_dir / "pdf.json").write_text(json.dumps(info) + "\n")
    warn = [l for l in (bdir / f"{slug}.log").read_text(errors="replace").splitlines()
            if "Overfull \\hbox" in l or "undefined" in l.lower() or "Missing character" in l]
    print(f"[articles] {slug}.pdf: {pages} pages, {info['size']}" + (f", {len(warn)} warnings" if warn else ""))
    for w in warn[:12]:
        print("   ", w)
    return info


def build_preprint(meta, blocks, art_dir: Path, template: Path) -> Path:
    """Single-column, unbranded copy for a preprint server: articles/<slug>/preprint/
    holds <slug>-preprint.pdf and <slug>-preprint-latex.zip (every file needed to rebuild it)."""
    slug = meta["slug"]
    name = f"{slug}-preprint"
    pdir = art_dir / "preprint"
    if pdir.exists():
        shutil.rmtree(pdir)
    (pdir / "figures").mkdir(parents=True)
    for f in ("pdarticle.cls", "numbered.bst"):
        shutil.copy2(template / f, pdir / f)
    used = [b["file"] for b in blocks if b["t"] == "figure"]
    for f in used:
        shutil.copy2(art_dir / "figures" / f, pdir / "figures" / f)
    (pdir / f"{name}.tex").write_text(tex_document(meta, render_blocks_tex(blocks, meta, "figures"), variant="preprint"))
    r = subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", f"{name}.tex"],
                       cwd=pdir, capture_output=True, text=True)
    if r.returncode != 0:
        log = (pdir / f"{name}.log").read_text(errors="replace")
        errs = [l for l in log.splitlines() if l.startswith("!") or "Error" in l][:20]
        raise SystemExit(f"[articles] LaTeX failed for the {slug} preprint:\n" + "\n".join(errs))
    import zipfile
    with zipfile.ZipFile(pdir / f"{name}-latex.zip", "w", zipfile.ZIP_DEFLATED) as z:
        for f in [f"{name}.tex", "pdarticle.cls", "numbered.bst", f"{name}.pdf"] + [f"figures/{x}" for x in used]:
            z.write(pdir / f, f"{name}/{f}")
    pages = re.search(r"Pages:\s+(\d+)", subprocess.run(["pdfinfo", str(pdir / f"{name}.pdf")], capture_output=True, text=True).stdout).group(1)
    print(f"[articles] {name}.pdf: {pages} pages, single column; LaTeX sources in {name}-latex.zip")
    return pdir


def web_figures(meta, art_dir: Path) -> callable:
    slug = meta["slug"]
    dest = OUT_ASSETS / slug
    dest.mkdir(parents=True, exist_ok=True)
    for src in sorted((art_dir / "figures").glob("*.jpg")):
        full, thumb = dest / src.name, dest / (src.stem + "-thumb.jpg")
        if not full.exists() or full.stat().st_mtime < src.stat().st_mtime:
            im = Image.open(src).convert("RGB")
            im.thumbnail((1800, 1800), Image.LANCZOS)
            im.save(full, quality=84, optimize=True, progressive=True)
            im.thumbnail((640, 640), Image.LANCZOS)
            im.save(thumb, quality=80, optimize=True, progressive=True)

    def url(name: str, prefix: str = "../") -> tuple[str, str]:
        stem = Path(name).stem
        return f"{prefix}assets/articles/{slug}/{stem}.jpg", f"{prefix}assets/articles/{slug}/{stem}-thumb.jpg"
    return url


def load(art_dir: Path) -> dict:
    meta = yaml.safe_load((art_dir / "meta.yaml").read_text())
    if isinstance(meta["published"], str):
        meta["published"] = dt.date.fromisoformat(meta["published"])
    meta["affiliations"] = {str(k): v for k, v in meta["affiliations"].items()}
    pre = meta.get("preprint")
    if pre and isinstance(pre.get("posted"), str):
        pre["posted"] = dt.date.fromisoformat(pre["posted"])
    return meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdf", action="store_true", help="typeset the PDF with the OmicVerse LaTeX template")
    ap.add_argument("--template", type=Path, default=ARTICLES / "_template", help="checkout of omicverse-latex-template (cloned if missing)")
    ap.add_argument("--only", help="build only this slug")
    ap.add_argument("--no-chrome", action="store_true", help="skip running sync_chrome.py")
    ap.add_argument("--preprint", action="store_true",
                    help="with --only: also build the single-column, unbranded preprint copy (PDF + LaTeX zip)")
    args = ap.parse_args()
    if args.preprint and not args.only:
        ap.error("--preprint needs --only <slug>")

    dirs = sorted(d for d in ARTICLES.iterdir() if d.is_dir() and not d.name.startswith("_") and (d / "meta.yaml").exists())
    entries = []
    OUT_PAGES.mkdir(parents=True, exist_ok=True)
    for d in dirs:
        meta = load(d)
        blocks = parse_article((d / "article.md").read_text())
        fig_url = web_figures(meta, d)
        if args.pdf and (not args.only or args.only == meta["slug"]):
            build_pdf(meta, blocks, d, ensure_template(args.template))
        if args.preprint and args.only == meta["slug"]:
            build_preprint(meta, blocks, d, ensure_template(args.template))
        pdf_info = json.loads((d / "pdf.json").read_text()) if (d / "pdf.json").exists() else None
        body, toc, figs = render_blocks_html(blocks, meta, fig_url)
        n_refs = sum(len(b["items"]) for b in blocks if b["t"] == "refs")
        (OUT_PAGES / f"{meta['slug']}.html").write_text(article_page(meta, body, toc, figs, pdf_info, n_refs))
        thumb_name = meta.get("thumbnail")
        entries.append({"meta": meta, "thumb": fig_url(f"{thumb_name}.jpg", prefix="")[1] if thumb_name else None})
        print(f"[articles] {meta['slug']}: {len(blocks)} blocks, {len(figs)} figures, {n_refs} references")
    entries.sort(key=lambda e: e["meta"]["published"], reverse=True)
    (SRC / "articles.html").write_text(index_page(entries))
    if not args.no_chrome:
        subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "sync_chrome.py")], check=True, stdout=subprocess.DEVNULL)
    return 0


if __name__ == "__main__":
    sys.exit(main())
