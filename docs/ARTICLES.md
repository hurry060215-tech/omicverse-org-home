# Articles (technical reports)

`omicverse.org/articles` hosts technical reports for OmicVerse methods and rebuilt R packages. Each report has a web page in the style of a journal article and a two-column PDF typeset with the [OmicVerse LaTeX template](https://github.com/omicverse/omicverse-latex-template). Both are generated from one folder, so the page and the PDF never drift apart.

## Folder layout

```
articles/<slug>/
  meta.yaml      title, type, published date, version, authors, affiliations, abstract, keywords, links, thumbnail
  article.md     the body (dialect below)
  figures/       fig1.jpg … figN.jpg, figS1.jpg …  (about 2400 px wide; used for the PDF and resized for the web)
  pdf.json       written by the build: page count and file size shown on the page
```

`articles/_template/` (a clone of the LaTeX template) and `articles/<slug>/build/` (LaTeX output) are ignored by git.

## Build

```bash
python3 scripts/build_articles.py --pdf   # pages, web figures and the PDF (needs latexmk + pdflatex)
python3 scripts/build_articles.py         # pages only, reusing the committed PDF
npm run build                             # copy src/ to dist/
```

The build writes `src/articles.html`, `src/articles/<slug>.html` and `src/assets/articles/<slug>/` (web figures, thumbnails and `<slug>.pdf`), then runs `scripts/sync_chrome.py`. Commit those outputs together with the source folder.

## `meta.yaml`

```yaml
slug: ft-kernel
type: Technical Report            # shown on the list page and used by its "Article Type" filter
title: "…"
short_title: FT-Kernel            # used in the BibTeX key
published: 2026-11-08
version: 1                        # bump when the report is revised
authors:
  - {name: First Author, affiliations: [1, 2]}
  - {name: Last Author, affiliations: [1], email: last@example.org}   # an email marks a corresponding author
affiliations:
  1: Department, Institution, City, Country
abstract: >-
  …
keywords: [one, two]
links: {code: https://github.com/…, api: ov.single.Fate, tutorial: https://omicverse.readthedocs.io/…}
thumbnail: fig1                   # figure shown on the list page
```

## `article.md` dialect

| Write | Becomes |
|---|---|
| `# Results` | section (`\section`) |
| `## Subheading` | subsection |
| `#* Funding` | back-matter section (unnumbered, not in the page's section list) |
| `# References` then `1. Author… (Year). Title. Journal. 10.xxxx/yyy` | numbered reference list; a trailing DOI becomes a link |
| `# Supplementary information` | supplementary part (one column in the PDF, after the references) |
| `[@3]`, `[@2,5]`, `[@7-9]` | superscript citations (`\cite{r3}` in LaTeX) |
| `Fig. 2c`, `Fig. S1a` | links to the figure |
| `*gene*`, `**bold**`, `` `code` ``, `[text](url)`, bare URLs | inline formatting and links |
| `$x_i$` and a line `$$ … \tag{4} $$` | inline and numbered display math (KaTeX on the page) |
| `1. item` with indented `- sub-item` lines | numbered list with bullets |
| `:::figure 3 fig3.jpg` / title line / legend line / `:::` | full-width figure with a bold title and legend |
| `:::table 1` / title line / pipe table / `:::` | full-width table; a row whose other cells are empty is a group heading |
