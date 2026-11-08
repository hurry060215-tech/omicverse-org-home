#!/usr/bin/env python3
import base64
import html
import io
import json
import keyword
import re
import shutil
import tokenize
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sync_chrome import ASSET_VERSION, footer_html, header_html  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
BLOG_DIR = ROOT / "src" / "blog"
BLOG_ASSET_DIR = ROOT / "src" / "assets" / "blog"


POSTS = [
    {
        "slug": "omicverse-latex-template",
        "title": "A LaTeX Template for OmicVerse Papers",
        "category": "Resources",
        "manual_markdown": r"""
# A LaTeX Template for OmicVerse Papers

Most OmicVerse projects end in a manuscript, and until now each group set its preprints in its own style. The new [omicverse-latex-template](https://github.com/omicverse/omicverse-latex-template) repository gives community papers, preprints and reports one shared look. It holds a single LaTeX class, `pdarticle`, a numbered reference style, vector logos and worked examples in English and Chinese.

{actions} [Download ZIP](https://github.com/omicverse/omicverse-latex-template/archive/refs/heads/main.zip) [View on GitHub](https://github.com/omicverse/omicverse-latex-template)

![Pages 1 and 2 of the English example in the default rail design](../assets/blog/latex-template/rail-pages.jpg)

## What is in the class

- **Two layouts.** A single-column manuscript by default, and a two-column journal layout with the `twocolumn` option. Both use the same type sizes: body 8.5/11 pt, title 20/21 pt.
- **Four page designs.** `rail` (the default), `masthead`, `panel` and `classic`.
- **Brand presets.** `brand=omicverse` and `brand=primordecode`, each with its own logo, colors and corner graphic. Every color can be overridden from the class options.
- **Chinese papers.** The `zh` option, compiled with XeLaTeX.
- **Review copies.** The `lineno` option adds line numbers.
- **Ready-made elements.** Clickable figure and table mentions, numbered boxes, tinted table heads and two-part figure captions.

## Quick start

Download the ZIP above, or clone the repository:

```bash
git clone https://github.com/omicverse/omicverse-latex-template.git
```

Copy `template.tex` and start writing. Keep `pdarticle.cls`, `numbered.bst` and the `logos/` folder next to your main `.tex` file, then compile:

```bash
# English papers: pdfLaTeX, XeLaTeX or LuaLaTeX
latexmk -pdf template.tex
# Chinese papers (the zh option) need XeLaTeX
latexmk -xelatex example-zh.tex
```

On Overleaf, upload `pdarticle.cls`, `numbered.bst` and `logos/` together with your `.tex` and `.bib` files. pdfLaTeX works for English papers; set the compiler to XeLaTeX for `zh`.

A minimal paper looks like this:

```latex
\documentclass[brand=omicverse]{pdarticle}

\journalname{OmicVerse}
\articletype{Research Article}
\journalinfo{OmicVerse Preprint, \today}
\copyrightholder{\the\year\ OmicVerse}

\title{Title of the paper}
\author{First Author,\affmark{1} and Last Author\affmark{1,*}}
\affiliation{1}{[Department, Institution, City, Country]}
\correspondence{author@example.com}

\summary{Summary text.}
\keywords{keyword one, keyword two}
\bibliographystyle{numbered}

\begin{document}
\maketitle
\section{Introduction}
Text.\cite{zeng2024omicverse}
\bibliography{example}
\end{document}
```

## Page designs

Pick a design with `design=`. Each image below is page 1 of the same example.

![Page 1 in the rail, masthead and panel designs](../assets/blog/latex-template/designs.jpg)

- `rail` puts a brand-color rail down the left edge of every page, with the journal line and article type set vertically in it. The summary sits between two rules.
- `masthead` runs a tinted band across the top of page 1 for the logo and the brand graphic, and puts the summary in a tinted box.
- `panel` sets a tinted panel behind the title block on page 1 and a thin tinted strip on later pages.
- `classic` keeps the earlier design for compatibility. It is not recommended for new papers.

Later pages in every design carry a small logo at the top left, the page number at the top right and a hairline rule.

## Two columns and line numbers

![Page 1 of the two-column layout, and page 2 of the single-column layout with line numbers](../assets/blog/latex-template/layouts.jpg)

`twocolumn` sets two columns with a 13 pt gap. A full-width `figure*` floats to the top of the next page, so define it one page before the place it should appear. `lineno` puts small gray line numbers in the right margin, or on the outer side of each column in the two-column layout.

## Your own colors and logo

![The OmicVerse and PrimorDecode brand presets](../assets/blog/latex-template/brands.jpg)

Each preset has three color roles: `tab` for the rail, `primary` for the summary, headings, figure titles and links, and `secondary` for the article type and the marker before each heading. Override any of them from the class options, for example `primary=1F6F8B, secondary=B23A48, tab=1F6F8B`. Roles you leave out keep the preset. `\brandlogo[height=0.4in]{path/to/logo}` swaps the logo, and `\brandgraphic{...}` replaces the TikZ corner graphic. For OmicVerse that graphic is a DNA helix that fades into a scatter of points, like a single-cell embedding.

## Chinese papers

![The English example and the Chinese example, page 1](../assets/blog/latex-template/chinese.jpg)

With `zh`, headings are never upper-cased, the labels become 摘要, 图, 表, 专栏, 通讯作者 and 参考文献, and line spacing widens by 1.08×. The CJK font is the first one installed of Source Han Sans SC, Noto Sans CJK SC, Hiragino Sans GB and PingFang SC.

## License

`pdarticle.cls` and the example and template files are distributed under the LaTeX Project Public License 1.3c or later. `numbered.bst` is based on `model6-num-names.bst` from Elsevier's elsarticle bundle, also under the LPPL. The OmicVerse and PrimorDecode logos belong to their owners and are not covered by the LPPL.

Questions and improvements are welcome as issues or pull requests on the [repository](https://github.com/omicverse/omicverse-latex-template).

{actions} [Download ZIP](https://github.com/omicverse/omicverse-latex-template/archive/refs/heads/main.zip) [Read the full README](https://github.com/omicverse/omicverse-latex-template#readme)
""",
        "excerpt": "pdarticle gives OmicVerse papers and preprints one shared look: single- or two-column layouts, four page designs, brand colors, line numbers and Chinese support. Download the ZIP and start from template.tex.",
    },
    {
        "slug": "nanostring-cosmx-workflow",
        "title": "Analyzing NanoString CosMx Data",
        "category": "Spatial",
        "source": PROJECT / "omicverse-tutorials/docs/Tutorials-space/t_nanostring_preprocess.ipynb",
        "wechat_source": PROJECT / "omicverse-tutorials/wechat_md/t_nanostring_preprocess_wechat.md",
        "excerpt": "Read CosMx data, validate FOV-level spatial structure, inspect segmentation quality, and prepare inputs for downstream spatial modeling.",
    },
    {
        "slug": "visium-hd-workflow",
        "title": "Analyzing Visium HD Data",
        "category": "Spatial",
        "source": PROJECT / "omicverse-tutorials/docs/Tutorials-space/t_visium_hd_preprocess.ipynb",
        "wechat_source": PROJECT / "omicverse-tutorials/wechat_md/t_visium_hd_preprocess_wechat.md",
        "excerpt": "Work with Visium HD bin-level output, tissue images, spatially variable genes, low-dimensional representations, and segmentation-derived views.",
    },
    {
        "slug": "xenium-workflow",
        "title": "Analyzing Xenium Data",
        "category": "Spatial",
        "source": PROJECT / "omicverse-tutorials/docs/Tutorials-space/t_xenium_preprocess.ipynb",
        "wechat_source": PROJECT / "omicverse-tutorials/wechat_md/t_xenium_preprocess_wechat.md",
        "excerpt": "Load a 10x Xenium outs folder into AnnData, keep cell polygon boundaries, run QC and Leiden clustering, then project clusters back onto tissue space.",
    },
    {
        "slug": "atera-workflow",
        "title": "Reading and Preprocessing Atera Data",
        "category": "Spatial",
        "source": PROJECT / "omicverse-core/omicverse_guide/docs/Tutorials-space/t_atera_preprocess.ipynb",
        "wechat_source": PROJECT / "omicverse-tutorials/wechat_md/t_atera_preprocess_wechat.md",
        "excerpt": "Use OmicVerse to read 10x Atera whole-transcriptome spatial data, align expression, centroids, boundaries, nuclei, morphology channels, and annotations.",
    },
    {
        "slug": "cellpose-visium-hd",
        "title": "Cell Segmentation for Visium HD with Cellpose",
        "category": "Segmentation",
        "source": PROJECT / "omicverse-tutorials/docs/Tutorials-space/t_cellpose.ipynb",
        "wechat_source": PROJECT / "omicverse-tutorials/wechat_md/t_cellpose_wechat.md",
        "excerpt": "Connect Visium HD 2 um bins, H&E images, Cellpose segmentation, and bin-to-cell aggregation in a pure Python workflow.",
    },
    {
        "slug": "plot1cell-python",
        "title": "Circular UMAPs with the Python Version of plot1cell",
        "category": "Visualization",
        "source": PROJECT / "omicverse-tutorials/docs/Tutorials-plotting/t_plot1cell.ipynb",
        "wechat_source": PROJECT / "omicverse-tutorials/wechat_md/t_plot1cell_wechat.md",
        "excerpt": "Use ov.pl.plot1cell to place UMAP or t-SNE cells inside a circular layout and add cluster-aware metadata rings around the embedding.",
    },
    {
        "slug": "leiden-resolution-selection",
        "title": "How to Choose Leiden Resolution",
        "category": "Single-cell",
        "wechat_source": PROJECT / "omicverse-tutorials/wechat_md/compare_resolution_methods_pbmc8k_wechat.md",
        "manual_markdown": """
# How to Choose Leiden Resolution: Comparing `leiden`, `auto_resolution`, and `champ`

When running single-cell clustering, almost everyone eventually asks the same question: what should the Leiden `resolution` be?

If the value is too low, distinct biological populations can be merged. If it is too high, a stable cell population can be fragmented into many small clusters. The numbers themselves, such as `0.5`, `1.0`, or `1.2`, do not have intrinsic biological meaning; they are controls for the granularity of community detection.

This tutorial asks a practical question: which OmicVerse tools can help choose a Leiden resolution?

We use `ov.datasets.pbmc8k()` because it includes curated `predicted_celltype` labels, making it possible to compare clustering results against a biological reference.

| Method | Function | Main idea |
|---|---|---|
| Manual Leiden | `ov.pp.leiden(adata, resolution=r)` | The user chooses the resolution directly |
| bootstrap-ARI | `ov.single.auto_resolution()` | Choose the resolution that is most stable under data perturbation |
| CHAMP | `ov.single.auto_resolution(method='champ')` | Choose the partition that is most stable across the modularity landscape |

The key point is that `auto_resolution` and `champ` answer different questions. The bootstrap-ARI workflow asks which resolution is stable under resampling. CHAMP asks which partition is stable as the gamma parameter changes. Neither method replaces biological interpretation; both provide stronger starting points than an arbitrary default.

## Environment setup

```python
import os, io, contextlib
os.environ.setdefault("OMICVERSE_DISABLE_LLM", "1")

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import omicverse as ov
ov.style()
```

## Load and preprocess PBMC8k

```python
adata = ov.datasets.pbmc8k()
adata.var_names_make_unique()
```

The dataset includes `predicted_celltype`, so later results can be scored against the curated cell type annotation.

```python
adata = ov.pp.qc(
    adata,
    tresh={"mito_perc": 0.2, "nUMIs": 500, "detected_genes": 250},
    doublets_method="scrublet",
)
ov.pp.preprocess(adata, mode="shiftlog|pearson", n_HVGs=2000)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
ov.pp.scale(adata)
ov.pp.pca(adata, layer="scaled", n_pcs=50)
ov.pp.neighbors(adata, n_neighbors=15, use_rep="scaled|original|X_pca")
ov.pp.umap(adata)
```

After QC, the example retains 7,677 cells and 2,000 highly variable genes. All methods use the same neighbor graph, so the comparison is controlled.

## Baseline: manual Leiden

```python
ov.pp.leiden(adata, resolution=0.5, key_added="leiden_r05")
ov.pp.leiden(adata, resolution=1.0, key_added="leiden_r10")
```

In this run, `resolution=0.5` produced 13 clusters and `resolution=1.0` produced 18 clusters. That behavior is expected: higher resolution usually means more clusters. The problem is that more clusters do not necessarily mean better biology.

## Bootstrap-ARI with `auto_resolution`

`ov.single.auto_resolution()` uses null-adjusted bootstrap-ARI stability by default. It repeatedly subsamples and reclusters at candidate resolutions, then subtracts a null background estimated from permuted data.

```python
with contextlib.redirect_stdout(io.StringIO()):
    _, autores_best, autores_scores = ov.single.auto_resolution(
        adata,
        resolutions=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5],
        n_subsamples=5,
        n_null_subsamples=3,
        random_state=0,
        key_added="leiden_auto",
    )
print(f"auto_resolution chose r = {autores_best}")
```

In the PBMC8k example, `auto_resolution` selected `r = 0.4`.

## CHAMP

CHAMP follows a different logic. Instead of asking which resolution survives bootstrap perturbation, it looks for stable partitions across the modularity landscape.

```python
_, champ_best_r, champ_df = ov.single.auto_resolution(
    adata,
    method="champ",
    key_added="leiden_champ",
)
```

In the same PBMC8k run, CHAMP returned an equivalent Leiden resolution around `r = 0.356`, with a stable gamma range near `[0.08, 0.21]`.

## What the comparison showed

On PBMC8k, the automatic methods performed better than manually choosing `resolution=0.5` or `resolution=1.0`. CHAMP had a slightly higher ARI, while bootstrap `auto_resolution` had a slightly higher NMI.

| Method | Clusters | ARI | NMI |
|---|---:|---:|---:|
| Manual Leiden, r=0.5 | 13 | lower | lower |
| Manual Leiden, r=1.0 | 18 | lower | lower |
| auto_resolution, r=0.4 | 11 | 0.822 | 0.807 |
| CHAMP, gamma in [0.08, 0.21] | 8 | 0.829 | 0.794 |

The practical reading is simple: automatic resolution selection is useful, but it should be treated as a reproducible starting point. Marker genes, known cell states, sample design, and downstream biological questions still decide whether a clustering is useful.
""",
        "excerpt": "Compare manual Leiden, bootstrap stability, and CHAMP on PBMC8k to choose a reproducible starting resolution for single-cell clustering.",
    },
    {
        "slug": "omicverse-2-1-update",
        "title": "What Changed in OmicVerse 2.1.x",
        "category": "Release",
        "manual_markdown": """
# What Changed in OmicVerse 2.1.x

OmicVerse `2.1.x` is not just a small maintenance update. From `v2.0.0` to the current development tree, the project accumulated hundreds of commits, hundreds of changed files, and close to one hundred thousand added lines of code. Two major branches, metabolomics and 16S / microbiome analysis, were also moving toward the merge window.

The update focuses on trajectory inference, spatial omics, preprocessing and I/O, plotting, the agent runtime, and new multi-omics modules.

To update:

```bash
pip install -U omicverse
```

## Update overview

### PP module

- Added a pure-Python `pydoubletfinder` backend for doublet detection.
- Automatically detects mitochondrial prefixes, avoiding manual `MT-` versus `mt-` handling.
- Upgraded Harmony to `harmonypy v0.2.0`, with GPU, CPU NumPy, and Apple Silicon MLX backends.
- Improved memory handling in `torch_pca` and `covariance_eigh`.
- Added the Rust-backed `anndataoom` out-of-memory AnnData reader.

### Single module

- Added `ov.single.Monocle`, a pure Python rebuild of Monocle 2 that avoids `rpy2` and a separate R environment.
- Continued work on dynamic trajectory tools, including lineage-aware trend plotting and improvements around Palantir and Slingshot visualization.

### Space module

- Added `ov.io.read_xenium` for full Xenium `outs/` ingestion.
- Converts Xenium `cell_boundaries` into WKT polygons that can be passed directly to `ov.pl.spatialseg`.
- Added `ov.io.write_visium_hd_cellseg` to export cell-level AnnData back into a Space Ranger v4-compatible directory structure.
- Added CosMx FOV-aware plotting, including multi-FOV layout and per-FOV overlays.
- Added CellSAM support and renamed `stardist()` to the more general `cellseg()`.
- Added `method="cellcharter"` to `ov.utils.cluster`.

### Plotting module

- Added Marsilea-based heatmap plotting.
- Added a set of cell-cell communication plotting APIs.
- Added `ov.pl.create_custom_colormap`.
- Added half-violin boxplots and started consolidating older plotting APIs.
- Improved legend behavior for subset plotting.

### Agent module

- Reworked `ov.Agent`, Jarvis, and OVAgent runtime internals.
- Added Codex OAuth support.
- Moved multiple message channels onto a unified `MessageRuntime`.
- Added sandboxing, timeouts, stdout guards, retries, and security hardening.

### New modules

- Added `ov.metabol` for metabolomics workflows, including ID mapping, MSEA, Mummichog, SERRF, DGCA, ASCA / MixedLM, ROC, and biomarker panels.
- Added `ov.micro` and `ov.alignment` for 16S / microbiome analysis from amplicons to microbiome AnnData objects.

## Why the update matters

The theme of this release is reducing dependency friction while making high-resolution spatial workflows more complete. Xenium, Visium HD, and CosMx are no longer only file formats that can be read; they now have fuller analysis paths that connect I/O, preprocessing, segmentation, plotting, and export.

The pure-Python Monocle 2 rebuild points in the same direction. When core methods no longer require an R bridge, they become easier to install, easier to run in notebooks and pipelines, and easier to expose to agents.

For day-to-day users, the headline is practical: update OmicVerse if you need modern spatial data support, more robust preprocessing, richer plotting, or a more capable agent runtime.
""",
        "excerpt": "A compact English release note for the 2.1.x line: pure-Python Monocle, Xenium and Visium HD workflows, plotting updates, agent runtime work, metabolomics, and microbiome support.",
    },
    {
        "slug": "anthropic-mythos-security",
        "title": "What Anthropic Mythos Preview Signals for Agent Security",
        "category": "Agents",
        "wechat_source": PROJECT / "omicverse-tutorials/wechat_md/anthropic_mythos_preview_wechat.md",
        "manual_markdown": """
# What Anthropic Mythos Preview Signals for Agent Security

The original WeChat note discussed Anthropic's security report on Claude Mythos Preview and Project Glasswing. The relevant takeaway for a scientific computing community is not that every model suddenly becomes an autonomous attacker. It is that capable agents are moving closer to workflows where they can inspect systems, reason about vulnerabilities, assemble tools, and act across a computer environment.

That matters for omics agents because the same capabilities that make an analysis agent useful also create risk:

- file-system access,
- package installation,
- shell execution,
- notebook execution,
- remote data fetching,
- credential handling,
- report generation.

## From code generation to operational access

Earlier AI coding tools mostly produced snippets. Modern agents can run commands, inspect outputs, retry after failures, and chain multiple tools. In bioinformatics, that is attractive because workflows are messy: data formats differ, dependencies break, and notebooks need iterative execution.

The security implication is direct: once an agent can operate a real environment, the boundary is no longer just the prompt. The boundary includes filesystem permissions, network access, secrets, sandbox policy, package provenance, and human review checkpoints.

## What agent-ready omics needs

OmicVerse agents should be useful without being opaque. A safe workflow needs:

- explicit input and output contracts,
- provenance for commands, packages, and datasets,
- bounded execution environments,
- guarded access to credentials,
- logs that a human can audit,
- recoverable failure states,
- reproducible reports rather than hidden notebook state.

## Practical direction

The long-term lesson is not to avoid agents. It is to design them like scientific instruments: powerful, constrained, inspectable, and calibrated. For omics, that means the agent layer should sit on top of stable APIs, documented workflows, and clear permission boundaries.
""",
        "excerpt": "A security-oriented note on why increasingly capable agents need bounded execution, provenance, credential discipline, and auditable scientific workflows.",
    },
]


def escape_attr(value):
    return html.escape(str(value), quote=True)


def inline_md(text):
    text = html.escape(text)
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", lambda m: f'<a href="{escape_attr(m.group(2))}" target="_blank" rel="noopener">{m.group(1)}</a>', text)
    return text


def highlight_python(code):
    out = []
    last = (1, 0)
    builtins = {
        "print", "len", "range", "str", "int", "float", "list", "dict", "set",
        "tuple", "zip", "enumerate", "open", "Path", "True", "False", "None",
    }
    try:
        tokens = tokenize.generate_tokens(io.StringIO(code).readline)
        for tok in tokens:
            token_type, token, start, end, _ = tok
            if token_type in {tokenize.ENCODING, tokenize.ENDMARKER}:
                continue
            line, col = start
            last_line, last_col = last
            if line > last_line:
                out.append("\n" * (line - last_line))
                out.append(" " * col)
            else:
                out.append(" " * max(0, col - last_col))

            cls = ""
            if token_type == tokenize.COMMENT:
                cls = "tok-comment"
            elif token_type == tokenize.STRING:
                cls = "tok-string"
            elif token_type == tokenize.NUMBER:
                cls = "tok-number"
            elif token_type == tokenize.NAME:
                if keyword.iskeyword(token):
                    cls = "tok-keyword"
                elif token in builtins:
                    cls = "tok-builtin"
            elif token_type == tokenize.OP:
                cls = "tok-op"

            escaped = html.escape(token)
            out.append(f'<span class="{cls}">{escaped}</span>' if cls else escaped)
            last = end
    except tokenize.TokenError:
        return html.escape(code)
    return "".join(out)


def highlight_code(code, lang):
    lang = (lang or "").lower()
    if lang in {"python", "py"}:
        return highlight_python(code)
    if lang in {"bash", "sh", "shell", "zsh"}:
        escaped = html.escape(code)
        escaped = re.sub(r"(^|\n)(\s*#.*)", r'\1<span class="tok-comment">\2</span>', escaped)
        escaped = re.sub(r"\b(pip|python|python3|conda|git|cd|mkdir|export)\b", r'<span class="tok-builtin">\1</span>', escaped)
        return escaped
    return html.escape(code)


def extract_wechat_figures(path, max_figures=12):
    if not path or not path.exists():
        return ""
    matches = re.findall(r"!\[([^\]]*)\]\(([^)]+)\)", path.read_text())
    if not matches:
        return ""
    lines = ["## Figures"]
    for idx, (_, src) in enumerate(matches[:max_figures], 1):
        lines.append(f"![Figure {idx}]({src})")
    return "\n\n".join(lines)


def markdown_to_html(markdown):
    lines = markdown.strip().splitlines()
    blocks = []
    para = []
    list_items = []
    table = []
    code = []
    in_code = False
    code_lang = ""

    def flush_para():
        nonlocal para
        if para:
            blocks.append(f"<p>{inline_md(' '.join(x.strip() for x in para))}</p>")
            para = []

    def flush_list():
        nonlocal list_items
        if list_items:
            blocks.append("<ul>" + "".join(f"<li>{inline_md(item)}</li>" for item in list_items) + "</ul>")
            list_items = []

    def flush_table():
        nonlocal table
        if table:
            rows = []
            for i, row in enumerate(table):
                cells = [c.strip() for c in row.strip("|").split("|")]
                tag = "th" if i == 0 else "td"
                if i == 1 and all(set(c.replace(":", "").replace("-", "")) == set() for c in cells):
                    continue
                rows.append("<tr>" + "".join(f"<{tag}>{inline_md(c)}</{tag}>" for c in cells) + "</tr>")
            blocks.append("<table>" + "".join(rows) + "</table>")
            table = []

    for raw in lines:
        line = raw.rstrip()
        if line.startswith("```"):
            if in_code:
                code_text = chr(10).join(code)
                blocks.append(f'<pre><code class="language-{escape_attr(code_lang)}">{highlight_code(code_text, code_lang)}</code></pre>')
                code = []
                code_lang = ""
                in_code = False
            else:
                flush_para(); flush_list(); flush_table()
                code_lang = line.strip("`").strip()
                in_code = True
            continue
        if in_code:
            code.append(line)
            continue

        if not line.strip():
            flush_para(); flush_list(); flush_table()
            continue
        if line.startswith("<!--"):
            continue
        if line.startswith("# "):
            flush_para(); flush_list(); flush_table()
            blocks.append(f"<h1>{inline_md(line[2:].strip())}</h1>")
            continue
        if line.startswith("## "):
            flush_para(); flush_list(); flush_table()
            blocks.append(f"<h2>{inline_md(line[3:].strip())}</h2>")
            continue
        if line.startswith("### "):
            flush_para(); flush_list(); flush_table()
            blocks.append(f"<h3>{inline_md(line[4:].strip())}</h3>")
            continue
        if line.lstrip().startswith("- "):
            flush_para(); flush_table()
            list_items.append(line.lstrip()[2:].strip())
            continue
        if re.match(r"^\s*\d+\.\s+", line):
            flush_para(); flush_table()
            list_items.append(re.sub(r"^\s*\d+\.\s+", "", line).strip())
            continue
        if line.startswith("|") and line.endswith("|"):
            flush_para(); flush_list()
            table.append(line)
            continue
        if line.strip().startswith("{actions}"):
            # a row of buttons: the first link is primary, the rest secondary; ZIP links download in place
            flush_para(); flush_list(); flush_table()
            buttons = []
            for i, (label, href) in enumerate(re.findall(r"\[([^\]]+)\]\(([^)]+)\)", line)):
                cls = "btn btn-primary" if i == 0 else "btn btn-ghost"
                target = "" if href.endswith(".zip") else ' target="_blank" rel="noopener"'
                buttons.append(f'<a class="{cls}" href="{escape_attr(href)}"{target}>{html.escape(label)}</a>')
            blocks.append('<div class="article-actions">' + "".join(buttons) + "</div>")
            continue
        if re.match(r"!\[[^\]]*\]\([^)]+\)", line.strip()):
            flush_para(); flush_list(); flush_table()
            alt, src = re.match(r"!\[([^\]]*)\]\(([^)]+)\)", line.strip()).groups()
            caption = f"<figcaption>{inline_md(alt)}</figcaption>" if alt else ""
            blocks.append(f'<figure><img src="{escape_attr(src)}" alt="{escape_attr(alt)}" loading="lazy" />{caption}</figure>')
            continue
        para.append(line)

    flush_para(); flush_list(); flush_table()
    return "\n".join(blocks)


def output_image_markdown(cell, slug, cell_index):
    lines = []
    figure_index = 0
    asset_dir = BLOG_ASSET_DIR / slug
    asset_dir.mkdir(parents=True, exist_ok=True)
    for output in cell.get("outputs", []):
        data = output.get("data") or {}
        for mime, ext in (("image/png", "png"), ("image/jpeg", "jpg")):
            payload = data.get(mime)
            if not payload:
                continue
            figure_index += 1
            encoded = "".join(payload) if isinstance(payload, list) else payload
            filename = f"cell-{cell_index:03d}-{figure_index}.{ext}"
            (asset_dir / filename).write_bytes(base64.b64decode(encoded))
            lines.append(f"![Output figure](../assets/blog/{slug}/{filename})")
            break
    return "\n\n".join(lines)


def notebook_to_markdown(path, slug):
    nb = json.loads(path.read_text())
    chunks = []
    asset_dir = BLOG_ASSET_DIR / slug
    if asset_dir.exists():
        shutil.rmtree(asset_dir)
    for idx, cell in enumerate(nb.get("cells", []), 1):
        source = "".join(cell.get("source", [])).strip()
        if not source:
            continue
        if cell.get("cell_type") == "markdown":
            chunks.append(source)
        elif cell.get("cell_type") == "code":
            chunks.append("```python\n" + source + "\n```")
            figures = output_image_markdown(cell, slug, idx)
            if figures:
                chunks.append(figures)
    return "\n\n".join(chunks)


def article_shell(post, article_html):
    title = escape_attr(post["title"])
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{title} — OmicVerse</title>
  <link rel="icon" type="image/png" sizes="32x32" href="../assets/brand/favicon-32.png" />
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="../styles.css?v={ASSET_VERSION}" />
  <script src="../script.js?v={ASSET_VERSION}" defer></script>
</head>
<body>
{header_html("blog", "../")}

  <main class="article page-width">
    <a class="article-back" href="../blog.html">← Back to blog</a>
    <p class="overline">{html.escape(post["category"])}</p>
    {article_html}
  </main>

{footer_html("../")}
</body>
</html>
"""


def post_thumb(post):
    """First figure of the post (local or remote), or None."""
    md = post.get("_markdown") or ""
    m = re.search(r"!\[[^\]]*\]\(([^)\s]+)\)", md)
    if not m:
        return None
    src = m.group(1)
    return src[3:] if src.startswith("../") else src


def post_card(post, prefix=""):
    thumb = post.get("_thumb")
    if thumb:
        src = thumb if thumb.startswith("http") else prefix + thumb
        fig = f'<span class="post-thumb"><img src="{escape_attr(src)}" alt="" loading="lazy" /></span>'
    else:
        fig = f'<span class="post-thumb is-placeholder">{html.escape(post["category"])}</span>'
    return f"""        <a class="post-card" href="{prefix}blog/{escape_attr(post["slug"])}.html">
          {fig}
          <span class="post-chip">{html.escape(post["category"])}</span>
          <span class="post-title">{html.escape(post["title"])}</span>
          <span class="post-desc">{html.escape(post["excerpt"])}</span>
        </a>"""


def blog_index(posts):
    cards = [post_card(p) for p in posts]
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Blog — OmicVerse</title>
  <link rel="icon" type="image/png" sizes="32x32" href="assets/brand/favicon-32.png" />
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="styles.css?v={ASSET_VERSION}" />
  <script src="script.js?v={ASSET_VERSION}" defer></script>
</head>
<body>
{header_html("blog")}

  <main>
    <section class="page-hero page-width">
      <p class="overline">Blog</p>
      <h1>Tutorials and technical notes.</h1>
    </section>

    <section class="page-width wide" style="padding-bottom:64px;">
      <div class="post-grid">
{chr(10).join(cards)}
      </div>
    </section>
  </main>

{footer_html()}
</body>
</html>
"""


def main():
    BLOG_DIR.mkdir(parents=True, exist_ok=True)
    for old in BLOG_DIR.glob("*.html"):
        old.unlink()
    for post in POSTS:
        if post.get("manual_markdown"):
            markdown = post["manual_markdown"]
            figures = extract_wechat_figures(post.get("wechat_source"))
        else:
            markdown = notebook_to_markdown(post["source"], post["slug"])
            figures = ""
        if figures:
            markdown = f"{markdown}\n\n{figures}"
        post["_markdown"] = markdown
        post["_thumb"] = post_thumb(post)
        article_html = markdown_to_html(markdown)
        (BLOG_DIR / f"{post['slug']}.html").write_text(article_shell(post, article_html))
    (ROOT / "src" / "blog.html").write_text(blog_index(POSTS))


if __name__ == "__main__":
    main()
