# Scoped rewrite validation reports

Four unpublished technical-report drafts, with Haoliang Zhu and Zehua Zeng as authors. Each has three multi-panel figures, measured results and explicit validation boundaries. The existing website builder excludes this `_drafts` directory; merging these sources does not publish pages.

- `splitst`: original-R score/balance/shift comparisons and repaired receiver/shift semantics. Weights are constructed; biological decontamination accuracy is not claimed.
- `spata2py`: original SPATA2 geometry comparison, including a deliberately documented convex/concave outline difference.
- `pycytotrace`: six real-data comparisons across two marrow datasets and three sizes, with input-gene sensitivity.
- `scissor-validation`: original and tighter Cox full-path comparisons. At solver tolerance 1e-9, zero entries exceed the unchanged comparison tolerance; tighter-threshold CV was not rerun.

`benchmark/` contains runnable scripts, metric tables, provenance and a checksummed evidence archive. It requires the documented source/input workspace; it does not install dependencies. Full article PDFs are not typeset. Figure PDFs are supplied. Browser layout validation was blocked by the local-file URL policy; static image links and exported figures were checked. Author/publication review remains pending.
