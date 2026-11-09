// ─── Articles: index filters, math, section tracker, copy buttons ─────────
(function () {
  "use strict";

  const form = document.querySelector("[data-articles-filter]");
  if (form) {
    const rows = Array.from(document.querySelectorAll("[data-articles-list] .ar-row"));
    const empty = document.querySelector("[data-articles-empty]");
    const apply = () => {
      const type = form.elements.type.value;
      const year = form.elements.year.value;
      let shown = 0;
      rows.forEach((r) => {
        const ok = (!type || r.dataset.type === type) && (!year || r.dataset.year === year);
        r.hidden = !ok;
        if (ok) shown += 1;
      });
      if (empty) empty.hidden = shown > 0;
    };
    form.addEventListener("change", apply);
    form.addEventListener("submit", (e) => e.preventDefault());
  }

  const main = document.querySelector(".paper-main");
  if (main && typeof window.renderMathInElement === "function") {
    window.renderMathInElement(main, {
      delimiters: [
        { left: "$$", right: "$$", display: true },
        { left: "$", right: "$", display: false },
      ],
      throwOnError: false,
    });
  }

  const tabs = document.querySelector("[data-side-tabs]");
  if (tabs) {
    tabs.addEventListener("click", (e) => {
      const btn = e.target.closest("[data-tab]");
      if (!btn) return;
      tabs.querySelectorAll("[data-tab]").forEach((b) => b.setAttribute("aria-selected", String(b === btn)));
      tabs.querySelectorAll("[data-panel]").forEach((p) => { p.hidden = p.dataset.panel !== btn.dataset.tab; });
    });
  }

  const tocLinks = Array.from(document.querySelectorAll(".side-toc a"));
  if (tocLinks.length && "IntersectionObserver" in window) {
    const byId = new Map(tocLinks.map((a) => [a.getAttribute("href").slice(1), a]));
    const observer = new IntersectionObserver((entries) => {
      entries.forEach((entry) => {
        if (!entry.isIntersecting) return;
        const link = byId.get(entry.target.id);
        if (!link) return;
        tocLinks.forEach((a) => a.classList.toggle("is-active", a === link));
      });
    }, { rootMargin: "-90px 0px -70% 0px" });
    byId.forEach((_, id) => {
      const el = document.getElementById(id);
      if (el) observer.observe(el);
    });
  }

  document.querySelectorAll("[data-copy-target]").forEach((btn) => {
    btn.addEventListener("click", async () => {
      const el = document.querySelector(btn.dataset.copyTarget);
      if (!el) return;
      try {
        await navigator.clipboard.writeText(el.textContent.trim());
        const label = btn.textContent;
        btn.textContent = "Copied";
        window.setTimeout(() => { btn.textContent = label; }, 1400);
      } catch (_) { /* clipboard unavailable */ }
    });
  });
})();
