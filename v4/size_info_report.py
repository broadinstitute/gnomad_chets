"""Build the variant size-info report — figures, stats.json, and
``report.md`` — from a per-variant size-info Hail Table.

This is the in-pipeline equivalent of ``/tmp/chr19_size_report_md.py``,
which read a TSV exported from an earlier run. Adapted to:

* Read the size-info HT directly (no manual TSV export). The HT is
  exported to a local TSV via ``ht.export`` once, then loaded into
  pandas for the existing pandas/matplotlib analysis. Keeps the report
  logic identical to what we validated for chr19; only the input
  reader changes.
* Apply ``heavy_contribution_cutoff`` at report time (the HT has no
  ``is_heavy`` column; the build step is pure-data per the size-info
  refactor). Heavy iff ``_contribution >= cutoff``.
* Upload all artifacts (figures + stats.json + report.md) to a GCS
  output directory via ``hl.hadoop_copy``.

Public entrypoint: :func:`build_report`.
"""
import json
import os
import shutil
import sys
import tempfile
import urllib.parse
import urllib.request
from typing import Optional

import hail as hl


TARGET_HEAVY_PARTITION_BYTES = 500 * 1024**2


def _fmt_b(b):
    if b is None:
        return "n/a"
    try:
        f = float(b)
    except (TypeError, ValueError):
        return "n/a"
    if not (f == f):  # NaN
        return "n/a"
    for u in ["B", "KB", "MB", "GB", "TB", "PB", "EB"]:
        if abs(f) < 1024 or u == "EB":
            return f"{f:.2f} {u}"
        f /= 1024


def _fetch_gene_symbols(ensembl_ids, cache_path):
    """Best-effort symbol lookup via mygene.info. Caches to a JSON
    file. Returns ``{ensembl_id: symbol}`` — missing entries map to
    themselves so callers always get a valid label."""
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            cache = json.load(f)
    todo = [i for i in set(ensembl_ids) if i not in cache]
    if not todo:
        for i in ensembl_ids:
            cache.setdefault(i, i)
        return cache
    print(f"Fetching {len(todo)} symbols from mygene.info...", file=sys.stderr)
    for i in range(0, len(todo), 500):
        chunk = todo[i:i + 500]
        data = urllib.parse.urlencode(
            {"ids": ",".join(chunk), "fields": "symbol"}
        ).encode()
        try:
            req = urllib.request.Request(
                "http://mygene.info/v3/gene", data=data,
            )
            resp = json.loads(
                urllib.request.urlopen(req, timeout=60).read()
            )
            for d in resp:
                cache[d["query"]] = d.get("symbol", d["query"])
        except Exception as e:
            print(f"  mygene.info batch failed: {e}", file=sys.stderr)
            for q in chunk:
                cache.setdefault(q, q)
    with open(cache_path, "w") as f:
        json.dump(cache, f, indent=2, sort_keys=True)
    for i in ensembl_ids:
        cache.setdefault(i, i)
    return cache


def _hist_log_x(data, ax, label, color, bins, alpha=0.65):
    import numpy as np
    arr = np.asarray(data, dtype=np.float64)
    arr = arr[arr > 0]
    if len(arr) == 0:
        return
    ax.hist(arr, bins=bins, label=label, color=color,
            alpha=alpha, edgecolor="black", linewidth=0.3)


def _save(fig, fig_dir, fname):
    import matplotlib.pyplot as plt
    path = f"{fig_dir}/{fname}"
    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return f"figures/{fname}"


def _quantiles(s, qs=(0.0, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1.0)):
    return {f"p{q*100:g}": int(s.quantile(q)) for q in qs}


def _generate(
    df,
    local_dir: str,
    heavy_contribution_cutoff: int,
    fetch_symbols: bool,
):
    """Pandas/matplotlib analysis. Writes figures + stats.json +
    report.md into ``local_dir``."""
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = f"{local_dir}/figures"
    os.makedirs(fig_dir, exist_ok=True)

    # Heavy decision applied here (not in the source HT).
    df["is_heavy"] = df["_contribution"] >= heavy_contribution_cutoff
    heavy = df[df.is_heavy]
    light = df[~df.is_heavy]
    heavy_s = df[df.is_heavy & (df.split_count > 1)]
    has_gene_id = "gene_id" in df.columns

    n_total = len(df)
    total_contrib = int(df._contribution.sum())
    heavy_contrib = int(heavy._contribution.sum())
    heavy_s_contrib = int(heavy_s._contribution.sum())
    total_bytes = int(df._bytes.sum())
    heavy_bytes = int(heavy._bytes.sum())
    heavy_s_bytes = int(heavy_s._bytes.sum())

    # Per-gene aggregates require gene_id.
    if has_gene_id:
        # gene_id is array<str>; explode to one row per (variant, gene)
        # so per-gene aggregates capture every gene a variant touches.
        df_x = df.explode("gene_id")
        df_x = df_x[df_x.gene_id.notna()]
        gene = df_x.groupby("gene_id").agg(
            n_variants=("v_idx", "count"),
            n_heavy_variants=("is_heavy", "sum"),
            sum_contribution=("_contribution", "sum"),
            sum_bytes=("_bytes", "sum"),
            max_contribution=("_contribution", "max"),
            max_bytes=("_bytes", "max"),
            max_split_count=("split_count", "max"),
        )
        gene_h = df_x[df_x.is_heavy].groupby("gene_id").agg(
            heavy_sum_contribution=("_contribution", "sum"),
            heavy_sum_bytes=("_bytes", "sum"),
        )
        gene_s = df_x[df_x.is_heavy & (df_x.split_count > 1)].groupby(
            "gene_id"
        ).agg(
            n_heavy_split_variants=("v_idx", "count"),
            heavy_split_sum_contribution=("_contribution", "sum"),
            heavy_split_sum_bytes=("_bytes", "sum"),
            heavy_split_max_split_count=("split_count", "max"),
        )
        gene = (
            gene.join(gene_h, how="left")
            .join(gene_s, how="left")
            .fillna({
                "heavy_sum_contribution": 0,
                "heavy_sum_bytes": 0,
                "n_heavy_split_variants": 0,
                "heavy_split_sum_contribution": 0,
                "heavy_split_sum_bytes": 0,
                "heavy_split_max_split_count": 0,
            })
            .reset_index()
        )
        gene["heavy_frac"] = gene.n_heavy_variants / gene.n_variants
        if fetch_symbols:
            symbols = _fetch_gene_symbols(
                gene.gene_id.tolist(),
                cache_path=f"{local_dir}/gene_symbols.json",
            )
            gene["symbol"] = gene.gene_id.map(lambda g: symbols.get(g, g))
            gene["display"] = gene.apply(
                lambda r: r.symbol if r.symbol and r.symbol != r.gene_id
                else r.gene_id,
                axis=1,
            )
        else:
            gene["symbol"] = gene["gene_id"]
            gene["display"] = gene["gene_id"]
    else:
        gene = None

    figs = {}

    # F1: split_count distribution
    sc = df.split_count.value_counts().sort_index()
    tail = sc[sc.index > 1]
    for yscale, suffix in [("log", ""), ("linear", "_linear")]:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        axes[0].bar(sc.index, sc.values, width=0.9, color="#4477aa")
        axes[0].set_yscale(yscale)
        if yscale == "linear":
            axes[0].set_xlim(1.5, 50)
            ymax = (
                int(tail[tail.index <= 50].max() * 1.05)
                if (tail.index <= 50).any() else 1
            )
            axes[0].set_ylim(0, ymax)
            axes[0].set_title("split_count distribution (2–50, linear y)")
        else:
            axes[0].set_xlim(0, 50)
            axes[0].set_title("split_count distribution (0–50, log y)")
        axes[0].set_xlabel("split_count")
        axes[0].set_ylabel("# variants" if yscale == "linear"
                          else "# variants (log)")
        axes[0].grid(True, alpha=0.3)
        if len(tail):
            axes[1].bar(tail.index, tail.values, width=1.0, color="#ee6677")
        axes[1].set_yscale(yscale)
        axes[1].set_xlabel("split_count")
        axes[1].set_ylabel("# variants" if yscale == "linear"
                          else "# variants (log)")
        axes[1].set_title(f"split_count > 1 (full tail, {yscale} y)")
        axes[1].grid(True, alpha=0.3)
        figs[f"split_count{suffix}"] = _save(
            fig, fig_dir, f"01_split_count{suffix}.png"
        )

    # F2: contribution histogram
    light_max_contrib = int(light._contribution.max()) if len(light) else 0
    bins = np.logspace(3, 12, 40)
    for yscale, suffix in [("log", ""), ("linear", "_linear")]:
        fig, ax = plt.subplots(figsize=(10, 5))
        _hist_log_x(df._contribution, ax, "all variants", "#4477aa", bins)
        _hist_log_x(heavy._contribution, ax, "heavy", "#ee6677", bins)
        _hist_log_x(heavy_s._contribution, ax,
                    "heavy & split>1", "#228833", bins, alpha=0.85)
        if light_max_contrib > 0:
            ax.axvline(light_max_contrib, color="#cc3311",
                       linestyle="--", linewidth=1, alpha=0.8,
                       label=f"max(light) = {_fmt_b(light_max_contrib).strip()}")
        ax.axvline(TARGET_HEAVY_PARTITION_BYTES, color="#222222",
                   linestyle=":", linewidth=1.2, alpha=0.8,
                   label="TARGET = 500 MB")
        ax.set_xscale("log"); ax.set_yscale(yscale)
        ax.set_xlabel("_contribution (bytes, log)")
        ax.set_ylabel("# variants" if yscale == "linear"
                      else "# variants (log)")
        ax.set_title(f"Per-variant _contribution ({yscale} y)")
        ax.legend(loc="upper right", fontsize=8); ax.grid(True, alpha=0.3)
        figs[f"contrib_hist{suffix}"] = _save(
            fig, fig_dir, f"02_contribution_hist{suffix}.png"
        )

    # F3: bytes histogram
    bins = np.logspace(6, 14, 40)
    for yscale, suffix in [("log", ""), ("linear", "_linear")]:
        fig, ax = plt.subplots(figsize=(10, 5))
        _hist_log_x(df._bytes, ax, "all variants", "#4477aa", bins)
        _hist_log_x(heavy._bytes, ax, "heavy", "#ee6677", bins)
        _hist_log_x(heavy_s._bytes, ax,
                    "heavy & split>1", "#228833", bins, alpha=0.85)
        ax.set_xscale("log"); ax.set_yscale(yscale)
        ax.set_xlabel("_bytes (log)")
        ax.set_ylabel("# variants" if yscale == "linear"
                      else "# variants (log)")
        ax.set_title(f"Per-variant _bytes ({yscale} y)")
        ax.legend(loc="upper left", fontsize=8); ax.grid(True, alpha=0.3)
        figs[f"bytes_hist{suffix}"] = _save(
            fig, fig_dir, f"03_bytes_hist{suffix}.png"
        )

    # F9: concentration
    sorted_contrib = np.sort(df._contribution.values)[::-1]
    cum = np.cumsum(sorted_contrib) / sorted_contrib.sum()
    xs = np.arange(1, len(cum) + 1) / len(cum)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xs * 100, cum * 100, color="#cc3311", linewidth=1.5)
    ax.plot([0, 100], [0, 100], "k--", alpha=0.5, label="uniform")
    ax.set_xlabel("Top-X% of variants")
    ax.set_ylabel("Cumulative % of total _contribution")
    ax.set_title("Concentration of _contribution")
    ax.grid(True, alpha=0.3); ax.legend()
    concentration_points = {}
    for pct in [0.1, 0.5, 1, 2, 5, 10, 25, 50]:
        idx = max(1, int(pct * len(cum) / 100)) - 1
        v = float(cum[idx] * 100)
        concentration_points[f"top_{pct}pct"] = v
        if pct in (1, 5, 10, 25):
            ax.plot(pct, v, "o", color="#222222", markersize=6)
            ax.annotate(f"top {pct}% → {v:.1f}%", (pct, v),
                        fontsize=8, xytext=(7, -3),
                        textcoords="offset points")
    figs["concentration"] = _save(fig, fig_dir, "09_concentration.png")

    # F11: cutoff sweep
    if has_gene_id and gene is not None:
        g_max = df_x.groupby("gene_id")._contribution.max().values
        sweep_xs = np.logspace(4, 12, 100)
        n_heavy_curve = np.array([(g_max >= x).sum() for x in sweep_xs])
        n_split_curve = np.array([
            (g_max > max(x, TARGET_HEAVY_PARTITION_BYTES)).sum()
            for x in sweep_xs
        ])
        fig, ax = plt.subplots(figsize=(11, 6))
        ax.plot(sweep_xs, n_heavy_curve, color="#ee6677",
                linewidth=2, label="genes with ≥1 heavy variant")
        ax.plot(sweep_xs, n_split_curve, color="#228833",
                linewidth=2, label="genes with ≥1 heavy split>1 variant")
        markers = [
            (f"cutoff in use ({_fmt_b(heavy_contribution_cutoff).strip()})",
             heavy_contribution_cutoff, "#cc3311"),
            ("TARGET (500 MB)", TARGET_HEAVY_PARTITION_BYTES, "#222222"),
            ("exclude 30 GB", 30 * 1024**3, "#aa3377"),
            ("exclude 50 GB", 50 * 1024**3, "#aa3377"),
        ]
        for txt, x, c in markers:
            ax.axvline(x, color=c, linestyle="--", linewidth=1, alpha=0.7)
            n_h = (g_max >= x).sum()
            ax.annotate(
                f"{txt}\n{n_h}/{len(gene)} "
                f"({100*n_h/len(gene):.0f}%)",
                (x, max(n_h, 1)), fontsize=8, color=c,
                xytext=(5, 5), textcoords="offset points",
            )
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("cutoff on max gene _contribution (log)")
        ax.set_ylabel("# genes (log)")
        ax.set_title("Genes with heavy variants vs cutoff")
        ax.legend(loc="lower left"); ax.grid(True, alpha=0.3)
        figs["cutoff_sweep"] = _save(
            fig, fig_dir, "11_cutoff_sweep.png"
        )

        # F8: top genes bar chart
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        top = gene.nlargest(15, "sum_contribution")
        axes[0].barh(range(len(top)), top.sum_contribution / (1024**3),
                     color="#4477aa")
        axes[0].set_yticks(range(len(top)))
        axes[0].set_yticklabels(top.display, fontsize=8)
        axes[0].invert_yaxis()
        axes[0].set_xlabel("Σ _contribution (GB)")
        axes[0].set_title("Top 15 genes by Σ _contribution")
        axes[0].set_xscale("log"); axes[0].grid(True, alpha=0.3, axis="x")
        top = gene.nlargest(15, "sum_bytes")
        axes[1].barh(range(len(top)), top.sum_bytes / (1024**4),
                     color="#ee6677")
        axes[1].set_yticks(range(len(top)))
        axes[1].set_yticklabels(top.display, fontsize=8)
        axes[1].invert_yaxis()
        axes[1].set_xlabel("Σ _bytes (TB)")
        axes[1].set_title("Top 15 genes by Σ _bytes")
        axes[1].set_xscale("log"); axes[1].grid(True, alpha=0.3, axis="x")
        figs["top_genes"] = _save(fig, fig_dir, "08_top_genes.png")

    # Stats dump
    stats = {
        "heavy_contribution_cutoff": int(heavy_contribution_cutoff),
        "TARGET_HEAVY_PARTITION_BYTES": TARGET_HEAVY_PARTITION_BYTES,
        "totals": {
            "n_variants": n_total,
            "n_heavy": int(len(heavy)),
            "n_light": int(len(light)),
            "n_heavy_split_gt1": int(len(heavy_s)),
            "sum_contribution": total_contrib,
            "sum_bytes": total_bytes,
            "heavy_sum_contribution": heavy_contrib,
            "heavy_sum_bytes": heavy_bytes,
            "heavy_split_sum_contribution": heavy_s_contrib,
            "heavy_split_sum_bytes": heavy_s_bytes,
            "heavy_share_contribution_pct":
                100 * heavy_contrib / total_contrib if total_contrib else 0,
            "heavy_share_bytes_pct":
                100 * heavy_bytes / total_bytes if total_bytes else 0,
        },
        "split_count": {
            "n_split_eq_1": int((df.split_count == 1).sum()),
            "n_split_gt_1": int((df.split_count > 1).sum()),
            "max_observed": int(df.split_count.max()),
        },
        "per_variant_contribution": {
            "all": _quantiles(df._contribution),
            "heavy": _quantiles(heavy._contribution) if len(heavy) else {},
            "heavy_split_gt1":
                _quantiles(heavy_s._contribution) if len(heavy_s) else {},
        },
        "per_variant_bytes": {
            "all": _quantiles(df._bytes),
            "heavy": _quantiles(heavy._bytes) if len(heavy) else {},
            "heavy_split_gt1":
                _quantiles(heavy_s._bytes) if len(heavy_s) else {},
        },
        "concentration": concentration_points,
        "figures": figs,
    }
    if gene is not None:
        stats["per_gene"] = {
            "n_genes": int(len(gene)),
            "n_genes_with_heavy": int((gene.n_heavy_variants > 0).sum()),
            "n_genes_with_heavy_split_gt1":
                int((gene.n_heavy_split_variants > 0).sum()),
        }
        stats["top_genes"] = {
            "by_sum_contribution": gene.nlargest(15, "sum_contribution")[
                ["gene_id", "symbol", "n_variants", "n_heavy_variants",
                 "sum_contribution", "sum_bytes", "max_contribution",
                 "max_split_count"]
            ].to_dict(orient="records"),
            "by_sum_bytes": gene.nlargest(15, "sum_bytes")[
                ["gene_id", "symbol", "n_variants", "n_heavy_variants",
                 "sum_contribution", "sum_bytes", "max_bytes",
                 "max_split_count"]
            ].to_dict(orient="records"),
            "by_heavy_split_sum_contribution": gene[
                gene.heavy_split_sum_contribution > 0
            ].nlargest(15, "heavy_split_sum_contribution")[
                ["gene_id", "symbol", "n_variants",
                 "n_heavy_split_variants",
                 "heavy_split_sum_contribution",
                 "heavy_split_max_split_count"]
            ].to_dict(orient="records"),
        }

    def _coerce(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, dict):
            return {k: _coerce(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_coerce(x) for x in o]
        return o

    with open(f"{local_dir}/stats.json", "w") as f:
        json.dump(_coerce(stats), f, indent=2, default=str)

    # Markdown report
    md = []
    md.append("# Variant Size-Info Report\n")
    md.append(
        f"Heavy cutoff used: **{_fmt_b(heavy_contribution_cutoff).strip()}**"
        f" (heavy iff `_contribution >= cutoff`).\n"
    )
    md.append("## Totals\n")
    md.append("| metric | value |\n|---|---|")
    md.append(f"| variants | {n_total:,} |")
    md.append(
        f"| heavy | {len(heavy):,} "
        f"({100*len(heavy)/n_total:.1f}%) |"
    )
    md.append(
        f"| heavy & split>1 | {len(heavy_s):,} "
        f"({100*len(heavy_s)/n_total:.2f}%) |"
    )
    md.append(f"| Σ `_contribution` | {_fmt_b(total_contrib)} |")
    md.append(f"| Σ `_bytes` | {_fmt_b(total_bytes)} |")
    if total_contrib:
        md.append(
            f"| heavy share of contribution | "
            f"{100*heavy_contrib/total_contrib:.2f}% |"
        )
    if total_bytes:
        md.append(
            f"| heavy share of bytes | "
            f"{100*heavy_bytes/total_bytes:.2f}% |"
        )
    md.append("")
    md.append("## Concentration\n")
    md.append(f"![concentration]({figs['concentration']})\n")
    for pct in [0.1, 1, 5, 10, 25]:
        md.append(
            f"- top **{pct}%** of variants → "
            f"**{concentration_points[f'top_{pct}pct']:.1f}%** "
            "of Σ `_contribution`"
        )
    md.append("")
    md.append("## split_count distribution\n")
    md.append(
        f"log y: ![split_count]({figs['split_count']})\n"
        f"linear y: ![split_count linear]({figs['split_count_linear']})\n"
    )
    md.append(
        f"- max observed: **{int(df.split_count.max())}**\n"
        f"- {int((df.split_count > 1).sum()):,} variants with "
        f"`split_count > 1`\n"
    )
    md.append("## Per-variant `_contribution`\n")
    md.append(f"log y: ![contrib]({figs['contrib_hist']})\n")
    md.append(f"linear y: ![contrib lin]({figs['contrib_hist_linear']})\n")
    md.append("## Per-variant `_bytes`\n")
    md.append(f"log y: ![bytes]({figs['bytes_hist']})\n")
    md.append(f"linear y: ![bytes lin]({figs['bytes_hist_linear']})\n")
    if gene is not None:
        md.append("## Top 15 genes by Σ `_contribution`\n")
        md.append(f"![top genes]({figs['top_genes']})\n")
        md.append(
            "| rank | symbol | n_var | n_heavy | Σ contrib | "
            "Σ bytes | max contrib | max split |"
        )
        md.append("|---|---|---|---|---|---|---|---|")
        for i, r in gene.nlargest(15, "sum_contribution").reset_index(
            drop=True
        ).iterrows():
            md.append(
                f"| {i+1} | **{r.display}** | "
                f"{int(r.n_variants):,} | "
                f"{int(r.n_heavy_variants):,} | "
                f"{_fmt_b(r.sum_contribution)} | "
                f"{_fmt_b(r.sum_bytes)} | "
                f"{_fmt_b(r.max_contribution)} | "
                f"{int(r.max_split_count)} |"
            )
        md.append("")
        md.append("## Gene-cutoff sweep\n")
        md.append(f"![sweep]({figs['cutoff_sweep']})\n")
    md.append("")
    md.append(
        "*Source: variant size-info HT. Figures + stats generated by "
        "`v4.size_info_report.build_report`.*"
    )
    with open(f"{local_dir}/report.md", "w") as f:
        f.write("\n".join(md))


def build_report(
    size_info_path: str,
    output_dir: str,
    heavy_contribution_cutoff: int,
    fetch_symbols: bool = True,
) -> str:
    """Read the variant size-info HT, build figures + stats.json +
    ``report.md`` under ``{output_dir}/size_info_report/``, and upload
    everything via :func:`hl.hadoop_copy`.

    :param size_info_path: Path to the variant size-info HT.
    :param output_dir: Directory under which to create
        ``size_info_report/``. Typically a GCS prefix.
    :param heavy_contribution_cutoff: Bytes. Heavy iff
        ``_contribution >= cutoff``.
    :param fetch_symbols: Try to populate gene symbols via mygene.info.
    :return: GCS path of the written ``report.md``.
    """
    import pandas as pd  # noqa: F401 — verify importable up-front

    out_prefix = f"{output_dir.rstrip('/')}/size_info_report"
    tmp_root = tempfile.mkdtemp(prefix="size_info_report_")
    try:
        # 1. Read HT and pull straight to pandas on the driver. The
        # size-info HT row schema is small (~tens of bytes per row); at
        # chr19 scale (~2M rows) this is well under 200 MB on the
        # driver. ht.export(...) on Dataproc writes via the cluster FS
        # (HDFS/GCS), not local /tmp, so to_pandas is also simpler.
        ht = hl.read_table(size_info_path)
        df = ht.to_pandas()
        # Hail's array<str> arrives as a Python list (or numpy array of
        # objects). Normalize to lists.
        if "gene_id" in df.columns:
            df["gene_id"] = df["gene_id"].apply(
                lambda v: list(v) if v is not None else []
            )
        # Drop locus/alleles — we don't use them in the report and they
        # carry Hail-specific Python objects.
        df = df.drop(columns=[c for c in ("locus", "alleles") if c in df.columns])
        # 2. Generate figures + stats + markdown locally.
        _generate(
            df, tmp_root,
            heavy_contribution_cutoff=heavy_contribution_cutoff,
            fetch_symbols=fetch_symbols,
        )
        # 3. Upload everything to the GCS prefix.
        for root, _, files in os.walk(tmp_root):
            for f in files:
                local = os.path.join(root, f)
                rel = os.path.relpath(local, tmp_root)
                remote = f"{out_prefix}/{rel}"
                hl.hadoop_copy(f"file://{local}", remote)
        return f"{out_prefix}/report.md"
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
