import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .config import PROJECT_ROOT, CANDIDATE_GENES

SEED = 20260908
PROLIFERATION_GENES = [
    "MKI67", "PCNA", "CCNB1", "CCNE1", "MCM2", "MCM4",
    "MCM6", "TOP2A", "BIRC5", "CDC20", "CDC6",
]
DDR_ANNOTATION = {
    "ATR", "ATRIP", "CHEK1", "CHEK2", "WEE1", "PARP1", "PARP2", "PARP3",
    "BRCA1", "BRCA2", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D",
    "POLQ", "ATM", "MRE11", "NBN", "RAD50", "RPA1", "RPA2", "RPA3",
    "FANCD2", "FANCI", "FANCA", "FANCC", "FANCL", "BRIP1", "BARD1",
    "TOPBP1", "CLSPN", "TIMELESS", "TIPIN", "RRM1", "RRM2", "PRKDC",
}


def _symbol(column):
    return re.sub(r"\s*\([^()]+\)\s*$", "", str(column)).strip()


def _bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    out[ok] = restored
    return out


def resource_status(root=PROJECT_ROOT):
    config = json.loads((root / "config/data_sources.yaml").read_text())
    rows = []
    needed = [
        "depmap_effect", "depmap_model", "depmap_expression", "depmap_default",
        "depmap_profiles", "depmap_condition", "depmap_readme", "depmap_probability",
        "prism_response", "prism_treatment", "prism_models", "prism_readme",
        "gdsc_response", "gdsc_models", "gdsc_compounds", "gpl96",
    ]
    for r in config["resources"]:
        if r["id"] not in needed:
            continue
        local = root / "data/external" / r["id"] / r["filename"]
        rows.append({
            "resource_id": r["id"],
            "file_id": r["url"].rsplit("/", 1)[-1],
            "filename": r["filename"],
            "size_bytes": r["expected_bytes"],
            "url": r["url"],
            "expected_local_path": str(local),
            "status": "present" if local.is_file() else "missing",
            "role": "optional probability sensitivity" if r["id"] == "depmap_probability" else "required within component or alternative drug resource",
        })
    return pd.DataFrame(rows)


def _row_id_to_model_ids(root, row_ids):
    ids = pd.Series(row_ids, dtype="string")
    if ids.str.startswith("ACH-").mean() > 0.9:
        return ids.astype(str)
    default_path = root / "data/external/depmap_default/OmicsDefaultModelProfiles.csv"
    profiles_path = root / "data/external/depmap_profiles/OmicsProfiles.csv"
    mapping = {}
    if default_path.is_file():
        d = pd.read_csv(default_path, usecols=["ModelID", "ProfileID"])
        mapping.update(dict(zip(d.ProfileID.astype(str), d.ModelID.astype(str))))
    if profiles_path.is_file():
        p = pd.read_csv(profiles_path, usecols=["ProfileID", "ModelID"])
        p = p.dropna().drop_duplicates("ProfileID")
        mapping.update(dict(zip(p.ProfileID.astype(str), p.ModelID.astype(str))))
    return ids.map(mapping)


def _curated_model_sets(model):
    text_cols = [
        c for c in [
            "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype",
            "PatientSubtypeFeatures", "ModelSubtypeFeatures", "PublicComments",
            "TissueOrigin", "CellLineName", "CCLEName",
        ] if c in model.columns
    ]
    joined = model[text_cols].fillna("").astype(str).agg(" | ".join, axis=1)
    breast = joined.str.contains(r"breast|mammary", case=False, regex=True, na=False)
    # Curated TNBC metadata
    tnbc = breast & joined.str.contains(r"\btnbc\b|triple[ -]?negative", case=False, regex=True, na=False)
    return model.loc[breast, "ModelID"].astype(str).tolist(), model.loc[tnbc, "ModelID"].astype(str).tolist(), text_cols


def _read_expression_predictors(root, eligible_ids):
    path = root / "data/external/depmap_expression/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    header = pd.read_csv(path, nrows=0).columns.tolist()
    id_col = header[0]
    wanted = set(CANDIDATE_GENES) | set(PROLIFERATION_GENES)
    gene_cols = [c for c in header[1:] if _symbol(c) in wanted]
    expr = pd.read_csv(path, usecols=[id_col] + gene_cols)
    expr["ModelID"] = _row_id_to_model_ids(root, expr[id_col])
    expr = expr[expr.ModelID.isin(set(eligible_ids))].drop_duplicates("ModelID").set_index("ModelID")

    by_symbol = {}
    for c in gene_cols:
        by_symbol.setdefault(_symbol(c), []).append(c)
    collapsed = pd.DataFrame(index=expr.index)
    for gene, cols in by_symbol.items():
        collapsed[gene] = expr[cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)

    cand_available = [g for g in CANDIDATE_GENES if g in collapsed.columns]
    cand_missing = [g for g in CANDIDATE_GENES if g not in collapsed.columns]
    if cand_missing:
        return collapsed, None, None, cand_available, cand_missing

    z = (collapsed[list(CANDIDATE_GENES)] - collapsed[list(CANDIDATE_GENES)].mean()) / collapsed[list(CANDIDATE_GENES)].std(ddof=0).replace(0, np.nan)
    signature = z.mean(axis=1)

    prolif_available = [g for g in PROLIFERATION_GENES if g in collapsed.columns]
    if len(prolif_available) >= 5:
        pz = (collapsed[prolif_available] - collapsed[prolif_available].mean()) / collapsed[prolif_available].std(ddof=0).replace(0, np.nan)
        proliferation = pz.mean(axis=1)
    else:
        proliferation = pd.Series(np.nan, index=collapsed.index, name="proliferation")
    return collapsed, signature, proliferation, cand_available, cand_missing


def _read_gene_effect(root, eligible_ids):
    path = root / "data/external/depmap_effect/CRISPRGeneEffect.csv"
    wanted = set(eligible_ids)
    pieces = []
    # Bounded matrix chunks
    for chunk in pd.read_csv(path, chunksize=64):
        id_col = chunk.columns[0]
        model_ids = _row_id_to_model_ids(root, chunk[id_col])
        keep = model_ids.isin(wanted)
        if keep.any():
            sub = chunk.loc[keep].copy()
            sub["ModelID"] = model_ids.loc[keep].values
            pieces.append(sub)
    if not pieces:
        return pd.DataFrame()
    ge = pd.concat(pieces, ignore_index=True).drop_duplicates("ModelID").set_index("ModelID")
    first = ge.columns[0] if len(ge.columns) else None
    if first is not None and first.startswith("Unnamed"):
        ge = ge.drop(columns=[first])
    return ge


def _fit_gene(y, x, cov=None):
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    if cov is None:
        ok = np.isfinite(y) & np.isfinite(x)
        X = np.column_stack([np.ones(ok.sum()), x[ok]])
    else:
        c = np.asarray(cov, float)
        ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(c)
        X = np.column_stack([np.ones(ok.sum()), x[ok], c[ok]])
    yy = y[ok]
    n, k = X.shape
    if n <= k + 2 or np.nanstd(x[ok]) == 0 or np.nanstd(yy) == 0:
        return n, np.nan, np.nan, np.nan, np.nan, np.nan
    try:
        inv = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return n, np.nan, np.nan, np.nan, np.nan, np.nan
    b = inv @ X.T @ yy
    resid = yy - X @ b
    df = n - k
    mse = float(resid @ resid / df)
    se = float(np.sqrt(max(inv[1, 1] * mse, 0)))
    beta = float(b[1])
    t = beta / se if se > 0 else np.nan
    p = float(2 * stats.t.sf(abs(t), df)) if np.isfinite(t) else np.nan
    ci = stats.t.ppf(0.975, df) * se if np.isfinite(se) else np.nan
    return n, beta, se, beta - ci, beta + ci, p


def _genome_wide_screen(gene_effect, signature, proliferation):
    ids = gene_effect.index.intersection(signature.dropna().index)
    ge = gene_effect.loc[ids]
    x = signature.loc[ids]
    pscore = proliferation.reindex(ids)
    records = []
    for col in ge.columns:
        y = pd.to_numeric(ge[col], errors="coerce").values
        n, beta, se, lo, hi, p = _fit_gene(y, x.values)
        na, ba, sea, loa, hia, pa = _fit_gene(y, x.values, pscore.values) if pscore.notna().sum() >= 8 else (n, np.nan, np.nan, np.nan, np.nan, np.nan)
        records.append({
            "dependency": _symbol(col), "depmap_column": col, "N": n,
            "beta": beta, "se": se, "ci_low": lo, "ci_high": hi, "p": p,
            "N_adjusted": na, "beta_adjusted": ba, "se_adjusted": sea,
            "ci_low_adjusted": loa, "ci_high_adjusted": hia, "p_adjusted": pa,
        })
    out = pd.DataFrame(records)
    out["fdr"] = _bh(out.p.values)
    out["fdr_adjusted"] = _bh(out.p_adjusted.values)
    out["stronger_dependency_high_signature"] = out.beta_adjusted.lt(0) if out.p_adjusted.notna().any() else out.beta.lt(0)
    out["annotation"] = np.where(out.dependency.isin(DDR_ANNOTATION), "DDR/replication-stress", "")
    return out.sort_values(["fdr_adjusted", "fdr", "p_adjusted", "p"], na_position="last").reset_index(drop=True)


def _make_depmap_figures(root, results, summary):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    figdir = root / "outputs/v2/figures/stage3"
    figdir.mkdir(parents=True, exist_ok=True)
    if results.empty:
        return
    pcol = "p_adjusted" if results.p_adjusted.notna().any() else "p"
    bcol = "beta_adjusted" if results.p_adjusted.notna().any() else "beta"
    d = results[np.isfinite(results[pcol]) & np.isfinite(results[bcol])].copy()
    if d.empty:
        return
    d["neglog10p"] = -np.log10(np.clip(d[pcol], 1e-300, 1))
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(d[bcol], d.neglog10p, s=10, alpha=.55)
    sig = d[(d.fdr_adjusted < .05) & (d[bcol] < 0)] if "fdr_adjusted" in d else d.iloc[0:0]
    for _, r in sig.head(12).iterrows():
        ax.annotate(r.dependency, (r[bcol], r.neglog10p), fontsize=7)
    ax.axvline(0, linewidth=.8)
    ax.set_xlabel("Gene-effect regression beta per 1 SD candidate signature")
    ax.set_ylabel("-log10(P)")
    ax.set_title("DepMap genome-wide CRISPR dependency screen")
    fig.tight_layout()
    fig.savefig(figdir / "depmap_volcano.pdf")
    fig.savefig(figdir / "depmap_volcano.png", dpi=300)
    plt.close(fig)

    top = d.sort_values("fdr_adjusted" if d.fdr_adjusted.notna().any() else "p").head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.errorbar(top[bcol], np.arange(len(top)), xerr=[top[bcol] - top["ci_low_adjusted" if bcol == "beta_adjusted" else "ci_low"], top["ci_high_adjusted" if bcol == "beta_adjusted" else "ci_high"] - top[bcol]], fmt="o")
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top.dependency)
    ax.axvline(0, linewidth=.8)
    ax.set_xlabel("Regression beta; negative = stronger dependency at higher signature")
    ax.set_title("Top DepMap dependency associations")
    fig.tight_layout()
    fig.savefig(figdir / "depmap_top_dependencies.pdf")
    fig.savefig(figdir / "depmap_top_dependencies.png", dpi=300)
    plt.close(fig)


def functional_results(root=PROJECT_ROOT):
    dest = root / "outputs/v2/stage3"
    dest.mkdir(parents=True, exist_ok=True)
    resources = resource_status(root)
    resources.to_csv(dest / "missing_resources.csv", index=False)

    dep_ids = ["depmap_effect", "depmap_model", "depmap_expression"]
    dep_missing = resources[resources.resource_id.isin(dep_ids) & (resources.status == "missing")]
    if not dep_missing.empty:
        dep = {
            "status": "unavailable", "reason": "Required DepMap files absent: " + ",".join(dep_missing.resource_id),
            "breast_model_n": None, "tnbc_model_n": None, "analysis_model_n": None,
            "FDR_significant_dependencies": None, "strongest_dependency": None,
        }
        pd.DataFrame([dep]).to_csv(dest / "depmap_results.csv", index=False)
    else:
        model = pd.read_csv(root / "data/external/depmap_model/Model.csv")
        breast_ids, tnbc_ids, metadata_fields = _curated_model_sets(model)
        # Curated TNBC eligibility
        primary_ids = tnbc_ids if len(tnbc_ids) >= 8 else breast_ids
        cohort_label = "curated_TNBC" if len(tnbc_ids) >= 8 else "breast_lineage_exploratory"
        expr, signature, proliferation, cand_avail, cand_missing = _read_expression_predictors(root, primary_ids)
        if signature is None:
            dep = {
                "status": "unavailable", "reason": "Candidate genes missing from DepMap expression: " + ",".join(cand_missing),
                "breast_model_n": len(breast_ids), "tnbc_model_n": len(tnbc_ids), "analysis_model_n": 0,
                "analysis_cohort": cohort_label, "candidate_genes_available": cand_avail,
                "candidate_genes_missing": cand_missing, "FDR_significant_dependencies": None,
                "strongest_dependency": None,
            }
            pd.DataFrame([dep]).to_csv(dest / "depmap_results.csv", index=False)
        else:
            ge = _read_gene_effect(root, signature.dropna().index.tolist())
            common = ge.index.intersection(signature.dropna().index)
            signature = signature.loc[common]
            proliferation = proliferation.reindex(common)
            ge = ge.loc[common]
            screen = _genome_wide_screen(ge, signature, proliferation)
            screen.to_csv(dest / "depmap_results.csv", index=False)
            adjusted_available = screen.p_adjusted.notna().any()
            if adjusted_available:
                sig_hits = screen[(screen.fdr_adjusted < .05) & (screen.beta_adjusted < 0)]
                ranking = screen.sort_values(["fdr_adjusted", "p_adjusted"], na_position="last")
            else:
                sig_hits = screen[(screen.fdr < .05) & (screen.beta < 0)]
                ranking = screen.sort_values(["fdr", "p"], na_position="last")
            strongest = None
            if len(ranking):
                r = ranking.iloc[0]
                strongest = {
                    "gene": r.dependency,
                    "beta": float(r.beta_adjusted if adjusted_available else r.beta),
                    "p": float(r.p_adjusted if adjusted_available else r.p),
                    "fdr": float(r.fdr_adjusted if adjusted_available else r.fdr),
                    "annotation": r.annotation,
                }
            dep = {
                "status": "estimated",
                "analysis_cohort": cohort_label,
                "breast_model_n": int(len(set(breast_ids))),
                "tnbc_model_n": int(len(set(tnbc_ids))),
                "analysis_model_n": int(len(common)),
                "metadata_fields_used_for_curated_eligibility": metadata_fields,
                "tnbc_definition": "literal TNBC/triple-negative annotation in curated non-expression model metadata only",
                "expression_receptor_thresholds_used": False,
                "candidate_genes_available": cand_avail,
                "candidate_genes_missing": cand_missing,
                "proliferation_genes_available": [g for g in PROLIFERATION_GENES if g in expr.columns],
                "FDR_significant_dependencies": int(len(sig_hits)),
                "strongest_dependency": strongest,
                "DDR_FDR_significant_dependencies": sig_hits[sig_hits.dependency.isin(DDR_ANNOTATION)].dependency.astype(str).tolist(),
                "effect_direction": "negative beta means stronger CRISPR dependency with higher candidate signature",
                "matched_random_signature_control": "not run in this implementation; genome-wide FDR screen retained as the prespecified functional test",
            }
            (dest / "depmap_summary.json").write_text(json.dumps(dep, indent=2), encoding="utf-8")
            _make_depmap_figures(root, screen, dep)

    prism = ["prism_response", "prism_models", "prism_treatment"]
    gdsc = ["gdsc_response", "gdsc_models", "gdsc_compounds"]
    prism_ready = resources[resources.resource_id.isin(prism)].status.eq("present").all()
    gdsc_ready = resources[resources.resource_id.isin(gdsc)].status.eq("present").all()
    if prism_ready or gdsc_ready:
        drug = {
            "status": "unavailable",
            "reason": "Drug files are present but compound/model crosswalk workflow has not been implemented; no drug association claimed",
            "FDR_significant_compounds": None,
        }
    else:
        drug = {
            "status": "unavailable",
            "reason": "No complete local PRISM or GDSC response/metadata bundle",
            "FDR_significant_compounds": None,
        }
    pd.DataFrame([drug]).to_csv(dest / "drug_results.csv", index=False)
    return dep, drug, resources


def stopped_results(root=PROJECT_ROOT):
    return functional_results(root)
