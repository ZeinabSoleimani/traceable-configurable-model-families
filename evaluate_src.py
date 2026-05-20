#!/usr/bin/env python3
"""
Evaluation runner for evidence_cpt / VACPD-style configurable process-tree discovery.

It evaluates the proposed package against the two essential discovery baselines
for variant-labeled logs:

  1. proposed_cpt      : evidence-guided configurable process tree, projected per variant
  2. pooled_im         : one process tree discovered from the union/pooled log
  3. separate_im       : one process tree discovered independently per variant
  4. separate_exact    : optional exact-log memorization baseline / upper-bound coverage check

Metrics written:
  - observed-trace coverage on train and test splits (fitness proxy)
  - bounded directly-follows precision/recall proxy
  - model-family size and simplicity metrics
  - compactness gain against separate per-variant models
  - shared/variant-specific node ratios for the proposed configurable tree
  - C-DFG evidence counts and residual-localization proxy where applicable
  - runtime

Run from the repository root, or set PYTHONPATH to the directory containing evidence_cpt.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import random
import statistics
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from evidence_cpt.cdfg import CDFG, SINK, SOURCE, directly_follows_counts
from evidence_cpt.data import Case, Trace, VariantLog, depression_example_log
from evidence_cpt.discovery import DiscoveryConfig, exact_log_tree, mine_process_tree
from evidence_cpt.enrichment import EnrichmentConfig, EnrichmentResult, enrich_process_tree
from evidence_cpt.metrics import residual_relation_localization
from evidence_cpt.ptree import ACT, LOOP, PAR, SEQ, TAU, XOR, ProcessTree
try:
    # Internal helper used only for bounded evaluation on large loop-heavy logs.
    from evidence_cpt.ptree import _accept_trace
except Exception:  # pragma: no cover
    _accept_trace = None


@dataclass(frozen=True)
class DatasetSpec:
    name: str
    path: Path | None
    fmt: str


@dataclass(frozen=True)
class MethodResult:
    method: str
    family: Mapping[str, ProcessTree]
    family_nodes: int
    family_visible_leaves: int
    family_operator_nodes: int
    family_max_depth: int
    runtime_sec: float
    proposed_result: EnrichmentResult | None = None
    combined_tree: ProcessTree | None = None


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    datasets = resolve_datasets(args)
    if not datasets:
        print("No input datasets found. Use --inputs or --example.", file=sys.stderr)
        return 2

    discovery_cfg = DiscoveryConfig(
        max_depth=args.max_depth,
        collapse_consecutive_repetitions=args.discovery_collapse_repetitions,
        max_exact_log_traces=args.max_exact_log_traces,
        max_exact_log_leaves=args.max_exact_log_leaves,
        use_flower_fallback=not args.no_flower_fallback,
    )
    enrichment_cfg = EnrichmentConfig(
        theta=args.theta,
        alpha=args.alpha,
        n_min=args.n_min,
        min_shared_activities=args.min_shared_activities,
        collapse_repeated_shared_activities=not args.no_collapse_repeated_shared_activities,
        discovery=discovery_cfg,
        compact_return_cycles=not args.no_compact_return_cycles,
    )

    per_variant_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    cdfg_rows: list[dict[str, object]] = []

    for spec in datasets:
        print(f"\n=== Dataset: {spec.name} ===")
        full_log = load_log(spec, args)
        print(f"cases={len(full_log.cases)}, variants={', '.join(full_log.variants)}, activities={len(full_log.activities)}")

        for seed in args.seeds:
            train_log, test_log = stratified_split(full_log, test_ratio=args.test_ratio, seed=seed)
            print(
                f"  seed={seed}: train={len(train_log.cases)}"
                + (f", test={len(test_log.cases)}" if test_log is not None else ", test=<none>")
            )

            method_results = run_methods(
                train_log,
                discovery_cfg=discovery_cfg,
                enrichment_cfg=enrichment_cfg,
                include_exact_baseline=args.include_exact_baseline,
            )

            # Used for compactness normalization.
            separate_nodes = method_results["separate_im"].family_nodes
            pooled_nodes = method_results["pooled_im"].family_nodes

            for method_name, method_result in method_results.items():
                train_eval = evaluate_family(
                    method_result.family,
                    train_log,
                    split="train",
                    max_loop_depth=args.max_loop_depth,
                    lang_loop_depth=args.lang_loop_depth,
                    max_lang_traces=args.max_lang_traces,
                    coverage_mode=args.coverage_mode,
                    coverage_max_cases_per_variant=args.coverage_max_cases_per_variant,
                    skip_df_proxy=args.skip_df_proxy,
                )
                test_eval = (
                    evaluate_family(
                        method_result.family,
                        test_log,
                        split="test",
                        max_loop_depth=args.max_loop_depth,
                        lang_loop_depth=args.lang_loop_depth,
                        max_lang_traces=args.max_lang_traces,
                        coverage_mode=args.coverage_mode,
                        coverage_max_cases_per_variant=args.coverage_max_cases_per_variant,
                        skip_df_proxy=args.skip_df_proxy,
                    )
                    if test_log is not None
                    else None
                )

                for row in train_eval["per_variant"]:
                    row.update(dataset=spec.name, seed=seed, method=method_name)
                    per_variant_rows.append(row)
                if test_eval is not None:
                    for row in test_eval["per_variant"]:
                        row.update(dataset=spec.name, seed=seed, method=method_name)
                        per_variant_rows.append(row)

                proposed_specific = proposed_tree_stats(method_result.proposed_result)
                loc_stats = localization_stats(method_result.proposed_result)

                summary_row = {
                    "dataset": spec.name,
                    "seed": seed,
                    "method": method_name,
                    "train_cases": len(train_log.cases),
                    "test_cases": len(test_log.cases) if test_log is not None else 0,
                    "coverage_mode": args.coverage_mode,
                    "coverage_max_cases_per_variant": args.coverage_max_cases_per_variant,
                    "max_loop_depth": args.max_loop_depth,
                    "skip_df_proxy": args.skip_df_proxy,
                    "train_coverage": train_eval["overall_coverage"],
                    "test_coverage": test_eval["overall_coverage"] if test_eval is not None else "",
                    "coverage_gap_train_minus_test": coverage_gap(train_eval["overall_coverage"], test_eval["overall_coverage"] if test_eval is not None else ""),
                    "mean_train_df_precision_proxy": train_eval["mean_df_precision"],
                    "mean_train_df_recall_proxy": train_eval["mean_df_recall"],
                    "mean_test_df_precision_proxy": test_eval["mean_df_precision"] if test_eval is not None else "",
                    "mean_test_df_recall_proxy": test_eval["mean_df_recall"] if test_eval is not None else "",
                    "family_nodes": method_result.family_nodes,
                    "family_visible_leaves": method_result.family_visible_leaves,
                    "family_operator_nodes": method_result.family_operator_nodes,
                    "family_max_depth": method_result.family_max_depth,
                    "compactness_gain_vs_separate_im": safe_ratio(separate_nodes - method_result.family_nodes, separate_nodes),
                    "size_ratio_vs_pooled_im": safe_ratio(method_result.family_nodes, pooled_nodes),
                    "runtime_sec": method_result.runtime_sec,
                    **proposed_specific,
                    **loc_stats,
                }
                summary_rows.append(summary_row)

            proposed = method_results["proposed_cpt"].proposed_result
            if proposed is not None:
                cdfg_rows.append({"dataset": spec.name, "seed": seed, **cdfg_stats(proposed.cdfg)})

    write_csv(out / "per_variant_metrics.csv", per_variant_rows)
    write_csv(out / "summary_metrics.csv", summary_rows)
    write_csv(out / "cdfg_metrics.csv", cdfg_rows)
    write_markdown_summary(out / "README_metrics.md", summary_rows, cdfg_rows)

    print(f"\nWrote evaluation files to: {out.resolve()}")
    print("  - summary_metrics.csv")
    print("  - per_variant_metrics.csv")
    print("  - cdfg_metrics.csv")
    print("  - README_metrics.md")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Evaluate evidence_cpt against discovery baselines.")
    p.add_argument("--inputs", nargs="*", default=[], help="CSV files or glob patterns, e.g. examples/_prepared_variant_runs/*.csv")
    p.add_argument("--example", action="store_true", help="also evaluate the built-in depression example")
    p.add_argument("--format", choices=["auto", "trace", "event"], default="auto")
    p.add_argument("--out", required=True, help="output directory for evaluation CSV files")

    p.add_argument("--case-col", default="case_id")
    p.add_argument("--variant-col", default="variant")
    p.add_argument("--trace-col", default="trace")
    p.add_argument("--trace-sep", default=",")
    p.add_argument("--activity-col", default="activity")
    p.add_argument("--timestamp-col", default="timestamp")

    p.add_argument("--test-ratio", type=float, default=0.20, help="stratified hold-out ratio; use 0 for train-only evaluation")
    p.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])

    p.add_argument("--theta", type=float, default=0.0)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--n-min", type=int, default=5)
    p.add_argument("--min-shared-activities", type=int, default=2)
    p.add_argument("--no-collapse-repeated-shared-activities", action="store_true")
    p.add_argument("--no-compact-return-cycles", action="store_true")

    p.add_argument("--max-depth", type=int, default=100)
    p.add_argument("--discovery-collapse-repetitions", action="store_true")
    p.add_argument("--max-exact-log-traces", type=int, default=100)
    p.add_argument("--max-exact-log-leaves", type=int, default=500)
    p.add_argument("--no-flower-fallback", action="store_true")

    p.add_argument("--include-exact-baseline", action="store_true", help="add separate_exact; can be large on real logs")
    p.add_argument("--coverage-mode", choices=["exact", "bounded", "skip"], default="exact",
                   help=("exact uses ProcessTree.accepts_observed_trace; bounded uses a fixed loop bound "
                         "and is recommended for large loop-heavy real logs; skip omits trace-coverage replay."))
    p.add_argument("--coverage-max-cases-per-variant", type=int, default=0,
                   help="0 means all cases; otherwise evaluate coverage on at most this many cases per variant/split")
    p.add_argument("--skip-df-proxy", action="store_true",
                   help="skip bounded language enumeration for the directly-follows precision/recall proxy")
    p.add_argument("--max-loop-depth", type=int, default=3, help="loop depth for observed-trace coverage; exact mode may internally require more")
    p.add_argument("--lang-loop-depth", type=int, default=2, help="loop depth for bounded language DF precision/recall proxy")
    p.add_argument("--max-lang-traces", type=int, default=5000, help="cap for bounded language enumeration")
    return p


def resolve_datasets(args: argparse.Namespace) -> list[DatasetSpec]:
    specs: list[DatasetSpec] = []
    seen: set[Path] = set()
    for pattern in args.inputs:
        matches = sorted(Path(p).resolve() for p in glob.glob(pattern))
        if not matches and Path(pattern).exists():
            matches = [Path(pattern).resolve()]
        for path in matches:
            if path.is_file() and path not in seen:
                seen.add(path)
                specs.append(DatasetSpec(path.stem, path, args.format))
    if args.example:
        specs.insert(0, DatasetSpec("depression_example", None, "trace"))
    return specs


def load_log(spec: DatasetSpec, args: argparse.Namespace) -> VariantLog:
    if spec.path is None:
        return depression_example_log()
    fmt = spec.fmt if spec.fmt != "auto" else detect_format(spec.path, args)
    if fmt == "trace":
        return VariantLog.from_trace_csv(
            spec.path,
            case_col=args.case_col,
            variant_col=args.variant_col,
            trace_col=args.trace_col,
            trace_sep=args.trace_sep,
        )
    return VariantLog.from_event_csv(
        spec.path,
        case_col=args.case_col,
        variant_col=args.variant_col,
        activity_col=args.activity_col,
        timestamp_col=args.timestamp_col,
    )


def detect_format(path: Path, args: argparse.Namespace) -> str:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = next(reader)
    cols = {c.strip() for c in header}
    if args.trace_col in cols:
        return "trace"
    if args.activity_col in cols:
        return "event"
    raise ValueError(
        f"Cannot auto-detect format for {path}. Expected column {args.trace_col!r} or {args.activity_col!r}."
    )


def stratified_split(log: VariantLog, *, test_ratio: float, seed: int) -> tuple[VariantLog, VariantLog | None]:
    if test_ratio <= 0:
        return log, None
    rng = random.Random(seed)
    train_cases: list[Case] = []
    test_cases: list[Case] = []
    for variant in log.variants:
        cases = [case for case in log.cases if case.variant == variant]
        rng.shuffle(cases)
        n = len(cases)
        if n < 3:
            n_test = 0
        else:
            n_test = max(1, round(n * test_ratio))
            n_test = min(n_test, n - 1)
        test_cases.extend(cases[:n_test])
        train_cases.extend(cases[n_test:])
    if not test_cases:
        return VariantLog(tuple(train_cases)), None
    return VariantLog(tuple(train_cases)), VariantLog(tuple(test_cases))


def run_methods(
    train_log: VariantLog,
    *,
    discovery_cfg: DiscoveryConfig,
    enrichment_cfg: EnrichmentConfig,
    include_exact_baseline: bool,
) -> dict[str, MethodResult]:
    results: dict[str, MethodResult] = {}

    t0 = time.perf_counter()
    proposed = enrich_process_tree(train_log, config=enrichment_cfg)
    runtime = time.perf_counter() - t0
    proposed_family = proposed.configured_trees()
    results["proposed_cpt"] = MethodResult(
        method="proposed_cpt",
        family=proposed_family,
        family_nodes=tree_size(proposed.tree),
        family_visible_leaves=visible_leaf_count(proposed.tree),
        family_operator_nodes=operator_node_count(proposed.tree),
        family_max_depth=tree_depth(proposed.tree),
        runtime_sec=runtime,
        proposed_result=proposed,
        combined_tree=proposed.tree,
    )

    t0 = time.perf_counter()
    pooled = mine_process_tree(train_log.all_traces(), config=discovery_cfg).simplified()
    runtime = time.perf_counter() - t0
    pooled_family = {variant: pooled for variant in train_log.variants}
    results["pooled_im"] = MethodResult(
        method="pooled_im",
        family=pooled_family,
        family_nodes=tree_size(pooled),
        family_visible_leaves=visible_leaf_count(pooled),
        family_operator_nodes=operator_node_count(pooled),
        family_max_depth=tree_depth(pooled),
        runtime_sec=runtime,
        combined_tree=pooled,
    )

    t0 = time.perf_counter()
    separate_family = {
        variant: mine_process_tree(train_log.sublog(variant), config=discovery_cfg).simplified()
        for variant in train_log.variants
    }
    runtime = time.perf_counter() - t0
    results["separate_im"] = family_method_result("separate_im", separate_family, runtime)

    if include_exact_baseline:
        t0 = time.perf_counter()
        exact_family = {
            variant: exact_log_tree(train_log.sublog(variant)).simplified()
            for variant in train_log.variants
        }
        runtime = time.perf_counter() - t0
        results["separate_exact"] = family_method_result("separate_exact", exact_family, runtime)

    return results


def family_method_result(method: str, family: Mapping[str, ProcessTree], runtime: float) -> MethodResult:
    return MethodResult(
        method=method,
        family=family,
        family_nodes=sum(tree_size(tree) for tree in family.values()),
        family_visible_leaves=sum(visible_leaf_count(tree) for tree in family.values()),
        family_operator_nodes=sum(operator_node_count(tree) for tree in family.values()),
        family_max_depth=max((tree_depth(tree) for tree in family.values()), default=0),
        runtime_sec=runtime,
    )


def evaluate_family(
    family: Mapping[str, ProcessTree],
    log: VariantLog,
    *,
    split: str,
    max_loop_depth: int,
    lang_loop_depth: int,
    max_lang_traces: int,
    coverage_mode: str,
    coverage_max_cases_per_variant: int,
    skip_df_proxy: bool,
) -> dict[str, object]:
    per_variant: list[dict[str, object]] = []
    total_covered = 0
    total_cases = 0
    precisions: list[float] = []
    recalls: list[float] = []

    for variant in log.variants:
        tree = family.get(variant)
        traces = log.sublog(variant)
        if tree is None:
            covered = 0
            coverage = 0.0 if traces else 1.0
            dfp = dfr = ""
            model_df_count = log_df_count = lang_size = ""
            lang_cap_hit = ""
        else:
            eval_traces = list(traces)
            coverage_eval_total = len(eval_traces)
            if coverage_max_cases_per_variant and coverage_max_cases_per_variant > 0:
                eval_traces = eval_traces[:coverage_max_cases_per_variant]
                coverage_eval_total = len(eval_traces)

            if coverage_mode == "skip":
                covered = ""
                coverage = ""
            elif coverage_mode == "bounded":
                if _accept_trace is None:
                    raise RuntimeError("coverage-mode=bounded requires evidence_cpt.ptree._accept_trace")
                covered = sum(
                    1 for trace in eval_traces
                    if _accept_trace(tree, tuple(trace), max_loop_depth=max_loop_depth)
                )
                coverage = covered / coverage_eval_total if coverage_eval_total else 1.0
            else:
                covered = sum(
                    1 for trace in eval_traces if tree.accepts_observed_trace(trace, max_loop_depth=max_loop_depth)
                )
                coverage = covered / coverage_eval_total if coverage_eval_total else 1.0

            if skip_df_proxy:
                dfp = dfr = ""
                model_df_count = log_df_count = lang_size = ""
                lang_cap_hit = ""
            else:
                df_proxy = bounded_df_precision_recall(
                    tree,
                    traces,
                    max_loop_depth=lang_loop_depth,
                    max_traces=max_lang_traces,
                )
                dfp = df_proxy["df_precision_proxy"]
                dfr = df_proxy["df_recall_proxy"]
                model_df_count = df_proxy["model_df_count"]
                log_df_count = df_proxy["log_df_count"]
                lang_size = df_proxy["bounded_language_size"]
                lang_cap_hit = df_proxy["bounded_language_cap_hit"]
                if isinstance(dfp, float):
                    precisions.append(dfp)
                if isinstance(dfr, float):
                    recalls.append(dfr)

        if isinstance(covered, int):
            total_covered += covered
            total_cases += coverage_eval_total if 'coverage_eval_total' in locals() else len(traces)
        per_variant.append(
            {
                "split": split,
                "variant": variant,
                "covered": covered,
                "total": len(traces),
                "coverage_evaluated_total": (coverage_eval_total if tree is not None and coverage_mode != "skip" else ""),
                "coverage": coverage,
                "df_precision_proxy": dfp,
                "df_recall_proxy": dfr,
                "model_df_count": model_df_count,
                "log_df_count": log_df_count,
                "bounded_language_size": lang_size,
                "bounded_language_cap_hit": lang_cap_hit,
                "variant_model_nodes": tree_size(tree) if tree is not None else "",
                "variant_model_visible_leaves": visible_leaf_count(tree) if tree is not None else "",
                "variant_model_operator_nodes": operator_node_count(tree) if tree is not None else "",
                "variant_model_depth": tree_depth(tree) if tree is not None else "",
            }
        )

    return {
        "per_variant": per_variant,
        "overall_coverage": total_covered / total_cases if total_cases else 1.0,
        "mean_df_precision": statistics.mean(precisions) if precisions else "",
        "mean_df_recall": statistics.mean(recalls) if recalls else "",
    }


def bounded_df_precision_recall(
    tree: ProcessTree,
    traces: Iterable[Trace],
    *,
    max_loop_depth: int,
    max_traces: int,
) -> dict[str, object]:
    observed_traces = [tuple(trace) for trace in traces]
    observed_counts, _ = directly_follows_counts(observed_traces)
    observed_df = set(observed_counts.keys())

    language = tree.enumerate_language(max_loop_depth=max_loop_depth, max_traces=max_traces)
    model_counts, _ = directly_follows_counts(list(language))
    model_df = set(model_counts.keys())

    if not model_df and not observed_df:
        precision = recall = 1.0
    else:
        precision = len(model_df & observed_df) / len(model_df) if model_df else 0.0
        recall = len(model_df & observed_df) / len(observed_df) if observed_df else 1.0

    return {
        "df_precision_proxy": precision,
        "df_recall_proxy": recall,
        "model_df_count": len(model_df),
        "log_df_count": len(observed_df),
        "bounded_language_size": len(language),
        "bounded_language_cap_hit": len(language) >= max_traces,
    }


def tree_size(tree: ProcessTree) -> int:
    return sum(1 for _ in tree.walk())


def visible_leaf_count(tree: ProcessTree) -> int:
    return sum(1 for node in tree.walk() if node.op == ACT)


def operator_node_count(tree: ProcessTree) -> int:
    return sum(1 for node in tree.walk() if node.op in {SEQ, XOR, PAR, LOOP})


def tree_depth(tree: ProcessTree) -> int:
    if not tree.children:
        return 1
    return 1 + max(tree_depth(child) for child in tree.children)


def proposed_tree_stats(result: EnrichmentResult | None) -> dict[str, object]:
    if result is None:
        return {
            "construction_mode": "",
            "construction_note": "",
            "well_formed_configurable_tree": "",
            "shared_node_ratio": "",
            "variant_specific_node_ratio": "",
            "nodes_active_all_variants": "",
            "nodes_active_one_variant": "",
        }
    variants = set(result.cdfg.variants)
    nodes = list(result.tree.walk())
    active_all = sum(1 for node in nodes if node.variants is not None and set(node.variants) == variants)
    active_one = sum(1 for node in nodes if node.variants is not None and len(node.variants) == 1)
    return {
        "construction_mode": result.construction_mode,
        "construction_note": result.construction_note,
        "well_formed_configurable_tree": result.tree.is_well_formed(),
        "shared_node_ratio": safe_ratio(active_all, len(nodes)),
        "variant_specific_node_ratio": safe_ratio(active_one, len(nodes)),
        "nodes_active_all_variants": active_all,
        "nodes_active_one_variant": active_one,
    }


def localization_stats(result: EnrichmentResult | None) -> dict[str, object]:
    if result is None:
        return {
            "residual_localization_precision": "",
            "residual_localization_recall": "",
            "residual_localization_f1": "",
            "residual_localization_tp": "",
            "residual_localization_fp": "",
            "residual_localization_fn": "",
        }
    if result.construction_mode != "anchor_localization" or not result.fragments:
        return {
            "residual_localization_precision": "n/a",
            "residual_localization_recall": "n/a",
            "residual_localization_f1": "n/a",
            "residual_localization_tp": "n/a",
            "residual_localization_fp": "n/a",
            "residual_localization_fn": "n/a",
        }
    loc = residual_relation_localization(result)
    return {
        "residual_localization_precision": loc.precision,
        "residual_localization_recall": loc.recall,
        "residual_localization_f1": loc.f1,
        "residual_localization_tp": loc.true_positive,
        "residual_localization_fp": loc.false_positive,
        "residual_localization_fn": loc.false_negative,
    }


def cdfg_stats(cdfg: CDFG) -> dict[str, object]:
    all_variants = set(cdfg.variants)
    business_edges = {
        edge: ev
        for edge, ev in cdfg.evidence.items()
        if edge[0] != SOURCE and edge[1] != SINK
    }
    return {
        "variants": ";".join(cdfg.variants),
        "activities": len(cdfg.activities),
        "relations_total": len(cdfg.evidence),
        "relations_business": len(business_edges),
        "relations_shared_business": sum(1 for ev in business_edges.values() if set(ev.availability) == all_variants),
        "relations_variant_specific_business": sum(1 for ev in business_edges.values() if len(ev.availability) == 1),
        "relations_partial_shared_business": sum(1 for ev in business_edges.values() if 1 < len(ev.availability) < len(all_variants)),
        "relations_dominant_business": sum(1 for ev in business_edges.values() if ev.dominance != "balanced"),
        "relations_balanced_business": sum(1 for ev in business_edges.values() if ev.dominance == "balanced"),
        "differential_relations_total": len(cdfg.differential_relations()),
    }


def coverage_gap(train_value: object, test_value: object) -> object:
    if isinstance(train_value, (int, float)) and isinstance(test_value, (int, float)):
        return train_value - test_value
    return ""


def safe_ratio(numer: int | float, denom: int | float) -> float:
    return float(numer) / float(denom) if denom else 0.0


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row.keys():
            if key not in seen:
                fieldnames.append(key)
                seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown_summary(path: Path, summary_rows: list[dict[str, object]], cdfg_rows: list[dict[str, object]]) -> None:
    lines = [
        "# Evidence-CPT evaluation outputs",
        "",
        "Recommended paper-level metrics:",
        "",
        "1. **Fitness / coverage**: `train_coverage`, `test_coverage`, and `coverage_gap_train_minus_test`.",
        "   If `coverage_mode=bounded`, coverage is bounded by `max_loop_depth` and should be reported as bounded replay coverage.",
        "2. **Precision proxy**: `mean_*_df_precision_proxy`; this is bounded directly-follows precision from enumerated model language.",
        "   For large logs, use `--skip-df-proxy` and report C-DFG relation statistics instead.",
        "3. **Generalization**: primarily `test_coverage`; report across repeated stratified splits.",
        "4. **Simplicity**: `family_nodes`, `family_visible_leaves`, `family_operator_nodes`, `family_max_depth`.",
        "5. **Consolidation**: `compactness_gain_vs_separate_im` and proposed `shared_node_ratio`.",
        "6. **Variant explainability**: C-DFG shared / variant-specific / dominant relation counts in `cdfg_metrics.csv`.",
        "7. **Differential localization**: residual-localization precision/recall/F1 when construction mode is `anchor_localization`.",
        "8. **Scalability**: `runtime_sec` and model size over increasing log sizes.",
        "",
        "Baselines included:",
        "",
        "- `pooled_im`: one model from the union log.",
        "- `separate_im`: one model per variant.",
        "- `separate_exact`: optional memorization upper-bound; use only for small logs or sanity checks.",
        "",
    ]
    if summary_rows:
        lines.append("## Mean summary by method")
        groups: dict[str, list[dict[str, object]]] = {}
        for row in summary_rows:
            groups.setdefault(str(row["method"]), []).append(row)
        lines.append("| method | mean train coverage | mean test coverage | mean nodes | mean runtime sec |")
        lines.append("|---|---:|---:|---:|---:|")
        for method, rows in sorted(groups.items()):
            train_cov = mean_numeric(rows, "train_coverage")
            test_cov = mean_numeric(rows, "test_coverage")
            nodes = mean_numeric(rows, "family_nodes")
            runtime = mean_numeric(rows, "runtime_sec")
            lines.append(f"| {method} | {train_cov:.3f} | {test_cov:.3f} | {nodes:.1f} | {runtime:.3f} |")
    if cdfg_rows:
        lines.extend(["", "## C-DFG rows", "", f"Total C-DFG result rows: {len(cdfg_rows)}"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def mean_numeric(rows: list[Mapping[str, object]], key: str) -> float:
    vals: list[float] = []
    for row in rows:
        value = row.get(key, "")
        if isinstance(value, (int, float)):
            vals.append(float(value))
        elif isinstance(value, str) and value.strip():
            try:
                vals.append(float(value))
            except ValueError:
                pass
    return statistics.mean(vals) if vals else 0.0


if __name__ == "__main__":
    raise SystemExit(main())
