"""Cross-solution validation of Python penalty functions against CONOP9 output."""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Any, cast

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from conop_py.cost import (  # noqa: E402
    ConopContext,
    build_section_observations,
    eventual_misfit,
    level_misfit,
    ordinal_misfit,
)
from conop_py.io import parse_loadfile, parse_solution, solution_to_sequence  # noqa: E402
from conop_py.plotting import PARAM_COLORS, PARAM_LABELS, PARAM_TAGS, init_plot, save_plot  # noqa: E402

PENALTY_RE = re.compile(r"(Level|Eventual|Ordinal) Penalty:\s*([0-9.]+)")
DEFAULT_OUT = ROOT / "results_py" / "cross_solution_validation"


def parse_conop9_penalties(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return {name.lower(): float(value) for name, value in PENALTY_RE.findall(text)}


def evaluate_solution(soln_path: Path, section_obs: dict) -> dict[str, float]:
    seq = solution_to_sequence(parse_solution(soln_path))
    ctx = ConopContext.build(seq, section_obs)
    return {
        "ordinal": ordinal_misfit(ctx),
        "level": level_misfit(ctx),
        "eventual": eventual_misfit(ctx),
    }


def collect_rows(results_dir: Path) -> list[dict[str, object]]:
    obs = parse_loadfile(ROOT / "CONOP-run" / "loadfile.dat")
    section_obs = build_section_observations(obs)
    rows: list[dict[str, object]] = []

    for run_dir in sorted(results_dir.glob("*/run_*")):
        soln = run_dir / "bestsoln.dat"
        outmain = run_dir / "outmain.txt"
        if not soln.exists() or not outmain.exists():
            continue

        py_scores = evaluate_solution(soln, section_obs)
        conop9_scores = parse_conop9_penalties(outmain)
        row: dict[str, object] = {
            "experiment": run_dir.parent.name,
            "run_id": run_dir.name,
        }
        for metric in ("ordinal", "level", "eventual"):
            py_score = py_scores[metric]
            ref_score = conop9_scores[metric]
            row[f"python_{metric}"] = py_score
            row[f"conop9_{metric}"] = ref_score
            row[f"delta_{metric}"] = py_score - ref_score
        rows.append(row)
    return rows


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    for metric in ("ordinal", "level", "eventual"):
        deltas = [float(row[f"delta_{metric}"]) for row in rows]
        abs_deltas = [abs(delta) for delta in deltas]
        summary.append({
            "metric": metric,
            "n": len(rows),
            "exact": sum(delta == 0 for delta in deltas),
            "mae": sum(abs_deltas) / len(abs_deltas) if abs_deltas else 0.0,
            "max_abs": max(abs_deltas) if abs_deltas else 0.0,
            "mean_delta": sum(deltas) / len(deltas) if deltas else 0.0,
        })
    return summary


def sequence_diversity(results_dir: Path) -> list[dict[str, object]]:
    sequences = [
        tuple(solution_to_sequence(parse_solution(path)))
        for path in sorted(results_dir.glob("*/run_*/bestsoln.dat"))
    ]
    return [{
        "archived_runs": len(sequences),
        "unique_permutations": len(set(sequences)),
    }]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_eventual_delta_plot(path: Path, rows: list[dict[str, object]]) -> None:
    """Plot Eventual deltas by parameter group for the paper."""
    init_plot()
    import matplotlib.pyplot as plt
    import numpy as np

    grouped = {
        tag: [
            float(row["delta_eventual"])
            for row in rows
            if row["experiment"] == tag
        ]
        for tag in PARAM_TAGS
    }
    tags = [tag for tag in PARAM_TAGS if grouped.get(tag)]
    x = np.arange(len(tags))
    offsets = np.array([-0.18, 0.0, 0.18])
    means = [sum(grouped[tag]) / len(grouped[tag]) for tag in tags]
    all_deltas = [float(row["delta_eventual"]) for row in rows]
    overall_mean = sum(all_deltas) / len(all_deltas)

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    for i, tag in enumerate(tags):
        vals = grouped[tag]
        color = PARAM_COLORS.get(tag, "#c44e52")
        for j, val in enumerate(vals):
            ax.scatter(
                x[i] + offsets[j % len(offsets)],
                val,
                s=64,
                color=color,
                edgecolor="black",
                linewidth=0.5,
                zorder=3,
            )
        ax.hlines(
            means[i],
            x[i] - 0.28,
            x[i] + 0.28,
            color="black",
            linewidth=2.0,
            zorder=4,
        )

    ax.axhline(0, color="black", linewidth=0.9)
    ax.axhline(
        overall_mean,
        color="#555555",
        linestyle="--",
        linewidth=1.0,
        label=f"overall mean = {overall_mean:.2f}",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([PARAM_LABELS.get(tag, tag) for tag in tags], fontsize=9)
    ax.set_ylabel("Python Eventual - CONOP9 Eventual", fontsize=10)
    ax.set_xlabel("CONOP9 parameter group", fontsize=10)
    ax.set_title("Eventual penalty bias by parameter group", fontsize=11)
    ax.tick_params(axis="y", labelsize=9)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    path.parent.mkdir(parents=True, exist_ok=True)
    save_plot(fig, path, dpi=220)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    rows = collect_rows(args.results_dir)
    if not rows:
        raise SystemExit(f"No archived CONOP9 runs found in {args.results_dir}")

    summary = summarize(rows)
    diversity = sequence_diversity(args.results_dir)
    write_csv(args.out_dir / "penalty_validation.csv", rows)
    write_csv(args.out_dir / "summary.csv", summary)
    write_csv(args.out_dir / "sequence_diversity.csv", diversity)
    write_eventual_delta_plot(args.out_dir / "eventual_delta_by_group.png", rows)

    print(f"Wrote {len(rows)} validation rows to {args.out_dir}")
    print(
        f"unique permutations: {diversity[0]['unique_permutations']}/"
        f"{diversity[0]['archived_runs']}"
    )
    for item in summary:
        print(
            f"{item['metric']:8s} n={item['n']:2d} exact={item['exact']:2d} "
            f"MAE={item['mae']:.2f} max_abs={item['max_abs']:.2f} "
            f"mean_delta={item['mean_delta']:.2f}"
        )

    if args.check:
        import math
        errors: list[str] = []
        runs = diversity[0]["archived_runs"]
        perms = diversity[0]["unique_permutations"]
        if runs != 21:
            errors.append(f"Expected 21 archived runs, found {runs}")
        if perms != 21:
            errors.append(f"Expected 21 unique archived permutations, found {perms}")

        baseline = {
            "ordinal": dict(exact=21, mae=0.0, max_abs=0.0, mean=0.0),
            "level": dict(exact=21, mae=0.0, max_abs=0.0, mean=0.0),
            "eventual": dict(exact=0, mae=33.76190476190476,
                             max_abs=83.0, mean=-33.76190476190476),
        }
        for item_raw in summary:
            item = cast(dict[str, Any], item_raw)
            metric = str(item["metric"])
            ref = baseline[metric]
            n = int(item["n"])
            exact = int(item["exact"])
            mae = float(item["mae"])
            max_abs = float(item["max_abs"])
            mean_delta = float(item["mean_delta"])
            if n != 21:
                errors.append(f"{metric}: expected n=21, found {n}")
            if exact != ref["exact"]:
                errors.append(
                    f"{metric}: expected exact={ref['exact']}, found {exact}")
            if not math.isclose(mae, ref["mae"], abs_tol=1e-6):
                errors.append(
                    f"{metric}: MAE changed from validated baseline "
                    f"{ref['mae']}, found {mae}")
            if not math.isclose(max_abs, ref["max_abs"], abs_tol=1e-6):
                errors.append(
                    f"{metric}: max_abs changed from validated baseline "
                    f"{ref['max_abs']}, found {max_abs}")
            if not math.isclose(mean_delta, ref["mean"], abs_tol=1e-6):
                errors.append(
                    f"{metric}: mean_delta changed from validated baseline "
                    f"{ref['mean']}, found {mean_delta}")

        if errors:
            raise SystemExit(
                "Cross-solution validation baseline changed:\n"
                + "\n".join(f"  - {e}" for e in errors))


if __name__ == "__main__":
    main()
