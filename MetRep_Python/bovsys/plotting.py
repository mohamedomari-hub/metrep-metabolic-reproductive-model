"""Reusable plotting helpers."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from bovsys.initial_conditions import STATE_INDEX, STATE_NAMES
from bovsys.simulate import SimulationResult


DEFAULT_OUTPUTS = ["Glucose", "Insulin", "IGF1", "P4", "E2", "Follicle", "CL"]

PLOT_LABELS = {
    "DMI": "Dry matter intake (g/day)",
    "Milk": "Milk yield (L/day)",
    "Glucose": "Glucose",
    "Insulin": "Insulin",
    "IGF1": "IGF-1",
    "P4": "Progesterone (P4)",
    "E2": "Estradiol (E2)",
    "Follicle": "Follicle",
    "CL": "Corpus luteum (CL)",
}

SCENARIO_LABELS = {
    "non_lactating_standard": "Non-lactating | standard",
    "non_lactating_acute": "Non-lactating | acute",
    "non_lactating_chronic": "Non-lactating | chronic",
    "lactating_c0_20": "Lactating | c0 = 20%",
    "lactating_c0_22_5": "Lactating | c0 = 22.5%",
    "lactating_c0_25": "Lactating | c0 = 25%",
    "lactating_c0_30": "Lactating | c0 = 30%",
}


def _apply_plot_style() -> None:
    """Apply readable defaults for exported analysis figures."""

    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 15,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
            "figure.titlesize": 17,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _label(name: str) -> str:
    return PLOT_LABELS.get(name, name.replace("_", " "))


def _scenario_label(name: str) -> str:
    return SCENARIO_LABELS.get(name, name.replace("_", " ").title())


def plot_selected_states(
    result: SimulationResult,
    output_path: Path,
    states: list[str] | None = None,
) -> Path:
    """Plot selected states from one simulation."""

    _apply_plot_style()
    states = DEFAULT_OUTPUTS if states is None else states
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(states),
        1,
        figsize=(11, max(4.5, 2.45 * len(states))),
        sharex=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)

    for ax, state in zip(axes, states):
        ax.plot(result.t, result.y[:, STATE_INDEX[state]], linewidth=2.0, color="tab:blue")
        ax.set_ylabel(_label(state))
        ax.grid(alpha=0.25, linewidth=0.8)
        ax.margins(x=0.01)

    axes[-1].set_xlabel("Time (days)")
    fig.suptitle(_scenario_label(result.scenario_name), y=1.01)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_all_states_grid(result: SimulationResult, output_path: Path) -> Path:
    """Plot all 22 model states/species in a grid for one simulation."""

    _apply_plot_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ncols = 4
    nrows = int(np.ceil(len(STATE_NAMES) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.4 * ncols, 2.9 * nrows),
        sharex=True,
        constrained_layout=True,
    )
    axes_flat = np.asarray(axes).ravel()

    for ax, state in zip(axes_flat, STATE_NAMES):
        ax.plot(result.t, result.y[:, STATE_INDEX[state]], linewidth=1.8, color="tab:blue")
        ax.set_title(_label(state), fontsize=12)
        ax.set_ylabel(_label(state), fontsize=10)
        ax.grid(alpha=0.22, linewidth=0.7)
        ax.margins(x=0.01)

    for ax in axes_flat[len(STATE_NAMES) :]:
        ax.axis("off")

    for ax in axes_flat[-ncols:]:
        if ax.has_data():
            ax.set_xlabel("Time (days)")

    fig.suptitle(f"{_scenario_label(result.scenario_name)}: all model species", y=1.01)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_scenario_comparison(
    results: list[SimulationResult],
    output_path: Path,
    states: list[str] | None = None,
    title: str | None = None,
) -> Path:
    """Plot selected states for multiple scenarios on shared axes."""

    _apply_plot_style()
    states = DEFAULT_OUTPUTS if states is None else states
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(states),
        1,
        figsize=(12, max(4.5, 2.55 * len(states))),
        sharex=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)

    for ax, state in zip(axes, states):
        idx = STATE_INDEX[state]
        for result in results:
            label = _scenario_label(result.scenario_name)
            ax.plot(result.t, result.y[:, idx], linewidth=1.8, label=label)
        ax.set_ylabel(_label(state))
        ax.grid(alpha=0.25, linewidth=0.8)
        ax.margins(x=0.01)

    axes[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    axes[-1].set_xlabel("Time (days)")
    if title:
        fig.suptitle(title, y=1.01)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_uncertainty_band(
    t: np.ndarray,
    samples: np.ndarray,
    state_name: str,
    output_path: Path,
    percentiles: tuple[float, float] = (5.0, 95.0),
) -> Path:
    """Plot median and percentile bands for one state."""

    _apply_plot_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lo, hi = np.percentile(samples, percentiles, axis=0)
    med = np.percentile(samples, 50.0, axis=0)

    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    ax.fill_between(t, lo, hi, alpha=0.25, color="tab:blue", label=f"{percentiles[0]:.0f}-{percentiles[1]:.0f}% band")
    ax.plot(t, med, linewidth=2.2, color="tab:blue", label="Median")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel(_label(state_name))
    ax.set_title(f"Uncertainty Propagation: {_label(state_name)}")
    ax.legend(loc="best", frameon=False)
    ax.grid(alpha=0.25, linewidth=0.8)
    ax.margins(x=0.01)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_compensation_network(
    edges: pd.DataFrame,
    nullspace: np.ndarray,
    parameter_names: list[str],
    output_path: Path,
    threshold_relative: float = 0.35,
    topk: int = 10,
    max_edges: int = 50,
    include_all_nodes: bool = False,
) -> Path | None:
    """Plot a parameter compensation network from nullspace edge weights.

    Nodes are parameters. Edges are dominant opposite-sign compensation pairs
    extracted from nullspace directions. Isolated orange nodes are parameters
    that dominate a nullspace direction alone.
    """

    if edges.empty and nullspace.size == 0:
        return None

    _apply_plot_style()
    try:
        import networkx as nx
    except Exception:
        return None

    output_path.parent.mkdir(parents=True, exist_ok=True)
    name_to_index = {name: index for index, name in enumerate(parameter_names)}

    def singleton_null_params() -> list[str]:
        singles: set[str] = set()
        if nullspace.size == 0:
            return []
        for vector in nullspace:
            values_abs = np.abs(vector.astype(float))
            max_abs = values_abs.max()
            if max_abs <= 0:
                continue
            keep = np.where(values_abs >= threshold_relative * max_abs)[0]
            if keep.size > topk:
                keep = np.argsort(values_abs)[::-1][:topk]
            if keep.size == 1:
                singles.add(parameter_names[int(keep[0])])
        return sorted(singles)

    def edge_sign_label(source: str, target: str) -> str:
        if nullspace.size == 0:
            return ""
        source_index = name_to_index.get(source)
        target_index = name_to_index.get(target)
        if source_index is None or target_index is None:
            return ""

        score = 0.0
        used = 0
        for vector in nullspace:
            values_abs = np.abs(vector)
            max_abs = values_abs.max()
            if max_abs <= 0:
                continue
            if (
                abs(vector[source_index]) >= threshold_relative * max_abs
                and abs(vector[target_index]) >= threshold_relative * max_abs
            ):
                score += np.sign(vector[source_index] * vector[target_index]) * (
                    abs(vector[source_index]) * abs(vector[target_index])
                )
                used += 1
        if used == 0:
            return ""
        return "+" if score > 0 else "-"

    graph = nx.DiGraph()
    if not edges.empty:
        if {"parameter_a", "parameter_b"}.issubset(edges.columns):
            source_column, target_column = "parameter_a", "parameter_b"
        elif {"increase_param", "decrease_param"}.issubset(edges.columns):
            source_column, target_column = "increase_param", "decrease_param"
        else:
            raise ValueError(
                "Compensation edge table must contain either parameter_a/parameter_b "
                "or increase_param/decrease_param columns."
            )
        top_edges = edges.sort_values("weight", ascending=False).head(max_edges)
        for row in top_edges.itertuples(index=False):
            source = getattr(row, source_column)
            target = getattr(row, target_column)
            graph.add_edge(source, target, weight=float(row.weight))

    core_nodes = {node for edge in graph.edges() for node in edge}
    singleton_only = sorted(set(singleton_null_params()) - core_nodes)
    for name in singleton_only:
        graph.add_node(name)
    if include_all_nodes:
        for name in parameter_names:
            graph.add_node(name)

    if graph.number_of_nodes() == 0:
        return None

    try:
        pos = nx.nx_agraph.graphviz_layout(graph, prog="neato")
        layout_name = "graphviz:neato"
    except Exception:
        try:
            pos = nx.nx_pydot.graphviz_layout(graph, prog="neato")
            layout_name = "pydot:neato"
        except Exception:
            pos = nx.kamada_kawai_layout(graph)
            layout_name = "kamada-kawai"

    node_colors = []
    node_sizes = []
    node_alphas = []
    for node in graph.nodes():
        if node in singleton_only:
            node_colors.append("tab:orange")
            node_sizes.append(3200)
            node_alphas.append(0.95)
        elif graph.degree(node) > 0:
            node_colors.append("tab:blue")
            node_sizes.append(3400)
            node_alphas.append(0.95)
        else:
            node_colors.append("lightgray")
            node_sizes.append(850 if include_all_nodes else 2400)
            node_alphas.append(0.28 if include_all_nodes else 0.75)
    if graph.number_of_edges() > 0:
        weights = np.array([graph[u][v]["weight"] for u, v in graph.edges()], dtype=float)
        widths = 1.5 + 7.0 * weights / weights.max()
    else:
        widths = []

    edge_labels = {
        (source, target): label
        for source, target in graph.edges()
        if (label := edge_sign_label(source, target))
    }

    fig, ax = plt.subplots(figsize=(22, 13))
    nx.draw_networkx_nodes(
        graph,
        pos,
        ax=ax,
        node_size=node_sizes,
        node_color=node_colors,
        alpha=node_alphas,
        linewidths=1.0,
        edgecolors="white",
    )
    nx.draw_networkx_edges(
        graph,
        pos,
        ax=ax,
        arrows=True,
        width=widths,
        alpha=0.8,
        arrowstyle="-|>",
        arrowsize=24,
        min_source_margin=14,
        min_target_margin=20,
        connectionstyle="arc3,rad=0.16",
    )
    labels = {
        node: node
        for node in graph.nodes()
        if graph.degree(node) > 0 or node in singleton_only
    }
    nx.draw_networkx_labels(
        graph,
        pos,
        labels=labels,
        ax=ax,
        font_size=11,
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "none", "alpha": 0.82},
    )
    if edge_labels:
        nx.draw_networkx_edge_labels(
            graph,
            pos,
            ax=ax,
            edge_labels=edge_labels,
            font_size=12,
            bbox={"boxstyle": "round,pad=0.15", "fc": "white", "ec": "none", "alpha": 0.85},
            label_pos=0.55,
            rotate=False,
        )

    legend_handles = [
        Patch(
            facecolor="tab:blue",
            edgecolor="white",
            label="Compensation parameter: appears in at least one nullspace pair.",
        ),
        Patch(
            facecolor="tab:orange",
            edgecolor="white",
            label="Singleton nullspace parameter: dominates a weakly identifiable direction alone.",
        ),
        Patch(
            facecolor="lightgray",
            edgecolor="white",
            alpha=0.35,
            label="Other analyzed parameter: included in SVD, not highlighted in this network.",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            linewidth=3.5,
            label="Edge: parameters can compensate for each other while preserving similar outputs.",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            linewidth=0,
            marker="$+$",
            markersize=12,
            label="+ label: parameters move in the same nullspace direction.",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            linewidth=0,
            marker="$-$",
            markersize=12,
            label="- label: parameters move in opposite nullspace directions.",
        ),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.01, 0.99),
        frameon=True,
        framealpha=0.94,
        facecolor="white",
        edgecolor="0.85",
        title="How to read this network",
        title_fontsize=12,
        fontsize=10,
    )
    explanation = (
        "All selected parameters are still used in the identifiability SVD.\n"
        "For readability, only blue/orange nodes are labeled; gray nodes are\n"
        "background parameters that were analyzed but not emphasized by the\n"
        "dominant nullspace compensation structure. Thicker edges indicate\n"
        "larger aggregate compensation weight."
    )
    ax.text(
        0.01,
        0.01,
        explanation,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        bbox={"boxstyle": "round,pad=0.45", "fc": "white", "ec": "0.85", "alpha": 0.94},
    )

    ax.set_title(
        "Parameter compensation network from nullspace\n"
        f"layout={layout_name}, edges={graph.number_of_edges()}, "
        f"singleton nodes={len(singleton_only)}, analyzed parameters={len(parameter_names)}"
    )
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path
