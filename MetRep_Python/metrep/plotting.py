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

from metrep.initial_conditions import STATE_INDEX, STATE_NAMES
from metrep.simulate import SimulationResult


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
    holistic_table: pd.DataFrame | None = None,
    layout_mode: str = "filtered",
    include_singletons: bool = True,
    title: str | None = None,
) -> Path | None:
    """Plot a parameter compensation network from nullspace edge weights.

    Nodes are parameters. Edges are dominant opposite-sign compensation pairs
    extracted from nullspace directions. When ``layout_mode="all-zones"`` and
    ``include_all_nodes=True``, all analyzed parameters are shown: compensated
    nodes are placed in the central spring-layout graph and no-edge parameters
    are placed in peripheral recommendation zones.
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
    if include_singletons:
        for name in singleton_only:
            graph.add_node(name)
    if include_all_nodes:
        for name in parameter_names:
            graph.add_node(name)

    if graph.number_of_nodes() == 0:
        return None

    class_column = None
    metadata = pd.DataFrame()
    if holistic_table is not None and not holistic_table.empty:
        parameter_column = "parameter" if "parameter" in holistic_table.columns else "param"
        if parameter_column in holistic_table.columns:
            metadata = holistic_table.set_index(parameter_column, drop=False)
            if "recommendation_3class" in metadata.columns:
                class_column = "recommendation_3class"
            elif "recommendation" in metadata.columns:
                class_column = "recommendation"

    recommendation_colors = {
        "Estimate": "#1f77b4",
        "Fix (anchor)": "#ff7f0e",
        "Fix (irrelevant)": "#8c8c8c",
    }

    def recommendation_for(node: str) -> str:
        if class_column and node in metadata.index:
            value = str(metadata.loc[node, class_column])
            if value in recommendation_colors:
                return value
        if node in singleton_only:
            return "Fix (anchor)"
        return "Estimate" if graph.degree(node) > 0 else "Fix (irrelevant)"

    def sensitivity_for(node: str) -> float:
        if "sensitivity_0to1" in metadata.columns and node in metadata.index:
            value = float(metadata.loc[node, "sensitivity_0to1"])
            return float(np.clip(value, 0.0, 1.0)) if np.isfinite(value) else 0.0
        return 0.65 if graph.degree(node) > 0 else 0.20

    def node_color(node: str) -> str:
        return recommendation_colors.get(recommendation_for(node), "#8c8c8c")

    def node_size(node: str) -> float:
        sensitivity = sensitivity_for(node)
        return 420.0 + 2600.0 * np.sqrt(max(sensitivity, 0.0))

    def circular_positions(nodes: list[str], center: tuple[float, float], radius: float) -> dict[str, tuple[float, float]]:
        if not nodes:
            return {}
        if len(nodes) == 1:
            return {nodes[0]: center}
        angles = np.linspace(0, 2 * np.pi, len(nodes), endpoint=False)
        return {
            node: (
                center[0] + radius * np.cos(angle),
                center[1] + 0.58 * radius * np.sin(angle),
            )
            for node, angle in zip(nodes, angles)
        }

    if layout_mode not in {"filtered", "all-zones"}:
        raise ValueError("layout_mode must be either 'filtered' or 'all-zones'.")

    if layout_mode == "all-zones" and include_all_nodes:
        edge_nodes = sorted({node for edge in graph.edges() for node in edge})
        isolated_nodes = sorted([node for node in graph.nodes() if graph.degree(node) == 0])
        pos: dict[str, tuple[float, float]] = {}
        if edge_nodes:
            central_graph = graph.subgraph(edge_nodes).copy()
            central_pos = nx.spring_layout(
                central_graph,
                seed=42,
                k=1.2 / np.sqrt(max(len(edge_nodes), 1)),
                iterations=300,
                weight="weight",
            )
            for node, (x_coord, y_coord) in central_pos.items():
                pos[node] = (4.2 * float(x_coord), 3.4 * float(y_coord))

        zone_estimate = [node for node in isolated_nodes if recommendation_for(node) == "Estimate"]
        zone_anchor = [node for node in isolated_nodes if recommendation_for(node) == "Fix (anchor)"]
        zone_irrelevant = [node for node in isolated_nodes if recommendation_for(node) == "Fix (irrelevant)"]
        zone_other = sorted(set(isolated_nodes) - set(zone_estimate) - set(zone_anchor) - set(zone_irrelevant))
        pos.update(circular_positions(zone_estimate, center=(-6.6, 3.8), radius=max(1.2, 0.18 * len(zone_estimate))))
        pos.update(circular_positions(zone_anchor, center=(6.6, 3.8), radius=max(1.2, 0.18 * len(zone_anchor))))
        pos.update(circular_positions(zone_irrelevant, center=(0.0, -5.2), radius=max(2.2, 0.13 * len(zone_irrelevant))))
        pos.update(circular_positions(zone_other, center=(0.0, 5.8), radius=max(1.2, 0.18 * len(zone_other))))
        layout_name = "all-zones"
    else:
        components = [sorted(component) for component in nx.weakly_connected_components(graph)]
        components = sorted(components, key=len, reverse=True)
        if len(components) > 1 and graph.number_of_edges() > 0:
            pos = {}
            ncols = min(3, len(components))
            x_spacing = 5.2
            y_spacing = 3.9
            for component_index, component in enumerate(components):
                subgraph = graph.subgraph(component).copy()
                local = nx.spring_layout(
                    subgraph,
                    seed=42 + component_index,
                    k=0.9 / np.sqrt(max(len(component), 1)),
                    iterations=300,
                    weight="weight",
                )
                row = component_index // ncols
                col = component_index % ncols
                x_offset = (col - (ncols - 1) / 2.0) * x_spacing
                y_offset = -row * y_spacing
                for node, (x_coord, y_coord) in local.items():
                    scale = 2.2 if len(component) > 3 else 1.4
                    pos[node] = (scale * float(x_coord) + x_offset, scale * float(y_coord) + y_offset)
            layout_name = "component-spring"
        else:
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

    node_colors = [node_color(node) for node in graph.nodes()]
    node_sizes = [node_size(node) for node in graph.nodes()]
    node_alphas = [0.95 if graph.degree(node) > 0 else (0.82 if layout_mode == "all-zones" else 0.42) for node in graph.nodes()]
    if graph.number_of_edges() > 0:
        weights = np.array([graph[u][v]["weight"] for u, v in graph.edges()], dtype=float)
        widths = 1.2 + 6.0 * weights / weights.max()
    else:
        widths = []

    edge_labels = {
        (source, target): label
        for source, target in graph.edges()
        if (label := edge_sign_label(source, target))
    }

    if layout_mode == "all-zones" and include_all_nodes:
        fig_width = max(22.0, min(34.0, 15.0 + 0.18 * graph.number_of_nodes()))
        fig_height = max(15.0, min(26.0, 10.0 + 0.11 * graph.number_of_nodes()))
    else:
        fig_width, fig_height = 22.0, 13.0
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
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
    labels = {node: node for node in graph.nodes() if graph.degree(node) > 0 or node in singleton_only}
    if layout_mode == "all-zones" and include_all_nodes:
        labels = {node: node for node in graph.nodes()}
    central_font = max(7, min(11, int(145 / max(len(labels), 1))))
    isolated_font = max(5, min(8, int(115 / max(len(labels), 1))))
    nx.draw_networkx_labels(
        graph,
        pos,
        labels=labels,
        ax=ax,
        font_size=central_font if layout_mode != "all-zones" else isolated_font,
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
        Patch(facecolor=recommendation_colors["Estimate"], edgecolor="white", label="Estimate"),
        Patch(facecolor=recommendation_colors["Fix (anchor)"], edgecolor="white", label="Fix (anchor)"),
        Patch(facecolor=recommendation_colors["Fix (irrelevant)"], edgecolor="white", label="Fix (irrelevant)"),
        Line2D([0], [0], color="black", linewidth=4.5, label="Thicker edge = stronger compensation"),
        Line2D([0], [0], color="black", linewidth=0, marker="o", markersize=12, label="Larger node = higher sensitivity"),
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
        bbox_to_anchor=(1.01, 0.99),
        frameon=True,
        framealpha=0.94,
        facecolor="white",
        edgecolor="0.85",
        title="How to read this network",
        title_fontsize=12,
        fontsize=10,
    )
    if layout_mode == "all-zones" and include_all_nodes:
        explanation = (
            "Node color = estimate/fix recommendation; node size = sensitivity;\n"
            "edge width = compensation strength. Isolated parameters are shown\n"
            "in peripheral zones."
        )
    else:
        explanation = (
            "Node color = estimate/fix recommendation; node size = sensitivity;\n"
            "edge width = compensation strength. This core view shows only\n"
            "parameters connected by the strongest compensation edges."
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

    title = "Parameter compensation network from nullspace" if title is None else title
    if layout_mode == "all-zones" and include_all_nodes:
        title = "Compensation network across all analyzed parameters"
    ax.set_title(
        f"{title}\n"
        f"layout={layout_name}, edges drawn={graph.number_of_edges()}, "
        f"singleton nodes={len(singleton_only) if include_singletons else 0}, "
        f"analyzed parameters={len(parameter_names)}"
    )
    if layout_mode == "all-zones" and include_all_nodes:
        ax.text(
            0.5,
            0.96,
            f"All {len(parameter_names)} parameters shown; only strongest compensation edges drawn.",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=11,
            bbox={"boxstyle": "round,pad=0.35", "fc": "white", "ec": "0.85", "alpha": 0.92},
        )
        zone_labels = [
            (-6.6, 5.9, "Estimated parameters\nwith weak/no compensation"),
            (6.6, 5.9, "Fixed anchor parameters\nwith weak/no compensation"),
            (0.0, -7.2, "Fixed irrelevant / low-sensitivity\nparameters"),
            (0.0, 3.2, "Central network:\nstrongest compensation edges"),
        ]
        for x_coord, y_coord, text in zone_labels:
            ax.text(
                x_coord,
                y_coord,
                text,
                ha="center",
                va="center",
                fontsize=10,
                color="0.25",
                bbox={"boxstyle": "round,pad=0.30", "fc": "white", "ec": "0.88", "alpha": 0.86},
            )
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def parameter_scenario_classes(
    holistic_table: pd.DataFrame,
    thresholds: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Return practical scenario labels for sensitivity/nullspace diagnostics."""

    thresholds = {} if thresholds is None else thresholds
    table = holistic_table.copy()
    class_column = "recommendation_3class" if "recommendation_3class" in table.columns else "recommendation"
    sens_hi = float(thresholds.get("sens_hi", table["sensitivity_0to1"].quantile(0.75)))
    sens_lo = float(thresholds.get("sens_lo", table["sensitivity_0to1"].quantile(0.25)))
    null_hi = float(thresholds.get("null_hi", table["nullspace_0to1"].quantile(0.75)))

    def classify(row: pd.Series) -> str:
        sensitivity = float(row["sensitivity_0to1"])
        nullspace = float(row["nullspace_0to1"])
        recommendation = str(row[class_column])
        if recommendation == "Fix (anchor)" and sensitivity >= sens_hi:
            return "Anchor candidate: high sensitivity but should be fixed"
        if sensitivity >= sens_hi and nullspace >= null_hi:
            return "Risky estimate: sensitive but compensated"
        if recommendation == "Estimate" and sensitivity >= sens_hi and nullspace < null_hi:
            return "Estimate candidate: sensitive and weakly compensated"
        if recommendation == "Fix (irrelevant)" or sensitivity <= sens_lo:
            return "Fix candidate: low sensitivity"
        return "Weak/isolated parameter"

    table["scenario_label"] = table.apply(classify, axis=1)
    columns = [
        "parameter",
        "sensitivity_0to1",
        "nullspace_0to1",
        "identifiable_0to1",
        class_column,
        "scenario_label",
    ]
    result = table[columns].rename(columns={class_column: "recommendation_3class"})
    return result.sort_values(["sensitivity_0to1", "nullspace_0to1"], ascending=[False, False])


def plot_parameter_scenario_map(
    holistic_table: pd.DataFrame,
    output_path: Path,
    thresholds: dict[str, float] | None = None,
) -> Path:
    """Plot sensitivity vs nullspace involvement for all analyzed parameters."""

    _apply_plot_style()
    thresholds = {} if thresholds is None else thresholds
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table = holistic_table.copy()
    class_column = "recommendation_3class" if "recommendation_3class" in table.columns else "recommendation"
    sens_hi = float(thresholds.get("sens_hi", table["sensitivity_0to1"].quantile(0.75)))
    null_hi = float(thresholds.get("null_hi", table["nullspace_0to1"].quantile(0.75)))

    color_map = {
        "Estimate": "#1f77b4",
        "Fix (anchor)": "#ff7f0e",
        "Fix (irrelevant)": "#8c8c8c",
    }
    fig, (ax, label_ax) = plt.subplots(
        1,
        2,
        figsize=(16, 9),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [3.2, 1.15]},
    )
    for recommendation, group in table.groupby(class_column):
        ax.scatter(
            group["nullspace_0to1"],
            group["sensitivity_0to1"],
            s=58,
            alpha=0.84,
            color=color_map.get(str(recommendation), "#9467bd"),
            edgecolor="white",
            linewidth=0.6,
            label=str(recommendation),
        )

    ax.axhline(sens_hi, color="0.25", linestyle="--", linewidth=1.2)
    ax.axvline(null_hi, color="0.25", linestyle="--", linewidth=1.2)
    left_x = max(0.03, null_hi * 0.45)
    right_x = min(0.98, null_hi + (1.0 - null_hi) * 0.50)
    low_y = max(0.05, sens_hi * 0.45)
    high_y = min(0.98, sens_hi + (1.0 - sens_hi) * 0.45)
    quadrant_style = {"fontsize": 10.5, "color": "0.25", "ha": "center", "va": "center"}
    ax.text(left_x, high_y, "Good estimate\ncandidates", **quadrant_style)
    ax.text(right_x, high_y, "Sensitive but\ncompensated", **quadrant_style)
    ax.text(left_x, low_y, "Low sensitivity /\nfix", **quadrant_style)
    ax.text(right_x, low_y, "Compensatory but\nweakly sensitive", **quadrant_style)

    important = set(table.nlargest(10, "sensitivity_0to1")["parameter"])
    important.update(table.nlargest(10, "nullspace_0to1")["parameter"])
    estimate_names = table.loc[table[class_column] == "Estimate", "parameter"].tolist()
    labels = table[table["parameter"].isin(important)].copy()
    labels = labels.sort_values(["sensitivity_0to1", "nullspace_0to1"], ascending=[False, False])
    for index, row in enumerate(labels.itertuples(index=False)):
        x_value = float(getattr(row, "nullspace_0to1"))
        y_value = float(getattr(row, "sensitivity_0to1"))
        offset_x = 0.010 if index % 2 == 0 else -0.010
        offset_y = 0.012 if index % 3 else -0.014
        ax.annotate(
            getattr(row, "parameter"),
            (x_value, y_value),
            xytext=(x_value + offset_x, y_value + offset_y),
            textcoords="data",
            fontsize=7.2,
            color="0.15",
            arrowprops={"arrowstyle": "-", "color": "0.70", "linewidth": 0.45},
            bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.72},
        )

    ax.set_xlim(-0.03, 1.04)
    ax.set_ylim(-0.03, 1.04)
    ax.set_xlabel("Nullspace involvement / compensation involvement (0..1)")
    ax.set_ylabel("Sensitivity (0..1)")
    ax.set_title("Parameter scenario map: sensitivity versus compensation involvement")
    ax.grid(alpha=0.22, linewidth=0.8)
    label_ax.axis("off")
    handles, labels_legend = ax.get_legend_handles_labels()
    label_ax.legend(
        handles,
        labels_legend,
        loc="lower left",
        frameon=True,
        framealpha=0.92,
        title="Recommendation",
    )
    wrapped = "\n".join(f"- {name}" for name in estimate_names)
    label_ax.text(
        0.0,
        1.0,
        "Estimate-class parameters\n(listed to avoid label overlap)\n\n" + wrapped,
        ha="left",
        va="top",
        fontsize=7.4,
        linespacing=1.18,
        bbox={"boxstyle": "round,pad=0.45", "fc": "white", "ec": "0.85", "alpha": 0.96},
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_sensitivity_ranked_by_class(holistic_table: pd.DataFrame, output_path: Path) -> Path:
    """Plot all parameters ranked by sensitivity and colored by recommendation."""

    _apply_plot_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table = holistic_table.sort_values("sensitivity_0to1", ascending=True).copy()
    class_column = "recommendation_3class" if "recommendation_3class" in table.columns else "recommendation"
    color_map = {
        "Estimate": "#1f77b4",
        "Fix (anchor)": "#ff7f0e",
        "Fix (irrelevant)": "#8c8c8c",
    }
    fig_height = max(13.0, 0.18 * len(table))
    fig, ax = plt.subplots(figsize=(13, fig_height), constrained_layout=True)
    ax.barh(
        table["parameter"],
        table["sensitivity_0to1"],
        color=[color_map.get(str(value), "#9467bd") for value in table[class_column]],
        height=0.76,
    )
    ax.set_xlabel("Sensitivity (0..1)")
    ax.set_ylabel("Parameter")
    ax.set_title("All analyzed parameters ranked by sensitivity")
    ax.grid(axis="x", alpha=0.25, linewidth=0.8)
    handles = [
        Patch(facecolor=color_map[label], label=label)
        for label in ["Estimate", "Fix (anchor)", "Fix (irrelevant)"]
        if label in set(table[class_column])
    ]
    ax.legend(handles=handles, loc="lower right", frameon=True, framealpha=0.92)
    ax.tick_params(axis="y", labelsize=7.5)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path
