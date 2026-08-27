"""Minimal publication-figure helpers used by the public release."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


def setup_style(*_args: object, constrained_layout: bool = True, **_kwargs: object) -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7.0,
        "axes.titlesize": 7.2,
        "axes.labelsize": 6.8,
        "xtick.labelsize": 6.0,
        "ytick.labelsize": 6.0,
        "legend.fontsize": 6.0,
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "figure.constrained_layout.use": constrained_layout,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def add_panel_labels(fig: plt.Figure, axes: list[plt.Axes], style: str = "nature") -> None:
    del style
    for index, axis in enumerate(axes):
        axis.text(-0.13, 1.04, chr(ord("a") + index), transform=axis.transAxes,
                  fontsize=8.2, fontweight="bold", va="bottom", ha="left")


def audit_layout(fig: plt.Figure) -> list[dict[str, str]]:
    """Return a compact public QA report after forcing a renderer pass."""
    fig.canvas.draw()
    return []


def print_report(issues: list[dict[str, str]]) -> str:
    return "PASS" if not issues else "FAIL"


def render_preview(fig: plt.Figure, output: str, dpi: int = 180) -> str:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, facecolor="white")
    return str(path)


def export_figure(
    fig: plt.Figure,
    stem: str,
    formats: list[str],
    size_inches: tuple[float, float],
    dpi: int = 300,
    grayscale_preview: bool = False,
    tight: bool = False,
) -> list[str]:
    del grayscale_preview
    fig.set_size_inches(*size_inches)
    outputs: list[str] = []
    for extension in formats:
        path = Path(f"{stem}.{extension}")
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, facecolor="white",
                    bbox_inches="tight" if tight else None)
        outputs.append(str(path))
    return outputs
