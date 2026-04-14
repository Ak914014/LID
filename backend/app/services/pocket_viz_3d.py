"""
Parametric 3D pocket + ligand + extended chain visualization for Plotly.

Not tied to real PDB coordinates — demonstrates layout, colors, and cavity mesh.
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np
import plotly.graph_objects as go


def _arc_pocket(
    n: int,
    center: np.ndarray,
    radius: float,
    t0: float,
    t1: float,
    z_jitter: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """C-shaped arc in xy (center below ligand); concave toward origin."""
    t = np.linspace(t0, t1, n)
    pts = np.zeros((n, 3))
    pts[:, 0] = center[0] + radius * np.cos(t)
    pts[:, 1] = center[1] + radius * np.sin(t)
    pts[:, 2] = center[2] + rng.uniform(-z_jitter, z_jitter, n)
    return pts


def _chain_from_point(
    start: np.ndarray,
    n_pts: int,
    tangent: np.ndarray,
    step: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Gentle arc extending from pocket; outward + oscillation so it reads as a chain."""
    tdir = tangent / np.linalg.norm(tangent)
    up = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(tdir, up)) > 0.95:
        up = np.array([0.0, 1.0, 0.0])
    b = np.cross(tdir, up)
    b = b / np.linalg.norm(b)
    n_pts = max(n_pts, 2)
    s = np.linspace(0.0, step * (n_pts - 1), n_pts)
    curve = 0.08 * np.sin(s * 0.35)
    wobble_z = 0.35 * np.sin(s * 0.5) + rng.uniform(-0.12, 0.12, n_pts)
    out = np.zeros((n_pts, 3))
    for i in range(n_pts):
        lateral = curve[i] * b
        out[i] = start + s[i] * tdir + lateral + np.array([0.0, 0.0, wobble_z[i]])
    return out


def _push_apart(
    pocket: np.ndarray,
    chain: np.ndarray,
    ligand: np.ndarray,
    min_d: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Iterative separation so residue spheres stay visually distinct."""
    pocket = pocket.copy()
    chain = chain.copy()
    for _ in range(28):
        for i in range(len(pocket)):
            for j in range(len(ligand)):
                d = pocket[i] - ligand[j]
                dn = np.linalg.norm(d)
                if dn < min_d and dn > 1e-9:
                    pocket[i] += 0.5 * (min_d - dn) * (d / dn)
        for i in range(len(chain)):
            for j in range(len(ligand)):
                d = chain[i] - ligand[j]
                dn = np.linalg.norm(d)
                if dn < min_d and dn > 1e-9:
                    chain[i] += 0.5 * (min_d - dn) * (d / dn)
        for i in range(len(chain)):
            for j in range(len(pocket)):
                d = chain[i] - pocket[j]
                dn = np.linalg.norm(d)
                if dn < min_d and dn > 1e-9:
                    chain[i] += 0.5 * (min_d - dn) * (d / dn)
        for i in range(len(pocket)):
            for j in range(i + 1, len(pocket)):
                d = pocket[i] - pocket[j]
                dn = np.linalg.norm(d)
                if dn < min_d * 0.99 and dn > 1e-9:
                    u = d / dn
                    pocket[i] += 0.5 * (min_d - dn) * u
                    pocket[j] -= 0.5 * (min_d - dn) * u
        for i in range(len(chain)):
            for j in range(i + 1, len(chain)):
                d = chain[i] - chain[j]
                dn = np.linalg.norm(d)
                if dn < min_d * 0.99 and dn > 1e-9:
                    u = d / dn
                    chain[i] += 0.5 * (min_d - dn) * u
                    chain[j] -= 0.5 * (min_d - dn) * u
    return pocket, chain


def _pocket_shell_mesh(
    center: np.ndarray,
    radius: float,
    t0: float,
    t1: float,
    thickness: float = 0.55,
) -> tuple[list[float], list[float], list[float], list[int], list[int], list[int]]:
    """Thin translucent shell hugging the concave pocket arc (grid → triangles)."""
    nu, nv = 28, 6
    u = np.linspace(float(t0), float(t1), nu)
    v = np.linspace(-thickness, thickness, nv)
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for ui in u:
        cu, su = np.cos(ui), np.sin(ui)
        base = center + radius * np.array([cu, su, 0.0])
        nvec = np.array([-su, cu, 0.0])
        for vi in v:
            p = base + vi * nvec * 0.35 + np.array([0.0, 0.0, vi * 0.9])
            xs.append(float(p[0]))
            ys.append(float(p[1]))
            zs.append(float(p[2]))

    def idx(iu: int, iv: int) -> int:
        return iu * nv + iv

    i_tri: list[int] = []
    j_tri: list[int] = []
    k_tri: list[int] = []
    for iu in range(nu - 1):
        for iv in range(nv - 1):
            a = idx(iu, iv)
            b = idx(iu + 1, iv)
            c = idx(iu + 1, iv + 1)
            d = idx(iu, iv + 1)
            i_tri.extend([a, a])
            j_tri.extend([b, c])
            k_tri.extend([c, d])
    return xs, ys, zs, i_tri, j_tri, k_tri


def build_pocket_demo_figure(seed: int = 42) -> go.Figure:
    rng = np.random.default_rng(seed)
    min_d = 2.05

    # Ligand — small cluster at origin
    n_lig = 14
    jitter = rng.normal(0, 0.22, size=(n_lig, 3))
    ligand = jitter - jitter.mean(axis=0, keepdims=True)

    # Pocket arc (residue numbers ~150–200): horseshoe around ligand
    pocket_start, pocket_end = 150, 167
    n_pocket = pocket_end - pocket_start + 1
    center = np.array([0.0, -3.55, 0.0])
    radius = 4.45
    t0, t1 = np.pi / 3.2, 2.12 * np.pi / 3  # ~ concave band
    pocket = _arc_pocket(n_pocket, center, radius, t0, t1, 0.28, rng)

    # Chain (300–340): starts near one arm, extends outward
    chain_start, chain_end = 300, 318
    n_chain = chain_end - chain_start + 1
    t_end = t1
    end_pt = center + radius * np.array([np.cos(t_end), np.sin(t_end), 0.0])
    end_pt[2] += rng.uniform(-0.15, 0.15)
    tangent = np.array([np.cos(t_end + 0.45), np.sin(t_end + 0.45), 0.25])
    chain = _chain_from_point(end_pt + tangent * 0.9, n_chain, tangent, 2.15, rng)

    pocket, chain = _push_apart(pocket, chain, ligand, min_d)

    pocket_nums = np.arange(pocket_start, pocket_end + 1)
    chain_nums = np.arange(chain_start, chain_end + 1)

    # Colors
    pocket_colors = [
        f"rgb({int(220 - 4 * i)}, {int(95 + 3 * i)}, {int(45 + 2 * i)})" for i in range(n_pocket)
    ]
    chain_colors = [
        f"rgb({int(45 + 5 * i)}, {int(110 + 4 * i)}, {int(195 - 3 * i)})" for i in range(n_chain)
    ]

    traces: list[Any] = []

    # Cavity shell (draw first = back)
    sx, sy, sz, si, sj, sk = _pocket_shell_mesh(center, radius * 1.02, t0, t1)
    traces.append(
        go.Mesh3d(
            x=sx,
            y=sy,
            z=sz,
            i=si,
            j=sj,
            k=sk,
            color="rgba(255,140,60,0.22)",
            flatshading=True,
            lighting=dict(ambient=0.45, diffuse=0.65),
            name="Pocket surface",
            showlegend=True,
            hoverinfo="skip",
        )
    )

    traces.append(
        go.Scatter3d(
            x=ligand[:, 0],
            y=ligand[:, 1],
            z=ligand[:, 2],
            mode="markers",
            marker=dict(size=7, color="#d946ef", opacity=0.95, line=dict(width=0.4, color="#701a75")),
            name="Ligand",
            hovertemplate="Ligand atom<br>%{x:.2f}, %{y:.2f}, %{z:.2f}<extra></extra>",
        )
    )

    traces.append(
        go.Scatter3d(
            x=pocket[:, 0],
            y=pocket[:, 1],
            z=pocket[:, 2],
            mode="markers+text",
            text=[str(n) for n in pocket_nums],
            textposition="top center",
            textfont=dict(size=9, color="rgba(120,30,10,0.95)"),
            marker=dict(
                size=11,
                color=pocket_colors,
                opacity=0.92,
                line=dict(width=0.5, color="rgba(100,40,20,0.8)"),
            ),
            name="Pocket residues",
            hovertemplate="Pocket %{text}<br>%{x:.2f}, %{y:.2f}, %{z:.2f}<extra></extra>",
        )
    )

    traces.append(
        go.Scatter3d(
            x=chain[:, 0],
            y=chain[:, 1],
            z=chain[:, 2],
            mode="markers+text",
            text=[str(n) for n in chain_nums],
            textposition="top center",
            textfont=dict(size=9, color="rgba(20,60,120,0.95)"),
            marker=dict(
                size=10,
                color=chain_colors,
                opacity=0.92,
                line=dict(width=0.45, color="rgba(30,80,140,0.75)"),
            ),
            name="Extended chain",
            hovertemplate="Chain %{text}<br>%{x:.2f}, %{y:.2f}, %{z:.2f}<extra></extra>",
        )
    )

    # Connector polyline pocket → chain (visual link)
    traces.append(
        go.Scatter3d(
            x=[pocket[-1, 0], chain[0, 0]],
            y=[pocket[-1, 1], chain[0, 1]],
            z=[pocket[-1, 2], chain[0, 2]],
            mode="lines",
            line=dict(color="rgba(100,100,120,0.45)", width=4),
            name="Pocket–chain",
            showlegend=False,
            hoverinfo="skip",
        )
    )

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(text="Binding pocket (demo geometry)", font=dict(size=16)),
        scene=dict(
            xaxis=dict(showbackground=False, gridcolor="rgba(0,0,0,0.08)", title="x (Å)"),
            yaxis=dict(showbackground=False, gridcolor="rgba(0,0,0,0.08)", title="y (Å)"),
            zaxis=dict(showbackground=False, gridcolor="rgba(0,0,0,0.08)", title="z (Å)"),
            aspectmode="data",
            bgcolor="rgba(248,250,252,0.95)",
            camera=dict(eye=dict(x=1.55, y=-1.45, z=0.85)),
        ),
        margin=dict(l=0, r=0, t=48, b=0),
        paper_bgcolor="white",
        legend=dict(x=0.02, y=0.98, bgcolor="rgba(255,255,255,0.8)"),
        height=640,
    )
    return fig


def figure_to_plotly_json(fig: go.Figure) -> dict[str, Any]:
    return json.loads(fig.to_json())


def write_standalone_html(path: str, seed: int = 42) -> None:
    fig = build_pocket_demo_figure(seed=seed)
    fig.write_html(path, include_plotlyjs="cdn", full_html=True)
