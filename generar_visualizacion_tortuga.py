#!/usr/bin/env python3
"""
Genera un visor 3D interactivo autocontenido para la tortuga con apoyo en la cabeza.

Salida:
  - tortuga_3d_interactiva.html
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import plotly.graph_objects as go


OUT = Path("tortuga_3d_interactiva.html")
FRONT_LOBE_SCALE = 2.555242
CM = np.array([0.507759, -0.001905, 0.138567], dtype=float)


def normalize(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("Vector nulo.")
    return v / n


def frame_from_direction(direction: np.ndarray) -> np.ndarray:
    ez = normalize(direction)
    helper = np.array([0.0, 0.0, 1.0], dtype=float)
    if abs(np.dot(ez, helper)) > 0.95:
        helper = np.array([0.0, 1.0, 0.0], dtype=float)
    ex = normalize(np.cross(helper, ez))
    ey = normalize(np.cross(ez, ex))
    return np.column_stack([ex, ey, ez])


def mesh_from_grid(x: np.ndarray, y: np.ndarray, z: np.ndarray, color: str, name: str, opacity: float = 1.0) -> go.Mesh3d:
    rows, cols = x.shape
    xx = x.ravel()
    yy = y.ravel()
    zz = z.ravel()
    i_idx: list[int] = []
    j_idx: list[int] = []
    k_idx: list[int] = []

    for r in range(rows - 1):
        for c in range(cols - 1):
            p = r * cols + c
            i_idx.extend([p, p + 1])
            j_idx.extend([p + cols, p + cols])
            k_idx.extend([p + cols + 1, p + 1])

    return go.Mesh3d(
        x=xx,
        y=yy,
        z=zz,
        i=i_idx,
        j=j_idx,
        k=k_idx,
        color=color,
        opacity=opacity,
        name=name,
        flatshading=False,
        lighting={"ambient": 0.45, "diffuse": 0.72, "specular": 0.18, "roughness": 0.92},
        lightposition={"x": 100, "y": 80, "z": 180},
        hoverinfo="skip",
        showscale=False,
    )


def ellipsoid_mesh(center: tuple[float, float, float], axes: tuple[float, float, float], color: str, name: str, opacity: float = 1.0) -> go.Mesh3d:
    u = np.linspace(0.0, 2.0 * math.pi, 42)
    v = np.linspace(0.0, math.pi, 24)
    uu, vv = np.meshgrid(u, v)
    ax, ay, az = axes
    cx, cy, cz = center
    x = cx + ax * np.cos(uu) * np.sin(vv)
    y = cy + ay * np.sin(uu) * np.sin(vv)
    z = cz + az * np.cos(vv)
    return mesh_from_grid(x, y, z, color=color, name=name, opacity=opacity)


def capsule_traces(p0: tuple[float, float, float], p1: tuple[float, float, float], radius: float, color: str, name: str) -> list[go.BaseTraceType]:
    p0a = np.array(p0, dtype=float)
    p1a = np.array(p1, dtype=float)
    direction = p1a - p0a
    length = np.linalg.norm(direction)
    rot = frame_from_direction(direction)

    theta = np.linspace(0.0, 2.0 * math.pi, 32)
    zz = np.linspace(0.0, length, 18)
    tt, ll = np.meshgrid(theta, zz)
    local = np.stack([radius * np.cos(tt), radius * np.sin(tt), ll], axis=-1)
    world = local @ rot.T + p0a

    return [
        mesh_from_grid(world[..., 0], world[..., 1], world[..., 2], color=color, name=name),
        sphere_mesh(tuple(p0a), radius, color=color, name=f"{name}-cap0"),
        sphere_mesh(tuple(p1a), radius, color=color, name=f"{name}-cap1"),
    ]


def sphere_mesh(center: tuple[float, float, float], radius: float, color: str, name: str, opacity: float = 1.0) -> go.Mesh3d:
    u = np.linspace(0.0, 2.0 * math.pi, 34)
    v = np.linspace(0.0, math.pi, 20)
    uu, vv = np.meshgrid(u, v)
    cx, cy, cz = center
    x = cx + radius * np.cos(uu) * np.sin(vv)
    y = cy + radius * np.sin(uu) * np.sin(vv)
    z = cz + radius * np.cos(vv)
    return mesh_from_grid(x, y, z, color=color, name=name, opacity=opacity)


def support_pad_trace() -> go.Surface:
    u = np.linspace(-1.0, 1.0, 42)
    v = np.linspace(-1.0, 1.0, 34)
    uu, vv = np.meshgrid(u, v)
    mask = uu * uu + vv * vv <= 1.0
    x = np.where(mask, 0.30 * uu, np.nan)
    y = np.where(mask, 0.18 * vv, np.nan)
    z = np.where(mask, -1.50, np.nan) # support_center_z is -1.44 relative to height
    return go.Surface(
        x=x + 0.50, # x center is 0.50
        y=y,
        z=z,
        surfacecolor=np.where(mask, 1.0, np.nan),
        colorscale=[[0.0, "#f59e0b"], [1.0, "#f59e0b"]],
        showscale=False,
        opacity=0.95,
        name="Base de apoyo",
        hoverinfo="skip",
    )


def support_axis_trace() -> go.Scatter3d:
    return go.Scatter3d(
        x=[0.50, 0.50],
        y=[0.0, 0.0],
        z=[-1.65, 1.25],
        mode="lines",
        line={"color": "#ef4444", "width": 6, "dash": "dash"},
        name="Eje de apoyo",
        hoverinfo="skip",
    )


def cm_trace() -> go.Scatter3d:
    return go.Scatter3d(
        x=[CM[0]],
        y=[CM[1]],
        z=[CM[2]],
        mode="markers",
        marker={"size": 7, "color": "#0f766e"},
        name="Centro de masa",
        hovertemplate="CM<br>x=%{x:.4f}<br>y=%{y:.4f}<br>z=%{z:.4f}<extra></extra>",
    )


def build_figure() -> go.Figure:
    traces: list[go.BaseTraceType] = []

    shell_color = "#4f6e3d"
    shell_high = "#688b4f"
    limb_color = "#6f9a5b"
    front_color = "#88ad66"

    # Shell shown as organic outer volume; the cavity is omitted in the viewer for clarity.
    traces.append(ellipsoid_mesh((-2.4, 0.0, 1.05), (1.95, 1.90, 0.95), shell_color, "Caparazon"))
    traces.append(ellipsoid_mesh((-2.55, 0.0, 1.48), (1.45, 1.50, 0.74), shell_high, "Boveda", opacity=0.95))

    # Cuello y cabeza
    traces.extend(capsule_traces((-1.5, 0.0, 0.6), (0.0, 0.0, -0.6), 0.35, limb_color, "Cuello gruesa"))
    traces.append(ellipsoid_mesh((0.2, 0.0, -0.9), (0.7, 0.5, 0.3), limb_color, "Cabeza realista"))
    traces.extend(capsule_traces((0.4, 0.0, -1.0), (0.5, 0.0, -1.4), 0.12, front_color, "Hocico inferior"))

    # Aletas delanteras
    traces.extend(capsule_traces((-1.0, 1.0, 0.4), (0.5, 1.8, 0.0), 0.25, front_color, "Hombro sup"))
    traces.extend(capsule_traces((-1.0, -1.0, 0.4), (0.5, -1.8, 0.0), 0.25, front_color, "Hombro inf"))
    
    s = FRONT_LOBE_SCALE
    traces.append(ellipsoid_mesh((1.2, 2.2, -0.1), (1.2 * s, 0.8 * s, 0.2 * s), "#95bb72", "Aleta principal sup"))
    traces.append(ellipsoid_mesh((1.2, -2.2, -0.1), (1.2 * s, 0.8 * s, 0.2 * s), "#95bb72", "Aleta principal inf"))

    # Patas traseras (conectadas)
    traces.extend(capsule_traces((-3.5, 1.0, 0.5), (-4.2, 1.5, 0.2), 0.15, "#799365", "Pata trasera sup"))
    traces.extend(capsule_traces((-3.5, -1.0, 0.5), (-4.2, -1.5, 0.2), 0.15, "#799365", "Pata trasera inf"))
    # Cola (conectada)
    traces.extend(capsule_traces((-4.0, 0.0, 0.8), (-5.0, 0.0, 0.2), 0.10, "#688259", "Cola adjunta"))

    traces.append(support_pad_trace())
    traces.append(support_axis_trace())
    traces.append(cm_trace())

    fig = go.Figure(data=traces)
    fig.update_layout(
        title={
            "text": "Tortuga 3D con apoyo en la cabeza y centro de masa",
            "x": 0.5,
            "font": {"size": 24},
        },
        scene={
            "xaxis": {"title": "x", "backgroundcolor": "rgba(255,255,255,0.05)", "gridcolor": "rgba(255,255,255,0.1)", "zerolinecolor": "rgba(255,255,255,0.2)"},
            "yaxis": {"title": "y", "backgroundcolor": "rgba(255,255,255,0.05)", "gridcolor": "rgba(255,255,255,0.1)", "zerolinecolor": "rgba(255,255,255,0.2)"},
            "zaxis": {"title": "z", "backgroundcolor": "rgba(255,255,255,0.05)", "gridcolor": "rgba(255,255,255,0.1)", "zerolinecolor": "rgba(255,255,255,0.2)"},
            "aspectmode": "data",
            "camera": {"eye": {"x": 1.95, "y": -1.75, "z": 0.95}},
        },
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font={"color": "#f8fafc"},
        margin={"l": 0, "r": 0, "b": 0, "t": 54},
        legend={"orientation": "h", "yanchor": "bottom", "y": 0.01, "xanchor": "center", "x": 0.5},
    )
    return fig


def main() -> None:
    fig = build_figure()
    fig.write_html(OUT, include_plotlyjs=True, full_html=True, config={"responsive": True, "displaylogo": False})
    print(f"Visor 3D guardado en: {OUT}")


if __name__ == "__main__":
    main()
