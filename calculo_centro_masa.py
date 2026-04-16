#!/usr/bin/env python3
"""
Calcula y calibra el centro de masa 2D del modelo de tortuga.

Uso rapido:
  python3 calculo_centro_masa.py --solve --plot vista_2d_tortuga.svg
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Geometry:
    # Caparazon
    body_center: tuple[float, float] = (-3.4, 0.0)
    body_a: float = 3.1
    body_b: float = 2.2

    # Cabeza
    head_center: tuple[float, float] = (-1.0, 0.0)
    head_radius: float = 0.95

    # Hocico (apoyo en (0,0))
    beak: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] = (
        (0.0, 0.0),
        (-0.8, 0.28),
        (-0.8, -0.28),
    )

    # Cola
    tail: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] = (
        (-6.5, 0.0),
        (-7.4, 0.65),
        (-7.4, -0.65),
    )

    # Aletas delanteras
    wing_up: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] = (
        (-1.6, 1.0),
        (5.0, 2.45),
        (1.8, 0.9),
    )
    wing_down: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] = (
        (-1.6, -1.0),
        (5.0, -2.45),
        (1.8, -0.9),
    )

    # Aletas traseras
    rear_up: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] = (
        (-4.2, 1.15),
        (-6.8, 2.2),
        (-5.9, 0.75),
    )
    rear_down: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] = (
        (-4.2, -1.15),
        (-6.8, -2.2),
        (-5.9, -0.75),
    )

    # Engrosamientos delanteros para calibrar masa
    weight_top_center: tuple[float, float] = (4.8, 2.1)
    weight_bottom_center: tuple[float, float] = (4.8, -2.1)

    x_min: float = -7.8
    x_max: float = 8.5
    y_min: float = -3.5
    y_max: float = 3.5

    @property
    def bbox_area(self) -> float:
        return (self.x_max - self.x_min) * (self.y_max - self.y_min)


GEOM = Geometry()


def tri_mask(
    x: np.ndarray,
    y: np.ndarray,
    p1: tuple[float, float],
    p2: tuple[float, float],
    p3: tuple[float, float],
) -> np.ndarray:
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    d1 = (x - x2) * (y1 - y2) - (x1 - x2) * (y - y2)
    d2 = (x - x3) * (y2 - y3) - (x2 - x3) * (y - y3)
    d3 = (x - x1) * (y3 - y1) - (x3 - x1) * (y - y1)

    has_neg = (d1 < 0) | (d2 < 0) | (d3 < 0)
    has_pos = (d1 > 0) | (d2 > 0) | (d3 > 0)
    return ~(has_neg & has_pos)


def inside_turtle_base(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    g = GEOM
    mask = ((x - g.body_center[0]) / g.body_a) ** 2 + ((y - g.body_center[1]) / g.body_b) ** 2 <= 1.0
    mask |= (x - g.head_center[0]) ** 2 + (y - g.head_center[1]) ** 2 <= g.head_radius**2
    mask |= tri_mask(x, y, *g.beak)
    mask |= tri_mask(x, y, *g.tail)
    mask |= tri_mask(x, y, *g.wing_up)
    mask |= tri_mask(x, y, *g.wing_down)
    mask |= tri_mask(x, y, *g.rear_up)
    mask |= tri_mask(x, y, *g.rear_down)
    return mask


def inside_turtle(x: np.ndarray, y: np.ndarray, weight_radius: float) -> np.ndarray:
    g = GEOM
    mask = inside_turtle_base(x, y)
    mask |= (x - g.weight_top_center[0]) ** 2 + (y - g.weight_top_center[1]) ** 2 <= weight_radius**2
    mask |= (x - g.weight_bottom_center[0]) ** 2 + (y - g.weight_bottom_center[1]) ** 2 <= weight_radius**2
    return mask


def build_grid(step: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    g = GEOM
    x = np.arange(g.x_min, g.x_max + step, step)
    y = np.arange(g.y_min, g.y_max + step, step)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    base = inside_turtle_base(xx, yy)
    return xx, yy, base, step * step


def centroid_grid(weight_radius: float, step: float = 0.005) -> tuple[float, float, float]:
    g = GEOM
    xx, yy, base, d_a = build_grid(step)
    mask = base.copy()
    mask |= (xx - g.weight_top_center[0]) ** 2 + (yy - g.weight_top_center[1]) ** 2 <= weight_radius**2
    mask |= (xx - g.weight_bottom_center[0]) ** 2 + (yy - g.weight_bottom_center[1]) ** 2 <= weight_radius**2

    area = float(np.sum(mask) * d_a)
    if area <= 0:
        raise RuntimeError("Area nula en la grilla. Revisa el modelo.")

    x_bar = float(np.sum(xx[mask]) * d_a / area)
    y_bar = float(np.sum(yy[mask]) * d_a / area)
    return area, x_bar, y_bar


def solve_weight_radius_grid(
    step: float = 0.005,
    lo: float = 0.4,
    hi: float = 2.8,
    iterations: int = 35,
) -> float:
    g = GEOM
    xx, yy, base, d_a = build_grid(step)

    def x_bar_for(r: float) -> float:
        mask = base.copy()
        mask |= (xx - g.weight_top_center[0]) ** 2 + (yy - g.weight_top_center[1]) ** 2 <= r**2
        mask |= (xx - g.weight_bottom_center[0]) ** 2 + (yy - g.weight_bottom_center[1]) ** 2 <= r**2
        area = float(np.sum(mask) * d_a)
        return float(np.sum(xx[mask]) * d_a / area)

    x_lo = x_bar_for(lo)
    x_hi = x_bar_for(hi)
    if x_lo * x_hi > 0:
        raise RuntimeError("No hubo cambio de signo en x_bar. Ajusta el intervalo [lo, hi].")

    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        x_mid = x_bar_for(mid)
        if x_lo * x_mid <= 0:
            hi = mid
        else:
            lo = mid
            x_lo = x_mid

    return 0.5 * (lo + hi)


def _svg_point(x: float, y: float, x_min: float, y_max: float, scale: float, pad: float) -> tuple[float, float]:
    sx = pad + (x - x_min) * scale
    sy = pad + (y_max - y) * scale
    return sx, sy


def _svg_polygon(points: tuple[tuple[float, float], ...], x_min: float, y_max: float, scale: float, pad: float) -> str:
    mapped = [_svg_point(px, py, x_min, y_max, scale, pad) for px, py in points]
    pts = " ".join(f"{px:.2f},{py:.2f}" for px, py in mapped)
    return pts


def save_svg(weight_radius: float, output: Path, grid_step: float) -> None:
    g = GEOM
    area, x_bar, y_bar = centroid_grid(weight_radius, step=grid_step)

    width = 1200
    height = 620
    pad = 40.0
    scale_x = (width - 2 * pad) / (g.x_max - g.x_min)
    scale_y = (height - 2 * pad) / (g.y_max - g.y_min)
    scale = min(scale_x, scale_y)

    body_cx, body_cy = _svg_point(g.body_center[0], g.body_center[1], g.x_min, g.y_max, scale, pad)
    head_cx, head_cy = _svg_point(g.head_center[0], g.head_center[1], g.x_min, g.y_max, scale, pad)
    wt_cx, wt_cy = _svg_point(g.weight_top_center[0], g.weight_top_center[1], g.x_min, g.y_max, scale, pad)
    wb_cx, wb_cy = _svg_point(g.weight_bottom_center[0], g.weight_bottom_center[1], g.x_min, g.y_max, scale, pad)
    cm_x, cm_y = _svg_point(x_bar, y_bar, g.x_min, g.y_max, scale, pad)
    support_x, support_y = _svg_point(0.0, 0.0, g.x_min, g.y_max, scale, pad)

    body_rx = g.body_a * scale
    body_ry = g.body_b * scale
    head_r = g.head_radius * scale
    weight_r = weight_radius * scale

    beak_poly = _svg_polygon(g.beak, g.x_min, g.y_max, scale, pad)
    tail_poly = _svg_polygon(g.tail, g.x_min, g.y_max, scale, pad)
    wing_up_poly = _svg_polygon(g.wing_up, g.x_min, g.y_max, scale, pad)
    wing_dn_poly = _svg_polygon(g.wing_down, g.x_min, g.y_max, scale, pad)
    rear_up_poly = _svg_polygon(g.rear_up, g.x_min, g.y_max, scale, pad)
    rear_dn_poly = _svg_polygon(g.rear_down, g.x_min, g.y_max, scale, pad)

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect x="0" y="0" width="{width}" height="{height}" fill="#f9fafb"/>
  <g fill="#ffd44d" stroke="#303030" stroke-width="2">
    <ellipse cx="{body_cx:.2f}" cy="{body_cy:.2f}" rx="{body_rx:.2f}" ry="{body_ry:.2f}"/>
    <circle cx="{head_cx:.2f}" cy="{head_cy:.2f}" r="{head_r:.2f}"/>
    <polygon points="{beak_poly}"/>
    <polygon points="{tail_poly}"/>
    <polygon points="{wing_up_poly}"/>
    <polygon points="{wing_dn_poly}"/>
    <polygon points="{rear_up_poly}"/>
    <polygon points="{rear_dn_poly}"/>
    <circle cx="{wt_cx:.2f}" cy="{wt_cy:.2f}" r="{weight_r:.2f}"/>
    <circle cx="{wb_cx:.2f}" cy="{wb_cy:.2f}" r="{weight_r:.2f}"/>
  </g>
  <g stroke-width="3" stroke-linecap="round">
    <line x1="{support_x - 8:.2f}" y1="{support_y - 8:.2f}" x2="{support_x + 8:.2f}" y2="{support_y + 8:.2f}" stroke="#d62828"/>
    <line x1="{support_x - 8:.2f}" y1="{support_y + 8:.2f}" x2="{support_x + 8:.2f}" y2="{support_y - 8:.2f}" stroke="#d62828"/>
  </g>
  <circle cx="{cm_x:.2f}" cy="{cm_y:.2f}" r="5" fill="#006d77"/>
  <text x="40" y="34" font-family="Arial, sans-serif" font-size="22" fill="#111">
    Tortuga 2D - r={weight_radius:.6f}, A={area:.4f}, x_bar={x_bar:.6f}, y_bar={y_bar:.6f}
  </text>
  <text x="40" y="62" font-family="Arial, sans-serif" font-size="17" fill="#333">
    X roja: hocico / punto de apoyo | Punto azul: centro de masa
  </text>
</svg>
"""

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Centro de masa de la tortuga en el hocico.")
    parser.add_argument("--radius", type=float, default=1.623037, help="Radio de engrosamiento en aletas delanteras.")
    parser.add_argument(
        "--solve",
        action="store_true",
        help="Calibra automaticamente el radio para forzar x_bar ~= 0.",
    )
    parser.add_argument(
        "--grid-step",
        type=float,
        default=0.005,
        help="Paso de integracion por grilla (mas pequeno = mas preciso).",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=Path("vista_2d_tortuga.svg"),
        help="Ruta de salida del esquema SVG (omitir con --no-plot).",
    )
    parser.add_argument("--no-plot", action="store_true", help="No generar esquema SVG.")
    args = parser.parse_args()

    radius = args.radius
    if args.solve:
        radius = solve_weight_radius_grid(step=args.grid_step)
        print(f"Radio calibrado (grilla): r = {radius:.6f}")

    area, x_bar, y_bar = centroid_grid(radius, step=args.grid_step)
    distance = float(np.hypot(x_bar, y_bar))

    print(f"Area aproximada: A = {area:.6f}")
    print(f"Centro de masa: x_bar = {x_bar:.6f}, y_bar = {y_bar:.6f}")
    print(f"Distancia al hocico (0,0): {distance:.6f}")

    if not args.no_plot:
        save_svg(radius, args.plot, grid_step=args.grid_step)
        print(f"Esquema guardado en: {args.plot}")


if __name__ == "__main__":
    main()
