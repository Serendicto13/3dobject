#!/usr/bin/env python3
"""
Verificacion 3D (Monte Carlo) del centro de masa de la tortuga realista.

Este script replica la geometria de `modelo_tortuga.scad` con primitivas analiticas
y calcula:
  - Volumen aproximado V
  - Centro de masa (x_bar, y_bar, z_bar)

Uso:
  python3 calculo_centro_masa_3d_tortuga.py
  python3 calculo_centro_masa_3d_tortuga.py --solve
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Params:
    shell_center_x: float = -2.0
    shell_a_x: float = 1.9
    front_lobe_scale: float = 1.383745

    x_min: float = -6.8
    x_max: float = 6.2
    y_min: float = -3.6
    y_max: float = 3.6
    z_min: float = -1.0
    z_max: float = 2.4

    @property
    def bbox_volume(self) -> float:
        return (self.x_max - self.x_min) * (self.y_max - self.y_min) * (self.z_max - self.z_min)


P = Params()


def ellipsoid_mask(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    c: tuple[float, float, float],
    a: tuple[float, float, float],
) -> np.ndarray:
    cx, cy, cz = c
    ax, ay, az = a
    return ((x - cx) / ax) ** 2 + ((y - cy) / ay) ** 2 + ((z - cz) / az) ** 2 <= 1.0


def sphere_mask(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    c: tuple[float, float, float],
    r: float,
) -> np.ndarray:
    cx, cy, cz = c
    return (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r**2


def capsule_mask(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    p0: tuple[float, float, float],
    p1: tuple[float, float, float],
    r: float,
) -> np.ndarray:
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    vx, vy, vz = x1 - x0, y1 - y0, z1 - z0
    wx, wy, wz = x - x0, y - y0, z - z0
    vv = vx * vx + vy * vy + vz * vz
    t = (wx * vx + wy * vy + wz * vz) / vv
    t = np.clip(t, 0.0, 1.0)
    qx = x0 + t * vx
    qy = y0 + t * vy
    qz = z0 + t * vz
    return (x - qx) ** 2 + (y - qy) ** 2 + (z - qz) ** 2 <= r**2


def inside_turtle_realistic(x: np.ndarray, y: np.ndarray, z: np.ndarray, front_lobe_scale: float) -> np.ndarray:
    # Caparazon ahuecado
    outer = ellipsoid_mask(x, y, z, (P.shell_center_x, 0.0, 0.25), (P.shell_a_x, 1.95, 1.25))
    outer |= ellipsoid_mask(x, y, z, (P.shell_center_x - 0.1, 0.0, 0.95), (P.shell_a_x * 0.8, 1.55, 1.0))
    inner = ellipsoid_mask(x, y, z, (P.shell_center_x - 0.05, 0.0, 0.35), (P.shell_a_x - 0.35, 1.55, 0.95))
    mask = outer & (~inner)

    # Cabeza + cuello + hocico
    mask |= ellipsoid_mask(x, y, z, (-0.78, 0.0, 0.0), (0.85, 0.65, 0.5))
    mask |= capsule_mask(x, y, z, (-1.6, 0.0, 0.03), (-1.1, 0.0, 0.0), 0.32)
    mask |= sphere_mask(x, y, z, (0.0, 0.0, 0.0), 0.14)

    # Aletas delanteras
    mask |= capsule_mask(x, y, z, (-1.5, 0.95, -0.06), (2.7, 2.25, -0.20), 0.24)
    mask |= capsule_mask(x, y, z, (-1.5, -0.95, -0.06), (2.7, -2.25, -0.20), 0.24)
    mask |= ellipsoid_mask(x, y, z, (3.1, 2.45, -0.20), (1.75, 0.75, 0.24))
    mask |= ellipsoid_mask(x, y, z, (3.1, -2.45, -0.20), (1.75, 0.75, 0.24))

    # Lobulos de masa frontal
    ax = 1.2 * front_lobe_scale
    ay = 0.82 * front_lobe_scale
    az = 0.36 * front_lobe_scale
    mask |= ellipsoid_mask(x, y, z, (3.55, 2.6, -0.22), (ax, ay, az))
    mask |= ellipsoid_mask(x, y, z, (3.55, -2.6, -0.22), (ax, ay, az))

    # Aletas traseras + cola
    mask |= capsule_mask(x, y, z, (-3.2, 1.02, -0.14), (-5.0, 1.9, -0.28), 0.22)
    mask |= capsule_mask(x, y, z, (-3.2, -1.02, -0.14), (-5.0, -1.9, -0.28), 0.22)
    mask |= capsule_mask(x, y, z, (-3.8, 0.0, -0.12), (-5.5, 0.0, -0.25), 0.13)

    return mask


def centroid_3d(samples: int, seed: int, front_lobe_scale: float) -> tuple[float, float, float, float]:
    rng = np.random.default_rng(seed)
    x = rng.uniform(P.x_min, P.x_max, samples)
    y = rng.uniform(P.y_min, P.y_max, samples)
    z = rng.uniform(P.z_min, P.z_max, samples)
    m = inside_turtle_realistic(x, y, z, front_lobe_scale=front_lobe_scale)
    if not np.any(m):
        raise RuntimeError("No se detecto volumen en la region de muestreo.")

    vol = float(np.mean(m) * P.bbox_volume)
    x_bar = float(np.mean(x[m]))
    y_bar = float(np.mean(y[m]))
    z_bar = float(np.mean(z[m]))
    return vol, x_bar, y_bar, z_bar


def solve_front_lobe_scale(samples: int, seed: int, lo: float = 0.8, hi: float = 2.4) -> float:
    rng = np.random.default_rng(seed)
    x = rng.uniform(P.x_min, P.x_max, samples)
    y = rng.uniform(P.y_min, P.y_max, samples)
    z = rng.uniform(P.z_min, P.z_max, samples)

    def xbar(scale: float) -> float:
        m = inside_turtle_realistic(x, y, z, front_lobe_scale=scale)
        return float(np.mean(x[m]))

    x_lo = xbar(lo)
    x_hi = xbar(hi)
    if x_lo * x_hi > 0:
        raise RuntimeError("No hay cambio de signo en x_bar para el intervalo dado.")

    for _ in range(34):
        mid = 0.5 * (lo + hi)
        x_mid = xbar(mid)
        if x_lo * x_mid <= 0:
            hi = mid
        else:
            lo = mid
            x_lo = x_mid

    return 0.5 * (lo + hi)


def main() -> None:
    parser = argparse.ArgumentParser(description="Centro de masa 3D de la tortuga realista.")
    parser.add_argument("--samples", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=20260414)
    parser.add_argument("--scale", type=float, default=P.front_lobe_scale, help="Escala de lobulo frontal.")
    parser.add_argument("--solve", action="store_true", help="Resuelve escala para x_bar ~= 0.")
    args = parser.parse_args()

    scale = args.scale
    if args.solve:
        scale = solve_front_lobe_scale(samples=args.samples, seed=args.seed)
        print(f"Escala calibrada: front_lobe_scale = {scale:.6f}")

    vol, x_bar, y_bar, z_bar = centroid_3d(samples=args.samples, seed=args.seed + 1, front_lobe_scale=scale)
    print(f"Volumen aproximado: V = {vol:.6f}")
    print(f"Centro de masa 3D: x_bar = {x_bar:.6f}, y_bar = {y_bar:.6f}, z_bar = {z_bar:.6f}")
    print(f"Distancia XY al hocico: {np.hypot(x_bar, y_bar):.6f}")


if __name__ == "__main__":
    main()
