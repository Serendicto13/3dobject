#!/usr/bin/env python3
"""
Verificacion 3D (Monte Carlo) del centro de masa para la tortuga con apoyo en la cabeza.

El modelo replica `modelo_tortuga_estable.scad` con:
  - caparazon hueco
  - cuello descendente
  - cabeza mas baja que el cuerpo
  - base eliptica de apoyo bajo la cabeza

Se calculan:
  - Volumen aproximado
  - Centro de masa 3D
  - Criterio de estabilidad sobre mesa plana o dedo ancho

Uso:
  python3 calculo_centro_masa_3d_tortuga.py
  python3 calculo_centro_masa_3d_tortuga.py --solve --samples 3000000
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace

import numpy as np


@dataclass(frozen=True)
class Params:
    x_min: float = -7.0
    x_max: float = 4.5
    y_min: float = -3.8
    y_max: float = 3.8
    z_min: float = -1.2
    z_max: float = 3.0

    @property
    def bbox_volume(self) -> float:
        return (self.x_max - self.x_min) * (self.y_max - self.y_min) * (self.z_max - self.z_min)


@dataclass(frozen=True)
class ModelConfig:
    front_lobe_scale: float = 1.0 # Will be calibrated
    scale_mm: float = 11.0
    support_ax: float = 0.30
    support_by: float = 0.18
    support_height: float = 0.12
    support_center_x: float = 0.50
    support_center_z: float = -1.44

    @property
    def support_plane_z(self) -> float:
        return self.support_center_z - 0.5 * self.support_height


P = Params()
DEFAULT_CONFIG = ModelConfig()


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


def elliptical_cylinder_mask(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    c: tuple[float, float, float],
    ax: float,
    by: float,
    h: float,
) -> np.ndarray:
    cx, cy, cz = c
    return (((x - cx) / ax) ** 2 + ((y - cy) / by) ** 2 <= 1.0) & (np.abs(z - cz) <= 0.5 * h)


def inside_turtle_head_support(x: np.ndarray, y: np.ndarray, z: np.ndarray, cfg: ModelConfig) -> np.ndarray:
    # Caparazon y boveda
    outer = ellipsoid_mask(x, y, z, (-2.4, 0.0, 1.05), (1.95, 1.90, 0.95))
    outer |= ellipsoid_mask(x, y, z, (-2.55, 0.0, 1.48), (1.45, 1.50, 0.74))
    
    # Cavidad movida fuertemente hacia la parte posterior para aligerar la cola
    inner = ellipsoid_mask(x, y, z, (-3.0, 0.0, 1.08), (2.0, 1.6, 0.8))
    mask = outer & (~inner)

    # Cuello gordo estirado hacia abajo y levemente adelante
    mask |= capsule_mask(x, y, z, (-1.5, 0.0, 0.6), (0.0, 0.0, -0.6), 0.35)
    # Cabeza ancha/plana (turtle head realista)
    mask |= ellipsoid_mask(x, y, z, (0.2, 0.0, -0.9), (0.7, 0.5, 0.3))
    # Hocico apoyado directo hacia la columna / el piso
    mask |= capsule_mask(x, y, z, (0.4, 0.0, -1.0), (0.5, 0.0, -1.4), 0.12)
    
    mask |= elliptical_cylinder_mask(
        x,
        y,
        z,
        (cfg.support_center_x, 0.0, cfg.support_center_z),
        cfg.support_ax,
        cfg.support_by,
        cfg.support_height,
    )

    # Aletas delanteras: Realistas, se extienden desde DE ADENTRO del caparazon hacia adelante
    # Base en (-1.0, 1.0, 0.4) garantizado adentro
    mask |= capsule_mask(x, y, z, (-1.0, 1.0, 0.4), (0.5, 1.8, 0.0), 0.25)
    mask |= capsule_mask(x, y, z, (-1.0, -1.0, 0.4), (0.5, -1.8, 0.0), 0.25)
    
    # Lobulo de aleta marina, posicionado en x=1.2 para dar suficiente contrapeso
    s = cfg.front_lobe_scale
    mask |= ellipsoid_mask(x, y, z, (1.2, 2.2, -0.1), (1.2 * s, 0.8 * s, 0.2 * s))
    mask |= ellipsoid_mask(x, y, z, (1.2, -2.2, -0.1), (1.2 * s, 0.8 * s, 0.2 * s))

    # Patas traseras cortas, conectadas adentro (-3.5, 1.0, 0.5)
    mask |= capsule_mask(x, y, z, (-3.5, 1.0, 0.5), (-4.2, 1.5, 0.2), 0.15)
    mask |= capsule_mask(x, y, z, (-3.5, -1.0, 0.5), (-4.2, -1.5, 0.2), 0.15)
    
    # Cola corta conectada adentro (-4.0, 0.0, 0.8)
    mask |= capsule_mask(x, y, z, (-4.0, 0.0, 0.8), (-5.0, 0.0, 0.2), 0.10)
    return mask


def centroid_3d(samples: int, seed: int, cfg: ModelConfig) -> tuple[float, float, float, float, float]:
    rng = np.random.default_rng(seed)
    x = rng.uniform(P.x_min, P.x_max, samples)
    y = rng.uniform(P.y_min, P.y_max, samples)
    z = rng.uniform(P.z_min, P.z_max, samples)
    mask = inside_turtle_head_support(x, y, z, cfg)
    if not np.any(mask):
        raise RuntimeError("No se detecto volumen en la region de muestreo.")

    vol = float(np.mean(mask) * P.bbox_volume)
    x_bar = float(np.mean(x[mask]))
    y_bar = float(np.mean(y[mask]))
    z_bar = float(np.mean(z[mask]))
    z_min = float(np.min(z[mask]))
    return vol, x_bar, y_bar, z_bar, z_min


def solve_front_lobe_scale(
    samples: int,
    seed: int,
    cfg: ModelConfig,
    lo: float = 1.0,
    hi: float = 4.0,
) -> float:
    rng = np.random.default_rng(seed)
    x = rng.uniform(P.x_min, P.x_max, samples)
    y = rng.uniform(P.y_min, P.y_max, samples)
    z = rng.uniform(P.z_min, P.z_max, samples)

    def xbar(scale: float) -> float:
        test_cfg = replace(cfg, front_lobe_scale=scale)
        mask = inside_turtle_head_support(x, y, z, test_cfg)
        return float(np.mean(x[mask])) - test_cfg.support_center_x

    x_lo = xbar(lo)
    x_hi = xbar(hi)
    if x_lo * x_hi > 0:
        raise RuntimeError("No hay cambio de signo en x_bar para el intervalo dado.")

    for _ in range(38):
        mid = 0.5 * (lo + hi)
        x_mid = xbar(mid)
        if x_lo * x_mid <= 0:
            hi = mid
        else:
            lo = mid
            x_lo = x_mid

    return 0.5 * (lo + hi)


def report_support(x_bar: float, y_bar: float, z_bar: float, z_min: float, cfg: ModelConfig) -> None:
    offset = float(np.hypot(x_bar - cfg.support_center_x, y_bar))
    offset_mm = offset * cfg.scale_mm
    support_test = ((x_bar - cfg.support_center_x) / cfg.support_ax) ** 2 + (y_bar / cfg.support_by) ** 2
    cm_height_mm = (z_bar - cfg.support_plane_z) * cfg.scale_mm

    print("Analisis del apoyo: base eliptica bajo la cabeza")
    print(f"Semieje mayor del apoyo: a = {cfg.support_ax:.6f} ({cfg.support_ax * cfg.scale_mm:.3f} mm)")
    print(f"Semieje menor del apoyo: b = {cfg.support_by:.6f} ({cfg.support_by * cfg.scale_mm:.3f} mm)")
    print(f"Plano de apoyo: z = {cfg.support_plane_z:.6f}")
    print(f"Punto mas bajo muestreado del solido: z_min = {z_min:.6f}")
    print(f"Offset horizontal del CM respecto al centro del apoyo: {offset:.6f} ({offset_mm:.3f} mm)")
    print(f"Altura del CM sobre el plano de apoyo: {cm_height_mm:.3f} mm")
    if support_test <= 1.0 and z_min >= cfg.support_plane_z - 0.01:
        print("Diagnostico: estable sobre mesa plana o dedo ancho")
        print("Motivo: la proyeccion del CM queda dentro de la elipse de apoyo y la cabeza es la parte mas baja.")
    else:
        print("Diagnostico: inestable o ambiguo")
        print("Motivo: la proyeccion del CM cae fuera del apoyo o existe otra zona mas baja que la cabeza.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Centro de masa 3D de la tortuga con apoyo en la cabeza.")
    parser.add_argument("--samples", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=20260416)
    parser.add_argument("--scale", type=float, default=DEFAULT_CONFIG.front_lobe_scale, help="Escala del lobulo frontal.")
    parser.add_argument("--solve", action="store_true", help="Recalibra la escala frontal para x_bar ~= 0.")
    parser.add_argument("--scale-mm", type=float, default=DEFAULT_CONFIG.scale_mm)
    parser.add_argument("--support-ax", type=float, default=DEFAULT_CONFIG.support_ax)
    parser.add_argument("--support-by", type=float, default=DEFAULT_CONFIG.support_by)
    parser.add_argument("--support-height", type=float, default=DEFAULT_CONFIG.support_height)
    parser.add_argument("--support-center-z", type=float, default=DEFAULT_CONFIG.support_center_z)
    args = parser.parse_args()

    cfg = ModelConfig(
        front_lobe_scale=args.scale,
        scale_mm=args.scale_mm,
        support_ax=args.support_ax,
        support_by=args.support_by,
        support_height=args.support_height,
        support_center_z=args.support_center_z,
    )

    if args.solve:
        scale = solve_front_lobe_scale(samples=args.samples, seed=args.seed, cfg=cfg)
        cfg = replace(cfg, front_lobe_scale=scale)
        print(f"Escala calibrada: front_lobe_scale = {scale:.6f}")

    vol, x_bar, y_bar, z_bar, z_min = centroid_3d(samples=args.samples, seed=args.seed + 1, cfg=cfg)
    print("Configuracion: apoyo en la cabeza")
    print(f"Volumen aproximado: V = {vol:.6f}")
    print(f"Centro de masa 3D: x_bar = {x_bar:.6f}, y_bar = {y_bar:.6f}, z_bar = {z_bar:.6f}")
    print(f"Distancia XY al centro del apoyo: {np.hypot(x_bar - cfg.support_center_x, y_bar):.6f}")
    report_support(x_bar, y_bar, z_bar, z_min, cfg)


if __name__ == "__main__":
    main()
