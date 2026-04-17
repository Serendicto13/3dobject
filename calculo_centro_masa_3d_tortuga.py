#!/usr/bin/env python3
"""
Verificacion 3D (Monte Carlo) del centro de masa de la tortuga.

La geometria replica EXACTAMENTE modelo_tortuga_estable.scad:
  - Caparazon hueco (dos elipsoides menos cavidad interior)
  - Cuello descendente y cabeza baja
  - Pad eliptico de apoyo bajo la cabeza
  - Aletas delanteras con lobulos variables (front_lobe_scale)
  - Aletas traseras y cola

Uso:
  python calculo_centro_masa_3d_tortuga.py
  python calculo_centro_masa_3d_tortuga.py --solve --samples 3000000
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace

import numpy as np


# ---------------------------------------------------------------------------
# Bounding box del modelo (unidades CAD, sin escalar)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Params:
    x_min: float = -6.0
    x_max: float = 4.8
    y_min: float = -4.0
    y_max: float = 4.0
    z_min: float = -0.7
    z_max: float = 2.5

    @property
    def bbox_volume(self) -> float:
        return (self.x_max - self.x_min) * (self.y_max - self.y_min) * (self.z_max - self.z_min)


# ---------------------------------------------------------------------------
# Parametros del modelo (deben coincidir con modelo_tortuga_estable.scad)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ModelConfig:
    front_lobe_scale: float = 1.685989   # front_lobe_scale del SCAD
    scale_mm: float = 11.0
    # Pad de apoyo (del SCAD: support_center_z = -0.47, translate([0, 0, -0.47]))
    support_ax: float = 0.30
    support_by: float = 0.18
    support_height: float = 0.12
    support_center_x: float = 0.0        # el pad esta en x=0, y=0
    support_center_z: float = -0.47      # support_center_z del SCAD

    @property
    def support_plane_z(self) -> float:
        """Plano inferior del pad de apoyo."""
        return self.support_center_z - 0.5 * self.support_height


P = Params()
DEFAULT_CONFIG = ModelConfig()


# ---------------------------------------------------------------------------
# Primitivas geometricas
# ---------------------------------------------------------------------------
def ellipsoid_mask(
    x: np.ndarray, y: np.ndarray, z: np.ndarray,
    c: tuple[float, float, float],
    a: tuple[float, float, float],
) -> np.ndarray:
    cx, cy, cz = c
    ax, ay, az = a
    return ((x - cx) / ax) ** 2 + ((y - cy) / ay) ** 2 + ((z - cz) / az) ** 2 <= 1.0


def capsule_mask(
    x: np.ndarray, y: np.ndarray, z: np.ndarray,
    p0: tuple[float, float, float],
    p1: tuple[float, float, float],
    r: float,
) -> np.ndarray:
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    vx, vy, vz = x1 - x0, y1 - y0, z1 - z0
    wx, wy, wz = x - x0, y - y0, z - z0
    vv = vx * vx + vy * vy + vz * vz
    t = np.clip((wx * vx + wy * vy + wz * vz) / vv, 0.0, 1.0)
    qx = x0 + t * vx
    qy = y0 + t * vy
    qz = z0 + t * vz
    return (x - qx) ** 2 + (y - qy) ** 2 + (z - qz) ** 2 <= r ** 2


def elliptical_cylinder_mask(
    x: np.ndarray, y: np.ndarray, z: np.ndarray,
    c: tuple[float, float, float],
    ax: float, by: float, h: float,
) -> np.ndarray:
    cx, cy, cz = c
    return (((x - cx) / ax) ** 2 + ((y - cy) / by) ** 2 <= 1.0) & (np.abs(z - cz) <= 0.5 * h)


# ---------------------------------------------------------------------------
# Modelo completo — replica exacta de modelo_tortuga_estable.scad
# ---------------------------------------------------------------------------
def inside_turtle(
    x: np.ndarray, y: np.ndarray, z: np.ndarray,
    cfg: ModelConfig,
) -> np.ndarray:
    """
    Replica de turtle_head_support() en modelo_tortuga_estable.scad.

    shell_hollow():
        outer = union(
            ellipsoid(c=[-2.4, 0, 1.05], a=[1.95, 1.90, 0.95]),
            ellipsoid(c=[-2.55, 0, 1.48], a=[1.45, 1.50, 0.74])
        )
        inner = ellipsoid(c=[-2.42, 0, 1.08], a=[1.48, 1.45, 0.62])
        shell = outer minus inner

    union con:
        capsule([-1.55,  0,    0.68], [-0.45,  0,    0.18], r=0.25)  cuello
        ellipsoid([-0.18, 0, -0.06], [0.62, 0.42, 0.30])              cabeza
        capsule([-0.02,  0,   -0.16], [ 0.18,  0,   -0.30], r=0.11)  hocico
        ellip_cylinder([0, 0, -0.47], ax=0.30, by=0.18, h=0.12)       pad apoyo
        capsule([-1.55,  0.95, 0.24], [1.45,  2.15, 0.08], r=0.20)   aleta izq
        capsule([-1.55, -0.95, 0.24], [1.45, -2.15, 0.08], r=0.20)   aleta der
        ellipsoid([1.85,  2.30, 0.02], [1.45, 0.78, 0.20])            extremo izq
        ellipsoid([1.85, -2.30, 0.02], [1.45, 0.78, 0.20])            extremo der
        ellipsoid([2.15,  2.25,-0.02], [1.25s, 0.82s, 0.30s])         lobulo izq
        ellipsoid([2.15, -2.25,-0.02], [1.25s, 0.82s, 0.30s])         lobulo der
        capsule([-3.00,  1.0,  0.16], [-4.45,  1.82, 0.02], r=0.16)  aleta trasera izq
        capsule([-3.00, -1.0,  0.16], [-4.45, -1.82, 0.02], r=0.16)  aleta trasera der
        capsule([-4.05,  0.0,  0.12], [-5.25,  0.0,  0.02], r=0.10)  cola
    """
    s = cfg.front_lobe_scale

    # -- Caparazon hueco --
    outer = ellipsoid_mask(x, y, z, (-2.4,  0.0, 1.05), (1.95, 1.90, 0.95))
    outer |= ellipsoid_mask(x, y, z, (-2.55, 0.0, 1.48), (1.45, 1.50, 0.74))
    inner  = ellipsoid_mask(x, y, z, (-2.42, 0.0, 1.08), (1.48, 1.45, 0.62))
    mask = outer & (~inner)

    # -- Cuello --
    mask |= capsule_mask(x, y, z, (-1.55, 0.0,  0.68), (-0.45, 0.0,  0.18), 0.25)
    # -- Cabeza --
    mask |= ellipsoid_mask(x, y, z, (-0.18, 0.0, -0.06), (0.62, 0.42, 0.30))
    # -- Hocico --
    mask |= capsule_mask(x, y, z, (-0.02, 0.0, -0.16), (0.18, 0.0, -0.30), 0.11)
    # -- Pad de apoyo --
    mask |= elliptical_cylinder_mask(
        x, y, z,
        (0.0, 0.0, cfg.support_center_z),
        cfg.support_ax, cfg.support_by, cfg.support_height,
    )
    # -- Aletas delanteras --
    mask |= capsule_mask(x, y, z, (-1.55,  0.95, 0.24), (1.45,  2.15, 0.08), 0.20)
    mask |= capsule_mask(x, y, z, (-1.55, -0.95, 0.24), (1.45, -2.15, 0.08), 0.20)
    mask |= ellipsoid_mask(x, y, z, (1.85,  2.30, 0.02), (1.45, 0.78, 0.20))
    mask |= ellipsoid_mask(x, y, z, (1.85, -2.30, 0.02), (1.45, 0.78, 0.20))
    # -- Lobulos delanteros (variable de diseno) --
    mask |= ellipsoid_mask(x, y, z, (2.15,  2.25, -0.02), (1.25*s, 0.82*s, 0.30*s))
    mask |= ellipsoid_mask(x, y, z, (2.15, -2.25, -0.02), (1.25*s, 0.82*s, 0.30*s))
    # -- Aletas traseras --
    mask |= capsule_mask(x, y, z, (-3.00,  1.0, 0.16), (-4.45,  1.82, 0.02), 0.16)
    mask |= capsule_mask(x, y, z, (-3.00, -1.0, 0.16), (-4.45, -1.82, 0.02), 0.16)
    # -- Cola --
    mask |= capsule_mask(x, y, z, (-4.05, 0.0, 0.12), (-5.25, 0.0, 0.02), 0.10)

    return mask


# ---------------------------------------------------------------------------
# Calculo del centroide por Monte Carlo
# ---------------------------------------------------------------------------
def centroid_3d(
    samples: int, seed: int, cfg: ModelConfig,
) -> tuple[float, float, float, float, float]:
    rng = np.random.default_rng(seed)
    x = rng.uniform(P.x_min, P.x_max, samples)
    y = rng.uniform(P.y_min, P.y_max, samples)
    z = rng.uniform(P.z_min, P.z_max, samples)
    mask = inside_turtle(x, y, z, cfg)
    if not np.any(mask):
        raise RuntimeError("No se detecto volumen en la region de muestreo.")
    vol   = float(np.mean(mask) * P.bbox_volume)
    x_bar = float(np.mean(x[mask]))
    y_bar = float(np.mean(y[mask]))
    z_bar = float(np.mean(z[mask]))
    z_min = float(np.min(z[mask]))
    return vol, x_bar, y_bar, z_bar, z_min


# ---------------------------------------------------------------------------
# Calibracion de front_lobe_scale
# ---------------------------------------------------------------------------
def solve_front_lobe_scale(
    samples: int, seed: int, cfg: ModelConfig,
    lo: float = 1.0, hi: float = 4.0,
) -> float:
    """Biseccion: encuentra front_lobe_scale tal que x_bar == support_center_x."""
    rng = np.random.default_rng(seed)
    x = rng.uniform(P.x_min, P.x_max, samples)
    y = rng.uniform(P.y_min, P.y_max, samples)
    z = rng.uniform(P.z_min, P.z_max, samples)

    def xbar_offset(scale: float) -> float:
        test_cfg = replace(cfg, front_lobe_scale=scale)
        mask = inside_turtle(x, y, z, test_cfg)
        return float(np.mean(x[mask])) - cfg.support_center_x

    x_lo = xbar_offset(lo)
    x_hi = xbar_offset(hi)
    if x_lo * x_hi > 0:
        raise RuntimeError(
            f"Sin cambio de signo en [{lo}, {hi}]: "
            f"x_bar(lo)={x_lo:.4f}, x_bar(hi)={x_hi:.4f}"
        )

    for _ in range(42):
        mid = 0.5 * (lo + hi)
        x_mid = xbar_offset(mid)
        if x_lo * x_mid <= 0:
            hi = mid
        else:
            lo = mid
            x_lo = x_mid

    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Reporte de estabilidad
# ---------------------------------------------------------------------------
def report_support(
    x_bar: float, y_bar: float, z_bar: float, z_min: float, cfg: ModelConfig,
) -> None:
    offset    = float(np.hypot(x_bar - cfg.support_center_x, y_bar))
    offset_mm = offset * cfg.scale_mm
    elipse_test = (
        ((x_bar - cfg.support_center_x) / cfg.support_ax) ** 2
        + (y_bar / cfg.support_by) ** 2
    )
    cm_height_mm = (z_bar - cfg.support_plane_z) * cfg.scale_mm

    print(f"\nAnalisis del apoyo — pad eliptico en (x={cfg.support_center_x}, y=0, z={cfg.support_center_z})")
    print(f"  Semiejes: a={cfg.support_ax} ({cfg.support_ax * cfg.scale_mm:.2f} mm), "
          f"b={cfg.support_by} ({cfg.support_by * cfg.scale_mm:.2f} mm)")
    print(f"  Plano inferior del pad: z = {cfg.support_plane_z:.4f}")
    print(f"  Punto mas bajo del solido: z_min = {z_min:.4f}")
    print(f"  Offset horizontal CM vs centro apoyo: {offset:.6f} ({offset_mm:.3f} mm)")
    print(f"  Altura del CM sobre el plano de apoyo: {cm_height_mm:.2f} mm")
    print(f"  Test elipse de apoyo: {elipse_test:.4f} (debe ser <= 1.0)")

    pad_is_lowest = z_min >= cfg.support_plane_z - 0.02
    if elipse_test <= 1.0 and pad_is_lowest:
        print("  Diagnostico: ESTABLE — la proyeccion del CM queda dentro del pad de apoyo")
        print("  y la cabeza es la parte mas baja del solido.")
    else:
        motivos = []
        if elipse_test > 1.0:
            motivos.append(f"proyeccion del CM fuera del pad (test={elipse_test:.4f})")
        if not pad_is_lowest:
            motivos.append(
                f"hay zona mas baja que el pad "
                f"(z_min={z_min:.4f} < {cfg.support_plane_z - 0.02:.4f})"
            )
        print(f"  Diagnostico: INESTABLE — {'; '.join(motivos)}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Centro de masa 3D (Monte Carlo) de la tortuga — replica de modelo_tortuga_estable.scad"
    )
    parser.add_argument("--samples",  type=int,   default=2_000_000,
                        help="Numero de muestras Monte Carlo (default: 2000000)")
    parser.add_argument("--seed",     type=int,   default=20260416)
    parser.add_argument("--scale",    type=float, default=DEFAULT_CONFIG.front_lobe_scale,
                        help="front_lobe_scale (default: valor del SCAD)")
    parser.add_argument("--solve",    action="store_true",
                        help="Calibra front_lobe_scale para que x_bar == 0")
    args = parser.parse_args()

    cfg = replace(DEFAULT_CONFIG, front_lobe_scale=args.scale)

    if args.solve:
        print(f"Calibrando front_lobe_scale con {args.samples:,} muestras ...")
        scale = solve_front_lobe_scale(samples=args.samples, seed=args.seed, cfg=cfg)
        cfg = replace(cfg, front_lobe_scale=scale)
        print(f"front_lobe_scale calibrado = {scale:.6f}")
        print(f"  >> Actualiza esta linea en modelo_tortuga_estable.scad:")
        print(f"  >>   front_lobe_scale = {scale:.6f};")

    vol, x_bar, y_bar, z_bar, z_min = centroid_3d(
        samples=args.samples, seed=args.seed + 1, cfg=cfg
    )

    print(f"\nfront_lobe_scale = {cfg.front_lobe_scale:.6f}")
    print(f"Volumen aproximado: V = {vol:.4f}  (escalado: {vol * cfg.scale_mm**3 / 1000:.2f} ml)")
    print(f"Centro de masa 3D:  x_bar = {x_bar:.6f},  y_bar = {y_bar:.6f},  z_bar = {z_bar:.6f}")
    print(f"Distancia XY al centro del apoyo: {np.hypot(x_bar - cfg.support_center_x, y_bar):.6f}")
    report_support(x_bar, y_bar, z_bar, z_min, cfg)


if __name__ == "__main__":
    main()
