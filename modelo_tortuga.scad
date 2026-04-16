/*
  Tortuga 3D realista (organica) equilibrada en el hocico.

  Idea:
  - Geometria suave por esferas, capsulas (hull de esferas) y elipsoides.
  - Caparazon ahuecado para reducir masa trasera.
  - Aletas delanteras anchas para desplazar masa hacia adelante.
  - Punto de apoyo en (0,0,0): punta del hocico.
*/

$fn = 120;

// Escala final de impresion
scale_mm = 11;

// Parametro calibrado por Monte Carlo 3D (ver calculo_centro_masa_3d_tortuga.py)
front_lobe_scale = 1.383745;

module ellipsoid(c = [0, 0, 0], a = [1, 1, 1]) {
    translate(c) scale(a) sphere(r = 1);
}

module capsule(p0 = [0, 0, 0], p1 = [1, 0, 0], r = 0.2) {
    hull() {
        translate(p0) sphere(r = r);
        translate(p1) sphere(r = r);
    }
}

module shell_hollow() {
    difference() {
        union() {
            // Caparazon base
            ellipsoid(c = [-2.0, 0, 0.25], a = [1.9, 1.95, 1.25]);
            // Boveda superior
            ellipsoid(c = [-2.1, 0, 0.95], a = [1.52, 1.55, 1.0]);
        }
        // Cavidad interna para quitar masa atras
        ellipsoid(c = [-2.05, 0, 0.35], a = [1.55, 1.55, 0.95]);
    }
}

module turtle_realista() {
    union() {
        shell_hollow();

        // Cuello + cabeza + hocico
        capsule(p0 = [-1.6, 0, 0.03], p1 = [-1.1, 0, 0.0], r = 0.32);
        ellipsoid(c = [-0.78, 0, 0.0], a = [0.85, 0.65, 0.5]);
        sphere(r = 0.14); // apoyo del hocico en (0,0,0)

        // Aletas delanteras (brazos)
        capsule(p0 = [-1.5, 0.95, -0.06], p1 = [2.7, 2.25, -0.20], r = 0.24);
        capsule(p0 = [-1.5, -0.95, -0.06], p1 = [2.7, -2.25, -0.20], r = 0.24);

        // Aletas delanteras (palas anchas)
        ellipsoid(c = [3.1, 2.45, -0.20], a = [1.75, 0.75, 0.24]);
        ellipsoid(c = [3.1, -2.45, -0.20], a = [1.75, 0.75, 0.24]);

        // Lobulos delanteros: masa util para equilibrio, integrados al estilo organico
        ellipsoid(c = [3.55, 2.6, -0.22], a = [1.2 * front_lobe_scale, 0.82 * front_lobe_scale, 0.36 * front_lobe_scale]);
        ellipsoid(c = [3.55, -2.6, -0.22], a = [1.2 * front_lobe_scale, 0.82 * front_lobe_scale, 0.36 * front_lobe_scale]);

        // Aletas traseras
        capsule(p0 = [-3.2, 1.02, -0.14], p1 = [-5.0, 1.9, -0.28], r = 0.22);
        capsule(p0 = [-3.2, -1.02, -0.14], p1 = [-5.0, -1.9, -0.28], r = 0.22);

        // Cola
        capsule(p0 = [-3.8, 0, -0.12], p1 = [-5.5, 0, -0.25], r = 0.13);
    }
}

scale([scale_mm, scale_mm, scale_mm]) turtle_realista();
