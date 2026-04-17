/*
  Tortuga 3D con apoyo en la cabeza.

  Idea fisica:
  - El caparazon queda levantado.
  - La cabeza sale hacia adelante y desciende.
  - La parte mas baja del objeto es una pequena base eliptica bajo la cabeza.

  Asi, sobre mesa plana o sobre un dedo ancho, la tortuga descansa visualmente en la cabeza
  y el resto del cuerpo permanece horizontal.
*/

$fn = 160;

scale_mm = 11;
front_lobe_scale = 1.685989;

support_ax = 0.30;
support_by = 0.18;
support_height = 0.12;
support_center_z = -0.47;

module ellipsoid(c = [0, 0, 0], a = [1, 1, 1]) {
    translate(c) scale(a) sphere(r = 1);
}

module capsule(p0 = [0, 0, 0], p1 = [1, 0, 0], r = 0.2) {
    hull() {
        translate(p0) sphere(r = r);
        translate(p1) sphere(r = r);
    }
}

module support_head_pad() {
    translate([0, 0, support_center_z])
        scale([support_ax, support_by, 1])
        cylinder(h = support_height, r = 1, center = true);
}

module shell_hollow() {
    difference() {
        union() {
            ellipsoid(c = [-2.4, 0, 1.05], a = [1.95, 1.90, 0.95]);
            ellipsoid(c = [-2.55, 0, 1.48], a = [1.45, 1.50, 0.74]);
        }
        ellipsoid(c = [-2.42, 0, 1.08], a = [1.48, 1.45, 0.62]);
    }
}

module turtle_head_support() {
    union() {
        shell_hollow();

        // Cuello descendente y cabeza mas baja que el caparazon
        capsule(p0 = [-1.55, 0, 0.68], p1 = [-0.45, 0, 0.18], r = 0.25);
        ellipsoid(c = [-0.18, 0, -0.06], a = [0.62, 0.42, 0.30]);
        capsule(p0 = [-0.02, 0, -0.16], p1 = [0.18, 0, -0.30], r = 0.11);

        // Base plana de apoyo bajo la cabeza
        support_head_pad();

        // Aletas delanteras elevadas
        capsule(p0 = [-1.55, 0.95, 0.24], p1 = [1.45, 2.15, 0.08], r = 0.20);
        capsule(p0 = [-1.55, -0.95, 0.24], p1 = [1.45, -2.15, 0.08], r = 0.20);
        ellipsoid(c = [1.85, 2.30, 0.02], a = [1.45, 0.78, 0.20]);
        ellipsoid(c = [1.85, -2.30, 0.02], a = [1.45, 0.78, 0.20]);

        // Engrosamientos delanteros para llevar el centro de masa hacia la cabeza
        ellipsoid(c = [2.15, 2.25, -0.02], a = [1.25 * front_lobe_scale, 0.82 * front_lobe_scale, 0.30 * front_lobe_scale]);
        ellipsoid(c = [2.15, -2.25, -0.02], a = [1.25 * front_lobe_scale, 0.82 * front_lobe_scale, 0.30 * front_lobe_scale]);

        // Aletas traseras y cola, mas altas que el plano de apoyo
        capsule(p0 = [-3.00, 1.0, 0.16], p1 = [-4.45, 1.82, 0.02], r = 0.16);
        capsule(p0 = [-3.00, -1.0, 0.16], p1 = [-4.45, -1.82, 0.02], r = 0.16);
        capsule(p0 = [-4.05, 0.0, 0.12], p1 = [-5.25, 0.0, 0.02], r = 0.10);
    }
}

scale([scale_mm, scale_mm, scale_mm]) turtle_head_support();
