/*
  Modelo 3D de tortuga equilibrada en el hocico.
  - El punto (0, 0) es el punto de apoyo en el hocico.
  - El espesor frontal de las aletas delanteras (circulos) se calibra
    para ubicar el centro de masa cerca de (0, 0).
*/

$fn = 96;

front_flipper_radius = 1.623037; // calibrado con calculo_centro_masa.py
scale_mm = 7;                    // 1 unidad CAD = 7 mm
thickness_mm = 8;                // espesor de extrusion

module turtle_2d(r = front_flipper_radius) {
    union() {
        // Caparazon (elipse)
        translate([-3.4, 0]) scale([3.1, 2.2]) circle(r = 1);

        // Cabeza
        translate([-1.0, 0]) circle(r = 0.95);

        // Hocico (punto de apoyo en (0,0))
        polygon(points = [
            [0.0, 0.0],
            [-0.8, 0.28],
            [-0.8, -0.28]
        ]);

        // Cola
        polygon(points = [
            [-6.5, 0.0],
            [-7.4, 0.65],
            [-7.4, -0.65]
        ]);

        // Aleta delantera superior
        polygon(points = [
            [-1.6, 1.0],
            [5.0, 2.45],
            [1.8, 0.9]
        ]);

        // Aleta delantera inferior
        polygon(points = [
            [-1.6, -1.0],
            [5.0, -2.45],
            [1.8, -0.9]
        ]);

        // Aleta trasera superior
        polygon(points = [
            [-4.2, 1.15],
            [-6.8, 2.2],
            [-5.9, 0.75]
        ]);

        // Aleta trasera inferior
        polygon(points = [
            [-4.2, -1.15],
            [-6.8, -2.2],
            [-5.9, -0.75]
        ]);

        // Engrosamiento frontal de aletas (ajuste de masa)
        translate([4.8, 2.1]) circle(r = r);
        translate([4.8, -2.1]) circle(r = r);
    }
}

linear_extrude(height = thickness_mm, center = true, convexity = 10)
    scale([scale_mm, scale_mm]) turtle_2d(front_flipper_radius);

// Marcador visual del apoyo
%translate([0, 0, thickness_mm / 2 + 0.5]) cylinder(h = 1, r = 0.35, center = true);
