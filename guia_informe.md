# Guia rapida: tortuga con centro de masa en el hocico

## 1) Modelo 2D usado para el calculo
Se uso una lamina de densidad uniforme definida por la union de regiones simples:

1. Caparazon: elipse.
2. Cabeza: circulo.
3. Hocico: triangulo.
4. Cola: triangulo.
5. Aletas delanteras: dos triangulos.
6. Aletas traseras: dos triangulos.
7. Engrosamientos delanteros: dos circulos de radio `r`.

El punto `(0,0)` se define como el punto de apoyo en el hocico.

## 2) Centro de masa en 2D (integrales dobles)
Con densidad superficial constante `rho`:

`M = rho * double_integral_R( dA )`

`Mx = rho * double_integral_R( x dA )`

`My = rho * double_integral_R( y dA )`

`x_bar = Mx/M = double_integral_R( x dA ) / double_integral_R( dA )`

`y_bar = My/M = double_integral_R( y dA ) / double_integral_R( dA )`

La condicion de equilibrio en el apoyo es:

`x_bar ~= 0` y `y_bar ~= 0`.

## 3) Extension a 3D (integrales triples y volumen)
El solido se construye por extrusion uniforme de la region 2D:

`S = { (x,y,z) : (x,y) in R, z in [-t/2, t/2] }`

`V = triple_integral_S( dV ) = t * A`

Con extrusion centrada:

`x_bar_3D = x_bar_2D`

`y_bar_3D = y_bar_2D`

`z_bar_3D = 0`

Por eso, si la tortuga 2D esta calibrada en el hocico, la pieza 3D tambien.

## 4) Resultado numerico actual
Al calibrar `r`:

- `r ~= 1.623037` (unidades CAD)
- `x_bar ~= 0.000031`
- `y_bar ~= 0.000007`
- Distancia al hocico: `~ 0.000032`

## 5) Como correr la verificacion
```bash
python3 calculo_centro_masa.py --solve --plot vista_2d_tortuga.svg
```

## 6) Modelo 3D para impresion
Archivo principal (realista):

- `modelo_tortuga.scad`
- `modelo_tortuga_plana.scad` queda como respaldo.

Pasos:

1. Abrir en OpenSCAD.
2. Renderizar (`F6`).
3. Exportar STL.
4. Imprimir.

## 7) Verificacion 3D (modelo organico)
Para sustentar integrales triples y volumen del modelo realista:

```bash
python3 calculo_centro_masa_3d_tortuga.py --samples 2000000
```

Se obtiene volumen aproximado y centro de masa 3D `(x_bar, y_bar, z_bar)`.
La distancia `sqrt(x_bar^2 + y_bar^2)` confirma equilibrio horizontal en el hocico.

## 8) Texto corto sugerido para el video
"Disenamos una tortuga 2D para el analisis con integrales dobles y ubicamos su centro de masa en el hocico. En paralelo construimos una tortuga 3D organica con caparazon abombado, cabeza y aletas suaves, y verificamos su volumen y centro de masa con un modelo numerico 3D basado en integrales triples."
