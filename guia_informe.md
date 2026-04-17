# Guia rapida: tortuga con centro de masa alineado al apoyo universal

## 1) Historia matematica recomendada
Presenta el proyecto en dos niveles:

1. Modelo 2D analitico para desarrollo a lapiz y papel.
2. Modelo 3D organico para la pieza real y la impresion.

La clave del discurso es:
las aletas y los lobulos delanteros son grandes porque agregan momento positivo en `x`
y empujan el centro de masa hacia la cabeza.

## 2) Modelo 2D usado para el calculo
La lamina se define como union de regiones simples:

1. Caparazon: elipse.
2. Cabeza: circulo.
3. Hocico: triangulo.
4. Cola: triangulo.
5. Aletas delanteras: dos triangulos.
6. Aletas traseras: dos triangulos.
7. Engrosamientos delanteros: dos circulos de radio `r`.

El punto `(0,0)` se define como el eje comun de apoyo.

## 3) Centro de masa en 2D (integrales dobles)
Con densidad superficial constante `rho`:

`A(r) = double_integral_R( dA )`

`Mx(r) = double_integral_R( x dA )`

`My(r) = double_integral_R( y dA )`

`x_bar(r) = Mx(r) / A(r)`

`y_bar(r) = My(r) / A(r)`

Por simetria respecto al eje `x`, se tiene `My(r) = 0`, luego `y_bar(r) = 0`.

En el modelo simplificado:

`A(r) = A0 + 2*pi*r^2`

`Mx(r) = Mx0 + 2*pi*r^2*(4.8)`

La condicion de equilibrio horizontal es:

`x_bar(r) ~= 0`

## 3) Extension a 3D (integrales triples y volumen)
El modelo real ya no es una extrusion plana. El solido se define como union y diferencia de:

1. Elipsoides para el caparazon y la cabeza.
2. Capsulas para cuello, aletas y cola.
3. Un anillo inferior para apoyo plano.
4. Una cavidad esferica para apoyo puntual en lapiz.

La formulacion rigurosa se hace con funcion caracteristica:

`V = triple_integral chi_T(x,y,z) dV`

`Mx = triple_integral x*chi_T(x,y,z) dV`

`My = triple_integral y*chi_T(x,y,z) dV`

`Mz = triple_integral z*chi_T(x,y,z) dV`

`x_bar = Mx/V`, `y_bar = My/V`, `z_bar = Mz/V`

## 4) Resultado numerico actual
Modelo 2D refinado:

- `r ~= 1.623037` (unidades CAD)
- `x_bar ~= 0.000031`
- `y_bar ~= 0.000007`
- Distancia al eje: `~ 0.000032`

Modelo 3D universal:

- `V ~= 35.793248`
- `x_bar ~= -0.000496`
- `y_bar ~= 0.004452`
- `z_bar ~= 0.143180`
- Offset horizontal a escala final: `~ 0.049 mm`

## 5) Como correr la verificacion
```bash
python3 calculo_centro_masa.py --solve --plot vista_2d_tortuga.svg
```

## 6) Modelo 3D para impresion
Archivo principal recomendado:

- `modelo_tortuga_estable.scad`
- Exportar como `STL` binario en `mm`.

Pasos:

1. Abrir en OpenSCAD.
2. Renderizar (`F6`).
3. Exportar STL.
4. Imprimir.

## 7) Verificacion 3D (modelo organico)
Para sustentar integrales triples, volumen y estabilidad:

```bash
python3 calculo_centro_masa_3d_tortuga.py --shell-mode solid --support-mode universal --samples 2000000
```

Se obtiene volumen aproximado y centro de masa 3D `(x_bar, y_bar, z_bar)`.
La distancia `sqrt(x_bar^2 + y_bar^2)` confirma alineacion horizontal del apoyo.

## 8) Texto corto sugerido para el video
"Disenamos primero un modelo 2D de la tortuga para trabajar con integrales dobles y calcular el centro de masa a mano. Luego construimos una version 3D organica definida por elipsoides, capsulas y una cavidad de apoyo, y calculamos su volumen y centro de masa con integrales triples usando una funcion caracteristica. Las aletas y los lobulos delanteros son mas grandes porque desplazan masa hacia la cabeza y hacen que el centro de masa coincida con el eje de apoyo."
