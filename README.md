# Proyecto: Tortuga 3D con centro de masa alineado al apoyo

## Archivos
- `modelo_tortuga.scad`: modelo 3D organico (realista, no extrusion plana).
- `modelo_tortuga_estable.scad`: variante solida con apoyo universal para mesa, dedo y punta de lapiz.
- `modelo_tortuga_plana.scad`: version plana anterior (respaldo).
- `calculo_centro_masa.py`: calculo 2D del centro de masa en la vista plana.
- `calculo_centro_masa_3d_tortuga.py`: verificacion 3D por Monte Carlo y diagnostico de estabilidad del apoyo.
- `generar_visualizacion_tortuga.py`: genera el visor 3D interactivo local.
- `tortuga_3d_interactiva.html`: visor 3D embebible de la tortuga.
- `guia_informe.md`: guia para documento y exposicion.
- `documentacion_proyecto_tortuga.html`: desarrollo matematico y visual completo.

## Verificar centro de masa
```bash
python3 calculo_centro_masa.py --solve --plot vista_2d_tortuga.svg
```

## Generar visor 3D interactivo
```bash
python3 generar_visualizacion_tortuga.py
```

## Verificar centro de masa 3D del modelo realista
```bash
python3 calculo_centro_masa_3d_tortuga.py --samples 2000000
```

## Verificar variante universal para mesa, dedo y lapiz
```bash
python3 calculo_centro_masa_3d_tortuga.py --shell-mode solid --support-mode universal --samples 2000000
```

## Abrir documentacion completa (HTML)
```bash
open documentacion_proyecto_tortuga.html
```

## Arranque en un solo comando
```bash
./iniciar_proyecto.sh
```

## Generar STL
1. Abrir `modelo_tortuga.scad` en OpenSCAD.
2. Render (`F6`).
3. Exportar STL (`File -> Export -> Export as STL`).

Para la variante recomendada de apoyo universal:
1. Abrir `modelo_tortuga_estable.scad` en OpenSCAD.
2. Render (`F6`).
3. Exportar STL (`File -> Export -> Export as STL`).
4. Mantener la escala en mm; el modelo ya esta definido con `scale_mm = 11`.

Para forzar STL binario desde terminal:
```bash
openscad --export-format binstl -o modelo_tortuga_estable.stl modelo_tortuga_estable.scad
```

Archivo binario ya generado en este proyecto:
`modelo_tortuga_estable_binario.stl`

Comando directo para abrir OpenSCAD con el modelo:
```bash
open -a /Applications/OpenSCAD-2021.01.app modelo_tortuga.scad
```
Para la variante estable:
```bash
open -a /Applications/OpenSCAD-2021.01.app modelo_tortuga_estable.scad
```
Si no funciona por nombre de version:
```bash
open -a /Applications/OpenSCAD.app modelo_tortuga.scad
```
Y para la variante universal:
```bash
open -a /Applications/OpenSCAD.app modelo_tortuga_estable.scad
```
