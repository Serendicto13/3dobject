# Proyecto: Tortuga 3D con centro de masa en el hocico

## Archivos
- `modelo_tortuga.scad`: modelo 3D organico (realista, no extrusion plana).
- `modelo_tortuga_plana.scad`: version plana anterior (respaldo).
- `calculo_centro_masa.py`: calculo 2D del centro de masa en la vista plana.
- `calculo_centro_masa_3d_tortuga.py`: verificacion 3D por Monte Carlo.
- `guia_informe.md`: guia para documento y exposicion.

## Verificar centro de masa
```bash
python3 calculo_centro_masa.py --solve --plot vista_2d_tortuga.svg
```

## Verificar centro de masa 3D del modelo realista
```bash
python3 calculo_centro_masa_3d_tortuga.py --samples 2000000
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

Comando directo para abrir OpenSCAD con el modelo:
```bash
open -a /Applications/OpenSCAD-2021.01.app modelo_tortuga.scad
```
Si no funciona por nombre de version:
```bash
open -a /Applications/OpenSCAD.app modelo_tortuga.scad
```
