#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

echo "1) Recalculando centro de masa..."
python3 calculo_centro_masa.py --solve --plot vista_2d_tortuga.svg

echo "1.1) Verificando centro de masa 3D (muestra rapida)..."
python3 calculo_centro_masa_3d_tortuga.py --samples 2000000

echo "2) Abriendo documentacion HTML..."
open documentacion_proyecto_tortuga.html

echo "3) Abriendo modelo 3D en OpenSCAD..."
if [ -d "/Applications/OpenSCAD-2021.01.app" ]; then
  open -a /Applications/OpenSCAD-2021.01.app modelo_tortuga.scad
elif [ -d "/Applications/OpenSCAD.app" ]; then
  open -a /Applications/OpenSCAD.app modelo_tortuga.scad
else
  echo "No se encontro OpenSCAD en /Applications."
  echo "Instala o reinstala con: brew reinstall --cask openscad"
  exit 1
fi

echo "Listo. Revisa navegador + OpenSCAD."
