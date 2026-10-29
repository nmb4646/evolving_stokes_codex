# Shaded curvature and curved velocity arrows

The current renderer uses PyVista/VTK with saturated diffuse color, restrained Phong highlights, a subtle material texture, and surface-conforming white arrows. The colormap is applied directly, with no extra gamma darkening; two lights provide depth without broad environment reflections washing out the colors. The MATLAB entry point `curvature_flow_panels` now calls this renderer.

## Run the sample

The local Python environment is already installed in `remesh/.venv-curvature-flow`.

```matlab
addpath('/home/nickbroussinos/Code/evolving_stokes_codex/remesh')
curvature_flow_panels( ...
    '/home/nickbroussinos/Code/evolving_stokes_codex/remesh/data/fs_batch_data/Sd_3.00ep00_Da_0.00ep00_gamy_p3.00ep01_v_8.01em01', ...
    '/home/nickbroussinos/Documents/ChatGPT/Membranes/curvature_flow_viz/shear_curvature_flow_panels.png', ...
    [35 85 135 185])
```

The same renderer can be called directly from the remesh directory:

```bash
.venv-curvature-flow/bin/python curvature_flow_render.py \
  data/fs_batch_data/Sd_3.00ep00_Da_0.00ep00_gamy_p3.00ep01_v_8.01em01 \
  --output /home/nickbroussinos/Documents/ChatGPT/Membranes/curvature_flow_viz \
  --frames 35 85 135 185 --resolution 1100
```

Use one frame for a close-up or four for a combined figure. Controls include `--smoothing`, `--subdivisions`, `--arrows`, `--resolution`, and `--view X Y Z`. The MATLAB wrapper also accepts an optional fourth argument with a Python executable, or reads `MEMBRANE_RENDER_PYTHON`.

## What the figure shows

- Color starts from the same signed edge-dihedral mean curvature as `Geometry.v_mean_curvature ./ Geometry.v_area`. The Python calculation was checked against MATLAB on geo85: maximum difference 1.78e-15.
- Velocity is the tangential component of the saved velocity. Both curvature and tangent velocity undergo four immediate-neighbor averaging passes by default, with self weight 2.
- Three Loop subdivision levels create a smooth display surface and interpolate the fields. This rounds the coarse simulation triangulation; the display mesh is for rendering, not numerical measurement. Source MAT files are read only.
- Each short arrow follows the averaged tangent velocity using midpoint integration. Every intermediate point is projected onto the refined surface. Arrowheads also curve along that surface; they are not flat triangles in a single tangent plane.
- Arrow arc length scales with seed speed on one shared scale across the selected frames, capped at 1.25 times the reference length. Seeds below 10% of the reference speed are omitted. The reference is the pooled 95th percentile.
- Camera, lighting, curvature limits, and velocity scale are shared across the four panels. The color limits use extrema of the averaged source-vertex curvature over the selected snapshots.
- The subtle grain is an illustrative material normal map affecting highlights only. It is not measured membrane roughness or a displacement of the mesh.
- Physical time labels are used only when the saved run is complete from geo0 and its saved dt is constant. Otherwise labels use the source frame number.

The sample is the shear run Sd=3, Da=0, gamy=30, v=0.801. It uses geo35, geo85, geo135, and geo185.

## Outputs and dependencies

Each output directory contains the combined PNG, an individual `geoN_curved.png` per snapshot, and `render_manifest.json` with source paths, smoothing settings, shared limits, camera direction, and arrow statistics.

To recreate the environment on Linux:

```bash
python3 -m venv --without-pip .venv-curvature-flow
python3 -m pip --python .venv-curvature-flow/bin/python install -r curvature_flow_requirements.txt
```

This setup uses PyVista 0.46.5 and VTK 9.5.2. Houdini is not required for this renderer.
