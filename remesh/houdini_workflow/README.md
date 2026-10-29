# Membrane render and animation workflow

One command reads MATLAB geometry, renders the approved green Houdini look,
writes cleaned image copies, and creates a GIF directly. Add --mp4 to also
create an H.264 video.

## Run the 20-frame shear sample

From this directory:

    python3 workflow.py render --data ../data/fs_batch_data/Sd_3.00ep00_Da_0.00ep00_gamy_p3.00ep01_v_8.01em01 --output outputs/shear_green_sample --frames 20 --resolution 512 --samples 3 --threads 4 --fps 5 --name shear_green --mp4

The sample uses source steps:
0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 99, 109, 119, 129, 139, 149,
159, 169, 179, 189.

This is 20 renders over simulation time 0 to 0.0567, with shear rate 30.
The source solver defines background flow u_x = gamy * z. The membrane
rotates by approximately 60 degrees during this interval.

The GIF plays forward at 5 fps (4 seconds), then restarts. Playback speed is
a presentation choice, independent of simulation time. No synthetic rotation,
per-frame alignment, reverse playback, or interpolated geometry is added.

## Outputs

Within the chosen output directory:

- raw/frame_0001.png through frame_0020.png: untouched Mantra RGBA renders.
- clean_frames/frame_000000.png through frame_000019.png: white background,
  covered bottom-right logo area, and time labels.
- shear_green.gif: looping animation with one palette shared by all frames.
- shear_green.mp4: optional H.264 video.
- animation.hipnc: editable scene; timeline frames 1 through 20 select the
  listed source timesteps.
- render_plan.json: exact source files and hashes, camera, rendering settings,
  physical timing information, and source mesh counts.
- manifest.json: render progress, per-frame import validation, output hashes,
  and final animation settings.

Change the MAT Data Directory or Source Geometry Indices on
/obj/membrane/CONTROLS in the saved Houdini scene. The render node is
/out/render, and the camera is /obj/cam1.

## Assemble an existing PNG folder

    python3 workflow.py assemble --input /path/to/render_frames --output /path/to/movie_output --fps 10 --name membrane --mp4

Both integer filenames (1.png, 2.png, 10.png) and prefixed filenames
(render001.png, render002.png) sort numerically. Conflicting frame IDs or mixed
filename families produce an error instead of silently changing frame order.
Images with different sizes are centered on a common even-sized white canvas.

The corner box automatically detects the detached logo band in transparent
Houdini PNGs. It checks for a gap from the subject and covers the logo's full
bounds with a small margin. This adapts to the logo's size at different render
resolutions. For opaque images or inconclusive detection it uses the legacy
273 x 79 pixel box. Set --mask WIDTH_FRACTION HEIGHT_FRACTION to override the
box, or --keep-logo to disable it. The box is applied before padding. Keep the
subject clear of this corner when using the legacy fallback or a custom box;
the sample camera leaves the detected logo area empty.

The GIF is encoded directly from cleaned PNGs; there is no lossy intermediate
AVI. All frames, including visually identical successive frames, are kept.
A common palette and disabled dithering reduce palette changes between frames.
GIF timing is rounded to the format's 10 ms increments and the actual rate
is recorded. MP4 uses the requested frame rate.

## Render controls and resume

Use --start / --end to limit the source step interval, or supply an explicit
increasing list with --indices "0,10,20,30". --frames evenly samples the
available states within the interval. Missing requested states cause an error.

The orthographic camera is fitted to the union of every selected geometry,
then held fixed. --view X Y Z sets its direction from the subject; the default
(-1.2, -1.4, 1.0) gives an oblique view of this run's XZ shear plane.
--padding controls the camera margin. The source simulation coordinates
and orientation are preserved throughout.

Increase --resolution (up to 960 for this noncommercial template) and --samples
for a finer render. --threads limits the Houdini and Mantra CPU thread count.

Rerun the identical render command with --resume after an interruption.
Existing frames are reused only when the plan, source files, template, renderer
code, and recorded PNG hashes match. Completed outputs can also be reassembled
with different timing via the assemble command and --force. Raw renders are
never changed by assembly. Use a new output directory for a different render plan.

Time labels are inferred only when all saved states from geo0 are present and
p.dt is constant throughout. Otherwise labels show the source step, unless
--dt explicitly supplies a known constant timestep. The solver's p.total_time
field records wall-clock computation time and is never used as physical time.
Use --no-labels to omit text from rendered animation frames.

## Green appearance and dependencies

assets/green_membrane.hipnc contains the approved original green material,
lighting, fixed teal-to-green curvature ramp, and both remesh stages. This
workflow's display geometry is smoothed and remeshed like the legacy scenes;
the MAT files are read-only. The green material uses its original gamma 1.0
output setting. Image texture and bump mapping are disabled in that material.

The packaged template was derived from evolving_stokes_2/remesh/green_membrane.hipnc.
The normal workflow does not need the old repository or the previous Codex
workspace. create_template.py can rebuild the template from a compatible
legacy scene when run with Houdini's hython.

Requirements:
- Python 3.10+ with Pillow for the CLI and animation assembly.
- Houdini with hython and Mantra, plus SciPy/NumPy in its Python environment.
- FFmpeg with libx264 or libopenh264 only when requesting MP4.

These are already available on this workstation. Houdini and FFmpeg are
detected automatically, including installations under /opt/hfs*. Override
with --hython or --ffmpeg if needed.

## Verification

    python3 -m unittest discover -s . -p test_animation.py -v

Tests cover numeric ordering, duplicate IDs, alpha compositing, mixed sizes,
corner masking, preservation of source images, exact GIF frame counts,
shared palettes, timing, and output protection. Every actual render also
checks imported coordinates and triangle connectivity against its source MAT.


## Smooth green previews from very coarse saved meshes

`green_surface_samples.py` provides a PyVista/VTK preview when a saved mesh has
become too coarse for a clean silhouette. It applies Taubin smoothing and Loop
subdivision, restores the source volume and volume centroid, and exports the
refined surface as PLY for a Houdini File SOP. Its uniform teal-green material
approximates the approved look; it does not run the Houdini shader.

From this directory:

    ../.venv-curvature-flow/bin/python green_surface_samples.py ../data/fs_batch_data/Sd_3.00ep00_Da_0.00ep00_gamy_p3.01ep00_v_8.02em01 --output /home/nickbroussinos/Documents/ChatGPT/Membranes/green_late_samples --frames 5000 9000

Outputs include individual renders, a two-sample figure, source/refined
comparisons, PLY meshes, and a manifest with source hashes and geometric checks.
Use --smoothing-iterations, --pass-band, and --subdivisions to adjust the display
reconstruction. This smoothing removes visual faceting; it cannot recover
geometry absent from the original coarse mesh.

### PLY normals and residual pinches

PLY exports from `green_surface_samples.py` contain geometry without imported
VTK normals. In this Houdini import path the supplied VTK normals opposed the
face normals and produced incorrect pointed specular highlights. Add a Normal
SOP after a File SOP, computing point normals with a 180-degree cusp angle.

Loop subdivision may still preserve local pinches in very coarse irregular
meshes. The corrected late-time example uses Houdini VDB Gaussian smoothing
(voxel size 0.02, radius 4 voxels, 4 iterations), then restores source volume and
centroid and recomputes normals. Its live editable network and reproduction
script are in:

`/home/nickbroussinos/Documents/ChatGPT/Membranes/green_late_samples/houdini_rounded/`

## Rounded 60-state animation with lighting B

`rounded_animation.py` uses an approved rounded Houdini scene as the lighting
and material template. It selects evenly spaced saved states, prepares the
initial smooth meshes, applies the corrected Houdini VDB reconstruction,
restores each state's source volume and centroid, recomputes normals, and
caches the final surfaces. A single camera is fitted to all final states.

Run from this directory with the project's PyVista environment:

```bash
../.venv-curvature-flow/bin/python rounded_animation.py \
  --data ../data/fs_batch_data/Sd_3.00ep00_Da_0.00ep00_gamy_p3.01ep00_v_8.02em01 \
  --template /home/nickbroussinos/Documents/ChatGPT/Membranes/green_late_samples/lighting_geo9000/B_front.hipnc \
  --output outputs/shear_B_60_full_run \
  --frames 60 --resolution 384 --samples 3 --threads 8 --fps 5 --name shear_B_60
```

The sample covers source steps 0 through 10000 (physical time 0 through
29.9003322) and plays forward at 5 fps for a 12-second loop. The restart is a
loop boundary; no interpolated states, synthetic rotation, or reverse playback
are added. Use `--start` and `--end` to sample a shorter part of a run.

Outputs include the GIF, clean and raw PNG frames, an editable `animation.hipnc`,
local initial meshes, final `.bgeo.sc` caches, source provenance, per-frame
geometry checks, and rendering progress. All final surfaces are checked for
closedness, a single connected component, conserved source volume, and
conserved source volume centroid. The display reconstruction changes local
shape and surface area; it does not add information to the source simulation.

The scene uses `/obj/membrane/DISPLAY_SOURCE` to select the verified caches
(input 0) or the live reconstruction (input 1). Input 1 allows editing VDB
rounding controls interactively; cached renders reproduce the checked output.

Add `--resume` after an interruption. Inputs, settings, and scripts must match
the saved plan; existing mesh, geometry, and PNG files are reused only if their
recorded hashes still match.


## Approved green animation backup (2026-10-28)

The exact user-approved camera and lighting B are frozen in [backups/2026-10-28_approved_green_membrane](backups/2026-10-28_approved_green_membrane/README.md). This snapshot includes all 10,001 source timesteps for future denser animations. Use its HIP file as the authoritative layout; the historical rounded renderer refits the camera.


## Permeable shear preview: approved green view (2026-10-29)

The completed 40-frame, 5-second preview of geo0–geo300 is in [outputs/permeable_Ca10p1_v0p9_steps0-300_40f](outputs/permeable_Ca10p1_v0p9_steps0-300_40f/README.md). Its scene retains the approved camera orientation, roll, lens, membrane transform, material, and lighting B. The camera distance and centering are fitted once to the entire new sequence. `approved_view_sequence.py` records and runs this specific 40-frame preview, including geometry validation and exact five-second MP4/GIF timing.
