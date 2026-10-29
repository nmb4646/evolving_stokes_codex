#!/usr/bin/env python3
"""Render a rounded membrane sequence with lighting from an approved Houdini scene.
Run with the project's .venv-curvature-flow Python (PyVista, SciPy, Pillow).
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import numpy as np
from green_surface_samples import load_mesh, refine, save_geometry_ply, volume_centroid, constant_time_step
from animation import assemble_frames
from workflow import find_hython

ROOT = Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save_json(path, data):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(data, indent=2) + '\n')
    temp.replace(path)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--template', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--frames', type=int, default=60)
    parser.add_argument('--start', type=int)
    parser.add_argument('--end', type=int)
    parser.add_argument('--resolution', type=int, default=384)
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--fps', type=float, default=5)
    parser.add_argument('--padding', type=float, default=1.45)
    parser.add_argument('--name', default='shear_B_60')
    parser.add_argument('--hython')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    data, template, out = (p.expanduser().resolve() for p in (args.data,args.template,args.output))
    if not 64 <= args.resolution <= 960 or not 1 <= args.samples <= 16 or not 1 <= args.threads <= 64:
        raise ValueError('Invalid render resolution, sample count, or thread count')
    if args.padding < 1.3 or not 0 < args.fps <= 100:
        raise ValueError('Use padding >= 1.3 and 0 < fps <= 100')
    available = {}
    for path in data.glob('geo*.mat'):
        match = re.fullmatch(r'geo(\d+)\.mat', path.name)
        if match:
            step = int(match.group(1))
            if step in available:
                raise ValueError('Duplicate geometry index')
            available[step] = path
    eligible = [i for i in sorted(available) if (args.start is None or i >= args.start) and (args.end is None or i <= args.end)]
    if not 2 <= args.frames <= len(eligible):
        raise ValueError('Frame count must be between 2 and the number of available states')
    indices = [eligible[i] for i in np.rint(np.linspace(0,len(eligible)-1,args.frames)).astype(int)]
    config = {'data':str(data), 'template':str(template), 'template_sha256':digest(template),
              'resolution':args.resolution, 'samples':args.samples, 'threads':args.threads,
              'fps':args.fps, 'padding':args.padding, 'name':args.name,
              'worker_sha256':digest(ROOT/'rounded_sequence.py'), 'preparer_sha256':digest(__file__),
              'refiner_sha256':digest(ROOT/'green_surface_samples.py'),
              'surface':{'taubin_iterations':20,'pass_band':.15,'loop_subdivisions':4,
                         'vdb_voxel_size':.02,'gaussian_radius_voxels':4,'gaussian_iterations':4},
              'sources':[{'step':i,'file':str(available[i]),'sha256':digest(available[i])} for i in indices]}
    plan_path = out/'plan.json'
    if plan_path.exists():
        if not args.resume or json.loads(plan_path.read_text()) != config:
            raise ValueError('Output already contains a plan; resume requires identical settings and inputs')
    elif out.exists() and any(out.iterdir()):
        raise ValueError('Choose an empty output folder')
    out.mkdir(parents=True,exist_ok=True)
    (out/'meshes').mkdir(exist_ok=True)
    save_json(plan_path,config)
    print('Verifying physical timestep across saved states...',flush=True)
    dt = constant_time_step(data)
    prep_path=out/'prepared.json'
    old = json.loads(prep_path.read_text()) if args.resume and prep_path.exists() else {}
    previous = {r['frame']:r for r in old.get('frames',[])}
    prep = {'constant_dt':dt,'physical_time_method':'step times p.dt, verified constant across complete saved run' if dt else 'unavailable; labels use source step', 'frames':[]}
    for frame, source in enumerate(config['sources'],1):
        mesh_path=out/'meshes'/('frame_%04d.ply'%frame)
        cached=previous.get(frame)
        if cached and mesh_path.exists() and digest(mesh_path)==cached['mesh_sha256']:
            prep['frames'].append(cached)
            continue
        mesh=load_mesh(Path(source['file']))
        volume,center=volume_centroid(mesh)
        display,report=refine(mesh,4,20,.15)
        save_geometry_ply(display,mesh_path)
        record={'frame':frame,'source_step':source['step'],'source_file':source['file'],
                'source_sha256':source['sha256'],'physical_time':source['step']*dt if dt else None,
                'source_volume':volume,'source_centroid':center.tolist(),
                'mesh_file':str(mesh_path.relative_to(out)),'mesh_sha256':digest(mesh_path),
                'initial_refinement':report}
        prep['frames'].append(record)
        save_json(prep_path,prep)
        if frame==1 or frame%5==0 or frame==args.frames:
            print('Prepared %d/%d: geo%d (%d source vertices)'%(frame,args.frames,source['step'],mesh.n_points),flush=True)
    save_json(prep_path,prep)
    env=os.environ.copy();env['HOUDINI_MAXTHREADS']=str(args.threads)
    subprocess.run([find_hython(args.hython),'-j',str(args.threads),str(ROOT/'rounded_sequence.py'),str(out)],env=env,check=True)
    labels=[('t = %.3f   |   step %d'%(r['physical_time'],r['source_step'])) if dt else ('step %d'%r['source_step']) for r in prep['frames']]
    result=assemble_frames(out/'raw',out,fps=args.fps,labels=labels,name=args.name,
                           cover_corner='auto',force=args.resume)
    manifest_path=out/'manifest.json'
    manifest=json.loads(manifest_path.read_text());manifest.update(status='complete',animation=result)
    save_json(manifest_path,manifest)
    print('GIF_COMPLETE',result['gif'], 'frames',result['frame_count'],'seconds',result['gif_duration_seconds'],flush=True)

if __name__=='__main__':
    main()
