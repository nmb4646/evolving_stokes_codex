#!/usr/bin/env python3
"""Refine coarse saved membranes and render green late-time samples."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.io import loadmat
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

# Reuse the same conservative check for physical-time labels.
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from curvature_flow_render import constant_time_step

GREEN=(0.080,0.380,0.330)

def unit(v):
    return v/np.linalg.norm(v)

def volume_centroid(mesh):
    faces=mesh.faces.reshape((-1,4))[:,1:]
    points=np.asarray(mesh.points,dtype=np.float64)
    a,b,c=(points[faces[:,j]] for j in range(3))
    volumes=np.einsum("ij,ij->i",a,np.cross(b,c))/6
    volume=volumes.sum()
    if abs(volume)<1e-12:
        raise ValueError("Mesh has no enclosed volume.")
    centroid=(volumes[:,None]*(a+b+c)/4).sum(axis=0)/volume
    return float(abs(volume)),centroid

def load_mesh(path):
    data=loadmat(path,variable_names=["P","M"])
    p=np.asarray(data["P"],dtype=float)
    f=np.asarray(data["M"],dtype=np.int64)-1
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():
        raise ValueError(f"Invalid points in {path}")
    if f.ndim!=2 or f.shape[1]!=3 or f.min()<0 or f.max()>=len(p):
        raise ValueError(f"Invalid faces in {path}")
    mesh=pv.PolyData(p,np.c_[np.full(len(f),3),f])
    if mesh.n_open_edges:
        raise ValueError("Expected a closed membrane.")
    return mesh

def refine(mesh,subdivisions,iterations,pass_band):
    target_volume,target_center=volume_centroid(mesh)
    display=mesh.copy()
    if iterations:
        display=display.smooth_taubin(n_iter=iterations,pass_band=pass_band,
            normalize_coordinates=True,boundary_smoothing=False,
            feature_smoothing=False)
    display=display.subdivide(subdivisions,subfilter="loop")
    smoothed_volume,smoothed_center=volume_centroid(display)
    factor=(target_volume/smoothed_volume)**(1/3)
    display.points=(np.asarray(display.points,dtype=np.float64)-smoothed_center)*factor+target_center
    display=display.compute_normals(point_normals=True,cell_normals=False,
        split_vertices=False,consistent_normals=True,auto_orient_normals=True)
    final_volume,final_center=volume_centroid(display)
    _,closest=display.find_closest_cell(mesh.points,return_closest_point=True)
    displacement=np.linalg.norm(closest-mesh.points,axis=1)
    report={"source_vertices":mesh.n_points,"source_triangles":mesh.n_cells,
        "display_vertices":display.n_points,"display_triangles":display.n_cells,
        "source_volume":target_volume,"display_volume":final_volume,
        "relative_volume_error":abs(final_volume-target_volume)/target_volume,
        "centroid_error":float(np.linalg.norm(final_center-target_center)),
        "source_area":float(mesh.area),"display_area":float(display.area),
        "volume_restoration_scale":factor,
        "source_vertex_to_display_rms":float(np.sqrt(np.mean(displacement**2))),
        "source_vertex_to_display_max":float(displacement.max()),
        "loop_subdivisions":subdivisions,"taubin_iterations":iterations,
        "taubin_pass_band":pass_band}
    return display,report

def save_geometry_ply(mesh, path):
    """Export geometry without VTK normals; Houdini must compute its own normals.

    On this import path the supplied VTK normals oppose Houdini's face normals,
    producing inverted specular shading. Keep normals on the preview mesh, but
    omit them from the exchange file so the receiving renderer can derive them.
    """
    exchange = mesh.copy()
    exchange.clear_data()
    exchange.save(path, recompute_normals=False)

def camera_for(meshes,direction):
    back=unit(np.asarray(direction,dtype=float))
    hint=np.array([0.,0.,1.])
    if abs(back@hint)>.98:hint=np.array([0.,1.,0.])
    right=unit(np.cross(hint,back));up=np.cross(back,right)
    axes=np.column_stack((right,up,back))
    p=np.vstack([m.points for m in meshes])@axes
    lo,hi=p.min(axis=0),p.max(axis=0)
    center=((lo+hi)/2)@axes.T
    scale=max(hi[0]-lo[0],hi[1]-lo[1])*.59
    return back,right,up,center,scale

def render(mesh,path,camera,resolution,raw=False):
    back,right,up,center,scale=camera
    size=scale*2
    plotter=pv.Plotter(off_screen=True,window_size=(resolution,resolution),
                       lighting="none")
    plotter.set_background("white")
    plotter.camera_position=[center+back*size*3,center,up]
    plotter.enable_parallel_projection();plotter.camera.parallel_scale=scale
    plotter.add_mesh(mesh,color=GREEN,smooth_shading=not raw,
        show_edges=raw,edge_color=(.09,.23,.20),line_width=.65,
        pbr=False,ambient=.40,diffuse=.80,specular=.14,specular_power=28)
    plotter.add_light(pv.Light(position=center+(-right+up+back)*size,
        focal_point=center,intensity=1.0,color=(1.,.97,.93)))
    plotter.add_light(pv.Light(position=center+(right*.6+back*1.1)*size,
        focal_point=center,intensity=.50,color=(.91,.97,1.)))
    plotter.enable_anti_aliasing("ssaa")
    plotter.enable_ssao(radius=size*.04,bias=size*.002,kernel_size=64)
    plotter.screenshot(str(path));plotter.close()

def font(size):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",size)

def assemble(records,out,resolution,comparison=False):
    n=len(records);header=100;label=55;footer=65;gap=30
    rows=2 if comparison else 1
    canvas=Image.new("RGB",(n*resolution+(n+1)*gap,
        header+rows*(resolution+label)+footer),"white")
    draw=ImageDraw.Draw(canvas)
    title="Coarse mesh / refined display surface" if comparison else "Green membrane - refined late-time samples"
    draw.text((gap,24),title,font=font(34),fill="#263a34")
    for j,r in enumerate(records):
        x=gap+j*(resolution+gap)
        time=f"t = {r['physical_time']:.3f}" if r["physical_time"] is not None else f"geo{r['frame']}"
        for row in range(rows):
            y=header+row*(resolution+label)
            raw=comparison and row==0
            name=f"geo{r['frame']}_{'raw' if raw else 'green_smooth'}.png"
            canvas.paste(Image.open(out/name).convert("RGB"),(x,y+label))
            desc=(f"  |  {r['source_vertices']} source vertices" if raw else "")
            if comparison and not raw:desc=f"  |  {r['display_vertices']:,} display vertices"
            draw.text((x+18,y+5),f"geo{r['frame']}  |  {time}"+desc,
                      font=font(25),fill="#35463e")
    note="Surface refinement for display; enclosed volume and center are restored."
    draw.text((gap,canvas.height-47),note,font=font(22),fill="#536158")
    filename="source_vs_refined.png" if comparison else "green_late_samples.png"
    canvas.save(out/filename)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("data",type=Path)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--frames",nargs="+",type=int,default=[5000,9000])
    p.add_argument("--resolution",type=int,default=1200)
    p.add_argument("--subdivisions",type=int,default=4)
    p.add_argument("--smoothing-iterations",type=int,default=20)
    p.add_argument("--pass-band",type=float,default=.15)
    p.add_argument("--view",nargs=3,type=float,default=[-1.2,-1.4,1.0])
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    sources=[];displays=[];records=[]
    print("Preparing smooth display meshes...",flush=True)
    for frame in args.frames:
        source=args.data/f"geo{frame}.mat"
        mesh=load_mesh(source)
        display,report=refine(mesh,args.subdivisions,args.smoothing_iterations,args.pass_band)
        save_geometry_ply(mesh,args.output/f"geo{frame}_source.ply")
        save_geometry_ply(display,args.output/f"geo{frame}_refined.ply")
        sources.append(mesh);displays.append(display)
        records.append({"frame":frame,"source":str(source.resolve()),
            "source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),**report})
    print("Checking saved timesteps...",flush=True)
    dt=constant_time_step(args.data)
    camera=camera_for(sources+displays,args.view)
    for source,display,record in zip(sources,displays,records):
        frame=record["frame"]
        record["physical_time"]=frame*dt if dt is not None else None
        render(source,args.output/f"geo{frame}_raw.png",camera,args.resolution,raw=True)
        render(display,args.output/f"geo{frame}_green_smooth.png",camera,args.resolution)
        print(f"geo{frame}: {record['source_vertices']} -> {record['display_vertices']} vertices; volume error {record['relative_volume_error']:.2e}",flush=True)
    assemble(records,args.output,args.resolution)
    assemble(records,args.output,args.resolution,comparison=True)
    metadata={"data":str(args.data.resolve()),"renderer":"PyVista/VTK",
        "look":"Uniform teal-green approximation to the approved Houdini material",
        "color_rgb":GREEN,"camera_direction":unit(np.asarray(args.view)).tolist(),
        "geometry":"Taubin smoothing, Loop subdivision, volume/centroid restoration",
        "purpose":"Display reconstruction of coarse geometry; no simulation data overwritten",
        "constant_dt":dt,"frames":records}
    (args.output/"manifest.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(f"FINAL {args.output/'green_late_samples.png'}",flush=True)

if __name__=="__main__":
    main()

