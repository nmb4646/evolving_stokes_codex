#!/usr/bin/env python3
"""Shaded membrane renders with short, surface-conforming velocity arrows."""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.sparse import csr_matrix
from scipy.ndimage import gaussian_filter
import pyvista as pv
import vtk
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFont

PALETTE = LinearSegmentedColormap.from_list("membrane", [
    (0.07,0.02,0.18), (0.26,0.06,0.40), (0.54,0.10,0.53),
    (0.79,0.19,0.43), (0.95,0.40,0.28), (0.99,0.69,0.31),
    (0.99,0.94,0.66)])

def unit(x):
    return x / np.maximum(np.linalg.norm(x, axis=-1, keepdims=True), 1e-30)

def native_fields(points, faces):
    """Same signed edge-dihedral curvature as MATLAB Geometry.m."""
    cross = np.cross(points[faces[:,1]]-points[faces[:,0]],
                     points[faces[:,2]]-points[faces[:,1]])
    area = np.linalg.norm(cross, axis=1)/2
    fn = unit(cross)
    vn = np.zeros_like(points)
    va = np.zeros(len(points))
    for j in range(3):
        np.add.at(vn, faces[:,j], cross/2)
        np.add.at(va, faces[:,j], area/3)
    vn = unit(vn)
    src = faces.T.reshape(-1)
    dst = faces[:,[1,2,0]].T.reshape(-1)
    fi = np.tile(np.arange(len(faces)), 3)
    keys = src * len(points) + dst
    order = np.argsort(keys)
    reverse = dst * len(points) + src
    flip = order[np.searchsorted(keys[order], reverse)]
    if not np.array_equal(keys[flip], reverse):
        raise ValueError("Expected a closed, consistently oriented triangle mesh.")
    edge = points[dst]-points[src]
    angle = np.arctan2(np.sum(np.cross(fn[fi],fn[fi[flip]])*unit(edge),axis=1),
                       np.sum(fn[fi]*fn[fi[flip]],axis=1))
    h = np.bincount(src, weights=.25*np.linalg.norm(edge,axis=1)*angle,
                    minlength=len(points))/va
    return h, vn

def neighbor_average(field, faces, passes):
    src = faces.T.reshape(-1)
    dst = faces[:,[1,2,0]].T.reshape(-1)
    adjacency = csr_matrix((np.ones(2*len(src)),
        (np.r_[src,dst], np.r_[dst,src])), shape=(len(field),len(field)))
    adjacency.data[:] = 1
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    denominator = (2+degree)[:,None] if field.ndim == 2 else 2+degree
    for _ in range(passes):
        field = (2*field+adjacency @ field)/denominator
    return field

def constant_time_step(folder):
    saved={}
    for path in folder.glob("geo*.mat"):
        match=re.fullmatch(r"geo(\d+)\.mat",path.name)
        if match:
            saved[int(match.group(1))]=path
    ids=sorted(saved)
    if not ids or ids!=list(range(ids[-1]+1)):
        return None
    dts=[]
    for frame in ids:
        d=loadmat(saved[frame],variable_names=["p"],squeeze_me=True,struct_as_record=False)
        dts.append(float(getattr(d.get("p"),"dt",float("nan"))))
    dts=np.asarray(dts)
    if not np.isfinite(dts).all() or np.any(dts<=0):
        return None
    return float(dts[0]) if np.allclose(dts,dts[0],rtol=1e-10,atol=0) else None

def load_state(folder, frame, passes, subdivisions):
    source = folder/f"geo{frame}.mat"
    d = loadmat(source, squeeze_me=True, struct_as_record=False)
    p = np.asarray(d["P"],dtype=float)
    f = np.asarray(d["M"],dtype=np.int64)-1
    if not np.isfinite(p).all() or f.min()<0 or f.max()>=len(p):
        raise ValueError(f"Invalid mesh in {source}")
    h, n = native_fields(p,f)
    velocity = np.asarray(d["velocity"]).reshape((-1,3),order="F")
    if velocity.shape != p.shape or not np.isfinite(velocity).all():
        raise ValueError(f"Invalid velocity in {source}")
    u = velocity-n*np.sum(velocity*n,axis=1)[:,None]
    h = neighbor_average(h,f,passes)
    u = neighbor_average(u,f,passes)
    u -= n*np.sum(u*n,axis=1)[:,None]
    base = pv.PolyData(p,np.c_[np.full(len(f),3),f])
    base["H"] = h
    base["U"] = u
    refined = base.subdivide(subdivisions,subfilter="loop")
    refined = refined.compute_normals(point_normals=True,cell_normals=False,
        split_vertices=False,consistent_normals=True,auto_orient_normals=True)
    n = np.asarray(refined.point_normals)
    refined["U"] -= n*np.sum(refined["U"]*n,axis=1)[:,None]
    dt = float(getattr(d.get("p"),"dt",float("nan")))
    return {"frame":frame,"source":str(source),"mesh":refined,"dt":dt,
            "source_vertices":len(p),"source_faces":len(f),
            "curvature_limits":[float(h.min()),float(h.max())],
            "speeds":np.linalg.norm(u,axis=1)}

class Surface:
    def __init__(self,mesh):
        self.points = np.asarray(mesh.points)
        self.faces = mesh.faces.reshape((-1,4))[:,1:]
        self.normals = np.asarray(mesh.point_normals)
        self.velocity = np.asarray(mesh["U"])
        self.locator = vtk.vtkStaticCellLocator()
        self.locator.SetDataSet(mesh)
        self.locator.BuildLocator()
        self.cell_id = vtk.reference(0)
        self.sub_id = vtk.reference(0)
        self.distance2 = vtk.reference(0.)
        self.closest = [0.,0.,0.]

    def sample(self,point):
        self.locator.FindClosestPoint(point,self.closest,self.cell_id,
                                      self.sub_id,self.distance2)
        p = np.array(self.closest)
        tri = self.faces[int(self.cell_id)]
        a,b,c = self.points[tri]
        x,y,q = b-a,c-a,p-a
        xx,xy,yy = x@x,x@y,y@y
        determinant = xx*yy-xy*xy
        bx = ((q@x)*yy-(q@y)*xy)/determinant
        by = ((q@y)*xx-(q@x)*xy)/determinant
        weights = np.array([1-bx-by,bx,by])
        normal = unit(weights @ self.normals[tri])
        velocity = weights @ self.velocity[tri]
        velocity -= (velocity@normal)*normal
        return p,normal,velocity

    def trace(self,seed,length,step,direction):
        p,n,u = self.sample(seed)
        points,normals = [p],[n]
        for _ in range(max(1,int(np.ceil(length/step)))):
            if np.linalg.norm(u) < 1e-10:
                break
            h = min(step,length)
            _,_,mid_u = self.sample(p+direction*.5*h*unit(u))
            q,qn,qu = self.sample(p+direction*h*unit(mid_u))
            if np.linalg.norm(q-p)<1e-8:
                break
            points.append(q); normals.append(qn)
            p,n,u = q,qn,qu
            length -= h
            if length <= 1e-10:
                break
        return np.asarray(points),np.asarray(normals)

def spaced_seeds(mesh,back,right,up,count):
    # Sample in screen coordinates to avoid clusters caused by foreshortening.
    normals = np.asarray(mesh.point_normals)
    ids = np.flatnonzero(normals @ back > .08)
    screen = np.c_[mesh.points[ids] @ right,mesh.points[ids] @ up]
    chosen = []
    index = np.argmin(np.linalg.norm(screen-screen.mean(axis=0),axis=1))
    distance = np.full(len(ids),np.inf)
    for _ in range(min(count,len(ids))):
        chosen.append(ids[index])
        distance = np.minimum(distance,np.sum((screen-screen[index])**2,axis=1))
        index = int(np.argmax(distance))
    return chosen

def curved_arrows(mesh,back,right,up,size,speed_reference,count):
    surface = Surface(mesh)
    tube_parts,head_points,head_faces = [],[],[]
    stats = []
    for seed_id in spaced_seeds(mesh,back,right,up,count):
        seed = mesh.points[seed_id]
        _,_,u = surface.sample(seed)
        relative_speed = np.linalg.norm(u)/max(speed_reference,1e-14)
        if relative_speed < .10:
            continue
        length = size*.16*min(relative_speed,1.25)
        # Centered glyph: integrate backward to its tail and forward to its tip.
        bp,bn = surface.trace(seed,length*.48,size*.0035,-1)
        fp,fn = surface.trace(seed,length*.52,size*.0035,1)
        points = np.vstack((bp[:0:-1],fp))
        normals = np.vstack((bn[:0:-1],fn))
        if len(points)<6:
            continue
        distance = np.r_[0,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))]
        distance /= distance[-1]
        # Hide glyphs crossing the far side naturally with the depth buffer.
        offset = size*.0015
        radius = size*.00095
        body_end = max(2,int(np.searchsorted(distance,.76)))
        centerline = points[:body_end+1]+normals[:body_end+1]*(offset+radius)
        line = pv.lines_from_points(centerline)
        tube_parts.append(line.tube(radius=radius,n_sides=8,capping=True))
        # A tapered, curved ribbon head is projected onto the same surface.
        head_ids = np.arange(body_end,len(points))
        for k in head_ids:
            tangent = points[min(k+1,len(points)-1)]-points[max(k-1,0)]
            side = unit(np.cross(normals[k],unit(tangent)))
            fraction = (distance[k]-distance[body_end])/(1-distance[body_end])
            width = min(size*.006,length*.085)*(1-fraction)
            for sign in (-1,1):
                vertex,normal,_ = surface.sample(points[k]+sign*width*side)
                head_points.append(vertex+normal*(offset+radius))
        base = len(head_points)-2*len(head_ids)
        for j in range(len(head_ids)-1):
            a=base+2*j
            head_faces.extend([[3,a,a+1,a+2],[3,a+1,a+3,a+2]])
        # Check that the path bends and remains on the refined membrane.
        chord=np.linalg.norm(points[-1]-points[0])
        arc=np.linalg.norm(np.diff(points,axis=0),axis=1).sum()
        stats.append({"arc_length":float(arc),"chord_length":float(chord),
                      "seed_speed":float(np.linalg.norm(u))})
    if not tube_parts:
        return None,None,stats
    tubes=pv.merge(tube_parts,merge_points=False)
    heads=pv.PolyData(np.asarray(head_points),np.asarray(head_faces).ravel())
    heads=heads.compute_normals(point_normals=True,split_vertices=False)
    return tubes,heads,stats

def surface_finish(mesh):
    # Material-only microtexture; it does not displace geometry or alter H/U.
    textured=mesh.texture_map_to_sphere(prevent_seam=True)
    tangents=vtk.vtkPolyDataTangents()
    tangents.SetInputData(textured)
    tangents.Update()
    textured=pv.wrap(tangents.GetOutput())
    grain=gaussian_filter(np.random.default_rng(7).standard_normal((512,1024)),.85,
                          mode="wrap")
    gy,gx=np.gradient(grain)
    slope=np.stack((-gx,-gy,np.ones_like(gx)),axis=-1)
    slope[...,:2] *= .08/max(float(np.std(gx)),1e-8)
    normal_map=((unit(slope)*.5+.5)*255).astype(np.uint8)
    texture=pv.Texture(normal_map)
    texture.InterpolateOn()
    texture.MipmapOn()
    texture.RepeatOn()
    return textured,texture

def render_state(state,destination,limits,camera,resolution,arrow_count):
    mesh=state["mesh"]
    back,right,up,center,scale=camera
    size=2*scale
    tubes,heads,stats=curved_arrows(mesh,back,right,up,size,
                                   state["speed_reference"],arrow_count)
    plotter=pv.Plotter(off_screen=True,window_size=(resolution,resolution),
                       lighting="none")
    plotter.set_background("white")
    plotter.camera_position=[center+back*size*3,center,up]
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale=scale
    shaded_mesh,normal_texture = surface_finish(mesh)
    # Apply the colormap directly. Diffuse lighting preserves saturation;
    # restrained specular highlights provide the glossy surface cues.
    actor=plotter.add_mesh(shaded_mesh,scalars="H",cmap=PALETTE,clim=limits,
        show_scalar_bar=False,smooth_shading=True,pbr=False,
        ambient=.28,diffuse=.78,specular=.30,specular_power=40)
    actor.GetProperty().SetNormalTexture(normal_texture)
    actor.GetProperty().SetNormalScale(.10)
    plotter.add_light(pv.Light(position=center+(-right+up+back)*size,
                              focal_point=center,intensity=.95))
    plotter.add_light(pv.Light(position=center+(right*.6+back*.9-up*.1)*size,
                              focal_point=center,intensity=.35))
    if tubes is not None:
        for item in (tubes,heads):
            plotter.add_mesh(item,color=(.96,.97,.99),smooth_shading=True,
                pbr=False,ambient=.65,diffuse=.55,specular=.12,specular_power=24)
    plotter.enable_anti_aliasing("ssaa")
    plotter.enable_ssao(radius=size*.025,bias=size*.0015,kernel_size=64)
    plotter.screenshot(str(destination))
    plotter.close()
    state["arrows"]=stats
    print(f"geo{state['frame']}: {len(stats)} curved arrows; wrote {destination}",flush=True)

def font(size):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",size)

def compose(states,output,resolution,limits,image_name):
    # Tight scientific layout with one common color scale.
    margin=55; gap=35; header=96; footer=96; bar_space=150
    width=2*resolution+2*margin+gap+bar_space
    height=2*resolution+header+footer+gap
    canvas=Image.new("RGB",(width,height),"white")
    draw=ImageDraw.Draw(canvas)
    draw.text((margin,26),"Mean curvature and tangential velocity",fill="#24232b",
              font=font(38))
    for j,state in enumerate(states):
        x=margin+(j%2)*(resolution+gap)
        y=header+(j//2)*(resolution+gap)
        canvas.paste(Image.open(output/f"geo{state['frame']}_curved.png").convert("RGB"),(x,y))
        t=state["physical_time"]
        label=(f"({chr(97+j)})   t = {t:.4g}" if t is not None else
               f"({chr(97+j)})   geo{state['frame']}")
        draw.text((x+25,y+15),label,fill="#33313b",font=font(30))
    bx=width-118; by=header+70; bh=2*resolution+gap-140; bw=25
    gradient=np.linspace(1,0,bh)
    color=(PALETTE(gradient)[:,:3]*255).astype(np.uint8)
    bar=Image.fromarray(np.repeat(color[:,None,:],bw,axis=1))
    canvas.paste(bar,(bx,by))
    draw.text((bx-7,by-48),"H",fill="#33313b",font=font(30))
    ticks=np.linspace(limits[0],limits[1],5)
    for value in ticks:
        y=by+int((limits[1]-value)/(limits[1]-limits[0])*(bh-1))
        draw.line((bx+bw,y,bx+bw+6,y),fill="#777777",width=1)
        draw.text((bx+bw+11,y-13),f"{value:.2g}",fill="#44414a",font=font(23))
    draw.text((margin,height-65),
        "White arrows follow the surface flow. Color: averaged mean curvature.",
        fill="#55515c",font=font(26))
    destination=output/image_name
    canvas.save(destination)
    return destination

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_folder",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--image-name",default="shear_curvature_flow_panels.png")
    parser.add_argument("--frames",type=int,nargs="+",default=[35,85,135,185])
    parser.add_argument("--smoothing",type=int,default=4)
    parser.add_argument("--subdivisions",type=int,default=3)
    parser.add_argument("--resolution",type=int,default=1100)
    parser.add_argument("--arrows",type=int,default=85)
    parser.add_argument("--view",type=float,nargs=3,default=[-.325,-.894,.309])
    args=parser.parse_args()
    if len(args.frames) not in (1,4):
        parser.error("Choose one frame for a close-up or four for a panel figure.")
    if args.smoothing < 0 or args.subdivisions < 0 or args.resolution < 256:
        parser.error("Smoothing/subdivisions must be nonnegative and resolution at least 256.")
    if args.arrows < 1 or np.linalg.norm(args.view) < 1e-12:
        parser.error("Use a positive arrow count and nonzero view direction.")
    if Path(args.image_name).name != args.image_name:
        parser.error("--image-name must be a filename without a directory.")
    args.output.mkdir(parents=True,exist_ok=True)
    states=[load_state(args.run_folder,n,args.smoothing,args.subdivisions) for n in args.frames]
    limits=[min(s["curvature_limits"][0] for s in states),
            max(s["curvature_limits"][1] for s in states)]
    if limits[1]<=limits[0]: limits[1]=limits[0]+1e-8
    speed_reference=float(np.percentile(np.concatenate([s["speeds"] for s in states]),95))
    back=unit(np.asarray(args.view))
    right=unit(np.cross([0.,0.,1.],back)); up=np.cross(back,right)
    points=np.vstack([s["mesh"].points for s in states])
    axes=np.column_stack((right,up,back))
    proj=points @ axes
    lower,upper=proj.min(axis=0),proj.max(axis=0)
    center=((lower+upper)/2) @ axes.T
    scale=max(upper[0]-lower[0],upper[1]-lower[1])*.57
    camera=back,right,up,center,scale
    dt=constant_time_step(args.run_folder)
    for state in states:
        state["physical_time"]=state["frame"]*dt if dt is not None else None
        state["speed_reference"]=speed_reference
        render_state(state,args.output/f"geo{state['frame']}_curved.png",
                     limits,camera,args.resolution,args.arrows)
    if len(states)==4:
        result=compose(states,args.output,args.resolution,limits,args.image_name)
    else:
        result=args.output/f"geo{states[0]['frame']}_curved.png"
    report={"curvature_limits":limits,"speed_reference":speed_reference,
            "smoothing_passes":args.smoothing,"loop_subdivisions":args.subdivisions,
            "camera_direction":back.tolist(),"curvature_definition":"Geometry.m signed edge-dihedral mean curvature",
            "display_geometry":"Loop subdivision; fields interpolated from source-vertex data",
            "material_texture":"Fixed procedural normal map; strength 0.10; no geometry displacement",
            "surface_shading":"Phong: ambient 0.28, diffuse 0.78, specular 0.30, exponent 40; direct colormap without extra gamma transform",
            "arrow_method":"RK2 tangent integration with projection to refined surface at every step",
            "frames":[{k:v for k,v in s.items() if k not in ("mesh","speeds")} for s in states]}
    (args.output/"render_manifest.json").write_text(json.dumps(report,indent=2)+"\n")
    print(f"FINAL {result}",flush=True)

if __name__=="__main__":
    main()

