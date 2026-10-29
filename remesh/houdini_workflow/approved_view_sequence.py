"""Reconstruct a new membrane sequence while retaining an approved Houdini view.

Run with the project PyVista Python; child Houdini phases run in hython.
"""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)

def prepare_one(task):
    import numpy as np
    from scipy.io import loadmat
    from green_surface_samples import load_mesh,refine,save_geometry_ply,volume_centroid
    frame,source,destination=task
    path=Path(source);destination=Path(destination)
    mesh=load_mesh(path);volume,center=volume_centroid(mesh)
    refined,report=refine(mesh,4,20,.15)
    save_geometry_ply(refined,destination)
    p=loadmat(path,variable_names=['p'],simplify_cells=True)['p']
    return {'frame':frame,'source_step':int(path.stem[3:]),'source_file':str(path),
            'source_sha256':digest(path),'source_volume':volume,'source_centroid':center.tolist(),
            'mesh_file':'meshes/'+destination.name,'mesh_sha256':digest(destination),'initial_refinement':report,
            'parameters':{name:float(p[name]) for name in ['Sd','Da','gamy','v','dt']}}

def load_live_scene(plan,out):
    import hou
    hou.hipFile.load(plan['template'],suppress_save_prompt=True,ignore_load_warnings=True)
    hou.setUpdateMode(hou.updateMode.AutoUpdate)
    # Read transforms before manual mode: an uncooked camera may report identity.
    camera=hou.node('/obj/cam1').worldTransform().asTuple()
    obj=hou.node('/obj/membrane');object_matrix=obj.worldTransform().asTuple()
    hou.setUpdateMode(hou.updateMode.Manual);hou.setFrame(1)
    obj.node('READ_REFINED_SURFACE').parm('file').set(str(out/'meshes/frame_$F4.ply'))
    records=json.loads((out/'prepared.json').read_text())['frames']
    targets=[{'volume':r['source_volume'],'centroid':r['source_centroid']} for r in records]
    restore=obj.node('RESTORE_VOLUME_AND_CENTER')
    lines=restore.parm('python').evalAsString().splitlines()
    assert sum(line.startswith('targets = ') for line in lines)==1
    restore.parm('python').set('\n'.join('targets = '+repr(targets) if line.startswith('targets = ') else line for line in lines))
    obj.node('SURFACE_DISTANCE').setParms({'voxelsize':.02,'exteriorbandvoxels':12,'interiorbandvoxels':12})
    obj.node('ROUND_PINCHES').setParms({'operation':'gaussian','radius':4,'iterations':4})
    obj.node('SMOOTH_NORMALS').setParms({'type':0,'cuspangle':180,'normalize':1})
    hou.setUpdateMode(hou.updateMode.AutoUpdate)
    return obj,camera,object_matrix,records

def geometry_worker(plan,out,start,end):
    import hou
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    obj,_,_,records=load_live_scene(plan,out)
    checkpoint=out/('geometry_%03d_%03d.json'%(start,end));checks=[]
    previous={r['frame']:r for r in json.loads(checkpoint.read_text())} if checkpoint.exists() else {}
    for r in records[start-1:end]:
        frame=r['frame'];cache=out/'cache'/('frame_%04d.bgeo.sc'%frame)
        if frame in previous and cache.exists() and digest(cache)==previous[frame]['cache_sha256']:
            checks.append(previous[frame]);save(checkpoint,checks);continue
        begin=time.monotonic();hou.setFrame(frame)
        geo=obj.node('Material_Membrane').geometry()
        if obj.node('Material_Membrane').errors():raise RuntimeError(obj.node('Material_Membrane').errors())
        points=np.array(geo.pointFloatAttribValues('P')).reshape(-1,3)
        faces=np.array([[v.point().number() for v in p.vertices()] for p in geo.prims()])
        assert faces.shape[1]==3 and np.isfinite(points).all()
        a,b,c=(points[faces[:,j]] for j in range(3))
        signed=np.einsum('ij,ij->i',a,np.cross(b,c))/6
        volume=abs(signed.sum());center=(signed[:,None]*(a+b+c)/4).sum(0)/signed.sum()
        verr=abs(volume-r['source_volume'])/r['source_volume'];cerr=np.linalg.norm(center-r['source_centroid'])
        edges=np.sort(np.vstack((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]])),axis=1)
        edges,counts=np.unique(edges,axis=0,return_counts=True)
        graph=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(points),len(points)))
        components=connected_components(graph,directed=False,return_labels=False)
        assert verr<1e-6 and cerr<1e-6 and np.all(counts==2) and components==1
        geo.saveToFile(str(cache))
        checks.append({'frame':frame,'source_step':r['source_step'],'vertices':len(points),'triangles':len(faces),
                       'relative_volume_error':float(verr),'centroid_error':float(cerr),
                       'closed_surface':True,'connected_components':int(components),'cache_sha256':digest(cache)})
        save(checkpoint,checks)
        print('GEOMETRY %d/%d geo%d %.1fs'%(frame,plan['frames'],r['source_step'],time.monotonic()-begin),flush=True)

def finalize_scene(plan,out):
    import hou
    import numpy as np
    obj,camera,object_matrix,records=load_live_scene(plan,out)
    original=np.array(camera).reshape(4,4);basis=original[:3,:3]
    transform=np.array(object_matrix).reshape(4,4)
    lo=np.full(3,np.inf);hi=-lo;all_bounds=[]
    for r in records:
        geo=hou.Geometry();geo.loadFromFile(str(out/'cache'/('frame_%04d.bgeo.sc'%r['frame'])))
        points=np.array(geo.pointFloatAttribValues('P')).reshape(-1,3)
        world=points@transform[:3,:3]+transform[3,:3]
        projected=world@basis.T;lo=np.minimum(lo,projected.min(0));hi=np.maximum(hi,projected.max(0))
        all_bounds.append((projected.min(0),projected.max(0)))
    target=((lo+hi)/2)@basis
    width=float(max((hi-lo)[:2])*1.35)
    cam=hou.node('/obj/cam1')
    focal=cam.parm('focal').eval();aperture=cam.parm('aperture').eval()
    assert cam.parm('projection').eval()==0,'Expected the approved perspective camera'
    distance=width*focal/aperture
    matrix=original.copy();matrix[3,:3]=target+basis[2]*distance
    cam.setWorldTransform(hou.Matrix4(matrix.ravel().tolist()))
    cam.setParms({'resx':plan['resolution'],'resy':plan['resolution'],'vm_background':'$HIP/white.png'})
    center_projected=target@basis.T
    occupancy=0.
    for low,high in all_bounds:
        depth=distance-(high[2]-center_projected[2])
        half=depth*aperture/focal/2
        occupancy=max(occupancy,float(np.max(np.abs(np.r_[low[:2]-center_projected[:2],high[:2]-center_projected[:2]]))/half))
    assert occupancy<.80,'Framing margin insufficient'
    obj.node('RENDER_CACHE').parm('file').set(str(out/'cache/frame_$F4.bgeo.sc'))
    obj.node('DISPLAY_SOURCE').parm('input').set(0)
    obj.node('DISPLAY_SOURCE').setDisplayFlag(True);obj.node('DISPLAY_SOURCE').setRenderFlag(True)
    hou.setFps(plan['fps']);hou.playbar.setFrameRange(1,plan['frames']);hou.playbar.setPlaybackRange(1,plan['frames'])
    rop=hou.node('/out/render')
    rop.setParms({'camera':'/obj/cam1','vm_picture':'$HIP/raw/frame_$F4.png','vm_samplesx':plan['samples'],
                  'vm_samplesy':plan['samples'],'vm_usemaxthreads':0,'vm_threadcount':plan['threads'],
                  'soho_pipecmd':'mantra -j %d'%plan['threads'],'vm_gamma':1.0})
    hou.setFrame(1);hou.hipFile.save(str(out/'animation.hipnc'))
    obj.node('RENDER_CACHE').parm('file').set('$HIP/cache/frame_$F4.bgeo.sc')
    obj.node('READ_REFINED_SURFACE').parm('file').set('$HIP/meshes/frame_$F4.ply')
    hou.hipFile.save(str(out/'animation.hipnc'))
    assert np.allclose(np.array(cam.worldTransform().asTuple()).reshape(4,4)[:3,:3],basis,atol=1e-8)
    save(out/'camera.json',{'world_matrix':matrix.tolist(),'template_world_matrix':original.tolist(),
                           'object_world_matrix':transform.tolist(),'view_angle_preserved':True,
                           'focal_length_preserved':focal,'target':target.tolist(),'width':width,
                           'maximum_frame_halfwidth_fraction':occupancy,
                           'lighting_and_material_preserved':True})
    print('SCENE_READY width %.4f occupancy %.3f'%(width,occupancy),flush=True)

def render_worker(plan,out,start,end):
    import hou
    import numpy as np
    from PIL import Image
    from animation import _auto_corner_mask
    scene=out/'animation.hipnc'
    hou.hipFile.load(str(scene),suppress_save_prompt=True,ignore_load_warnings=True)
    hou.setUpdateMode(hou.updateMode.AutoUpdate)
    expected=np.array(json.loads((out/'camera.json').read_text())['world_matrix']).ravel()
    assert np.allclose(hou.node('/obj/cam1').worldTransform().asTuple(),expected,atol=1e-8)
    checkpoint=out/('render_%03d_%03d.json'%(start,end));records=[]
    previous={r['frame']:r for r in json.loads(checkpoint.read_text())} if checkpoint.exists() else {}
    tracked={path:hou.node(path).worldTransform().asTuple() for path in
             ['/obj/cam1','/obj/membrane','/obj/KEY_FRONT','/obj/SOFT_FRONT_FILL']}
    rop=hou.node('/out/render')
    for frame in range(start,end+1):
        target=out/'raw'/('frame_%04d.png'%frame)
        if frame in previous and target.exists() and digest(target)==previous[frame]['sha256']:
            records.append(previous[frame]);save(checkpoint,records);continue
        begin=time.monotonic();hou.setFrame(frame)
        for path,transform in tracked.items():assert np.allclose(hou.node(path).worldTransform().asTuple(),transform,atol=1e-8)
        temporary=target.with_suffix('.tmp.png');rop.parm('vm_picture').set(str(temporary))
        rop.render(frame_range=(frame,frame),verbose=False,output_progress=False)
        with Image.open(temporary) as im:
            im.load();assert im.size==(plan['resolution'],)*2
            rectangle,mask=_auto_corner_mask(im.convert('RGBA'))
        assert mask['method']=='alpha_corner_components' and not mask.get('protected_foreground_overlap'),mask
        temporary.replace(target)
        records.append({'frame':frame,'sha256':digest(target),'seconds':time.monotonic()-begin,'corner_mask':mask})
        save(checkpoint,records)
        print('RENDER %d/%d %.1fs'%(frame,plan['frames'],records[-1]['seconds']),flush=True)

def assemble(plan,out):
    import numpy as np
    from PIL import Image,ImageDraw,ImageFont,GifImagePlugin
    from animation import assemble_frames,_palette
    frames=json.loads((out/'prepared.json').read_text())['frames']
    result=assemble_frames(out/'raw',out,fps=plan['fps'],labels=['t = %.2f'%r['physical_time'] for r in frames],
                           make_mp4=True,ffmpeg='/opt/hfs20.0.724/bin/hffmpeg',name=plan['name'],force=True)
    # Distribute centisecond GIF delays to preserve the requested total duration.
    paths=sorted((out/'clean_frames').glob('frame_*.png'))
    ticks=np.rint(np.linspace(0,plan['duration']*100,plan['frames']+1)).astype(int)
    delays=(np.diff(ticks)*10).tolist()
    palette=_palette(paths);colors=palette.getpalette()
    background=min(range(len(colors)//3),key=lambda i:sum((255-v)**2 for v in colors[3*i:3*i+3]))
    destination=Path(result['gif']);temporary=destination.with_suffix('.tmp.gif')
    with temporary.open('wb') as stream:
        for i,path in enumerate(paths):
            with Image.open(path) as source:frame=source.convert('RGB').quantize(palette=palette,dither=Image.Dither.NONE)
            if i==0:
                header,_=GifImagePlugin.getheader(frame,info={'duration':delays[0],'background':background,'optimize':False})
                for block in header:stream.write(block)
                stream.write(b'\x21\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00')
            for block in GifImagePlugin.getdata(frame,duration=delays[i],disposal=2,optimize=False,include_color_table=False):stream.write(block)
        stream.write(b';')
    with Image.open(temporary) as gif:
        assert gif.n_frames==plan['frames'];total=0
        for i in range(gif.n_frames):gif.seek(i);total+=gif.info['duration']
        assert total==round(plan['duration']*1000)
    temporary.replace(destination)
    result.update(gif_frame_durations_ms=delays,gif_frame_duration_ms=None,gif_duration_seconds=total/1000,gif_fps=plan['fps'])
    decoded=subprocess.run([result['ffmpeg'],'-hide_banner','-i',result['mp4'],'-f','rawvideo','-pix_fmt','rgb24','-'],
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
    size=plan['resolution'];unit=size*size*3
    assert len(decoded.stdout)==plan['frames']*unit
    duration_match=re.search(r'Duration: (\d+):(\d+):([\d.]+)',decoded.stderr.decode())
    assert duration_match,decoded.stderr.decode()
    hours,minutes,seconds=map(float,duration_match.groups())
    video_duration=hours*3600+minutes*60+seconds
    assert abs(video_duration-plan['duration'])<=.02
    indices=np.rint(np.linspace(0,plan['frames']-1,6)).astype(int)
    canvas=Image.new('RGB',(1152,832),'white');draw=ImageDraw.Draw(canvas)
    font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',16)
    for i,frame in enumerate(indices):
        image=Image.frombytes('RGB',(size,size),decoded.stdout[frame*unit:(frame+1)*unit]);image.thumbnail((384,384))
        x=i%3*384;y=i//3*416;canvas.paste(image,(x,y+32))
        draw.text((x+12,y+7),'geo%d'%frames[frame]['source_step'],font=font,fill='#264239')
    canvas.save(out/'storyboard.png')
    checks=[]
    for path in sorted(out.glob('geometry_*.json')):checks.extend(json.loads(path.read_text()))
    assert len(checks)==plan['frames']
    verification={'decoded_frames':plan['frames'],'unique_frames':len({hashlib.sha256(decoded.stdout[i*unit:(i+1)*unit]).hexdigest() for i in range(plan['frames'])}),
                  'mp4_duration_seconds':video_duration,'gif_duration_seconds':total/1000,
                  'max_relative_volume_error':max(r['relative_volume_error'] for r in checks),
                  'all_surfaces_closed_and_connected':True,'mp4_sha256':digest(result['mp4']),'gif_sha256':digest(result['gif'])}
    save(out/'manifest.json',{'status':'complete','plan':plan,'camera':json.loads((out/'camera.json').read_text()),
                             'scene':str(out/'animation.hipnc'),'source_frames':frames,'animation':result,'verification':verification})
    (out/'README.md').write_text('# Permeable shear flow: approved green Houdini view\n\n'
        f"{plan['frames']} actual saved timesteps from geo{plan['source_steps'][0]} through geo{plan['source_steps'][-1]}, "
        f"{size} × {size}, {plan['fps']:g} fps, {plan['duration']:g} seconds.\n\n"
        'The camera orientation, roll, perspective lens, membrane object transform, material, lighting B, and reconstruction settings '
        'come from the frozen approved scene. Camera centering and distance are fitted once to all selected surfaces. '
        'Each display surface has its original enclosed volume and centroid restored.\n\n'
        f"Files: animation.hipnc, {plan['name']}.mp4, {plan['name']}.gif, manifest.json. "
        'Source mapping, input hashes, settings, and verification are in the JSON records.\n')
    print('VIDEO_COMPLETE',result['mp4'],json.dumps(verification),flush=True)

def launch(plan,out,phase,ranges):
    env=dict(os.environ,HOUDINI_MAXTHREADS=str(plan['threads']),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS=str(plan['threads']))
    processes=[subprocess.Popen(['/opt/hfs20.0.724/bin/hython',str(Path(__file__).resolve()),'--output',str(out),
                                 '--phase',phase,'--range',str(start),str(end)],env=env) for start,end in ranges]
    codes=[p.wait() for p in processes]
    if any(codes):raise RuntimeError('%s workers failed: %s'%(phase,codes))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--phase',choices=['geometry','scene','render']);parser.add_argument('--range',nargs=2,type=int)
    parser.add_argument('--start',type=int,default=0);parser.add_argument('--end',type=int,default=300)
    selection=parser.add_mutually_exclusive_group()
    selection.add_argument('--frames',type=int);selection.add_argument('--stride',type=int)
    parser.add_argument('--duration',type=float,default=5.);parser.add_argument('--name')
    parser.add_argument('--resolution',type=int,default=512);parser.add_argument('--samples',type=int,default=3)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--data',type=Path);parser.add_argument('--template',type=Path)
    args=parser.parse_args();out=args.output.resolve()
    if args.phase:
        plan=json.loads((out/'plan.json').read_text())
        if args.range is None:args.range=[1,plan['frames']]
        if args.phase=='geometry':geometry_worker(plan,out,*args.range)
        elif args.phase=='scene':finalize_scene(plan,out)
        else:render_worker(plan,out,*args.range)
        return
    import numpy as np
    from scipy.io import loadmat
    if args.data is None or args.template is None:raise ValueError('Data and template are required')
    if not 0<=args.start<=args.end or args.duration<=0:raise ValueError('Invalid interval or duration')
    if args.stride is not None:
        if args.stride<1:raise ValueError('Stride must be positive')
        indices=list(range(args.start,args.end+1,args.stride))
    else:
        count=args.frames if args.frames is not None else 40
        if not 2<=count<=args.end-args.start+1:raise ValueError('Invalid frame count')
        indices=np.rint(np.linspace(args.start,args.end,count)).astype(int).tolist()
    if len(indices)<2 or not .01<=len(indices)/args.duration<=100:raise ValueError('Invalid playback rate')
    if abs(args.duration*100-round(args.duration*100))>1e-7:raise ValueError('GIF duration must be a whole number of centiseconds')
    if not 64<=args.resolution<=960 or not 1<=args.samples<=16 or not 1<=args.threads<=64:raise ValueError('Invalid rendering settings')
    name=args.name or 'permeable_shear_%d-%d'%(args.start,args.end)
    if Path(name).name!=name or name in ('','.','..'):raise ValueError('Invalid output name')
    plan={'data':str(args.data.resolve()),'template':str(args.template.resolve()),'template_sha256':digest(args.template),
          'source_steps':indices,'frames':len(indices),'duration':args.duration,'fps':len(indices)/args.duration,
          'resolution':args.resolution,'samples':args.samples,'threads':args.threads,'name':name,
          'worker_script_sha256':digest(__file__),'surface_settings':{'taubin_iterations':20,'pass_band':.15,'loop_subdivisions':4,
          'vdb_voxel_size':.02,'gaussian_radius':4,'gaussian_iterations':4},'source_hashes':{str(i):digest(args.data/('geo%d.mat'%i)) for i in indices}}
    out.mkdir(parents=True,exist_ok=True)
    for name in ['meshes','cache','raw']:(out/name).mkdir(exist_ok=True)
    if (out/'plan.json').exists():assert json.loads((out/'plan.json').read_text())==plan
    else:save(out/'plan.json',plan)
    shutil.copy2(args.template.parent/'white.png',out/'white.png')
    # Verify actual timestep values over the entire requested interval.
    dt=[float(loadmat(args.data/('geo%d.mat'%i),variable_names=['p'],simplify_cells=True)['p']['dt']) for i in range(args.start,args.end+1)]
    assert np.allclose(dt,dt[0],rtol=1e-12,atol=0)
    if not (out/'prepared.json').exists():
        tasks=[(n,str(args.data/('geo%d.mat'%step)),str(out/'meshes'/('frame_%04d.ply'%n))) for n,step in enumerate(indices,1)]
        records=[]
        with concurrent.futures.ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
            for r in pool.map(prepare_one,tasks):
                r['physical_time']=r['source_step']*dt[0];records.append(r)
                print('PREPARED %d/%d geo%d'%(r['frame'],len(indices),r['source_step']),flush=True)
        save(out/'prepared.json',{'constant_dt':dt[0],'physical_time_method':'step times p.dt, verified across steps %d–%d'%(args.start,args.end),'frames':records})
    split=(len(indices)+1)//2
    ranges=[(1,split),(split+1,len(indices))]
    launch(plan,out,'geometry',ranges)
    launch(plan,out,'scene',[(1,len(indices))])
    launch(plan,out,'render',ranges)
    assemble(plan,out)

if __name__=='__main__':main()
