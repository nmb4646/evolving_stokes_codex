"""Houdini worker: cache rounded surfaces, fit one camera, and render every state."""
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import sys
import time
import hou
import numpy as np
from PIL import Image
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from render_sequence import camera_for_sequence, digest, save_json

ROOT=Path(__file__).resolve().parent
OUT=Path(sys.argv[1]).resolve()
plan=json.loads((OUT/'plan.json').read_text())
prep=json.loads((OUT/'prepared.json').read_text())
(OUT/'cache').mkdir(exist_ok=True);(OUT/'raw').mkdir(exist_ok=True)
shutil.copyfile(ROOT/'assets/white.png',OUT/'white.png')
hou.hipFile.load(plan['template'],suppress_save_prompt=True,ignore_load_warnings=True)
hou.setUpdateMode(hou.updateMode.Manual)
obj=hou.node('/obj/membrane')
reader=obj.node('READ_REFINED_SURFACE')
reader.parm('file').set(str(OUT/'meshes/frame_$F4.ply'))
restore=obj.node('RESTORE_VOLUME_AND_CENTER')
targets=[{'volume':r['source_volume'],'centroid':r['source_centroid']} for r in prep['frames']]
lines=restore.parm('python').evalAsString().splitlines()
if sum(line.startswith('targets = ') for line in lines)!=1:
    raise ValueError('Template needs the rounded-surface volume restoration SOP')
restore.parm('python').set('\n'.join('targets = '+repr(targets) if line.startswith('targets = ') else line for line in lines))
obj.node('SURFACE_DISTANCE').setParms({'voxelsize':.02,'exteriorbandvoxels':12,'interiorbandvoxels':12})
obj.node('ROUND_PINCHES').setParms({'operation':'gaussian','radius':4,'iterations':4})
obj.node('SMOOTH_NORMALS').setParms({'type':0,'cuspangle':180,'normalize':1})
# Evaluate the loaded camera before reading its transform; manual mode can
# return an uncooked identity matrix for a newly loaded camera.
hou.setUpdateMode(hou.updateMode.AutoUpdate)
cam=hou.node('/obj/cam1')
axes=np.asarray(cam.worldTransform().asTuple()).reshape(4,4)[:3,:3].T
back=axes[:,2]
final_node=obj.node('Material_Membrane')
hou.setUpdateMode(hou.updateMode.AutoUpdate)
check_path=OUT/'geometry_checks.json'
checks=json.loads(check_path.read_text()) if check_path.exists() else []
previous={r['frame']:r for r in checks}
checks=[];bounds=[]
for r in prep['frames']:
    started=time.monotonic();frame=r['frame']
    cache=OUT/'cache'/('frame_%04d.bgeo.sc'%frame)
    old=previous.get(frame)
    if old and cache.exists() and digest(cache)==old['cache_sha256']:
        geo=hou.Geometry();geo.loadFromFile(str(cache));record=old
    else:
        hou.setFrame(frame)
        geo=final_node.geometry()
        errors={n.path():n.errors() for n in obj.children() if n.errors() and n.name() not in ('READ_MAT','remesh1','remesh2')}
        if errors:raise RuntimeError(str(errors))
        points=np.asarray(geo.pointFloatAttribValues('P')).reshape(-1,3)
        faces=np.asarray([[v.point().number() for v in p.vertices()] for p in geo.prims()],dtype=np.int64)
        if faces.ndim!=2 or faces.shape[1]!=3 or not np.isfinite(points).all():
            raise ValueError('Invalid final display surface at frame %d'%frame)
        a,b,c=(points[faces[:,j]] for j in range(3))
        signed=np.einsum('ij,ij->i',a,np.cross(b,c))/6
        volume=float(abs(signed.sum()))
        center=(signed[:,None]*(a+b+c)/4).sum(axis=0)/signed.sum()
        volume_error=abs(volume-r['source_volume'])/r['source_volume']
        center_error=float(np.linalg.norm(center-np.asarray(r['source_centroid'])))
        edges=np.sort(np.vstack((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]])),axis=1)
        edges,counts=np.unique(edges,axis=0,return_counts=True)
        closed=bool(np.all(counts==2))
        graph=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(points),len(points)))
        components=int(connected_components(graph,directed=False,return_labels=False))
        if volume_error>1e-6 or center_error>1e-6 or not closed or components!=1:
            raise RuntimeError('Surface verification failed at frame %d: V %.3g, center %.3g, closed %s, components %d'%(frame,volume_error,center_error,closed,components))
        geo.saveToFile(str(cache))
        record={'frame':frame,'source_step':r['source_step'],'vertices':len(points),'triangles':len(faces),
                'volume':volume,'relative_volume_error':volume_error,'centroid_error':center_error,
                'closed_surface':closed,'connected_components':components,'cache_file':str(cache.relative_to(OUT)),
                'cache_sha256':digest(cache),'seconds':round(time.monotonic()-started,3)}
    points=np.asarray(geo.pointFloatAttribValues('P')).reshape(-1,3)
    projected=points@axes
    lo,hi=projected.min(axis=0),projected.max(axis=0)
    bounds.append(np.asarray(list(itertools.product(*zip(lo,hi))))@axes.T)
    checks.append(record);save_json(check_path,checks)
    if frame==1 or frame%5==0 or frame==len(prep['frames']):
        print('Rounded and verified %d/%d: geo%d'%(frame,len(prep['frames']),r['source_step']),flush=True)
camera=camera_for_sequence(bounds,back,plan['padding'])
if not np.allclose(np.asarray(camera['world_matrix']).reshape(4,4)[:3,:3], axes.T, atol=1e-8):
    raise RuntimeError('Fitted camera must retain the approved template viewing orientation')
print('Approved camera orientation preserved', flush=True)
cam.setWorldTransform(hou.Matrix4(camera['world_matrix']))
cam.setParms({'projection':1,'orthowidth':camera['ortho_width'],'resx':plan['resolution'],'resy':plan['resolution'],
              'near':.001,'far':10000.,'vm_background':'$HIP/white.png'})
cache_reader=obj.createNode('file','RENDER_CACHE')
cache_reader.parm('file').set(str(OUT/'cache/frame_$F4.bgeo.sc'))
cache_reader.setComment('Cached verified display surfaces. Switch DISPLAY_SOURCE to input 1 to edit the live reconstruction.')
cache_reader.setGenericFlag(hou.nodeFlag.DisplayComment,True)
switch=obj.createNode('switch','DISPLAY_SOURCE')
switch.setInput(0,cache_reader);switch.setInput(1,final_node);switch.parm('input').set(0)
switch.setDisplayFlag(True);switch.setRenderFlag(True);obj.layoutChildren()
hou.setFps(plan['fps']);hou.playbar.setFrameRange(1,len(prep['frames']));hou.playbar.setPlaybackRange(1,len(prep['frames']))
rop=hou.node('/out/render')
rop.setParms({'vm_picture':'$HIP/raw/frame_$F4.png','vm_samplesx':plan['samples'],'vm_samplesy':plan['samples'],
              'vm_usemaxthreads':0,'vm_threadcount':plan['threads'],'soho_pipecmd':'mantra -j %d'%plan['threads'],'vm_gamma':1.0})
scene=OUT/'animation.hipnc';hou.setFrame(1);hou.hipFile.save(str(scene))
reader.parm('file').set('$HIP/meshes/frame_$F4.ply');cache_reader.parm('file').set('$HIP/cache/frame_$F4.bgeo.sc')
hou.hipFile.save(str(scene))
manifest_path=OUT/'manifest.json'
old_manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
old_frames={r['frame']:r for r in old_manifest.get('completed_frames',[])}
manifest={'status':'rendering','renderer':'Houdini %s / Mantra'%hou.applicationVersionString(),
          'lighting':'B: front key at 15 degrees elevation; fixed frontal fill and sky from supplied template',
          'scene':str(scene),'camera':camera,'plan':plan,'physical_time_method':prep['physical_time_method'],
          'constant_dt':prep['constant_dt'],'source_frames':prep['frames'],'completed_frames':[]}
save_json(manifest_path,manifest)
for r,check in zip(prep['frames'],checks):
    started=time.monotonic();frame=r['frame'];destination=OUT/'raw'/('frame_%04d.png'%frame)
    old=old_frames.get(frame)
    if old and destination.exists() and digest(destination)==old['png_sha256']:
        with Image.open(destination) as im:
            if im.size!=(plan['resolution'],plan['resolution']):raise RuntimeError('Invalid cached render size')
            im.verify()
        manifest['completed_frames'].append(old);save_json(manifest_path,manifest)
        print('[%d/%d] reused geo%d'%(frame,len(checks),r['source_step']),flush=True)
        continue
    hou.setFrame(frame)
    temp=destination.with_suffix('.tmp.png');rop.parm('vm_picture').set(str(temp))
    print('[%d/%d] rendering geo%d'%(frame,len(checks),r['source_step']),flush=True)
    rop.render(frame_range=(frame,frame),verbose=False,output_progress=False)
    with Image.open(temp) as im:
        if im.size!=(plan['resolution'],plan['resolution']):raise RuntimeError('Unexpected render dimensions')
        im.verify()
    temp.replace(destination)
    record={'frame':frame,'source_step':r['source_step'],'physical_time':r['physical_time'],
            'png_sha256':digest(destination),'geometry_cache_sha256':check['cache_sha256'],
            'seconds':round(time.monotonic()-started,3)}
    manifest['completed_frames'].append(record);save_json(manifest_path,manifest)
    print('[%d/%d] complete in %.1fs'%(frame,len(checks),record['seconds']),flush=True)
rop.parm('vm_picture').set('$HIP/raw/frame_$F4.png');hou.setFrame(1);hou.hipFile.save(str(scene))
manifest['status']='rendered';save_json(manifest_path,manifest)
print('RENDER_COMPLETE',OUT,flush=True)
