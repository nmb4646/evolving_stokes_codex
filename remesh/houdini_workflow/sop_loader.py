"""Embedded in the scene: read the source MAT selected by the timeline."""
from pathlib import Path
import hou
import numpy as np
from scipy.io import loadmat

node = hou.pwd()
geo = node.geometry()
geo.clear()
controls = node.parent().node("CONTROLS")
directory = controls.evalParm("data_directory")
if directory:
    indices = [int(v) for v in controls.evalParm("source_indices").replace(",", " ").split()]
    output_frame = int(round(hou.frame()))
    if not 1 <= output_frame <= len(indices):
        raise hou.NodeError("Timeline frame must be within the configured source sequence.")
    step = indices[output_frame - 1]
    path = Path(directory) / ("geo%d.mat" % step)
    data = loadmat(str(path), variable_names=["P", "M"])
    positions = np.asarray(data["P"], dtype=np.float64)
    faces = np.asarray(data["M"])
    if positions.ndim != 2 or positions.shape[1] != 3 or not np.isfinite(positions).all():
        raise hou.NodeError("P must be a finite N x 3 array.")
    if faces.ndim != 2 or faces.shape[1] != 3 or not np.isfinite(faces).all():
        raise hou.NodeError("M must be a finite F x 3 array.")
    if not np.equal(faces, np.floor(faces)).all():
        raise hou.NodeError("M contains noninteger indices.")
    faces = faces.astype(np.int64) - 1
    if faces.min() < 0 or faces.max() >= len(positions):
        raise hou.NodeError("M must use MATLAB one-based indices into P.")
    points = geo.createPoints(positions.tolist())
    for face in faces:
        polygon = geo.createPolygon()
        for index in face:
            polygon.addVertex(points[int(index)])
    geo.addAttrib(hou.attribType.Global, "source_mat", "")
    geo.setGlobalAttribValue("source_mat", str(path))
    geo.addAttrib(hou.attribType.Global, "source_step", step)
