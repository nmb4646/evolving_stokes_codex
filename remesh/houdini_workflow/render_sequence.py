"""Houdini worker for workflow.py; all frame transforms and color limits are fixed."""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import sys
import time
from pathlib import Path

import hou
import numpy as np
from PIL import Image
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def select_files(directory, count, start=None, end=None, explicit=None):
    files = {}
    for path in directory.glob("geo*.mat"):
        match = re.fullmatch(r"geo(\d+)\.mat", path.name)
        if match:
            index = int(match.group(1))
            if index in files:
                raise ValueError("Duplicate source index: %d" % index)
            files[index] = path
    if not files:
        raise ValueError("No geo<number>.mat files found in %s" % directory)
    available = sorted(files)
    if explicit:
        chosen = [int(v) for v in explicit.replace(",", " ").split()]
        if chosen != sorted(set(chosen)):
            raise ValueError("--indices must be unique and increasing.")
        missing = set(chosen) - files.keys()
        if missing:
            raise ValueError("Missing source steps: %s" % sorted(missing))
    else:
        eligible = [v for v in available if (start is None or v >= start) and (end is None or v <= end)]
        if count < 2 or count > len(eligible):
            raise ValueError("--frames must be between 2 and the number of available states.")
        chosen = [eligible[int(i)] for i in np.rint(np.linspace(0, len(eligible) - 1, count))]
    if len(chosen) < 2:
        raise ValueError("An animation needs at least two source states.")
    return files, chosen


def timing(files, explicit_dt=None):
    # In this solver p.total_time accumulates tic/toc: it is wall-clock time.
    # Automatically infer physical time only for a complete, constant-dt run.
    steps = sorted(files)
    dts, gammas = [], []
    for index in steps:
        p = loadmat(str(files[index]), variable_names=["p"], squeeze_me=True, struct_as_record=False).get("p")
        dts.append(float(getattr(p, "dt", float("nan"))))
        gammas.append(float(getattr(p, "gamy", float("nan"))))
    dt = None
    method = "unknown: labels use source step"
    if explicit_dt is not None:
        dt = float(explicit_dt)
        if not math.isfinite(dt) or dt <= 0:
            raise ValueError("--dt must be finite and positive.")
        method = "source step times user-specified constant dt"
    elif steps == list(range(steps[-1] + 1)) and np.isfinite(dts).all() and dts[0] > 0 and np.allclose(dts, dts[0], rtol=1e-10, atol=0):
        dt = dts[0]
        method = "source step times p.dt; constant dt verified across every saved state from geo0"
    gamma = gammas[0] if np.isfinite(gammas).all() and np.allclose(gammas, gammas[0]) else None
    return dt, gamma, method


def camera_for_sequence(points_list, direction, padding):
    back = np.asarray(direction, dtype=float)
    if not np.isfinite(back).all() or np.linalg.norm(back) < 1e-8:
        raise ValueError("Camera direction must be a finite nonzero vector.")
    back /= np.linalg.norm(back)
    up_hint = np.array([0., 0., 1.])
    if abs(np.dot(back, up_hint)) > 0.98:
        up_hint = np.array([0., 1., 0.])
    right = np.cross(up_hint, back)
    right /= np.linalg.norm(right)
    up = np.cross(back, right)
    axes = np.column_stack((right, up, back))
    projected = np.concatenate(points_list) @ axes
    lower, upper = projected.min(axis=0), projected.max(axis=0)
    center = ((lower + upper) / 2) @ axes.T
    width = float(max(upper[0] - lower[0], upper[1] - lower[1]) * padding)
    eye = center + back * max(width * 3, 8)
    matrix = np.eye(4)
    matrix[:3, :3] = axes.T
    matrix[3, :3] = eye
    return {"world_matrix": matrix.ravel().tolist(), "target": center.tolist(),
            "direction_from_target": back.tolist(), "ortho_width": width,
            "padding": padding, "fixed_for_entire_sequence": True}


def render(config):
    directory = Path(config["data"]).resolve()
    output = Path(config["output"]).resolve()
    template = Path(config["template"]).resolve()
    files, chosen = select_files(directory, config["frames"], config.get("start"), config.get("end"), config.get("indices"))
    dt, gamma, time_method = timing(files, config.get("dt"))
    points_list, frames = [], []
    for ordinal, step in enumerate(chosen, 1):
        path = files[step]
        data = loadmat(str(path), variable_names=["P", "M"])
        points = np.asarray(data["P"], dtype=float)
        faces = np.asarray(data["M"])
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("Invalid P in %s" % path)
        if faces.ndim != 2 or faces.shape[1] != 3 or not np.isfinite(faces).all() or not np.equal(faces, np.floor(faces)).all():
            raise ValueError("Invalid M in %s" % path)
        if faces.min() < 1 or faces.max() > len(points):
            raise ValueError("Invalid triangle indices in %s" % path)
        points_list.append(points)
        physical_time = step * dt if dt is not None else None
        frames.append({"frame": ordinal, "source_step": step, "source_file": str(path),
                       "source_sha256": digest(path), "vertices": len(points), "triangles": len(faces),
                       "physical_time": physical_time,
                       "shear_strain": gamma * physical_time if gamma is not None and physical_time is not None else None,
                       "raw_file": "raw/frame_%04d.png" % ordinal})
    camera = camera_for_sequence(points_list, config["view"], config["padding"])
    plan = {"schema": 1, "data_directory": str(directory), "template_sha256": digest(template),
            "worker_sha256": digest(__file__), "importer_sha256": digest(ROOT / "sop_loader.py"),
            "resolution": config["resolution"], "samples": config["samples"], "threads": config["threads"],
            "camera": camera, "physical_time_method": time_method, "constant_dt": dt,
            "shear_rate": gamma, "frames": frames,
            "look": "approved legacy green material, lighting, fixed color ramp, and both remesh stages",
            "simulation_orientation_preserved": True}
    plan_path = output / "render_plan.json"
    if plan_path.exists():
        if not config["resume"]:
            raise FileExistsError("Output already has a render plan. Use --resume for the same settings or a new output folder.")
        if json.loads(plan_path.read_text()) != plan:
            raise ValueError("Resume settings/source/template differ from the existing render plan. Choose a new output folder.")
    elif output.exists() and any(output.iterdir()):
        raise FileExistsError("Output folder is not empty and has no matching render plan.")
    output.mkdir(parents=True, exist_ok=True)
    raw_dir = output / "raw"
    raw_dir.mkdir(exist_ok=True)
    save_json(plan_path, plan)
    manifest_path = output / "manifest.json"
    manifest = {"plan": plan, "status": "rendering", "completed_frames": []}
    old_records = {}
    if config["resume"] and manifest_path.exists():
        old_records = {r["frame"]: r for r in json.loads(manifest_path.read_text()).get("completed_frames", [])}

    hou.hipFile.load(str(template), suppress_save_prompt=True, ignore_load_warnings=True)
    obj = hou.node("/obj/membrane")
    obj.setWorldTransform(hou.Matrix4(1))
    controls = obj.node("CONTROLS")
    controls.setParms({"data_directory": str(directory), "source_indices": " ".join(map(str, chosen))})
    obj.node("READ_MAT").parm("python").set((ROOT / "sop_loader.py").read_text())
    cam = hou.node("/obj/cam1")
    cam.setWorldTransform(hou.Matrix4(camera["world_matrix"]))
    cam.setParms({"projection": 1, "orthowidth": camera["ortho_width"],
                  "resx": config["resolution"], "resy": config["resolution"],
                  "near": 0.001, "far": 10000., "vm_background": "$HIP/white.png"})
    shutil.copyfile(ROOT / "assets/white.png", output / "white.png")
    rop = hou.node("/out/render")
    picture = "$HIP/raw/frame_$F4.png"
    rop.setParms({"vm_picture": picture, "vm_samplesx": config["samples"], "vm_samplesy": config["samples"],
                  "vm_usemaxthreads": 0, "vm_threadcount": config["threads"],
                  "soho_pipecmd": "mantra -j %d" % config["threads"], "vm_gamma": 1.0})
    hou.setFps(config["fps"])
    hou.playbar.setFrameRange(1, len(chosen))
    hou.playbar.setPlaybackRange(1, len(chosen))
    hou.setFrame(1)
    hou.setUpdateMode(hou.updateMode.AutoUpdate)
    scene = output / "animation.hipnc"
    hou.hipFile.save(str(scene))
    save_json(manifest_path, manifest)
    print("Plan: %d states, fixed camera, %dx%d, %d CPU threads" %
          (len(frames), config["resolution"], config["resolution"], config["threads"]), flush=True)
    for record, expected_points in zip(frames, points_list):
        frame = record["frame"]
        destination = output / record["raw_file"]
        old = old_records.get(frame)
        if config["resume"] and old and destination.exists() and old.get("png_sha256") == digest(destination):
            with Image.open(destination) as image:
                valid = image.size == (config["resolution"], config["resolution"])
                image.verify()
            if valid:
                manifest["completed_frames"].append(old)
                save_json(manifest_path, manifest)
                print("[%d/%d] reused geo%d" % (frame, len(frames), record["source_step"]), flush=True)
                continue
        started = time.monotonic()
        hou.setFrame(frame)
        raw = obj.node("READ_MAT").geometry()
        imported = np.asarray(raw.pointFloatAttribValues("P")).reshape(-1, 3)
        source_faces = loadmat(record["source_file"], variable_names=["M"])["M"].astype(int) - 1
        imported_faces = np.asarray([[p.number() for p in prim.points()] for prim in raw.prims()])
        if not np.allclose(imported, expected_points, rtol=1e-6, atol=1e-7) or not np.array_equal(imported_faces, source_faces):
            raise RuntimeError("Geometry import validation failed for %s" % record["source_file"])
        final_node = obj.node("Material_Membrane")
        final_geometry = final_node.geometry()
        errors = {n.path(): n.errors() for n in obj.children() if n.errors()}
        if errors:
            raise RuntimeError(str(errors))
        temp_png = destination.with_name(destination.stem + ".tmp.png")
        rop.parm("vm_picture").set(str(temp_png))
        print("[%d/%d] rendering geo%d" % (frame, len(frames), record["source_step"]), flush=True)
        rop.render(frame_range=(frame, frame), verbose=False, output_progress=False)
        with Image.open(temp_png) as image:
            if image.size != (config["resolution"], config["resolution"]):
                raise RuntimeError("Unexpected render dimensions.")
            image.verify()
        os.replace(temp_png, destination)
        done = {"frame": frame, "source_step": record["source_step"], "png_sha256": digest(destination),
                "seconds": round(time.monotonic() - started, 3),
                "maximum_import_position_error": float(np.max(np.abs(imported - expected_points))),
                "triangle_connectivity_exact": True, "display_vertices": len(final_geometry.points()),
                "display_triangles": len(final_geometry.prims())}
        manifest["completed_frames"].append(done)
        save_json(manifest_path, manifest)
        print("[%d/%d] complete in %.1fs" % (frame, len(frames), done["seconds"]), flush=True)
    rop.parm("vm_picture").set(picture)
    hou.setFrame(1)
    hou.hipFile.save(str(scene))
    manifest["status"] = "rendered"
    manifest["scene"] = str(scene)
    save_json(manifest_path, manifest)
    print("RENDER_COMPLETE", str(output), flush=True)


if __name__ == "__main__":
    render(json.loads(sys.argv[1]))
