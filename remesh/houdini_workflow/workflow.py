#!/usr/bin/env python3
"""Render MATLAB membrane states with Houdini, or assemble existing PNG frames."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from animation import assemble_frames

ROOT = Path(__file__).resolve().parent


def find_hython(explicit=None):
    if explicit:
        path = Path(explicit).expanduser()
        if not path.is_file():
            raise FileNotFoundError("Hython executable does not exist: %s" % path)
        return str(path)
    found = shutil.which("hython")
    if found:
        return found
    hfs = os.environ.get("HFS")
    if hfs and (Path(hfs) / "bin/hython").is_file():
        return str(Path(hfs) / "bin/hython")
    candidates = sorted(Path("/opt").glob("hfs*/bin/hython"), reverse=True)
    if candidates:
        return str(candidates[0])
    raise FileNotFoundError("Houdini hython was not found. Supply --hython /path/to/hython.")


def parser():
    cli = argparse.ArgumentParser(description=__doc__)
    commands = cli.add_subparsers(dest="command", required=True)
    render = commands.add_parser("render", help="MAT files -> green Houdini scene + raw PNGs + cleaned GIF")
    render.add_argument("--data", type=Path, required=True, help="Folder containing geo<number>.mat files")
    render.add_argument("--output", type=Path, required=True)
    render.add_argument("--frames", type=int, default=20)
    render.add_argument("--start", type=int)
    render.add_argument("--end", type=int)
    render.add_argument("--indices", help="Explicit increasing indices, comma or space separated")
    render.add_argument("--resolution", type=int, default=512, help="Square image size")
    render.add_argument("--samples", type=int, default=3, help="Mantra pixel samples per dimension")
    render.add_argument("--threads", type=int, default=4)
    render.add_argument("--view", nargs=3, type=float, default=[-1.2, -1.4, 1.0],
                        metavar=("X", "Y", "Z"), help="Fixed camera direction from the object")
    render.add_argument("--padding", type=float, default=1.4, help="Fixed camera margin, fitted to all selected states")
    render.add_argument("--dt", type=float, help="Explicit constant simulation timestep for time labels")
    render.add_argument("--template", type=Path, default=ROOT / "assets/green_membrane.hipnc")
    render.add_argument("--hython", help="Override detected Houdini Python executable")
    render.add_argument("--resume", action="store_true", help="Resume the identical plan, verifying source/output hashes")
    render.add_argument("--no-labels", action="store_true")
    assemble = commands.add_parser("assemble", help="Existing PNG folder -> cleaned GIF, optional MP4")
    assemble.add_argument("--input", type=Path, required=True)
    assemble.add_argument("--output", type=Path, required=True)
    assemble.add_argument("--force", action="store_true", help="Replace generated animation outputs")
    for command in (render, assemble):
        command.add_argument("--fps", type=float, default=5)
        command.add_argument("--name", default="animation")
        command.add_argument("--keep-logo", action="store_true", help="Disable the white corner box")
        command.add_argument("--mask", nargs=2, type=float, default=None,
                             metavar=("WIDTH_FRACTION", "HEIGHT_FRACTION"),
                             help="Override automatic logo box with fractions of each source image")
        command.add_argument("--mp4", action="store_true", help="Also encode H.264 MP4")
        command.add_argument("--ffmpeg", help="Override detected FFmpeg/hffmpeg executable")
    return cli


def main():
    args = parser().parse_args()
    if args.fps <= 0:
        raise ValueError("--fps must be positive.")
    if args.mask is not None and any(not 0 < v <= 1 for v in args.mask):
        raise ValueError("--mask fractions must be greater than zero and at most one.")
    labels = None
    output = args.output.expanduser().resolve()
    if args.command == "render":
        if not 64 <= args.resolution <= 960:
            raise ValueError("--resolution must be between 64 and 960 for this noncommercial template.")
        if not 1 <= args.samples <= 16 or not 1 <= args.threads <= 64:
            raise ValueError("Use 1..16 samples and 1..64 threads.")
        if args.padding < 1.3:
            raise ValueError("--padding must be at least 1.3 to leave room around the object.")
        config = {name: getattr(args, name) for name in
                  ("frames", "start", "end", "indices", "resolution", "samples", "threads",
                   "view", "padding", "dt", "resume", "fps")}
        config.update(data=str(args.data.expanduser().resolve()), output=str(output),
                      template=str(args.template.expanduser().resolve()))
        environment = os.environ.copy()
        environment["HOUDINI_MAXTHREADS"] = str(args.threads)
        command = [find_hython(args.hython), "-j", str(args.threads), str(ROOT / "render_sequence.py"), json.dumps(config)]
        subprocess.run(command, env=environment, check=True)
        manifest_path = output / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        frames = manifest["plan"]["frames"]
        if not args.no_labels:
            labels = [
                ("t = %.4f   |   step %d" % (f["physical_time"], f["source_step"]))
                if f["physical_time"] is not None else ("step %d" % f["source_step"])
                for f in frames
            ]
        input_dir = output / "raw"
        force = args.resume
    else:
        input_dir = args.input.expanduser().resolve()
        force = args.force
    result = assemble_frames(input_dir, output, fps=args.fps,
                             cover_corner=None if args.keep_logo else (tuple(args.mask) if args.mask else "auto"),
                             labels=labels, make_mp4=args.mp4, ffmpeg=args.ffmpeg,
                             name=args.name, force=force)
    if args.command == "render":
        manifest["status"] = "complete"
        manifest["animation"] = result
        temp = manifest_path.with_suffix(".json.tmp")
        temp.write_text(json.dumps(manifest, indent=2) + "\n")
        temp.replace(manifest_path)
    print(json.dumps({key: result[key] for key in ("frame_count", "dimensions", "gif_fps", "gif_duration_seconds", "gif", "mp4")}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, FileExistsError, RuntimeError, subprocess.CalledProcessError) as exc:
        sys.exit("Error: %s" % exc)
