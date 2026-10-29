"""Assemble a numbered PNG sequence into a GIF and optionally an H.264 MP4.

Only processed copies are masked/padded; source images are always read-only.
The GIF uses one palette for the entire sequence and retains duplicate frames.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from PIL import GifImagePlugin, Image, ImageDraw, ImageFont

# Pillow 9.0 ships with Houdini/system Python; newer versions use enum groups.
_RESAMPLE = getattr(Image, "Resampling", Image)
_QUANTIZE = getattr(Image, "Quantize", Image)
_DITHER = getattr(Image, "Dither", Image)


def numbered_pngs(directory: Path) -> list[Path]:
    """Return a single filename family's PNGs in numeric frame order.

    Leading zeros do not distinguish frame IDs: 1.png and 001.png conflict.
    Fixed digits in the prefix are allowed, e.g. camera2_frame001.png.
    Different prefixes/suffixes require separate input directories.
    """
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise ValueError(f"Input directory does not exist: {directory}")
    files = [p for p in directory.iterdir() if p.is_file() and p.suffix.lower() == ".png"]
    if not files:
        raise ValueError(f"No PNG frames found in {directory}")
    family = None
    indexed = {}
    for path in files:
        match = re.fullmatch(r"(.*?)([0-9]+)([^0-9]*)", path.stem)
        if not match:
            raise ValueError(f"PNG filename needs a numeric frame ID: {path.name}")
        prefix, digits, suffix = match.groups()
        current = (prefix.casefold(), suffix.casefold())
        if family is None:
            family = current
        elif current != family:
            raise ValueError("Ambiguous PNG sequence: multiple filename families; "
                             "put one numbered sequence in the input directory")
        index = int(digits)
        if index in indexed:
            raise ValueError(f"Duplicate frame ID {index}: {indexed[index].name} and {path.name}")
        indexed[index] = path
    return [indexed[i] for i in sorted(indexed)]



def _auto_corner_mask(rgba: Image.Image) -> tuple[list[int], dict]:
    """Find isolated alpha components in the bottom-right corner.

    Looking at components instead of whole-image row bands permits a membrane
    tip to reach the logo's rows while remaining safely to its left. A candidate
    must be wholly inside the bottom 20% and right half, have a broad horizontal
    extent, and extend toward the right/bottom edges. The final rectangle is
    checked pixel-by-pixel: it may contain only these components and transparency.
    Opaque or unconfident images use the legacy 273 x 79 pixel corner mask.
    """
    width, height = rgba.size
    fallback = [max(0, width - 273), max(0, height - 79), width, height]
    alpha = rgba.getchannel("A")
    if alpha.getextrema()[0] == 255:
        return fallback, {"method": "legacy_pixel_fallback", "reason": "opaque image",
                          "fallback_size": [273, 79]}
    values = alpha.tobytes()
    visited = bytearray(len(values))
    selected = bytearray(len(values))
    corner_x, corner_y = math.ceil(width * 0.5), math.ceil(height * 0.8)
    candidates = []
    for y in range(corner_y, height):
        for x in range(corner_x, width):
            seed = y * width + x
            if not values[seed] or visited[seed]:
                continue
            stack = [seed]
            visited[seed] = 1
            pixels = []
            x0, x1, y0, y1 = x, x + 1, y, y + 1
            while stack:
                index = stack.pop()
                py, px = divmod(index, width)
                pixels.append(index)
                x0, x1 = min(x0, px), max(x1, px + 1)
                y0, y1 = min(y0, py), max(y1, py + 1)
                for ny in range(max(0, py - 1), min(height, py + 2)):
                    for nx in range(max(0, px - 1), min(width, px + 2)):
                        neighbor = ny * width + nx
                        if values[neighbor] and not visited[neighbor]:
                            visited[neighbor] = 1
                            stack.append(neighbor)
            if x0 >= corner_x and y0 >= corner_y:
                candidates.append([x0, y0, x1, y1])
                for index in pixels:
                    selected[index] = 1
    reason = "no isolated corner components"
    if candidates:
        bbox = [min(b[0] for b in candidates), min(b[1] for b in candidates),
                max(b[2] for b in candidates), max(b[3] for b in candidates)]
        # The antialiased logo may be one connected component or many letters.
        if bbox[2] >= width * 0.75 and bbox[3] >= height * 0.85 and bbox[2] - bbox[0] >= max(width * 0.1, 2 * (bbox[3] - bbox[1])):
            rectangle = [max(0, bbox[0] - 2), max(0, bbox[1] - 2), width, height]
            protected = any(values[y * width + x] and not selected[y * width + x]
                            for y in range(rectangle[1], height)
                            for x in range(rectangle[0], width))
            if not protected:
                return rectangle, {"method": "alpha_corner_components", "detected_bbox": bbox,
                                   "component_count": len(candidates), "margin_pixels": 2,
                                   "protected_foreground_overlap": False}
            reason = "corner rectangle would intersect other foreground"
        else:
            reason = "isolated components do not match the expected corner location/extent"
    return fallback, {"method": "legacy_pixel_fallback", "reason": reason,
                      "fallback_size": [273, 79]}


def _font(size: int):
    for font in ("DejaVuSans.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(font, size)
        except OSError:
            pass
    return ImageFont.load_default()


def _palette(paths: list[Path]) -> Image.Image:
    # Include every frame while bounding the palette-training mosaic's memory.
    side = min(128, max(1, int(math.sqrt(1_000_000 / len(paths)))))
    columns = min(32, len(paths))
    mosaic = Image.new("RGB", (columns * side, math.ceil(len(paths) / columns) * side), "white")
    for i, path in enumerate(paths):
        with Image.open(path) as frame:
            thumb = frame.convert("RGB")
            thumb.thumbnail((side, side), _RESAMPLE.LANCZOS)
            mosaic.paste(thumb, ((i % columns) * side, (i // columns) * side))
    return mosaic.quantize(colors=256, method=_QUANTIZE.MEDIANCUT,
                           dither=_DITHER.NONE)


def _write_gif(paths: list[Path], destination: Path, duration_ms: int) -> None:
    palette = _palette(paths)
    colors = palette.getpalette()
    background = min(range(len(colors) // 3),
                     key=lambda i: sum((255 - v) ** 2 for v in colors[3 * i:3 * i + 3]))
    with destination.open("wb") as stream:
        for i, path in enumerate(paths):
            with Image.open(path) as source:
                frame = source.convert("RGB").quantize(palette=palette, dither=_DITHER.NONE)
            if i == 0:
                header, _ = GifImagePlugin.getheader(
                    frame, info={"duration": duration_ms, "background": background, "optimize": False})
                for block in header:
                    stream.write(block)
                # Explicit infinite-loop extension also works with Pillow 9.0,
                # whose getheader emits only the global color table/header.
                stream.write(b"\x21\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00")
            # Pillow's save_all coalesces identical frames. Explicit image blocks
            # retain every sampled timestep, all using the global palette.
            for block in GifImagePlugin.getdata(
                    frame, duration=duration_ms, disposal=2, optimize=False,
                    include_color_table=False):
                stream.write(block)
        stream.write(b";")
    with Image.open(destination) as check:
        if check.n_frames != len(paths):
            raise RuntimeError("GIF verification failed: output lost sequence frames")


def _find_ffmpeg(explicit: str | None) -> str:
    if explicit:
        executable = shutil.which(explicit)
        if executable:
            return executable
        raise ValueError(f"FFmpeg executable is unavailable: {explicit}")
    executable = shutil.which("ffmpeg")
    if executable:
        return executable
    houdini = sorted(glob.glob("/opt/hfs*/bin/hffmpeg"), reverse=True)
    if houdini:
        return houdini[0]
    raise RuntimeError("MP4 requested but FFmpeg was not found. Supply --ffmpeg /path/to/ffmpeg.")


def _write_mp4(directory: Path, destination: Path, fps: float, executable: str) -> str:
    result = subprocess.run([executable, "-hide_banner", "-encoders"],
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise RuntimeError(f"Cannot query FFmpeg encoders: {result.stderr[-2000:]}")
    available = {line.split()[1] for line in result.stdout.splitlines()
                 if len(line.split()) >= 2 and line.lstrip().startswith("V")}
    encoder = next((v for v in ("libx264", "libopenh264") if v in available), None)
    if encoder is None:
        raise RuntimeError("This FFmpeg build has no supported software H.264 encoder "
                           "(libx264 or libopenh264). Supply another executable with --ffmpeg.")
    command = [executable, "-hide_banner", "-loglevel", "error", "-y", "-threads", "4",
               "-framerate", str(fps), "-start_number", "0", "-i",
               str(directory / "frame_%06d.png"), "-an", "-c:v", encoder,
               "-threads", "4", "-pix_fmt", "yuv420p", "-movflags", "+faststart"]
    if encoder == "libx264":
        command += ["-crf", "18", "-preset", "medium"]
    else:
        command += ["-b:v", "8M"]
    command += [str(destination)]
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode or not destination.is_file() or destination.stat().st_size == 0:
        raise RuntimeError(f"FFmpeg MP4 encoding failed ({encoder}): {result.stderr[-4000:]}")
    return encoder


def assemble_frames(input_dir: Path, output_dir: Path, *, fps: float = 5,
                    cover_corner: tuple[float, float] | str | None = "auto",
                    labels: list[str] | None = None, make_mp4: bool = False,
                    ffmpeg: str | None = None, name: str = "animation",
                    force: bool = False) -> dict:
    """Write clean_frames/frame_000000.png etc. and name.gif atomically.

    cover_corner="auto" detects detached bottom-right alpha components; opaque or
    unconfident images fall back to a 273 x 79 pixel mask (clamped to frame size).
    An explicit (width_fraction, height_fraction) masks each ORIGINAL frame's
    bottom-right corner. None disables masking. Masking precedes centered padding.
    GIF frame duration rounds to the nearest 10 ms; the return value reports the
    actual GIF rate. MP4, when requested, uses the requested fps.
    Existing generated outputs require force=True. Unrecognized clean_frames
    contents are never removed, even with force=True.
    """
    fps = float(fps)
    if not math.isfinite(fps) or not 0 < fps <= 100:
        raise ValueError("fps must be finite, greater than zero, and at most 100")
    duration_ms = int(math.floor(100 / fps + 0.5)) * 10
    if duration_ms > 655350:
        raise ValueError("Requested frame duration exceeds GIF's 655.35 second limit")
    if not name or name in {".", ".."} or Path(name).name != name or "\\" in name:
        raise ValueError("name must be a filename stem without directory separators")
    if isinstance(cover_corner, str) and cover_corner != "auto":
        raise ValueError('cover_corner string must be "auto"')
    if cover_corner is not None and cover_corner != "auto":
        if len(cover_corner) != 2 or not all(math.isfinite(v) and 0 <= v <= 1 for v in cover_corner):
            raise ValueError("cover_corner must contain width and height fractions between 0 and 1")
        cover_corner = tuple(float(v) for v in cover_corner)
    paths = numbered_pngs(input_dir)
    if labels is not None and (len(labels) != len(paths) or not all(isinstance(v, str) for v in labels)):
        raise ValueError("labels must contain exactly one string per sorted input frame")
    output_dir = Path(output_dir).expanduser().resolve()
    clean_target = output_dir / "clean_frames"
    if clean_target == Path(input_dir).expanduser().resolve() or clean_target in paths[0].parents:
        raise ValueError("Output clean_frames must not contain the source images")
    gif_target = output_dir / f"{name}.gif"
    mp4_target = output_dir / f"{name}.mp4"
    targets = [clean_target, gif_target] + ([mp4_target] if make_mp4 else [])
    for target in targets:
        if target.is_symlink():
            raise ValueError(f"Refusing to replace a symlink output: {target}")
        if target.exists() and not force:
            raise FileExistsError(f"Output already exists: {target}; choose a new directory or use force=True")
    for target in targets[1:]:
        if target.exists() and not target.is_file():
            raise ValueError(f"Animation output exists and is not a regular file: {target}")
    if clean_target.exists():
        if not clean_target.is_dir():
            raise ValueError(f"clean_frames exists and is not a directory: {clean_target}")
        if any(not p.is_file() or p.is_symlink() or not re.fullmatch(r"frame_[0-9]{6,}\.png", p.name)
               for p in clean_target.iterdir()):
            raise ValueError("Existing clean_frames contains unrecognized files; choose a new output directory")
    dimensions = []
    for path in paths:
        with Image.open(path) as source:
            if source.format != "PNG" or getattr(source, "n_frames", 1) != 1:
                raise ValueError(f"Expected a static PNG: {path}")
            dimensions.append(source.size)
            source.verify()
    width = 2 * math.ceil(max(w for w, _ in dimensions) / 2)
    height = 2 * math.ceil(max(h for _, h in dimensions) / 2)
    executable = _find_ffmpeg(ffmpeg) if make_mp4 else None
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".animation-", dir=output_dir) as temporary:
        stage = Path(temporary)
        clean = stage / "clean_frames"
        clean.mkdir()
        processed = []
        mask_rectangles = []
        mask_detection = []
        font = _font(max(11, round(min(width, height) * 0.025)))
        margin = max(4, round(min(width, height) * 0.015))
        for i, path in enumerate(paths):
            with Image.open(path) as source:
                source.load()
                rgba = source.convert("RGBA")
            frame = Image.alpha_composite(Image.new("RGBA", rgba.size, "white"), rgba).convert("RGB")
            rectangle = None
            detection = {"method": "disabled"}
            if cover_corner == "auto":
                rectangle, detection = _auto_corner_mask(rgba)
            elif cover_corner is not None:
                mask_width = math.ceil(frame.width * cover_corner[0])
                mask_height = math.ceil(frame.height * cover_corner[1])
                if mask_width and mask_height:
                    rectangle = [frame.width - mask_width, frame.height - mask_height,
                                 frame.width, frame.height]
                    detection = {"method": "explicit_fractions"}
            if rectangle is not None:
                frame.paste("white", tuple(rectangle))
            mask_rectangles.append(rectangle)
            mask_detection.append(detection)
            padded = Image.new("RGB", (width, height), "white")
            padded.paste(frame, ((width - frame.width) // 2, (height - frame.height) // 2))
            if labels is not None and labels[i]:
                draw = ImageDraw.Draw(padded)
                bounds = draw.textbbox((margin, margin), labels[i], font=font)
                draw.rectangle((bounds[0] - 3, bounds[1] - 3, bounds[2] + 3, bounds[3] + 3), fill="white")
                draw.text((margin, margin), labels[i], fill=(35, 35, 35), font=font)
            destination = clean / f"frame_{i:06d}.png"
            padded.save(destination)
            processed.append(destination)
        staged_gif = stage / gif_target.name
        _write_gif(processed, staged_gif, duration_ms)
        encoder = _write_mp4(clean, stage / mp4_target.name, fps, executable) if make_mp4 else None
        # Every output is complete before publishing. Keep old outputs in the
        # private staging directory until all replacements succeed, for rollback.
        staged_targets = [(clean, clean_target), (staged_gif, gif_target)]
        if make_mp4:
            staged_targets.append((stage / mp4_target.name, mp4_target))
        backups = []
        published = []
        try:
            for i, (source, target) in enumerate(staged_targets):
                if target.exists():
                    if not force:
                        raise FileExistsError(f"Output appeared during assembly: {target}")
                    backup = stage / f"previous_{i}"
                    target.rename(backup)
                    backups.append((backup, target))
                source.rename(target)
                published.append((target, source))
        except Exception:
            for target, source in reversed(published):
                target.rename(source)
            for backup, target in reversed(backups):
                backup.rename(target)
            raise
    return {
        "frame_count": len(paths), "dimensions": [width, height],
        "source_frames": [str(p) for p in paths], "source_dimensions": [list(s) for s in dimensions],
        "requested_fps": fps, "gif_fps": 1000 / duration_ms,
        "gif_frame_duration_ms": duration_ms,
        "gif_duration_seconds": len(paths) * duration_ms / 1000,
        "gif": str(gif_target), "clean_frames": str(clean_target),
        "cover_corner": list(cover_corner) if isinstance(cover_corner, tuple) else cover_corner,
        "mask_rectangles": mask_rectangles, "mask_detection": mask_detection, "labels": labels,
        "mp4": str(mp4_target) if make_mp4 else None,
        "mp4_encoder": encoder, "ffmpeg": executable,
        "mp4_duration_seconds": len(paths) / fps if make_mp4 else None,
        "palette": "one shared 256-color palette; dithering disabled", "loop": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--fps", type=float, default=5)
    parser.add_argument("--cover-corner", nargs=2, type=float, metavar=("WIDTH", "HEIGHT"), default=None,
                        help="Override automatic alpha detection with width/height fractions")
    parser.add_argument("--no-cover-corner", action="store_true")
    parser.add_argument("--label", action="append", dest="labels", help="One label per frame, in numeric order; repeat for each frame")
    parser.add_argument("--mp4", action="store_true")
    parser.add_argument("--ffmpeg")
    parser.add_argument("--name", default="animation")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    try:
        result = assemble_frames(args.input_dir, args.output_dir, fps=args.fps,
                                 cover_corner=None if args.no_cover_corner else (tuple(args.cover_corner) if args.cover_corner else "auto"),
                                 labels=args.labels, make_mp4=args.mp4, ffmpeg=args.ffmpeg,
                                 name=args.name, force=args.force)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Animation assembly failed: {exc}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
