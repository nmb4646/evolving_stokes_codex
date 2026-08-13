#!/usr/bin/env python3
"""Create truncated Motion JPEG AVI files from annotated collapse frames."""

from __future__ import annotations

import argparse
import os
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path


END_FRAMES = {
    "collapse_0.20": 17000,
    "collapse_2.10": 14995,
    "collapse_3.00": 26900,
}


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_input = script_dir / "data" / "collapse_videos_annotated"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=default_input)
    parser.add_argument("--output", type=Path, default=default_input)
    parser.add_argument("--fps", type=int, default=14)
    parser.add_argument("--quality", type=int, default=95)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace AVI files that already exist in the output directory.",
    )
    return parser.parse_args()


def numbered_pngs(series_dir: Path, end_frame: int) -> list[tuple[int, Path]]:
    frames: list[tuple[int, Path]] = []
    for path in series_dir.glob("*.png"):
        try:
            frame_number = int(path.stem)
        except ValueError as exc:
            raise ValueError(f"PNG filename is not an integer frame: {path}") from exc
        if frame_number <= end_frame:
            frames.append((frame_number, path))

    frames.sort()
    if not frames or frames[-1][0] != end_frame:
        raise ValueError(f"Required ending frame {end_frame}.png is missing in {series_dir}")
    return frames


def avi_metadata(path: Path) -> tuple[int, int, int, int]:
    with path.open("rb") as handle:
        header = handle.read(1024 * 1024)

    chunk = header.find(b"avih")
    if chunk < 0 or chunk + 48 > len(header):
        raise ValueError(f"AVI main header not found in {path}")
    payload = chunk + 8
    microseconds_per_frame = struct.unpack_from("<I", header, payload)[0]
    total_frames = struct.unpack_from("<I", header, payload + 16)[0]
    width = struct.unpack_from("<I", header, payload + 32)[0]
    height = struct.unpack_from("<I", header, payload + 36)[0]
    return total_frames, width, height, microseconds_per_frame


def encode_series(
    series_dir: Path,
    output_path: Path,
    end_frame: int,
    fps: int,
    quality: int,
    force: bool,
) -> tuple[int, int, int, int]:
    frames = numbered_pngs(series_dir, end_frame)
    if output_path.exists() and not force:
        raise FileExistsError(
            f"Output already exists: {output_path} (use --force to replace it)"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_path.with_name(f".{output_path.name}.tmp")
    temporary_output.unlink(missing_ok=True)

    try:
        with tempfile.TemporaryDirectory(prefix=f"{series_dir.name}-avi-") as temp_name:
            temp_dir = Path(temp_name)
            for index, (_, source) in enumerate(frames):
                (temp_dir / f"{index:06d}.png").symlink_to(source.resolve())

            command = [
                "gst-launch-1.0",
                "-e",
                "-q",
                "multifilesrc",
                f"location={temp_dir / '%06d.png'}",
                "start-index=0",
                f"stop-index={len(frames) - 1}",
                f"caps=image/png,framerate={fps}/1",
                "!",
                "pngdec",
                "!",
                "videoconvert",
                "!",
                "jpegenc",
                f"quality={quality}",
                "!",
                "avimux",
                "!",
                "filesink",
                f"location={temporary_output}",
            ]
            subprocess.run(command, check=True)

        metadata = avi_metadata(temporary_output)
        total_frames, width, height, microseconds_per_frame = metadata
        if total_frames != len(frames):
            raise ValueError(
                f"Encoded {total_frames} frames for {series_dir.name}; "
                f"expected {len(frames)}"
            )
        if (width, height) != (810, 810):
            raise ValueError(f"Unexpected AVI dimensions: {width}x{height}")

        os.replace(temporary_output, output_path)
        return metadata
    finally:
        temporary_output.unlink(missing_ok=True)


def main() -> None:
    args = parse_args()
    input_root = args.input.resolve()
    output_root = args.output.resolve()

    if not input_root.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_root}")
    if not 1 <= args.fps <= 120:
        raise ValueError("--fps must be between 1 and 120")
    if not 0 <= args.quality <= 100:
        raise ValueError("--quality must be between 0 and 100")
    if shutil.which("gst-launch-1.0") is None:
        raise FileNotFoundError("gst-launch-1.0 is required to encode the AVI files")

    for series_name, end_frame in END_FRAMES.items():
        series_dir = input_root / series_name
        output_path = output_root / f"{series_name}_to_{end_frame}.avi"
        total, width, height, microseconds = encode_series(
            series_dir,
            output_path,
            end_frame,
            args.fps,
            args.quality,
            args.force,
        )
        print(
            f"{output_path.name}: {total} frames, {width}x{height}, "
            f"{1_000_000 / microseconds:.6g} fps"
        )


if __name__ == "__main__":
    main()
