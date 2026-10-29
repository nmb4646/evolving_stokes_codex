#!/usr/bin/env python3
"""Create truncated H.264 MP4 files from annotated collapse frames."""

from __future__ import annotations

import argparse
import os
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path

from create_collapse_avis import END_FRAMES, numbered_pngs


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_input = script_dir / "data" / "collapse_videos_annotated"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=default_input)
    parser.add_argument("--output", type=Path, default=default_input)
    parser.add_argument("--fps", type=int, default=28)
    parser.add_argument("--bitrate", type=int, default=8_000_000)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace MP4 files that already exist in the output directory.",
    )
    return parser.parse_args()


def mp4_sample_count(path: Path) -> int:
    """Read the video sample count from the MP4 sample-size atom."""
    data = path.read_bytes()
    atom_type = data.find(b"stsz")
    if atom_type < 4 or atom_type + 16 > len(data):
        raise ValueError(f"MP4 sample-size atom not found in {path}")

    atom_size = struct.unpack_from(">I", data, atom_type - 4)[0]
    if atom_size < 20 or atom_type - 4 + atom_size > len(data):
        raise ValueError(f"Invalid MP4 sample-size atom in {path}")
    return struct.unpack_from(">I", data, atom_type + 12)[0]


def require_gstreamer_component(component: str) -> None:
    result = subprocess.run(
        ["gst-inspect-1.0", component],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if result.returncode:
        raise RuntimeError(f"Required GStreamer component is unavailable: {component}")


def encode_series(
    series_dir: Path,
    output_path: Path,
    end_frame: int,
    fps: int,
    bitrate: int,
    force: bool,
) -> int:
    frames = numbered_pngs(series_dir, end_frame)
    if output_path.exists() and not force:
        raise FileExistsError(
            f"Output already exists: {output_path} (use --force to replace it)"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_path.with_name(f".{output_path.stem}.tmp.mp4")
    temporary_output.unlink(missing_ok=True)

    try:
        with tempfile.TemporaryDirectory(prefix=f"{series_dir.name}-mp4-") as temp_name:
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
                "video/x-raw,format=I420",
                "!",
                "openh264enc",
                f"bitrate={bitrate}",
                "rate-control=quality",
                "enable-frame-skip=false",
                "complexity=high",
                "usage-type=screen",
                f"gop-size={2 * fps}",
                "!",
                "h264parse",
                "!",
                "video/x-h264,stream-format=avc,alignment=au",
                "!",
                "mp4mux",
                "faststart=true",
                "!",
                "filesink",
                f"location={temporary_output}",
            ]
            result = subprocess.run(command, text=True, capture_output=True)
            if result.returncode:
                raise RuntimeError(
                    f"MP4 encoding failed for {series_dir.name}:\n{result.stderr.strip()}"
                )

        sample_count = mp4_sample_count(temporary_output)
        if sample_count != len(frames):
            raise ValueError(
                f"Encoded {sample_count} frames for {series_dir.name}; "
                f"expected {len(frames)}"
            )

        os.replace(temporary_output, output_path)
        return sample_count
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
    if args.bitrate <= 0:
        raise ValueError("--bitrate must be positive")
    if shutil.which("gst-launch-1.0") is None or shutil.which("gst-inspect-1.0") is None:
        raise FileNotFoundError("GStreamer is required to encode the MP4 files")
    for component in ("pngdec", "videoconvert", "openh264enc", "h264parse", "mp4mux"):
        require_gstreamer_component(component)

    for series_name, end_frame in END_FRAMES.items():
        series_dir = input_root / series_name
        output_path = output_root / f"{series_name}_to_{end_frame}.mp4"
        total = encode_series(
            series_dir,
            output_path,
            end_frame,
            args.fps,
            args.bitrate,
            args.force,
        )
        print(
            f"{output_path.name}: {total} frames, 810x810, "
            f"{args.fps} fps, H.264 MP4"
        )


if __name__ == "__main__":
    main()
