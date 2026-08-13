#!/usr/bin/env python3
"""Add Gamma and simulation-time labels to membrane-collapse PNG frames."""

from __future__ import annotations

import argparse
import os
import re
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-annotate-collapse")
)

import matplotlib as mpl
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.mathtext import MathTextParser
from PIL import Image


BASE_TIMESTEP = 1.2e-4
SERIES_PATTERN = re.compile(r"^collapse_(?P<gamma>[0-9]+(?:\.[0-9]+)?)$")
mpl.rcParams["mathtext.fontset"] = "cm"
MATH_TEXT_PARSER = MathTextParser("agg")


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_input = script_dir / "data" / "collapse_videos"
    default_output = script_dir / "data" / "collapse_videos_annotated"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=default_input)
    parser.add_argument("--output", type=Path, default=default_output)
    parser.add_argument(
        "--font-size",
        type=int,
        default=None,
        help="Computer Modern font size in pixels (default: 3.8%% of frame size).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace PNGs that already exist in the output directory.",
    )
    return parser.parse_args()


def numbered_pngs(series_dir: Path) -> list[tuple[int, Path]]:
    frames: list[tuple[int, Path]] = []
    for path in series_dir.glob("*.png"):
        try:
            frame_number = int(path.stem)
        except ValueError as exc:
            raise ValueError(f"PNG filename is not an integer frame: {path}") from exc
        frames.append((frame_number, path))

    if not frames:
        raise ValueError(f"No PNG frames found in {series_dir}")
    return sorted(frames)


def annotate_frame(
    source: Path,
    destination: Path,
    label: str,
    requested_font_size: int | None,
    force: bool,
) -> None:
    if destination.exists() and not force:
        raise FileExistsError(
            f"Output already exists: {destination} (use --force to replace it)"
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as opened:
        image = opened.convert("RGB")
        font_size = requested_font_size or max(20, round(min(image.size) * 0.038))
        parsed = MATH_TEXT_PARSER.parse(
            label,
            dpi=72,
            prop=FontProperties(size=font_size),
        )
        alpha = Image.fromarray(np.asarray(parsed[5]))
        rendered_label = Image.new("RGBA", alpha.size, (0, 0, 0, 0))
        rendered_label.putalpha(alpha)
        margin = round(min(image.size) * 0.03)
        position = (margin, image.height - margin - rendered_label.height)
        image.paste(rendered_label, position, rendered_label)

        temporary = destination.with_name(f".{destination.name}.tmp")
        image.save(temporary, format="PNG", compress_level=6)
        os.replace(temporary, destination)


def main() -> None:
    args = parse_args()
    input_root = args.input.resolve()
    output_root = args.output.resolve()

    if not input_root.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_root}")
    if args.font_size is not None and args.font_size <= 0:
        raise ValueError("--font-size must be positive")
    if output_root == input_root or input_root in output_root.parents:
        raise ValueError("Output must not be the input directory or one of its children")

    series: list[tuple[float, Path]] = []
    for directory in input_root.iterdir():
        if not directory.is_dir():
            continue
        match = SERIES_PATTERN.fullmatch(directory.name)
        if match:
            series.append((float(match.group("gamma")), directory))

    if not series:
        raise ValueError(f"No collapse_<Gamma> series found in {input_root}")

    total = 0
    for gamma, series_dir in sorted(series):
        frames = numbered_pngs(series_dir)
        frame_zero = frames[0][0]
        timestep = BASE_TIMESTEP * gamma**-2.0
        gamma_text = f"{gamma:g}"
        destination_dir = output_root / series_dir.name

        for frame_number, source in frames:
            time = (frame_number - frame_zero) * timestep
            label = rf"$\Gamma = {gamma_text},\ t = {time:.6f}$"
            annotate_frame(
                source,
                destination_dir / source.name,
                label,
                args.font_size,
                args.force,
            )
            total += 1

        print(
            f"{series_dir.name}: {len(frames)} frames, frame0={frame_zero}, "
            f"dt={timestep:.12g}"
        )

    print(f"Wrote {total} annotated frames to {output_root}")


if __name__ == "__main__":
    main()
