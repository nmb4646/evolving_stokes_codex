"""Small regression tests for ordering, image fidelity, and GIF sequence integrity."""
import hashlib
from pathlib import Path
import tempfile
import unittest

from PIL import Image, ImageDraw

from animation import assemble_frames, numbered_pngs, _auto_corner_mask


def gif_image_tables(data):
    """Read image descriptor flags without relying on Pillow's RGB conversion."""
    cursor = 13
    flags = data[10]
    if flags & 128:
        cursor += 3 * (2 ** ((flags & 7) + 1))
    local_tables = []
    while data[cursor] != 0x3B:
        kind = data[cursor]
        cursor += 1
        if kind == 0x21:  # extension label, then data subblocks
            cursor += 1
        elif kind == 0x2C:
            flags = data[cursor + 8]
            local_tables.append(bool(flags & 128))
            cursor += 9
            if flags & 128:
                cursor += 3 * (2 ** ((flags & 7) + 1))
            cursor += 1  # LZW minimum code size
        else:
            raise AssertionError(f"Unexpected GIF block {kind:#x}")
        while data[cursor]:
            cursor += data[cursor] + 1
        cursor += 1
    return local_tables


class AnimationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.raw = self.root / "raw"
        self.raw.mkdir()
        self.out = self.root / "output"

    def tearDown(self):
        self.temporary.cleanup()

    def frame(self, name, size=(16, 16), color="red", mode="RGB"):
        path = self.raw / name
        Image.new(mode, size, color).save(path)
        return path

    def test_numeric_sort_and_duplicate_ids(self):
        for name in ("render010.png", "render002.png", "render1.png"):
            self.frame(name)
        self.assertEqual([p.name for p in numbered_pngs(self.raw)],
                         ["render1.png", "render002.png", "render010.png"])
        self.frame("render01.png")
        with self.assertRaisesRegex(ValueError, "Duplicate frame ID"):
            numbered_pngs(self.raw)

    def test_ambiguous_filename_families(self):
        self.frame("frame1.png")
        self.frame("render2.png")
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            numbered_pngs(self.raw)

    def test_rgba_compositing_original_mask_and_even_centered_padding(self):
        self.frame("1.png", (5, 3), (255, 0, 0, 128), "RGBA")
        self.frame("2.png", (8, 6), (0, 0, 255, 0), "RGBA")
        originals = {p.name: hashlib.sha256(p.read_bytes()).digest() for p in self.raw.iterdir()}
        result = assemble_frames(self.raw, self.out, cover_corner=(0.4, 1 / 3))
        self.assertEqual(result["dimensions"], [8, 6])
        self.assertEqual(result["mask_rectangles"], [[3, 2, 5, 3], [4, 4, 8, 6]])
        with Image.open(self.out / "clean_frames/frame_000000.png") as image:
            self.assertEqual(image.mode, "RGB")
            self.assertEqual(image.getpixel((0, 0)), (255, 255, 255))
            self.assertEqual(image.getpixel((1, 1)), (255, 127, 127))
            self.assertEqual(image.getpixel((4, 3)), (255, 255, 255))
        with Image.open(self.out / "clean_frames/frame_000001.png") as image:
            self.assertEqual(image.getextrema(), ((255, 255), (255, 255), (255, 255)))
        self.assertEqual(originals, {p.name: hashlib.sha256(p.read_bytes()).digest() for p in self.raw.iterdir()})

    def test_gif_preserves_identical_frames_durations_and_global_palette(self):
        for index, color in [(1, "red"), (2, "red"), (10, "blue")]:
            self.frame(f"{index}.png", color=color)
        result = assemble_frames(self.raw, self.out, fps=24, cover_corner=None)
        self.assertEqual(result["gif_frame_duration_ms"], 40)
        self.assertEqual(result["gif_fps"], 25)
        self.assertAlmostEqual(result["gif_duration_seconds"], 0.12)
        gif = Path(result["gif"])
        with Image.open(gif) as animation:
            self.assertEqual(animation.n_frames, 3)
            self.assertEqual(animation.info["loop"], 0)
            for index, expected in enumerate([(255, 0, 0), (255, 0, 0), (0, 0, 255)]):
                animation.seek(index)
                self.assertEqual(animation.info["duration"], 40)
                self.assertEqual(animation.convert("RGB").getpixel((0, 0)), expected)
        self.assertEqual(gif_image_tables(gif.read_bytes()), [False, False, False])

    def test_auto_components_keep_membrane_tip_on_same_rows(self):
        frame = Image.new("RGBA", (100, 100), (0, 0, 0, 0))
        draw = ImageDraw.Draw(frame)
        draw.rectangle((20, 10, 65, 70), fill=(0, 180, 90, 255))
        draw.rectangle((40, 65, 51, 90), fill=(0, 180, 90, 255))
        for left in (75, 84, 93):
            draw.rectangle((left, 86, left + 3, 94), fill=(100, 100, 100, 255))
        frame.save(self.raw / "1.png")
        result = assemble_frames(self.raw, self.out)
        self.assertEqual(result["cover_corner"], "auto")
        self.assertEqual(result["mask_rectangles"], [[73, 84, 100, 100]])
        self.assertEqual(result["mask_detection"][0]["method"], "alpha_corner_components")
        self.assertFalse(result["mask_detection"][0]["protected_foreground_overlap"])
        with Image.open(self.out / "clean_frames/frame_000000.png") as cleaned:
            self.assertEqual(cleaned.getpixel((45, 88)), (0, 180, 90))
            self.assertEqual(cleaned.getpixel((76, 88)), (255, 255, 255))

    def test_auto_opaque_fallback_uses_legacy_pixels(self):
        self.frame("1.png", size=(320, 160))
        result = assemble_frames(self.raw, self.out)
        self.assertEqual(result["mask_rectangles"], [[47, 81, 320, 160]])
        self.assertEqual(result["mask_detection"][0]["method"], "legacy_pixel_fallback")
        with Image.open(self.out / "clean_frames/frame_000000.png") as cleaned:
            self.assertEqual(cleaned.getpixel((46, 100)), (255, 0, 0))
            self.assertEqual(cleaned.getpixel((47, 100)), (255, 255, 255))

    def test_actual_sample_alpha_masks_when_available(self):
        raw = Path(__file__).parent / "outputs/shear_green_sample/raw"
        paths = sorted(raw.glob("frame_[0-9][0-9][0-9][0-9].png"))
        if not paths:
            self.skipTest("Rendered membrane sample is not present")
        for path in paths:
            with self.subTest(frame=path.name), Image.open(path) as frame:
                rectangle, info = _auto_corner_mask(frame.convert("RGBA"))
                self.assertEqual(info["method"], "alpha_corner_components")
                self.assertEqual(rectangle, [276, 437, 512, 512])
                self.assertFalse(info["protected_foreground_overlap"])

    def test_output_protection_and_force(self):
        self.frame("1.png")
        first = assemble_frames(self.raw, self.out)
        with self.assertRaises(FileExistsError):
            assemble_frames(self.raw, self.out)
        self.frame("2.png", color="blue")
        second = assemble_frames(self.raw, self.out, force=True)
        self.assertEqual(second["frame_count"], 2)
        self.assertEqual(first["gif"], second["gif"])
        note = self.out / "clean_frames/user_notes.txt"
        note.write_text("Keep me")
        with self.assertRaisesRegex(ValueError, "unrecognized"):
            assemble_frames(self.raw, self.out, force=True)
        self.assertEqual(note.read_text(), "Keep me")

    def test_label_validation_and_source_folder_guard(self):
        self.frame("1.png")
        with self.assertRaisesRegex(ValueError, "labels"):
            assemble_frames(self.raw, self.out, labels=[])
        self.raw.rename(self.root / "clean_frames")
        with self.assertRaisesRegex(ValueError, "source"):
            assemble_frames(self.root / "clean_frames", self.root, force=True)


if __name__ == "__main__":
    unittest.main()
