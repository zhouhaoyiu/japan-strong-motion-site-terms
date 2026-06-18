#!/usr/bin/env python3
"""Create a preview-safe image PDF from a text PDF.

Some PDF preview components fail to render embedded CJK CID fonts and show only
Latin text. This helper renders each page with Poppler and wraps the rendered
pages into a simple image-only PDF, preserving visual layout for review.
"""

from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess
import tempfile
import zlib

from PIL import Image


def _pdf_object(num: int, body: bytes | str) -> tuple[int, bytes]:
    if isinstance(body, str):
        body = body.encode("latin1")
    return num, body


def build_image_pdf(image_paths: list[pathlib.Path], out_path: pathlib.Path, dpi: int) -> None:
    objects: list[tuple[int, bytes]] = []
    obj_id = 1
    catalog_id = obj_id
    obj_id += 1
    pages_id = obj_id
    obj_id += 1
    page_ids: list[int] = []

    for image_path in image_paths:
        image = Image.open(image_path).convert("RGB")
        width_px, height_px = image.size
        width_pt = width_px * 72.0 / dpi
        height_pt = height_px * 72.0 / dpi
        compressed = zlib.compress(image.tobytes(), 6)

        image_id = obj_id
        obj_id += 1
        content_id = obj_id
        obj_id += 1
        page_id = obj_id
        obj_id += 1

        image_stream = (
            f"<< /Type /XObject /Subtype /Image /Width {width_px} /Height {height_px} "
            f"/ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode "
            f"/Length {len(compressed)} >>\nstream\n"
        ).encode("latin1") + compressed + b"\nendstream"

        content = f"q {width_pt:.6f} 0 0 {height_pt:.6f} 0 0 cm /Im0 Do Q".encode(
            "latin1"
        )
        content_stream = (
            f"<< /Length {len(content)} >>\nstream\n".encode("latin1")
            + content
            + b"\nendstream"
        )
        page_body = (
            f"<< /Type /Page /Parent {pages_id} 0 R "
            f"/MediaBox [0 0 {width_pt:.6f} {height_pt:.6f}] "
            f"/Resources << /XObject << /Im0 {image_id} 0 R >> >> "
            f"/Contents {content_id} 0 R >>"
        )

        objects.append(_pdf_object(image_id, image_stream))
        objects.append(_pdf_object(content_id, content_stream))
        objects.append(_pdf_object(page_id, page_body))
        page_ids.append(page_id)

    kids = " ".join(f"{page_id} 0 R" for page_id in page_ids)
    objects.append(_pdf_object(pages_id, f"<< /Type /Pages /Count {len(page_ids)} /Kids [{kids}] >>"))
    objects.append(_pdf_object(catalog_id, f"<< /Type /Catalog /Pages {pages_id} 0 R >>"))
    objects.sort(key=lambda item: item[0])

    pdf = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offsets = [0] * (max(num for num, _ in objects) + 1)
    for num, body in objects:
        offsets[num] = len(pdf)
        pdf.extend(f"{num} 0 obj\n".encode("latin1"))
        pdf.extend(body)
        pdf.extend(b"\nendobj\n")

    xref_offset = len(pdf)
    pdf.extend(f"xref\n0 {len(offsets)}\n".encode("latin1"))
    pdf.extend(b"0000000000 65535 f \n")
    for num in range(1, len(offsets)):
        pdf.extend(f"{offsets[num]:010d} 00000 n \n".encode("latin1"))
    pdf.extend(
        (
            f"trailer << /Size {len(offsets)} /Root {catalog_id} 0 R >>\n"
            f"startxref\n{xref_offset}\n%%EOF\n"
        ).encode("latin1")
    )
    out_path.write_bytes(pdf)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_pdf", type=pathlib.Path)
    parser.add_argument("output_pdf", type=pathlib.Path)
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()

    temp_dir = pathlib.Path(tempfile.mkdtemp(prefix="pdf_preview_flatten_"))
    try:
        prefix = temp_dir / "page"
        subprocess.run(
            ["pdftoppm", "-png", "-r", str(args.dpi), str(args.source_pdf), str(prefix)],
            check=True,
        )
        image_paths = sorted(temp_dir.glob("page-*.png"))
        if not image_paths:
            raise RuntimeError("pdftoppm produced no rendered pages")
        build_image_pdf(image_paths, args.output_pdf, args.dpi)
        print(
            f"Wrote {args.output_pdf} with {len(image_paths)} pages, "
            f"size={args.output_pdf.stat().st_size} bytes"
        )
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
