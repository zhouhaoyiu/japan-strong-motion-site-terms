#!/usr/bin/env python3
"""Download public J-SHIS/NIED ground-motion flatfile archives with resume.

The official server supports HTTP range requests. This downloader stores chunks
separately, verifies sizes, assembles the zip, and tests the archive.
"""

from __future__ import annotations

import argparse
import concurrent.futures as futures
import math
import time
import zipfile
from pathlib import Path

import requests


URLS = {
    "v2024": "https://www.j-shis.bosai.go.jp/labs/ground-motion-flatfile/data/v2024/flatfile-v2024.zip",
    "sub1-v2024": "https://www.j-shis.bosai.go.jp/labs/ground-motion-flatfile/data/v2024/flatfile_sub1-v2024.zip",
    "sub2-v2024": "https://www.j-shis.bosai.go.jp/labs/ground-motion-flatfile/data/v2024/flatfile_sub2-v2024.zip",
}
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "work" / "external_data" / "jshis_gmf"
USER_AGENT = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Safari/537.36"
)


def remote_size(url: str) -> int:
    response = requests.head(url, allow_redirects=True, timeout=30)
    response.raise_for_status()
    size = response.headers.get("Content-Length")
    if not size:
        raise RuntimeError(f"No Content-Length for {url}")
    return int(size)


def download_range(
    url: str,
    part_path: Path,
    start: int,
    end: int,
    retries: int = 5,
) -> int:
    expected = end - start + 1
    if part_path.exists() and part_path.stat().st_size == expected:
        return expected
    tmp_path = part_path.with_suffix(part_path.suffix + ".tmp")
    headers = {
        "Range": f"bytes={start}-{end}",
        "User-Agent": USER_AGENT,
    }
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, headers=headers, stream=True, timeout=60) as response:
                if response.status_code != 206:
                    raise RuntimeError(f"Expected 206 range response, got {response.status_code}")
                with tmp_path.open("wb") as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            actual = tmp_path.stat().st_size
            if actual != expected:
                raise RuntimeError(f"Chunk size mismatch: expected {expected}, got {actual}")
            tmp_path.replace(part_path)
            return expected
        except Exception:
            if tmp_path.exists():
                tmp_path.unlink()
            if attempt == retries:
                raise
            time.sleep(2 * attempt)
    raise RuntimeError("unreachable")


def assemble(parts: list[Path], output_path: Path, expected_size: int) -> None:
    tmp_path = output_path.with_suffix(output_path.suffix + ".assembling")
    with tmp_path.open("wb") as out:
        for part in parts:
            with part.open("rb") as handle:
                while True:
                    chunk = handle.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
    actual = tmp_path.stat().st_size
    if actual != expected_size:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"Assembled size mismatch: expected {expected_size}, got {actual}")
    tmp_path.replace(output_path)


def test_zip(path: Path) -> None:
    if not zipfile.is_zipfile(path):
        raise RuntimeError(f"Not a valid zip file: {path}")
    with zipfile.ZipFile(path) as zf:
        bad = zf.testzip()
        if bad:
            raise RuntimeError(f"Zip test failed at member: {bad}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", choices=sorted(URLS), nargs="?", default="sub1-v2024")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--chunk-mb", type=int, default=16)
    parser.add_argument("--keep-parts", action="store_true")
    args = parser.parse_args()

    url = URLS[args.name]
    args.data_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.data_dir / Path(url).name
    if output_path.exists():
        test_zip(output_path)
        print(f"Archive already complete: {output_path}")
        return

    size = remote_size(url)
    chunk_size = args.chunk_mb * 1024 * 1024
    n_chunks = math.ceil(size / chunk_size)
    part_dir = args.data_dir / f"{output_path.name}.parts"
    part_dir.mkdir(parents=True, exist_ok=True)
    ranges = []
    for idx in range(n_chunks):
        start = idx * chunk_size
        end = min(size - 1, start + chunk_size - 1)
        part_path = part_dir / f"part-{idx:04d}"
        ranges.append((idx, part_path, start, end))

    done_bytes = sum(
        part.stat().st_size
        for _, part, start, end in ranges
        if part.exists() and part.stat().st_size == end - start + 1
    )
    print(
        f"Downloading {output_path.name}: {size / 1024 / 1024:.1f} MB, "
        f"{n_chunks} chunks, already {done_bytes / 1024 / 1024:.1f} MB",
        flush=True,
    )

    completed = 0
    completed_bytes = done_bytes
    with futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_map = {
            executor.submit(download_range, url, part, start, end): (idx, part, start, end)
            for idx, part, start, end in ranges
            if not (part.exists() and part.stat().st_size == end - start + 1)
        }
        for future in futures.as_completed(future_map):
            idx, _, start, end = future_map[future]
            bytes_done = future.result()
            completed += 1
            completed_bytes += bytes_done
            print(
                f"Chunk {idx + 1}/{n_chunks} complete; "
                f"{completed_bytes / size * 100:.1f}% assembled-ready",
                flush=True,
            )

    parts = [part for _, part, _, _ in ranges]
    assemble(parts, output_path, size)
    test_zip(output_path)
    print(f"Downloaded and verified: {output_path}", flush=True)

    if not args.keep_parts:
        for part in parts:
            part.unlink(missing_ok=True)
        part_dir.rmdir()


if __name__ == "__main__":
    main()
