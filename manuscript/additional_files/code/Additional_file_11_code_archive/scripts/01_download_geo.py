#!/usr/bin/env python
"""Download GEO processed files for the IPF oligonucleotide ML project.

The script downloads GEO series matrix files, family SOFT files, and
supplementary processed files from the NCBI GEO FTP HTTPS mirror.
It does not download SRA/FASTQ raw sequencing files.
"""

from __future__ import annotations

import argparse
import csv
import html.parser
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path


GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
GEO_QUERY_BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        for key, value in attrs:
            if key == "href" and value:
                self.links.append(value)


@dataclass
class DownloadRecord:
    dataset_id: str
    file_type: str
    url: str
    local_path: str
    status: str
    size_bytes: int | None = None
    message: str = ""


def load_dotenv(path: Path) -> None:
    if not path.exists():
        return
    for raw_line in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        os.environ.setdefault(key.strip(), value.strip().strip('"').strip("'"))


def gse_group(accession: str) -> str:
    match = re.fullmatch(r"GSE(\d+)", accession.upper())
    if not match:
        raise ValueError(f"Not a valid GSE accession: {accession}")
    digits = match.group(1)
    if len(digits) <= 3:
        return "GSEnnn"
    return f"GSE{digits[:-3]}nnn"


def gse_root(accession: str) -> str:
    accession = accession.upper()
    return f"{GEO_BASE}/{gse_group(accession)}/{accession}"


def open_url(url: str, timeout: int = 60, headers: dict[str, str] | None = None) -> urllib.response.addinfourl:
    request_headers = {"User-Agent": "ipf-oligo-ml/0.1"}
    if headers:
        request_headers.update(headers)
    request = urllib.request.Request(url, headers=request_headers)
    return urllib.request.urlopen(request, timeout=timeout)


def list_links(url: str) -> list[str]:
    with open_url(url) as response:
        text = response.read().decode("utf-8", errors="replace")
    parser = LinkParser()
    parser.feed(text)
    links = []
    for link in parser.links:
        if link.startswith("?") or link.startswith("/"):
            continue
        if link in {"../", "./"}:
            continue
        links.append(urllib.parse.urljoin(url, link))
    return links


def remote_size(url: str) -> int | None:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "ipf-oligo-ml/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            length = response.headers.get("Content-Length")
            return int(length) if length else None
    except Exception:
        return None


def download_file(url: str, dest: Path, retries: int = 3) -> tuple[str, int | None, str]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    size = remote_size(url)
    if dest.exists() and size is not None and dest.stat().st_size == size:
        return "already_exists", size, ""
    if dest.exists() and size is None and dest.stat().st_size > 0:
        return "already_exists", dest.stat().st_size, "remote size unavailable"

    temp = dest.with_suffix(dest.suffix + ".part")
    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            resume_from = temp.stat().st_size if temp.exists() else 0
            headers = {"Range": f"bytes={resume_from}-"} if resume_from > 0 else None
            with open_url(url, timeout=120, headers=headers) as response:
                mode = "ab" if resume_from > 0 and getattr(response, "status", None) == 206 else "wb"
                if mode == "wb" and resume_from > 0:
                    resume_from = 0
                with temp.open(mode) as handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
            if size is not None and temp.stat().st_size != size:
                if temp.stat().st_size < size:
                    last_error = f"incomplete download: {temp.stat().st_size} of {size} bytes"
                    time.sleep(2 * attempt)
                    continue
                return "failed", size, f"downloaded size {temp.stat().st_size} exceeds expected {size}"
            temp.replace(dest)
            final_size = dest.stat().st_size
            return "downloaded", final_size, ""
        except Exception as exc:
            last_error = f"attempt {attempt}/{retries}: {exc}"
            time.sleep(2 * attempt)
    return "failed", size, last_error


def query_soft_url(accession: str) -> str:
    query = urllib.parse.urlencode(
        {
            "acc": accession.upper(),
            "targ": "self",
            "form": "text",
            "view": "full",
        }
    )
    return f"{GEO_QUERY_BASE}?{query}"


def is_default_supplementary_file(filename: str, include_raw: bool) -> bool:
    lower = filename.lower()
    if lower in {"index.html", "filelist.txt"}:
        return False
    if not include_raw and ("_raw.tar" in lower or lower.endswith("_raw.tar")):
        return False
    return True


def discover_files(accession: str, include_raw: bool = False) -> list[tuple[str, str]]:
    root = gse_root(accession)
    wanted: list[tuple[str, str]] = []

    directories = {
        "matrix": f"{root}/matrix/",
        "soft": f"{root}/soft/",
        "supplementary": f"{root}/suppl/",
    }
    for file_type, directory in directories.items():
        try:
            links = list_links(directory)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                continue
            raise
        for link in links:
            filename = Path(urllib.parse.urlparse(link).path).name
            if not filename or filename.endswith("/"):
                continue
            if file_type == "matrix" and not filename.endswith(".gz"):
                continue
            if file_type == "soft" and not filename.endswith(".gz"):
                continue
            if file_type == "supplementary" and not is_default_supplementary_file(filename, include_raw):
                continue
            wanted.append((file_type, link))
    return wanted


def read_accessions(dataset_table: Path) -> list[str]:
    with dataset_table.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle)
        accessions = []
        for row in rows:
            dataset_id = (row.get("dataset_id") or "").strip().upper()
            if re.fullmatch(r"GSE\d+", dataset_id):
                accessions.append(dataset_id)
    return sorted(set(accessions), key=lambda x: int(x[3:]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--dataset-table", type=Path, default=None)
    parser.add_argument("--only", nargs="*", default=None, help="Optional list of GSE accessions to download")
    parser.add_argument("--no-supplementary", action="store_true", help="Skip GEO supplementary processed files")
    parser.add_argument(
        "--include-raw-supplementary",
        action="store_true",
        help="Also download GEO *_RAW.tar supplementary files. These can be very large.",
    )
    parser.add_argument(
        "--no-query-soft-fallback",
        action="store_true",
        help="Disable fallback download from www.ncbi.nlm.nih.gov/geo/query/acc.cgi.",
    )
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    load_dotenv(Path.home() / ".codex" / ".env")

    dataset_table = args.dataset_table or project_root / "metadata" / "dataset_table.csv"
    accessions = [x.upper() for x in args.only] if args.only else read_accessions(dataset_table)
    out_root = project_root / "data_raw" / "GEO"
    log_root = project_root / "logs"
    log_root.mkdir(parents=True, exist_ok=True)

    records: list[DownloadRecord] = []
    for accession in accessions:
        print(f"[{accession}] discovering files", flush=True)
        dataset_dir = out_root / accession
        try:
            files = discover_files(accession, include_raw=args.include_raw_supplementary)
        except Exception as exc:
            records.append(
                DownloadRecord(accession, "discovery", gse_root(accession), "", "failed", message=str(exc))
            )
            print(f"[{accession}] discovery failed: {exc}", flush=True)
            files = []

        if args.no_supplementary:
            files = [(kind, url) for kind, url in files if kind != "supplementary"]

        index_path = dataset_dir / "geo_file_index.json"
        dataset_dir.mkdir(parents=True, exist_ok=True)
        index_path.write_text(
            json.dumps([{"file_type": kind, "url": url, "size_bytes": remote_size(url)} for kind, url in files], indent=2),
            encoding="utf-8",
        )

        if not files:
            if args.no_query_soft_fallback:
                records.append(DownloadRecord(accession, "discovery", gse_root(accession), "", "no_files"))
                print(f"[{accession}] no downloadable GEO processed files found", flush=True)
                continue
            fallback_url = query_soft_url(accession)
            fallback_dest = dataset_dir / "query_soft" / f"{accession}_full.soft.txt"
            print(f"[{accession}] falling back to GEO query SOFT", flush=True)
            status, size, message = download_file(fallback_url, fallback_dest)
            records.append(
                DownloadRecord(
                    accession,
                    "query_soft",
                    fallback_url,
                    str(fallback_dest),
                    status,
                    size_bytes=size,
                    message=message,
                )
            )
            continue

        for file_type, url in files:
            filename = Path(urllib.parse.urlparse(url).path).name
            dest = dataset_dir / file_type / filename
            print(f"[{accession}] downloading {file_type}/{filename}", flush=True)
            status, size, message = download_file(url, dest)
            records.append(
                DownloadRecord(accession, file_type, url, str(dest), status, size_bytes=size, message=message)
            )
            if status == "failed":
                print(f"[{accession}] failed {filename}: {message}", flush=True)

        if not args.no_query_soft_fallback:
            fallback_url = query_soft_url(accession)
            fallback_dest = dataset_dir / "query_soft" / f"{accession}_full.soft.txt"
            print(f"[{accession}] downloading GEO query SOFT fallback copy", flush=True)
            status, size, message = download_file(fallback_url, fallback_dest)
            records.append(
                DownloadRecord(
                    accession,
                    "query_soft",
                    fallback_url,
                    str(fallback_dest),
                    status,
                    size_bytes=size,
                    message=message,
                )
            )

    manifest_csv = log_root / "geo_download_manifest.csv"
    with manifest_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(records[0]).keys()) if records else [])
        if records:
            writer.writeheader()
            for record in records:
                writer.writerow(asdict(record))

    manifest_json = log_root / "geo_download_manifest.json"
    manifest_json.write_text(json.dumps([asdict(r) for r in records], indent=2), encoding="utf-8")
    print(f"Manifest written to {manifest_csv}", flush=True)
    return 0 if all(r.status != "failed" for r in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
