"""Stage 10 local inputs, namespace checks and sequence alignment. No outcome inputs."""
import gzip
import hashlib
import json
import re
import shutil
import string
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
from .paths import get_project_root, load_dataset_config

ROOT = get_project_root()
CFG = load_dataset_config()
ANN = ROOT / "data/annotations"
OUT = ROOT / "outputs/tables/stage10"
PROC = ROOT / "data/processed/stage10"
FIG = ROOT / "outputs/figures/stage10"
BLOOD = ("GSE63061", "GSE63060", "GSE140829")
BRAIN = "GSE125583"
STATES = ("VALIDATED_UNIQUE_LNCRNA", "PROVISIONAL", "AMBIGUOUS", "NON_LNCRNA", "UNMAPPED")

def prepare():
    for folder in (ANN, OUT, PROC, FIG):
        folder.mkdir(parents=True, exist_ok=True)

def annotation(ds):
    return ROOT / CFG["datasets"][ds]["annotation_path"]

def source(key):
    return ROOT / CFG["annotations"][key]["path"]

def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()

def table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, quoting=3,
                       comment="#" if path.suffix == ".txt" else None)

def matrix_path(ds):
    directory = (ROOT / CFG["datasets"][ds]["expression_output"]).parent
    prefix = "raw_counts" if ds == BRAIN else "expression"
    return directory / f"{prefix}_samples_x_features.parquet"

def features(ds):
    schema = pq.ParquetFile(matrix_path(ds)).schema_arrow
    if json.loads(schema.metadata[b"pandas"])["index_columns"] != ["sample_id"]:
        raise ValueError(f"Unexpected sample index for {ds}")
    ids = [name for name in schema.names if name != "sample_id"]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate measured IDs in {ds}")
    return ids

def decode_nuid(probe):
    """Du et al. nuID encoding, including checksum and zero padding."""
    alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits + "_."
    if len(probe) < 2 or any(c not in alphabet for c in probe):
        return ""
    check, *payload = map(alphabet.index, probe)
    residue, pad = divmod(check, 3)
    if residue != sum(payload) % 21 or payload[-1] % (4**pad):
        return ""
    bases = "".join("ACGT"[(n >> shift) & 3] for n in payload for shift in (4, 2, 0))
    return bases[:len(bases)-pad]

def stable(identifier):
    return re.sub(r"\.\d+(?=_PAR_Y$|$)", "", identifier)

def align(sequences, fasta, output):
    """Compile the local scanner, then pipe decompressed transcripts through it."""
    compiler = shutil.which("g++")
    if not compiler:
        raise RuntimeError("Sequence verification requires a C++11 g++ compiler on PATH")
    with tempfile.TemporaryDirectory() as directory:
        tmp = Path(directory)
        exe = tmp / "match.exe"
        subprocess.run([compiler, "-std=c++11", "-O3", str(ROOT/"scripts/stage10_match.cpp"),
                        "-o", str(exe)], check=True)
        query = tmp / "queries.txt"
        query.write_text("\n".join(sorted(sequences)), encoding="ascii")
        with gzip.open(fasta, "rb") as src, (tmp/"hits.tsv").open("wb") as dst:
            with subprocess.Popen([str(exe), str(query)], stdin=subprocess.PIPE, stdout=dst) as proc:
                shutil.copyfileobj(src, proc.stdin, 1024*1024)
                proc.stdin.close()
                if proc.wait():
                    raise RuntimeError("Transcript alignment failed")
        shutil.copyfile(tmp/"hits.tsv", output)

def measured_platform(ds):
    platform = table(annotation(ds))
    platform["sequence"] = (platform.ID.map(decode_nuid) if ds == BLOOD[2] else platform.SEQUENCE.str.upper())
    ids = features(ds)
    return pd.DataFrame({"ID": ids}).merge(platform, how="left", on="ID", indicator=True).fillna("")

def subset(ds, ids, suffix):
    """Column projection in source order, retaining original GSM index and dtype."""
    frame = pd.read_parquet(matrix_path(ds), columns=[*ids,"sample_id"])
    if frame.index.has_duplicates or not frame.index.str.fullmatch("GSM[0-9]+").all():
        raise ValueError(f"Invalid sample identifiers in {ds}")
    frame.columns = frame.columns.astype(str)
    if list(frame.columns) != ids:
        raise ValueError("Projection changed feature order")
    frame.to_parquet(PROC / f"{ds}_{suffix}.parquet")

