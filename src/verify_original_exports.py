"""Verify recovered exports against the historical paper data without altering either."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRESSURES = [15, 16, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 33]
VARIABLES = ["X", "Y", "Density", "QX", "QY", "T", "U", "V", "Txy", "Mach", "Pressure", "Kundsen"]


def parse_zones(text):
    """Read these two ordered ASCII POINT formats; fail on missing/nonfinite data."""
    header = re.split(r"^\s*ZONE\b", text, maxsplit=1, flags=re.I | re.M)[0]
    variable_match = re.search(r"\bVARIABLES\s*=\s*(.*)", header, re.I | re.S)
    if variable_match is None:
        raise ValueError("Missing variable header")
    variable_text = variable_match[1].strip()
    names = re.findall(r'"([^"]+)"', variable_text)
    if not names:
        names = [name.strip() for name in variable_text.split(",")]
    if names != VARIABLES:
        raise ValueError(f"Unexpected variables: {names}")
    zones = []
    for line in text.splitlines():
        line = line.strip()
        if re.match(r"ZONE\b", line, re.I):
            zones.append({"rows": [], "packing": None})
        if not zones:
            continue
        zone = zones[-1]
        for axis in ("I", "J", "K"):
            match = re.search(rf"\b{axis}\s*=\s*(\d+)", line, re.I)
            if match:
                zone[axis] = int(match[1])
        packing = re.search(r"\bDATAPACKING\s*=\s*(\w+)", line, re.I)
        if packing:
            zone["packing"] = packing[1].upper()
        if re.match(r"^[+-]?(?:\d|\.\d)", line):
            values = [float(x.replace("D", "E").replace("d", "e")) for x in line.split()]
            if len(values) != len(VARIABLES):
                raise ValueError("Expected twelve values per POINT row")
            zone["rows"].append(values)
    if not zones:
        raise ValueError("No zones")
    arrays = []
    for zone in zones:
        if zone["packing"] != "POINT" or "I" not in zone or "J" not in zone:
            raise ValueError("Missing ordered POINT-zone dimensions/packing")
        array = np.asarray(zone["rows"], dtype=np.float64)
        expected = zone["I"] * zone["J"] * zone.get("K", 1)
        if array.shape != (expected, len(VARIABLES)) or not np.isfinite(array).all():
            raise ValueError("Invalid, incomplete or nonfinite zone")
        arrays.append(array)
    return arrays


def verify(root=ROOT):
    manifest = json.loads((root / "data/original_exports/manifest.json").read_text(encoding="utf-8"))
    cases = manifest["cases"]
    if [case["back_pressure_kpa"] for case in cases] != PRESSURES:
        raise ValueError("Expected exactly the fifteen ordered paper cases")
    report = {"historical_data_commit": manifest["historical_data_commit"], "cases": []}
    for case in cases:
        original = (root / case["original_path"]).read_bytes()
        published = (root / case["published_path"]).read_bytes()
        if hashlib.sha256(original).hexdigest() != case["original_sha256"]:
            raise ValueError(f"Original bytes changed: {case['original_path']}")
        # Historical Git data are LF; Windows working checkouts may use CRLF.
        published_lf = published.replace(b"\r\n", b"\n")
        if hashlib.sha256(published_lf).hexdigest() != case["published_lf_sha256"]:
            raise ValueError(f"Published data changed: {case['published_path']}")
        a = parse_zones(original.decode("utf-8-sig"))
        b = parse_zones(published.decode("utf-8-sig"))
        expected_shapes = [(3131, 12), (1271, 12)]
        if [x.shape for x in a] != expected_shapes or [x.shape for x in b] != expected_shapes:
            raise ValueError("Expected main 101x31 and buffer 31x41 zones")
        mismatch_counts = [int(np.count_nonzero(x.astype(np.float32) != y.astype(np.float32)))
                           for x, y in zip(a, b)]
        if any(mismatch_counts):
            raise ValueError(f"Values differ at published SINGLE precision: {case['original_path']}")
        report["cases"].append({"back_pressure_kpa": case["back_pressure_kpa"],
                                "sha256_verified": True,
                                "float32_mismatches_by_zone": mismatch_counts})
    report["verified_cases"] = len(cases)
    report["verified_values"] = len(cases) * (3131 + 1271) * 12
    report["all_match_at_published_precision"] = True
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional new report path (never overwrites)")
    args = parser.parse_args()
    report = verify()
    rendered = json.dumps(report, indent=2) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
