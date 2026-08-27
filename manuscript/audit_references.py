#!/usr/bin/env python3
"""Machine-check reference numbering, use, PDF rendering and LaTeX logs.

This is a publication QA check only.  It does not run any scientific analysis.
Semantic claim-to-source scope is reviewed separately in
``V3_SENTENCE_CLAIM_EVIDENCE_AUDIT_CN.md``.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "MANUSCRIPT_DRAFT_EN.md"
PDF = HERE / "NATURE_COMMUNICATIONS_SUBMISSION.pdf"
LOGS = [
    HERE / "NATURE_COMMUNICATIONS_SUBMISSION.log",
    HERE / "NATURE_COMMUNICATIONS_SUPPLEMENTARY_INFORMATION.log",
    HERE / "DESCRIPTION_OF_ADDITIONAL_SUPPLEMENTARY_FILES.log",
]
OUTPUT = HERE / "logs" / "v3_reference_integrity_audit.json"


def expand_cluster(value: str) -> list[int]:
    numbers: list[int] = []
    for raw in value.replace("–", "-").replace("—", "-").split(","):
        part = raw.strip()
        if not part:
            continue
        if "-" in part:
            start, end = (int(item) for item in part.split("-", 1))
            numbers.extend(range(start, end + 1))
        else:
            numbers.append(int(part))
    return numbers


text = SOURCE.read_text(encoding="utf-8")
body, reference_block = text.split("## References", 1)
reference_numbers = [int(item) for item in re.findall(r"(?m)^(\d+)\.\s", reference_block)]
clusters = re.findall(r"<sup>(.*?)</sup>", body)
cited_numbers = sorted({number for cluster in clusters for number in expand_cluster(cluster)})

pdf_text = subprocess.run(
    ["pdftotext", "-layout", str(PDF), "-"],
    check=True,
    text=True,
    capture_output=True,
).stdout
rendered_reference_numbers = sorted(
    {int(item) for item in re.findall(r"(?m)^\s*\[(\d+)\]\s", pdf_text)}
)

warning_pattern = re.compile(
    r"Overfull|Underfull|Undefined|undefined|LaTeX Warning|Citation.*undefined"
)
log_warnings = {
    path.name: warning_pattern.findall(path.read_text(encoding="utf-8", errors="replace"))
    for path in LOGS
}

expected = list(range(1, max(reference_numbers) + 1))
report = {
    "reference_count": len(reference_numbers),
    "citation_cluster_count": len(clusters),
    "references_consecutive": reference_numbers == expected,
    "unused_references": sorted(set(reference_numbers) - set(cited_numbers)),
    "unknown_citations": sorted(set(cited_numbers) - set(reference_numbers)),
    "pdf_rendered_references_complete": rendered_reference_numbers == reference_numbers,
    "pdf_unresolved_marker_count": pdf_text.count("[?]"),
    "paper1_public_citation_placeholders": body.count(
        "Paper 1 preprint citation to be inserted after public release"
    ),
    "latex_log_warning_matches": log_warnings,
}
report["pass"] = all(
    [
        report["references_consecutive"],
        not report["unused_references"],
        not report["unknown_citations"],
        report["pdf_rendered_references_complete"],
        report["pdf_unresolved_marker_count"] == 0,
        not any(log_warnings.values()),
    ]
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2, ensure_ascii=False))
raise SystemExit(0 if report["pass"] else 1)
