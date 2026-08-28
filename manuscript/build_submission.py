#!/usr/bin/env python3
"""Build Paper 2 V3 main and Supplementary Information LaTeX files.

The Markdown files remain the editable sources. The generated submission uses
the standard article class, numerical citations and embedded references. Main
and supplementary figures are read from the already audited ``figures``
directory; this builder performs no scientific analysis.
"""

from __future__ import annotations

import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
MAIN_SOURCE = HERE / "MANUSCRIPT_DRAFT_EN.md"
SI_SOURCE = HERE / "SUPPLEMENTARY_INFORMATION_DRAFT_EN.md"
ADDITIONAL_SOURCE = HERE / "DESCRIPTION_OF_ADDITIONAL_SUPPLEMENTARY_FILES.md"
MAIN_TEX = HERE / "NATURE_COMMUNICATIONS_SUBMISSION.tex"
SI_TEX = HERE / "NATURE_COMMUNICATIONS_SUPPLEMENTARY_INFORMATION.tex"
ADDITIONAL_TEX = HERE / "DESCRIPTION_OF_ADDITIONAL_SUPPLEMENTARY_FILES.tex"


def escape_latex(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
        "–": "--", "—": "---", "−": "$-$", "±": r"$\pm$",
        "×": r"$\times$", "≥": r"$\geq$", "≤": r"$\leq$",
        "é": r"\'e", "ö": r'\"o', "“": "``", "”": "''",
        "‘": "`", "’": "'", "‑": "-",
    }
    return "".join(replacements.get(char, char) for char in value)


def citation_keys(label: str) -> list[str]:
    keys: list[str] = []
    for item in label.replace("–", "-").replace("—", "-").split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            start, end = [int(part) for part in item.split("-", 1)]
            keys.extend(f"ref{number}" for number in range(start, end + 1))
        else:
            keys.append(f"ref{int(item)}")
    return keys


def inline_latex(value: str) -> str:
    tokens: dict[str, str] = {}

    def token(content: str) -> str:
        key = f"ZZTOKEN{len(tokens)}ZZ"
        tokens[key] = content
        return key

    value = re.sub(
        r"<sup>(.*?)</sup>",
        # A citation cluster is an unbreakable superscript box.  Permit a line
        # break immediately before it so long clusters do not protrude beyond
        # the text block while preserving the journal-style superscript form.
        lambda match: token(r"\allowbreak\textsuperscript{\cite{" + ",".join(citation_keys(match.group(1))) + "}}"),
        value,
    )
    value = re.sub(r"\\\((.+?)\\\)", lambda match: token("$" + match.group(1) + "$"), value)
    # ``\path`` retains monospaced semantics while allowing long product and
    # data-file identifiers to break at underscores and punctuation.
    value = re.sub(r"`([^`]+)`", lambda match: token(r"\path{" + match.group(1) + "}"), value)

    def replace_url(match: re.Match[str]) -> str:
        raw = match.group(0)
        trailing = ""
        while raw and raw[-1] in ".,;":
            trailing = raw[-1] + trailing
            raw = raw[:-1]
        return token(r"\url{" + raw + "}") + trailing

    value = re.sub(r"https?://[^\s)>]+", replace_url, value)
    value = re.sub(r"\*\*(.+?)\*\*", lambda match: token(r"\textbf{" + escape_latex(match.group(1)) + "}"), value)
    value = re.sub(r"\*([^*]+?)\*", lambda match: token(r"\textit{" + escape_latex(match.group(1)) + "}"), value)
    rendered = escape_latex(value)
    for key, content in tokens.items():
        rendered = rendered.replace(key, content)
    return rendered


def split_sections(lines: list[str]) -> tuple[str, dict[str, list[str]]]:
    title = lines[0].removeprefix("# ").strip()
    sections: dict[str, list[str]] = {}
    current: str | None = None
    for line in lines[1:]:
        if line.startswith("## "):
            current = line[3:].strip()
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    return title, sections


def render_body(lines: list[str], subsection_command: str = "subsection") -> str:
    output: list[str] = []
    paragraph: list[str] = []
    math: list[str] = []
    in_math = False

    def flush_paragraph() -> None:
        if paragraph:
            output.append(inline_latex(" ".join(item.strip() for item in paragraph)))
            output.append("")
            paragraph.clear()

    for line in lines:
        stripped = line.strip()
        if in_math:
            if stripped == r"\]":
                output.extend([r"\[", "\n".join(math), r"\]", ""])
                math.clear()
                in_math = False
            else:
                math.append(line)
            continue
        if stripped == r"\[":
            flush_paragraph()
            in_math = True
            continue
        if not stripped:
            flush_paragraph()
            continue
        if line.startswith("### "):
            flush_paragraph()
            output.append(rf"\{subsection_command}*{{{inline_latex(line[4:].strip())}}}")
            output.append("")
            continue
        paragraph.append(line)
    flush_paragraph()
    if in_math:
        raise RuntimeError("Unclosed display-math block in Markdown source")
    return "\n".join(output).rstrip()


def parse_legends(lines: list[str], label: str) -> list[tuple[int, str, str]]:
    pattern = re.compile(rf"^### {re.escape(label)} (\d+)\s*\|\s*(.*)$")
    blocks: list[tuple[int, str, str]] = []
    number: int | None = None
    title = ""
    body: list[str] = []

    def flush() -> None:
        nonlocal number, title, body
        if number is not None:
            blocks.append((number, title, " ".join(v.strip() for v in body if v.strip())))
        number, title, body = None, "", []

    for line in lines:
        match = pattern.match(line)
        if match:
            flush()
            number, title = int(match.group(1)), match.group(2).strip()
        elif number is not None:
            body.append(line)
    flush()
    return blocks


def render_references(lines: list[str]) -> str:
    output = [r"\begin{thebibliography}{99}"]
    for line in lines:
        match = re.match(r"^(\d+)\.\s+(.*)$", line)
        if match:
            output.extend([rf"\bibitem{{ref{match.group(1)}}} {inline_latex(match.group(2))}", ""])
    output.append(r"\end{thebibliography}")
    return "\n".join(output)


def render_results(lines: list[str], legend_blocks: list[tuple[int, str, str]]) -> str:
    figures = [
        "fig1_last_passage_decomposition.pdf",
        "fig2_cross_system_context.pdf",
        "fig3_dual_reference.pdf",
        "fig4_boundary_counterfactual.pdf",
    ]
    legends = {number: (title, body) for number, title, body in legend_blocks}
    blocks: list[tuple[str, list[str]]] = []
    prelude: list[str] = []
    heading: str | None = None
    body: list[str] = []
    for line in lines:
        if line.startswith("### "):
            if heading is not None:
                blocks.append((heading, body))
            heading, body = line[4:].strip(), []
        elif heading is not None:
            body.append(line)
        else:
            prelude.append(line)
    if heading is not None:
        blocks.append((heading, body))
    if len(blocks) != 4 or sorted(legends) != [1, 2, 3, 4]:
        raise RuntimeError("Main Results must contain four subsections and four legends")
    output: list[str] = [render_body(prelude), ""] if any(line.strip() for line in prelude) else []
    for number, ((heading, lines_for_block), figure) in enumerate(zip(blocks, figures), start=1):
        output.extend([r"\subsection*{" + inline_latex(heading) + "}", "", render_body(lines_for_block), ""])
        title, caption = legends[number]
        output.extend([
            r"\clearpage", r"\begin{center}",
            rf"\includegraphics[width=\textwidth,height=0.84\textheight,keepaspectratio]{{figures/{figure}}}",
            r"\end{center}",
            rf"\noindent\textbf{{Figure {number} | {inline_latex(title)}}} {inline_latex(caption)}\par", "",
        ])
    return "\n".join(output)


def render_additional_body(lines: list[str]) -> str:
    """Render the simple headings, paragraphs and bullet lists in the file description."""
    output: list[str] = []
    paragraph: list[str] = []
    bullets: list[str] = []

    def flush_paragraph() -> None:
        if paragraph:
            output.extend([inline_latex(" ".join(item.strip() for item in paragraph)), ""])
            paragraph.clear()

    def flush_bullets() -> None:
        if bullets:
            output.append(r"\begin{itemize}")
            output.extend(r"\item " + inline_latex(item) for item in bullets)
            output.extend([r"\end{itemize}", ""])
            bullets.clear()

    for line in lines:
        stripped = line.strip()
        if not stripped:
            flush_paragraph()
            flush_bullets()
        elif line.startswith("### "):
            flush_paragraph()
            flush_bullets()
            output.extend([r"\subsection*{" + inline_latex(line[4:].strip()) + "}", ""])
        elif line.startswith("## "):
            flush_paragraph()
            flush_bullets()
            output.extend([r"\section*{" + inline_latex(line[3:].strip()) + "}", ""])
        elif line.startswith("- "):
            flush_paragraph()
            bullets.append(line[2:].strip())
        else:
            flush_bullets()
            paragraph.append(line)
    flush_paragraph()
    flush_bullets()
    return "\n".join(output).rstrip()


def build_main() -> None:
    title, sections = split_sections(MAIN_SOURCE.read_text(encoding="utf-8").splitlines())
    required = ["Abstract", "Introduction", "Results", "Discussion", "Methods",
                "Data availability", "Code availability", "Figure legends", "References"]
    missing = [item for item in required if item not in sections]
    if missing:
        raise RuntimeError(f"Missing main sections: {missing}")
    legends = parse_legends(sections["Figure legends"], "Figure")
    preamble = rf"""\documentclass[12pt]{{article}}
\usepackage[utf8]{{inputenc}}
\usepackage[T1]{{fontenc}}
\usepackage{{lmodern}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{cite}}
\usepackage{{url}}
\usepackage[margin=1in]{{geometry}}
\emergencystretch=3em
\hbadness=10000
\urlstyle{{same}}
\title{{{inline_latex(title)}}}
\author{{Haobo Liu$^{{1}}$ and Ke Li$^{{1,*}}$\\[0.5em]
\small $^{{1}}$School of Aeronautic Science and Engineering, Beihang University,\\
\small 37 Xueyuan Road, Beijing 100191, China\\[0.35em]
\small $^{{*}}$Correspondence: like@buaa.edu.cn\\[0.2em]
\small ORCID: Haobo Liu, 0000-0001-9231-8961; Ke Li, 0000-0002-3694-1772}}
\date{{}}
\begin{{document}}
\maketitle
\begin{{abstract}}
{render_body(sections['Abstract'])}
\end{{abstract}}
"""
    output = [preamble]
    for section in ["Introduction", "Results", "Discussion", "Methods",
                    "Data availability", "Code availability"]:
        output.append(r"\section*{" + inline_latex(section) + "}")
        output.append(render_results(sections[section], legends) if section == "Results" else render_body(sections[section]))
        output.append("")
    output.extend([
        r"\clearpage", render_references(sections["References"]), "",
        r"\section*{Acknowledgements}",
        "This work was supported by the National Natural Science Foundation of China (NSFC; grant 61773039). We thank the original investigators and repositories for making the tracking, behavioural and environmental datasets available.", "",
        r"\section*{Author contributions}",
        "Conceptualization, H.L.; methodology, H.L.; software, H.L.; validation, H.L.; formal analysis, H.L.; investigation, H.L.; data curation, H.L.; writing---original draft, H.L.; writing---review and editing, H.L.; visualization, H.L.; supervision, K.L.; funding acquisition, K.L. Both authors read and approved the final manuscript.", "",
        r"\section*{Competing interests}",
        "The authors declare no competing interests.", "",
        r"\section*{Use of artificial intelligence-assisted technologies}",
        "A large language model was used for language drafting and manuscript organization under author supervision. It did not generate raw data or make final scientific decisions. The authors verified all numerical claims against archived analysis tables and take full responsibility for the manuscript.",
        r"\end{document}",
    ])
    MAIN_TEX.write_text("\n".join(output) + "\n", encoding="utf-8")


def split_si(lines: list[str]) -> tuple[str, str, list[str], dict[str, list[str]]]:
    document_title = lines[0].removeprefix("# ").strip()
    index = next(i for i, line in enumerate(lines[1:], start=1) if line.startswith("## "))
    article_title = lines[index].removeprefix("## ").strip()
    intro: list[str] = []
    sections: dict[str, list[str]] = {}
    current: str | None = None
    for line in lines[index + 1:]:
        if line.startswith("## "):
            current = line[3:].strip()
            sections[current] = []
        elif current is None:
            intro.append(line)
        else:
            sections[current].append(line)
    return document_title, article_title, intro, sections


def markdown_table(lines: list[str]) -> tuple[list[str], list[list[str]], list[str]]:
    table_indices = [i for i, line in enumerate(lines) if line.lstrip().startswith("|")]
    if len(table_indices) < 2:
        raise RuntimeError("Supplementary table block has no Markdown table")
    start = table_indices[0]
    table_lines = [lines[i] for i in table_indices]
    cells = lambda line: [value.strip() for value in line.strip().strip("|").split("|")]
    header = cells(table_lines[0])
    rows = [cells(line) for line in table_lines[2:]]
    if any(len(row) != len(header) for row in rows):
        raise RuntimeError("Inconsistent Supplementary table width")
    return header, rows, lines[:start]


def render_si_tables(lines: list[str]) -> str:
    pattern = re.compile(r"^### Supplementary Table (\d+)\s*\|\s*(.*)$")
    blocks: list[tuple[int, str, list[str]]] = []
    number: int | None = None
    title = ""
    body: list[str] = []
    for line in lines + ["### Supplementary Table 999 | END"]:
        match = pattern.match(line)
        if match:
            if number is not None:
                blocks.append((number, title, body))
            number, title, body = int(match.group(1)), match.group(2).strip(), []
        elif number is not None:
            body.append(line)
    blocks = [block for block in blocks if block[0] != 999]
    if [block[0] for block in blocks] != list(range(1, 8)):
        raise RuntimeError("Supplementary Tables must be numbered 1--7")
    output: list[str] = []
    for number, title, body in blocks:
        header, rows, intro = markdown_table(body)
        width = 0.92 / len(header)
        spec = "@{}" + "".join(r">{\raggedright\arraybackslash}p{" + f"{width:.3f}" + r"\linewidth}" for _ in header) + "@{}"
        output.extend([
            r"\clearpage", r"\begin{landscape}",
            rf"\section*{{Supplementary Table {number} | {inline_latex(title)}}}",
            render_body(intro), r"\begingroup", r"\scriptsize",
            r"\setlength{\tabcolsep}{2pt}", r"\renewcommand{\arraystretch}{1.15}",
            r"\begin{longtable}{" + spec + "}", r"\toprule",
            " & ".join(r"\textbf{" + inline_latex(cell) + "}" for cell in header) + r" \\",
            r"\midrule", r"\endfirsthead", r"\toprule",
            " & ".join(r"\textbf{" + inline_latex(cell) + "}" for cell in header) + r" \\",
            r"\midrule", r"\endhead",
        ])
        output.extend(" & ".join(inline_latex(cell) for cell in row) + r" \\" for row in rows)
        output.extend([r"\bottomrule", r"\end{longtable}", r"\endgroup", r"\end{landscape}", ""])
    return "\n".join(output)


def build_si() -> None:
    document_title, article_title, intro, sections = split_si(SI_SOURCE.read_text(encoding="utf-8").splitlines())
    required = ["Supplementary Methods", "Supplementary Results", "Supplementary Figures",
                "Supplementary Tables", "Supplementary Data manifest"]
    missing = [item for item in required if item not in sections]
    if missing:
        raise RuntimeError(f"Missing SI sections: {missing}")
    legends = parse_legends(sections["Supplementary Figures"], "Supplementary Figure")
    if [item[0] for item in legends] != [1, 2, 3, 4, 5]:
        raise RuntimeError("Supplementary Figures must be numbered 1--5")
    figures = [
        "figS1_observation_support.pdf", "figS2_complete_last_passage.pdf",
        "figS3_laysan_conditional.pdf", "figS4_behaviour_context.pdf",
        "figS5_boundary_counterfactual.pdf",
    ]
    preamble = rf"""\documentclass[11pt]{{article}}
\usepackage[utf8]{{inputenc}}
\usepackage[T1]{{fontenc}}
\usepackage{{lmodern}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{url}}
\usepackage{{array}}
\usepackage{{booktabs}}
\usepackage{{longtable}}
\usepackage{{pdflscape}}
\usepackage[margin=0.8in]{{geometry}}
\emergencystretch=3em
\hbadness=10000
\urlstyle{{same}}
\title{{{inline_latex(document_title)}\\[0.6em]\large {inline_latex(article_title)}}}
\author{{Haobo Liu$^{{1}}$ and Ke Li$^{{1,*}}$\\[0.5em]
\small $^{{1}}$School of Aeronautic Science and Engineering, Beihang University,\\
\small 37 Xueyuan Road, Beijing 100191, China\\[0.35em]
\small $^{{*}}$Correspondence: like@buaa.edu.cn\\[0.2em]
\small ORCID: Haobo Liu, 0000-0001-9231-8961; Ke Li, 0000-0002-3694-1772}}
\date{{}}
\begin{{document}}
\maketitle
{render_body(intro)}
"""
    output = [preamble]
    for section in ["Supplementary Methods", "Supplementary Results"]:
        output.extend([r"\section*{" + inline_latex(section) + "}", render_body(sections[section]), ""])
    output.append(r"\section*{Supplementary Figures}")
    for (number, title, caption), figure in zip(legends, figures):
        output.extend([
            r"\clearpage", r"\begin{center}",
            rf"\includegraphics[width=\textwidth,height=0.78\textheight,keepaspectratio]{{figures/{figure}}}",
            r"\end{center}",
            rf"\noindent\textbf{{Supplementary Figure {number} | {inline_latex(title)}}} {inline_latex(caption)}\par", "",
        ])
    output.extend([render_si_tables(sections["Supplementary Tables"]), r"\clearpage",
                   r"\section*{Supplementary Data manifest}",
                   render_body(sections["Supplementary Data manifest"]), r"\end{document}"])
    SI_TEX.write_text("\n".join(output) + "\n", encoding="utf-8")


def build_additional_description() -> None:
    lines = ADDITIONAL_SOURCE.read_text(encoding="utf-8").splitlines()
    title = lines[0].removeprefix("# ").strip()
    preamble = rf"""\documentclass[11pt]{{article}}
\usepackage[utf8]{{inputenc}}
\usepackage[T1]{{fontenc}}
\usepackage{{lmodern}}
\usepackage{{url}}
\usepackage[margin=0.8in]{{geometry}}
\emergencystretch=3em
\hbadness=10000
\urlstyle{{same}}
\title{{{inline_latex(title)}}}
\author{{}}
\date{{}}
\begin{{document}}
\maketitle
"""
    ADDITIONAL_TEX.write_text(
        preamble + render_additional_body(lines[1:]) + "\n\\end{document}\n",
        encoding="utf-8",
    )


def main() -> int:
    build_main()
    build_si()
    build_additional_description()
    print(f"wrote {MAIN_TEX}")
    print(f"wrote {SI_TEX}")
    print(f"wrote {ADDITIONAL_TEX}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
