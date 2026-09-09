#!/usr/bin/env python3
"""Convert paper/k-casimir.md to LaTeX and compile a PDF."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MD = ROOT / "k-casimir.md"
TEX = ROOT / "k-casimir.tex"
PDF = ROOT / "k-casimir.pdf"

PREAMBLE = r"""
\documentclass[11pt,letterpaper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{textcomp}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{tabularx}
\usepackage{graphicx}
\usepackage{caption}
\usepackage[margin=1in]{geometry}
\usepackage{setspace}
\usepackage{microtype}
\usepackage[hidelinks,pdfusetitle]{hyperref}
\graphicspath{{../figures/}}
\allowdisplaybreaks
\setstretch{1.08}
\captionsetup{font=small,labelfont=bf,skip=6pt}
\setlength{\parskip}{0.35em}
\setlength{\parindent}{1.2em}
""".strip()


def escape_text(s: str) -> str:
    s = s.replace("\\", r"\textbackslash{}")
    s = s.replace("&", r"\&")
    s = s.replace("%", r"\%")
    s = s.replace("#", r"\#")
    s = s.replace("_", r"\_")
    s = s.replace("~", r"\textasciitilde{}")
    return s


def _next_unescaped_span(s: str, i: int) -> tuple[str, str, int]:
    """Return (kind, value, new_index) for the next span starting at i."""
    n = len(s)
    if s.startswith("\\(", i):
        j = s.find("\\)", i + 2)
        if j < 0:
            raise ValueError(f"unclosed inline math: {s[i:i+40]!r}")
        return "math", s[i : j + 2], j + 2
    if s.startswith("`", i):
        j = s.find("`", i + 1)
        if j < 0:
            raise ValueError(f"unclosed code: {s[i:i+40]!r}")
        return "code", s[i + 1 : j], j + 1
    if s.startswith("**", i):
        j = i + 2
        while j < n:
            if s.startswith("\\(", j):
                j = s.find("\\)", j + 2)
                if j < 0:
                    raise ValueError("unclosed math inside bold")
                j += 2
                continue
            if s.startswith("**", j):
                inner = convert_inline(s[i + 2 : j])
                return "wrap", r"\textbf{" + inner + "}", j + 2
            j += 1
        return "text", s[i], i + 1
    if s.startswith("*", i):
        j = i + 1
        while j < n:
            if s.startswith("\\(", j):
                j = s.find("\\)", j + 2)
                if j < 0:
                    raise ValueError("unclosed math inside italic")
                j += 2
                continue
            if s.startswith("*", j) and not s.startswith("**", j):
                inner = convert_inline(s[i + 1 : j])
                return "wrap", r"\emph{" + inner + "}", j + 1
            j += 1
        return "text", s[i], i + 1
    nxt = n
    for mark in ("\\(", "`", "**", "*"):
        k = s.find(mark, i)
        if k != -1:
            nxt = min(nxt, k)
    return "text", s[i:nxt], nxt


def convert_inline(s: str) -> str:
    out: list[str] = []
    i, n = 0, len(s)
    while i < n:
        kind, val, i = _next_unescaped_span(s, i)
        if kind == "math":
            out.append(val)
        elif kind == "code":
            out.append(r"\texttt{" + escape_text(val) + "}")
        elif kind == "wrap":
            out.append(val)
        else:
            out.append(escape_text(val))
    return "".join(out)


def split_table_row(line: str) -> list[str]:
    line = line.strip().strip("|")
    return [c.strip() for c in line.split("|")]


def is_table_sep(line: str) -> bool:
    s = line.strip().strip("|").replace(" ", "")
    return bool(s) and all(set(cell) <= set("-:") and "-" in cell for cell in s.split("|"))


def emit_table(header: list[str], rows: list[list[str]]) -> list[str]:
    ncols = len(header)
    if ncols == 3:
        spec = r">{\raggedright\arraybackslash}p{0.34\textwidth}"
        spec += r">{\raggedright\arraybackslash}p{0.34\textwidth}"
        spec += r">{\raggedright\arraybackslash}p{0.20\textwidth}"
    elif ncols == 4:
        spec = r">{\raggedright\arraybackslash}p{0.22\textwidth}" * 4
    else:
        spec = "l" * ncols
    lines = [
        r"\begin{center}",
        r"{\small",
        rf"\begin{{tabular}}{{{spec}}}",
        r"\toprule",
        " & ".join(r"\textbf{" + convert_inline(c) + "}" for c in header) + r" \\",
        r"\midrule",
    ]
    for row in rows:
        row = (row + [""] * ncols)[:ncols]
        lines.append(" & ".join(convert_inline(c) for c in row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", "}", r"\end{center}", ""]
    return lines


def heading_to_tex(level: int, title: str) -> str:
    title = title.strip()
    numbered = re.match(r"^(\d+)\.\s+(.*)$", title)
    if numbered:
        title = numbered.group(2)
    appendix = re.match(r"^Appendix [A-Z]\.\s+(.*)$", title)
    if appendix:
        title = appendix.group(1)
    cmd = {1: "section", 2: "subsection", 3: "subsubsection"}[level]
    star = ""
    if re.match(r"^(Abstract|Key decisions|Open questions|PR plan|References)$", title):
        star = "*"
    converted = convert_inline(title)
    if r"\(" in title:
        plain = re.sub(r"\\\((.+?)\\\)", r"\1", title)
        plain = plain.replace("\\", "")
        converted = rf"\texorpdfstring{{{converted}}}{{{plain}}}"
    return f"\\{cmd}{star}{{{converted}}}"


def convert(md: str) -> str:
    lines = md.splitlines()
    body: list[str] = []
    i = 0
    title = ""
    blurb = ""
    abstract: list[str] = []
    in_abstract = False
    abstract_display = False
    in_display = False
    display_buf: list[str] = []
    in_list = False
    list_kind = ""
    appendix_started = False

    def close_list() -> None:
        nonlocal in_list, list_kind
        if in_list:
            body.append(r"\end{" + list_kind + "}")
            body.append("")
            in_list = False
            list_kind = ""

    while i < len(lines):
        line = lines[i]

        if in_display:
            display_buf.append(line)
            if r"\]" in line:
                body.extend(display_buf)
                body.append("")
                display_buf = []
                in_display = False
            i += 1
            continue

        if line.strip() == "---":
            close_list()
            i += 1
            continue

        if line.startswith("# ") and not title:
            title = line[2:].strip()
            i += 1
            continue

        if not title:
            i += 1
            continue

        if not blurb and line.startswith("**Kasimir theory.**"):
            blurb = convert_inline(line)
            i += 1
            continue

        if line.startswith("## Abstract"):
            in_abstract = True
            i += 1
            continue

        if in_abstract:
            if line.startswith("## "):
                in_abstract = False
                abstract_display = False
                # fall through to handle this heading
            else:
                stripped = line.strip()
                if r"\[" in stripped or abstract_display:
                    if stripped:
                        abstract.append(stripped)
                    if r"\[" in stripped and r"\]" not in stripped:
                        abstract_display = True
                    if r"\]" in stripped:
                        abstract_display = False
                        abstract.append("")
                elif stripped:
                    abstract.append(convert_inline(stripped))
                    abstract.append("")
                i += 1
                continue

        if line.startswith("## "):
            close_list()
            raw_title = line[3:].strip()
            if raw_title.startswith("Appendix ") and not appendix_started:
                body.append(r"\appendix")
                body.append("")
                appendix_started = True
            if raw_title.startswith("References"):
                body.append(r"\begin{sloppypar}")
            body.append(heading_to_tex(1, raw_title))
            body.append("")
            i += 1
            continue

        if line.startswith("### "):
            close_list()
            body.append(heading_to_tex(2, line[4:]))
            body.append("")
            i += 1
            continue

        if line.startswith("!["):
            close_list()
            m = re.match(r"!\[(.*)\]\((.*)\)", line.strip())
            if not m:
                raise ValueError(f"bad image: {line}")
            cap, path = m.group(1), m.group(2)
            fname = Path(path).name
            body += [
                r"\begin{figure}[htbp]",
                r"\centering",
                rf"\includegraphics[width=\textwidth]{{{fname}}}",
                rf"\caption{{{convert_inline(cap)}}}",
                r"\end{figure}",
                "",
            ]
            i += 1
            continue

        if "|" in line and i + 1 < len(lines) and is_table_sep(lines[i + 1]):
            close_list()
            header = split_table_row(line)
            i += 2
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                rows.append(split_table_row(lines[i]))
                i += 1
            body.extend(emit_table(header, rows))
            continue

        bullet = re.match(r"^- (.*)$", line)
        numbered = re.match(r"^(\d+)\.\s+(.*)$", line)
        if bullet or numbered:
            kind = "itemize" if bullet else "enumerate"
            text = bullet.group(1) if bullet else numbered.group(2)
            if not in_list or list_kind != kind:
                close_list()
                body.append(r"\begin{" + kind + "}")
                in_list = True
                list_kind = kind
            body.append(r"\item " + convert_inline(text))
            i += 1
            continue

        if not line.strip():
            close_list()
            body.append("")
            i += 1
            continue

        close_list()
        if r"\[" in line:
            in_display = r"\]" not in line
            if in_display:
                display_buf = [line]
            else:
                body.append(line)
                body.append("")
            i += 1
            continue

        body.append(convert_inline(line))
        body.append("")
        i += 1

    close_list()

    tex = [
        PREAMBLE,
        r"\begin{document}",
        rf"\title{{{convert_inline(title)}}}",
        r"\author{Kasimir}",
        r"\date{}",
        r"\maketitle",
        r"\begin{quote}",
        r"\small " + blurb,
        r"\end{quote}",
        r"\begin{abstract}",
        *abstract,
        r"\end{abstract}",
        "",
        *body,
        r"\end{sloppypar}",
        r"\end{document}",
        "",
    ]
    return "\n".join(tex)


def compile_pdf() -> None:
    md = MD.read_text(encoding="utf-8")
    tex = convert(md)
    TEX.write_text(tex, encoding="utf-8")
    cmd = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", TEX.name]
    for _ in range(2):
        proc = subprocess.run(
            cmd,
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout[-4000:])
            sys.stderr.write(proc.stderr[-2000:])
            raise SystemExit(f"pdflatex failed with code {proc.returncode}")
    print(f"wrote {PDF}")


if __name__ == "__main__":
    compile_pdf()
