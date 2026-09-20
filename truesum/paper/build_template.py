"""Builds the TOMS submission manuscript on ACM's own Word template.

The template supplies the styles, page geometry, headers and theme; this
script clears its sample body and emits the manuscript using its named
styles: Title_document, Short Title, Authors, Affiliation, Abstract,
CCSDescription, KeyWords, ACMRefHead/ACMRef, Head1/Head2, PostHeadPara/Para,
Image/FigureCaption, TableCaption, Algorithm/AlgorithmCaption, AckHead/
AckPara, ReferenceHead/Bib_entry, and the In-text code character style.

Citations are author-year, which the TOMS author guidelines specify and the
template permits; reference entries follow the ACM Reference Format the
template demonstrates, alphabetically by first author.
"""
import pathlib
import re

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

HERE = pathlib.Path(__file__).parent
IMG = HERE / "img"


def _find(name):
    """The ACM template is not redistributable, so it may sit beside the
    repository rather than inside it."""
    for d in (HERE, HERE.parent):
        if (d / name).exists():
            return d / name
    raise SystemExit(f"{name} not found in {HERE} or {HERE.parent}")


TEMPLATE = _find("acm_submission_template.docx")
OUT = HERE / "truesum-toms.docx"

REPO_URL = "https://github.com/ArchaeaSoftware/parallelprogrammer/tree/main/truesum"

TOKEN = re.compile(r"\{\{([cib]):(.*?)\}\}", re.S)


def emit(par, text, **_):
    """Inline markup: {{c:...}} code, {{i:...}} italic, {{b:...}} bold."""
    pos = 0
    for m in TOKEN.finditer(text):
        if m.start() > pos:
            par.add_run(text[pos:m.start()])
        kind, inner = m.group(1), m.group(2)
        r = par.add_run(inner)
        if kind == "c":
            try:
                r.style = doc.styles["In-text code"]
            except KeyError:
                r.font.name = "Courier New"
        elif kind == "i":
            r.italic = True
        else:
            r.bold = True
        pos = m.end()
    if pos < len(text):
        par.add_run(text[pos:])


def unnumber(p):
    """Turn off a style's automatic numbering for this paragraph."""
    pr = p._p.get_or_add_pPr()
    numpr = OxmlElement("w:numPr")
    for tag, val in (("w:ilvl", "0"), ("w:numId", "0")):
        e = OxmlElement(tag)
        e.set(qn("w:val"), val)
        numpr.append(e)
    pr.append(numpr)
    return p


def styled(style, text=None):
    p = doc.add_paragraph(style=doc.styles[style])
    if text is not None:
        emit(p, text)
    return p


# The body script calls these with LaTeX-era keywords; they are ignored here
# because the template's styles carry the formatting.
IN_ACKS = []


def para(_doc, text, first_indent=0.25, **_):
    # A paragraph directly after a heading takes the template's PostHeadPara,
    # which is the unindented variant. The body marks those with indent 0.
    if IN_ACKS:
        IN_ACKS.clear()
        return styled("AckPara", text)
    return styled("PostHeadPara" if first_indent == 0 else "Para", text)


def heading(_doc, text, size=12, **_):
    # The template's Head1 and Head2 number themselves, so drop the manual
    # section numbers the body carries.
    text = re.sub(r"^\d+(\.\d+)*\.?\s+", "", text)
    if text.strip().upper() == "ACKNOWLEDGMENTS":
        IN_ACKS.append(True)
        return styled("AckHead", text)
    return styled("Head1" if size >= 12 else "Head2", text)


def bullet(_doc, text, **_):
    # Not the template's List Paragraph, which carries automatic numbering.
    p = doc.add_paragraph(style=doc.styles["Para"])
    pf = p.paragraph_format
    pf.left_indent, pf.first_line_indent = Inches(0.45), Inches(-0.2)
    pf.space_after = Pt(3)
    emit(p, "\u2022\u2003" + text)
    return p


def figure(_doc, name, caption, width_in):
    p = styled("Image")
    p.add_run().add_picture(str(IMG / name), width=Inches(min(width_in, 5.9)))
    styled("AlgorithmCaption" if name.startswith("alg") else "FigureCaption",
           caption)


def table(_doc, caption, header, rows, **_):
    styled("TableCaption", caption)
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    for i, h in enumerate(header):
        cell = t.rows[0].cells[i]
        cell.paragraphs[0].style = doc.styles["TableCell"]
        if i:
            cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.RIGHT
        r = cell.paragraphs[0].add_run(h)
        r.bold = True
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].paragraphs[0].style = doc.styles["TableCell"]
            if i:
                cells[i].paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.RIGHT
            emit(cells[i].paragraphs[0], v)
    doc.add_paragraph(style=doc.styles["Para"])
    return t


# --------------------------------------------------------------------------
doc = Document(str(TEMPLATE))

# Clear the template's sample content, keeping styles, geometry and sectPr.
body_el = doc.element.body
for child in list(body_el):
    if not child.tag.endswith("}sectPr"):
        body_el.remove(child)

TITLE = ("Algorithm XXXX: truesum — Exact Elementwise Summation of "
         "Floating-Point Matrices on CPUs and GPUs")
styled("Title_document", TITLE)
styled("Short Title", "truesum: Exact Elementwise Summation of "
                      "Floating-Point Matrices")
styled("Authors", "Nicholas Wilt")
styled("Affiliation", "Archaea Software, LLC, United States, "
                      "nicholas@archaeasoftware.com")

ABSTRACT = (
    "Floating-point addition is not associative, so a sum of "
    "double-precision values depends on the order of evaluation, and "
    "parallel or distributed summations may vary from run to run. truesum "
    "is a C++17 library that sums floating-point matrices elementwise and "
    "exactly. Every entry of an accumulation matrix is a two's complement "
    "fixed-point integer wide enough that no rounding occurs, making the "
    "result independent of submission order, thread count, and device. "
    "Each column shares one exponent, placed at the least significant bit "
    "of any of its values rather than the most significant, and storage "
    "is limb-major, so that one entry's carry chain runs within a vector "
    "lane while rows fill the lanes. Unlike accumulators whose width is a "
    "property of the format, a column is sized from the dynamic range of "
    "its data. A twelve-byte survey per column, performed on input "
    "matrices, records the lowest true unit in the last place and the "
    "highest occupied bit and allows an accumulation matrix to be sized "
    "before any input matrices arrive. Readback is correctly rounded, "
    "optionally with an exact residual, and a mean over a caller-supplied "
    "count is rounded once rather than twice.")

styled("Abstract", ABSTRACT)

styled("CCSDescription",
       "CCS CONCEPTS • Mathematics of computing → Mathematical "
       "software • Mathematics of computing → Arbitrary-precision "
       "arithmetic • Computing methodologies → Parallel algorithms")
styled("KeyWords",
       "Additional Keywords and Phrases: Exact summation, reproducibility, "
       "block floating point, superaccumulator, fixed-point accumulation, "
       "correct rounding, AVX-512, CUDA, limb-major storage")
styled("ACMRefHead", "ACM Reference Format:")
styled("ACMRef",
       "Nicholas Wilt. 2026. " + TITLE + ". ACM Trans. Math. Softw. 0, 0, "
       "Article 0 (2026), 14 pages. "
       "https://doi.org/10.1145/nnnnnnn.nnnnnnn")
styled("AuthNotes",
       "Author's address: Nicholas Wilt, Archaea Software, LLC, United "
       "States; email: nicholas@archaeasoftware.com; ORCID: [to be "
       "supplied]. This work received no external funding. An earlier and "
       "substantially different draft of this material was posted by the "
       "author to a personal newsletter; it has not been published or "
       "submitted elsewhere.")

exec(open(HERE / "_body.py").read())

styled("ReferenceHead", "REFERENCES")
for entry in [
    "Willow Ahrens, James Demmel, and Hong Diep Nguyen. 2020. Algorithms for "
    "efficient reproducible floating point summation. {{i:ACM Transactions on "
    "Mathematical Software}} 46, 3, Article 22 (2020), 49 pages.",

    "Caroline Collange, David Defour, Stef Graillat, and Roman Iakymchuk. "
    "2015. Numerical reproducibility for the parallel reduction on multi- and "
    "many-core architectures. {{i:Parallel Computing}} 49 (2015), 83–97.",

    "James Demmel and Hong Diep Nguyen. 2013. Fast reproducible "
    "floating-point summation. In {{i:Proceedings of the 21st IEEE Symposium "
    "on Computer Arithmetic (ARITH-21)}}. IEEE, Austin, TX, 163–172.",

    "James Demmel and Hong Diep Nguyen. 2015. Parallel reproducible "
    "summation. {{i:IEEE Transactions on Computers}} 64, 7 (2015), "
    "2060–2070.",

    "James Demmel, Willow Ahrens, and Hong Diep Nguyen. 2016. {{i:Efficient "
    "Reproducible Floating Point Summation and BLAS}}. Technical Report "
    "UCB/EECS-2016-121. EECS Department, University of California, Berkeley, "
    "CA.",

    "Torbjörn Granlund and the GMP Development Team. 2023. {{i:GNU MP: "
    "The GNU Multiple Precision Arithmetic Library}} (6.3.0 ed.). Section "
    "3.2, Nomenclature and Types. Retrieved from "
    "https://gmplib.org/manual/Nomenclature-and-Types.",

    "Roman Iakymchuk, Caroline Collange, David Defour, and Stef Graillat. "
    "2015. ExBLAS: Reproducible and accurate BLAS library. In {{i:Proceedings "
    "of the Workshop on Numerical Reproducibility at Exascale (NRE at SC15)}}. "
    "ACM, Austin, TX.",

    "David S. Johnson. 2002. A theoretician's guide to the experimental "
    "analysis of algorithms. In {{i:Data Structures, Near Neighbor Searches, "
    "and Methodology: Fifth and Sixth DIMACS Implementation Challenges}}, "
    "Michael H. Goldwasser, David S. Johnson, and Catherine C. McGeoch "
    "(Eds.). DIMACS Series in Discrete Mathematics and Theoretical Computer "
    "Science, Vol. 59. American Mathematical Society, Providence, RI, "
    "215–250.",

    "Ulrich Kulisch. 2008. {{i:Computer Arithmetic and Validity: Theory, "
    "Implementation, and Applications}}. de Gruyter, Berlin, Germany.",

    "Ulrich Kulisch and Van Snyder. 2011a. The exact dot product as basic "
    "tool for long interval arithmetic. {{i:Computing}} 91, 3 (2011), "
    "307–313.",

    "Ulrich Kulisch and Van Snyder. 2011b. {{i:The Exact Dot Product}}. "
    "Proposal to the IEEE P1788 working group on interval arithmetic. "
    "Retrieved from "
    "https://grouper.ieee.org/groups/1788/email/pdf7OYfgSd9H4.pdf.",

    "Xiaojun Lei, Tongxiang Gu, Xiaowen Xu, and Stef Graillat. 2025. A "
    "general framework for reproducible parallel preconditioned Krylov "
    "methods using three BLAS variants. In {{i:Proceedings of the 9th "
    "International Conference on High Performance Compilation, Computing and "
    "Communications (HP3C 2025)}}. ACM, Jinan, China, 63–68.",

    "Radford M. Neal. 2015. {{i:Fast Exact Summation Using Small and Large "
    "Superaccumulators}}. arXiv:1505.05571. University of Toronto, Toronto, "
    "Canada.",

    "Katsuhisa Ozaki, Takeshi Ogita, Shin'ichi Oishi, and Siegfried M. Rump. "
    "2012. Error-free transformations of matrix multiplication by using fast "
    "routines of matrix multiplication and its applications. {{i:Numerical "
    "Algorithms}} 59, 1 (2012), 95–118.",

    "Nicholas Wilt. 2026. {{i:truesum}}: source repository. " + REPO_URL + ".",
]:
    unnumber(styled("Bib_entry", entry))

doc.save(str(OUT))
n = len(ABSTRACT.split())
print(f"wrote {OUT.name} on the ACM template")
print(f"abstract: {n} words -> {'OK' if 150 <= n <= 200 else 'OUT OF RANGE'}")
