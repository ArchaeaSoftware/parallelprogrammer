"""Renders paper/body.md through the document builders in build_template.py.

The Markdown is the source of the manuscript's prose. This module only turns it
back into the calls the builders already expect, so the Word output is the same
as it was when the body lived in Python.

Inline:  `code` -> {{c:}},  **bold** -> {{b:}},  *italic* -> {{i:}}
Blocks:  # H1, ## H2, paragraphs, "- " bullets,
         :::figure <png> <width-in> / caption / :::
         :::table / caption / GFM table / :::
         :::equation / one line / :::

A paragraph carries the template's indented style only when it follows another
paragraph; after a heading, figure, table, bullet run or equation it takes the
unindented variant. That was the rule the Python body followed by hand, and it
is inferred here rather than annotated.
"""
import re

_CODE = re.compile(r'`([^`]+)`')
_BOLD = re.compile(r'\*\*([^*]+)\*\*')
_ITAL = re.compile(r'(?<!\*)\*([^*]+)\*(?!\*)')


def remark(s):
    """Markdown emphasis back into the {{x:...}} tokens emit() understands."""
    s = _CODE.sub(r'{{c:\1}}', s)
    s = _BOLD.sub(r'{{b:\1}}', s)
    s = _ITAL.sub(r'{{i:\1}}', s)
    return s


def _blocks(text):
    """Split into blocks, keeping ::: fences whole."""
    out, buf, fence = [], [], None
    for line in text.split('\n'):
        if fence is None and line.startswith(':::'):
            if buf:
                out.append('\n'.join(buf)); buf = []
            fence = [line]
        elif fence is not None:
            if line.strip() == ':::':
                out.append('\n'.join(fence)); fence = None
            else:
                fence.append(line)
        elif not line.strip():
            if buf:
                out.append('\n'.join(buf)); buf = []
        else:
            buf.append(line)
    if buf:
        out.append('\n'.join(buf))
    if fence is not None:
        raise SystemExit('body.md: unclosed ::: fence')
    return out


def render(doc, text, *, heading, para, bullet, figure, table, equation):
    after_para = False
    blocks = [b for b in _blocks(text) if not b.lstrip().startswith('<!--')]
    for i, b in enumerate(blocks):
        first = b.split('\n', 1)[0]
        nxt = blocks[i + 1] if i + 1 < len(blocks) else ''

        if first.startswith('## '):
            heading(doc, remark(first[3:].strip()), size=11, before=12)
            after_para = False
        elif first.startswith('# '):
            heading(doc, remark(first[2:].strip()))
            after_para = False
        elif first.startswith(':::figure'):
            _, name, width = first.split()
            figure(doc, name, remark(' '.join(b.split('\n')[1:])), float(width))
            after_para = False
        elif first.startswith(':::equation'):
            equation(doc, remark(' '.join(b.split('\n')[1:])))
            after_para = False
        elif first.startswith(':::table'):
            lines = [l for l in b.split('\n')[1:] if l.strip()]
            rows = [l for l in lines if l.lstrip().startswith('|')]
            cap = ' '.join(l for l in lines if not l.lstrip().startswith('|'))
            cells = lambda r: [remark(c.strip()) for c in r.strip().strip('|').split('|')]
            table(doc, remark(cap), cells(rows[0]), [cells(r) for r in rows[2:]])
            after_para = False
        elif first.lstrip().startswith('- '):
            for item in re.split(r'\n(?=- )', b):
                bullet(doc, remark(' '.join(item.split())[2:]))
            after_para = False
        else:
            # space_after is tightened when an equation follows, as the Python
            # body did for the one equation it introduced.
            extra = {'space_after': 2} if nxt.startswith(':::equation') else {}
            para(doc, remark(' '.join(b.split())),
                 first_indent=0.25 if after_para else 0, **extra)
            after_para = True
