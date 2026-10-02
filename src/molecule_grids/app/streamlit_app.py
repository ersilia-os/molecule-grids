"""Molecule grids app: paste SMILES, get a slide- or print-ready figure.

Run with `molecule-grids app` or `streamlit run src/molecule_grids/app/streamlit_app.py`.
"""

import base64
import io
import random
import re
import zipfile
from pathlib import Path

import streamlit as st

from molecule_grids import FORMATS, draw_grid
from molecule_grids.data import examples
from molecule_grids.io import export
from molecule_grids.io.parse import parse_text

ASSETS = Path(__file__).parent / "assets"
FAVICON = str(ASSETS / "favicon.png")
WORDMARK = str(ASSETS / "ersilia_brand.png")

EXAMPLE = examples.as_text(examples.DEFAULT)


@st.cache_data(show_spinner="Drawing…")
def figure(smiles, names, fmt, style, width, size, number, group, frame, squeeze):
    """Draw the grid; returns (pages, summary)."""
    grid = draw_grid(
        smiles,
        names,
        format=fmt,
        style=style,
        width=width,
        size=size,
        number=number,
        group=group,
        frame=frame,
        squeeze=squeeze,
    )
    summary = {
        "columns": grid.columns,
        "capacity": grid.capacity,
        "bond": grid.bond_pt,
        "width": grid.width_mm,
        "height": grid.height_mm,
    }
    return grid.pages, summary


@st.cache_data(show_spinner="Packing pages…")
def all_pages_zip(pages, stem, dpi):
    """Every page as numbered SVG and PNG files in one ZIP archive."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for k, svg in enumerate(pages, 1):
            z.writestr(f"{stem}_{k}.svg", svg)
            z.writestr(f"{stem}_{k}.png", export.to_png(svg, dpi))
    return buf.getvalue()


def get_examples():
    st.session_state.molecules = examples.as_text(examples.shuffled())


def shuffle_order():
    lines = [ln for ln in st.session_state.molecules.splitlines() if ln.strip()]
    random.shuffle(lines)
    st.session_state.molecules = "\n".join(lines) + "\n"


# ---------------------------------------------------------------- page
# Theme lives in .streamlit/config.toml (Ersilia house style, as in rafiki-workshop-2026).
# The rules below cover what the theme cannot reach; keyed hooks (st-key-*) are stable,
# data-testid selectors are the fragile part - check them after a Streamlit upgrade.
st.set_page_config(page_title="Molecule grids", page_icon=FAVICON, layout="wide")
st.logo(WORDMARK, size="large", link="https://ersilia.io")
st.html("""
<style>
  /* No page title, so the card starts level with the wordmark, not under a 6rem pad. */
  [data-testid="stMainBlockContainer"] { padding: 2rem 2rem !important; }
  /* Sidebar collapsed: no logo in the header, just clear the toggle. */
  [data-testid="stHeader"] [data-testid="stLogoLink"] { display: none !important; }
  .stApp:has([data-testid="stSidebar"][aria-expanded="false"]) [data-testid="stMainBlockContainer"] {
      padding-top: 3.5rem !important; }

  /* Sidebar toggles: a hairline chevron drawn in CSS, instead of the heavy icon-font glyph
     (Streamlit's bundled Material Symbols has a single, fixed weight). */
  [data-testid="stSidebarCollapseButton"] button, [data-testid="stExpandSidebarButton"] {
      width: 26px; height: 26px; border-radius: 0.5rem; color: #9A93A6;
      display: flex; align-items: center; justify-content: center; }
  [data-testid="stSidebarCollapseButton"] button:hover, [data-testid="stExpandSidebarButton"]:hover {
      background: rgba(108, 92, 231, 0.07); color: #6C5CE7; }
  [data-testid="stSidebarCollapseButton"] button > span, [data-testid="stExpandSidebarButton"] > span {
      display: none !important; }
  [data-testid="stSidebarCollapseButton"] button::before, [data-testid="stExpandSidebarButton"]::before {
      content: ""; width: 6px; height: 6px; border-left: 1.25px solid currentColor;
      border-bottom: 1.25px solid currentColor; }
  [data-testid="stSidebarCollapseButton"] button::before { transform: translateX(1.5px) rotate(45deg); }
  [data-testid="stExpandSidebarButton"]::before { transform: translateX(-1.5px) rotate(-135deg); }

  /* Room for "SMILES name" on one line; min-width keeps the sidebar drag-resizable. */
  [data-testid="stSidebar"][aria-expanded="true"] { min-width: 440px !important; }

  /* SMILES are data: mono, one molecule per visual line, scroll sideways instead of wrapping. */
  .st-key-molecules textarea { white-space: pre; overflow-x: auto; font-size: 0.8rem; line-height: 1.6;
                               font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }

  /* The figure card: white surface on the near-white ground, plum-tinted hairline shadow. */
  .st-key-card-figure { background: #FFFFFF; box-shadow: 0 1px 3px rgba(80, 40, 90, .07); }
  .figure { overflow-x: auto; padding: 0.5rem 0 1.25rem; }  /* room under the figure, inside the card */
  .figure img { max-width: 100%; height: auto; display: block; margin: 0 auto; }

  /* Callouts: recessed panel, hairline, the hue only on the left edge (never a saturated fill). */
  [class*="st-key-callout-"] { background: #F4F4F8 !important; border: 1px solid #E6E6EE !important;
                               border-left: 3px solid var(--hue) !important; border-radius: 0.75rem !important;
                               padding: 0.8rem 1rem 0.65rem !important; }
  .st-key-callout-skipped { --hue: #C98A1E; }
  .st-key-callout-empty   { --hue: #6C5CE7; }
  .st-key-callout-pages { --hue: #C98A1E; }
  [class*="st-key-callout-"] [data-testid="stCaptionContainer"] p { color: #2C3E50; font-weight: 500; }
  [class*="st-key-callout-"] [data-testid="stMarkdownContainer"],
  [class*="st-key-callout-"] [data-testid="stCaptionContainer"] { margin-bottom: 0 !important; }
  [class*="st-key-callout-"] p, [class*="st-key-callout-"] ul { margin: 0 !important; }
  [class*="st-key-callout-"] ul { padding-left: 1.25rem !important; }
  [class*="st-key-callout-"] li { margin: 0 !important; font-size: 0.85rem; color: #6B6675; }
  [class*="st-key-callout-"] code { color: #2C3E50; background: #FFFFFF; }

  /* The shuffle action sits under the text box as a quiet text link. */
  .st-key-input-actions { margin-top: -0.5rem; }
  .st-key-input-actions button { min-height: 0; padding: 0.1rem 0; color: #6B6675; }
  .st-key-input-actions button p { font-size: 0.85rem; font-weight: 450; }
  .st-key-input-actions button:hover { color: #6C5CE7; }

  /* Buttons are furniture: quiet weight, colour carries state. */
  .st-key-card-figure button p { font-weight: 450; }
</style>""")

with st.sidebar:
    st.caption("Paste SMILES and get a grid of chemical structures, ready for slides or for a manuscript.")
    st.session_state.setdefault("molecules", EXAMPLE)
    text = st.text_area(
        "Molecules",
        height=240,
        key="molecules",
        help="One molecule per line: SMILES, then an optional name. "
        "CXSMILES extensions such as |&1:1,20| (racemate) are kept. Lines starting with # are ignored.",
    )
    smiles, names, errors = parse_text(text)
    with st.container(horizontal=True, vertical_alignment="center", key="input-actions"):
        st.button(
            "Get examples",
            icon=":material/science:",
            type="tertiary",
            on_click=get_examples,
            help=f"Replace the list with {examples.SHUFFLE_SIZE} random global-health drugs, "
            "including analogue families that share a scaffold.",
        )
        st.button(
            "Shuffle order",
            icon=":material/shuffle:",
            type="tertiary",
            on_click=shuffle_order,
            help="Put the molecules in a random order. Numbering follows the new order.",
        )
        if errors:  # details are in the callout beside the figure
            st.caption(
                f":orange[:material/warning: {len(errors)} line{'s' if len(errors) != 1 else ''} skipped]",
                width="content",
            )

    mode = st.segmented_control(
        "Format",
        ["slide", "print"],
        default="slide",
        required=True,
        width="stretch",
        format_func=lambda m: {"slide": ":material/desktop_windows: Slide", "print": ":material/print: Print"}[m],
        help="The page. Slide: 330 mm wide at full width, up to 186 mm tall (16:9), 20 pt bonds at "
        "Medium, plum frame. Print: 180 mm (two journal columns), up to 247 mm tall, 14.4 pt bonds "
        "(ACS), black frame; insert the SVG at 100 % in Word or Illustrator to keep the dimensions.",
    )
    style = st.segmented_control(
        "Drawing style",
        ["medicinal", "computational"],
        default="medicinal",
        required=True,
        width="stretch",
        format_func=lambda s: {"medicinal": "Medicinal chemist", "computational": "Computational"}[s],
        help="How molecules are drawn. Medicinal chemist: ChemDraw ACS 1996 style, black and white. "
        "Computational: RDKit's atom colours.",
    )
    full_mm = FORMATS[mode].size_in * 25.4
    share = st.slider(
        f"Proportion of full width · {st.session_state.get('share', 100) / 100 * full_mm:.0f} mm",
        25,
        100,
        100,
        5,
        key="share",
        format="%d%%",
        help=f"How much of the full width the figure takes. Full width is {FORMATS['print'].size_in * 25.4:.0f} mm "
        f"in print (two journal columns; half is about one column) and {FORMATS['slide'].size_in * 25.4:.0f} mm "
        "on slides (a wide slide), as in stylia. The width sets how many molecules fit per row; "
        "the height follows, up to one page.",
    )
    width = share / 100
    size = st.segmented_control(
        "Molecule size",
        ["small", "medium", "large"],
        default="medium",
        required=True,
        width="stretch",
        format_func=str.capitalize,
        help="70, 100 or 140 % of the standard bond: 14.4 pt in print (ACS) and 20 pt on slides at Medium. "
        "Bigger molecules mean fewer columns.",
    )
    squeeze = st.segmented_control(
        "Squeeze",
        ["rows", "grid", "free", False],
        default="rows",
        required=True,
        width="stretch",
        format_func=lambda s: {False: "Off", "grid": "Grid", "rows": "Rows", "free": "Free"}[s],
        help="Reorder and turn molecules to fill the figure, at the same bond length and width; "
        "only the height shrinks. Rows: rows flow like text. Grid: aligned columns fitted to their "
        "molecules. Free: molecules interlock by their outlines. Off: square cells in input order. "
        "Numbers follow the new order.",
    )
    group = st.toggle(
        "Group by scaffold",
        value=False,
        disabled=bool(squeeze),
        help="Keeps analogues with the same scaffold together. Square cells only: squeezing places "
        "molecules itself. Analogues share one orientation either way.",
    )
    number = st.toggle("Number compounds 1, 2, 3…", help="Replaces names with bold numbers in reading order.")
    frame = st.toggle(
        "Frame", value=True, help="A thin outline around the whole figure: plum on slides, black in print."
    )

    st.space("large")
    st.caption(
        "Brought to you by the [Ersilia Open Source Initiative](https://ersilia.io), "
        "a tech-nonprofit fueling sustainable research in the Global South."
    )

if errors:
    # A callout, not st.warning: recessed panel, hue on the left edge only (as in rafiki-workshop-2026).
    n = len(errors)
    with st.container(key="callout-skipped", gap="small"):
        st.caption(
            f":orange[:material/warning:] Skipped {n} line{'s' if n != 1 else ''}. "
            f"Fix or remove {'them' if n != 1 else 'it'} to include {'those molecules' if n != 1 else 'that molecule'}."
        )
        st.markdown(
            "\n".join(f"- Line {no}: `{line[:60]}` is {reason}." for no, line, reason in errors[:10])
            + (f"\n- …and {n - 10} more." if n > 10 else "")
        )

if not smiles:
    with st.container(key="callout-empty"):
        st.caption(
            ":violet[:material/arrow_back:] Paste SMILES in the sidebar, one per line, "
            "each optionally followed by a name."
        )
    st.stop()

pages, info = figure(tuple(smiles), tuple(names), mode, style, width, size, number, group, frame, squeeze)
n_pages = len(pages)
if n_pages > 1:
    with st.container(key="callout-pages"):
        st.caption(
            f":orange[:material/warning:] {len(smiles)} molecules need {n_pages} pages "
            f"({info['capacity']} fit on one). Pick a page below, download them all, or try a smaller "
            "size or another squeeze."
        )

stem = f"molecules_{mode}"
dpi = FORMATS[mode].png_dpi

with st.container(border=True, key="card-figure"):
    with st.container(horizontal=True, vertical_alignment="center"):
        page = 1
        if n_pages > 1:
            page = st.number_input(
                "Page",
                1,
                n_pages,
                1,
                key=f"page-{n_pages}",
                width=110,
                label_visibility="collapsed",
                help=f"Page 1 to {n_pages}.",
            )
        svg = pages[page - 1]
        height_mm = float(re.search(r"height='([\d.]+)pt'", svg).group(1)) * 25.4 / 72
        st.caption(
            (f"Page {page} of {n_pages} · " if n_pages > 1 else "")
            + f"{len(smiles)} structure{'s' if len(smiles) != 1 else ''} · "
            f"{mode.capitalize()} · {'medicinal chemist' if style == 'medicinal' else 'computational'} · "
            f"{info['width']:.0f} × {height_mm:.0f} mm · bond {info['bond']:.1f} pt"
        )
        name = f"{stem}_{page}" if n_pages > 1 else stem
        st.download_button(
            "SVG",
            svg,
            file_name=f"{name}.svg",
            mime="image/svg+xml",
            icon=":material/download:",
            help="Vector, for Word, Illustrator or slides." + (" This page." if n_pages > 1 else ""),
        )
        st.download_button(
            "PNG",
            export.to_png(svg, dpi),
            file_name=f"{name}.png",
            mime="image/png",
            icon=":material/download:",
            help=f"{dpi} dpi raster." + (" This page." if n_pages > 1 else ""),
        )
        if n_pages > 1:
            st.download_button(
                "All pages",
                all_pages_zip(tuple(pages), stem, dpi),
                file_name=f"{stem}.zip",
                mime="application/zip",
                icon=":material/folder_zip:",
                help="Every page as SVG and PNG, numbered in reading order.",
            )
    b64 = base64.b64encode(svg.encode()).decode()
    st.markdown(
        f"<div class='figure'><img src='data:image/svg+xml;base64,{b64}' "
        f"alt='Grid of chemical structures, page {page} of {n_pages}'/></div>",
        unsafe_allow_html=True,
    )
