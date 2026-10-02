"""Molecule grids app: paste SMILES, get a slide- or print-ready figure.

Run with `molecule-grids app` or `streamlit run src/molecule_grids/app/streamlit_app.py`.
"""

import base64
import random
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
def figure(smiles, names, fmt, width, zoom, number, group, frame, squeeze):
    """Draw the grid; returns (pages, summary)."""
    grid = draw_grid(
        smiles,
        names,
        format=fmt,
        width=width,
        zoom=zoom,
        number=number,
        group=group,
        frame=frame,
        squeeze=squeeze,
    )
    summary = {
        "columns": grid.columns,
        "capacity": grid.capacity,
        "bond": grid.bond_pt,
        "width": grid.width_in,
        "height": grid.height_in,
    }
    return grid.pages, summary


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
  .figure { overflow-x: auto; padding: 0.5rem 0 0.25rem; }
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
        help="Slide: RDKit colours, 13 in wide at full width. "
        "Print: ChemDraw ACS 1996 style, 7.09 in (two journal columns) at full width, 14.4 pt bonds "
        "at 100 %; insert the SVG at 100 % in Word or Illustrator to keep the dimensions.",
    )
    width = st.slider(
        "Width",
        0.25,
        1.0,
        1.0,
        0.05,
        format="%.2f",
        help=f"Fraction of the full width: {FORMATS['slide'].size_in:g} in for slides, "
        f"{FORMATS['print'].size_in:g} in for print (as in stylia). "
        "As many columns as fit; each molecule sits in a square cell, so the height follows.",
    )
    zoom = st.slider(
        "Magnifier",
        50,
        200,
        100,
        10,
        format="%d%%",
        help="Molecule size. 100 % is a 14.4 pt bond in print (ACS) and 20 pt on slides. "
        "Bigger molecules mean fewer columns.",
    )
    number = st.toggle("Number compounds 1, 2, 3…", help="Replaces names with bold numbers in input order.")
    group = st.toggle(
        "Group by scaffold",
        value=False,
        help="Keeps analogues with the same scaffold together, in the same orientation.",
    )
    squeeze = st.segmented_control(
        "Squeeze",
        [False, "grid", "rows", "free"],
        default=False,
        required=True,
        width="stretch",
        format_func=lambda s: {False: "Off", "grid": "Grid", "rows": "Rows", "free": "Free"}[s],
        help="Reorder and turn molecules to fill the figure, at the same bond length; the figure "
        "gets smaller. Grid: aligned columns fitted to their molecules. Rows: rows flow like text. "
        "Free: molecules interlock by their outlines. Numbers follow the new order.",
    )
    frame = st.toggle("Frame", value=False, help="A thin black outline around the whole figure.")

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

pages, info = figure(tuple(smiles), tuple(names), mode, width, zoom / 100, number, group, frame, squeeze)
if len(pages) > 1:
    with st.container(key="callout-pages"):
        st.caption(
            f":orange[:material/warning:] Showing the first {info['capacity']} of {len(smiles)} molecules. "
            f"Download the GIF for all of them ({len(pages)} frames)."
        )

svg = pages[0]
stem = f"molecules_{mode}"
dpi = FORMATS[mode].png_dpi
b64 = base64.b64encode(svg.encode()).decode()

with st.container(border=True, key="card-figure"):
    with st.container(horizontal=True, vertical_alignment="center"):
        st.caption(
            f"{len(smiles)} structure{'s' if len(smiles) != 1 else ''} · "
            f"{'ACS 1996 print style' if mode == 'print' else 'Slide style'} · "
            f"{info['width']:.2f} × {info['height']:.2f} in · bond {info['bond']:.1f} pt"
        )
        if len(pages) > 1:
            st.download_button(
                "GIF",
                export.to_gif(pages, dpi),
                file_name=f"{stem}.gif",
                mime="image/gif",
                icon=":material/download:",
                help=f"All molecules, {info['capacity']} per frame.",
            )
        else:
            st.download_button(
                "SVG",
                svg,
                file_name=f"{stem}.svg",
                mime="image/svg+xml",
                icon=":material/download:",
                help="Vector, for Word, Illustrator or slides.",
            )
            st.download_button(
                "PNG",
                export.to_png(svg, dpi),
                file_name=f"{stem}.png",
                mime="image/png",
                icon=":material/download:",
                help=f"{dpi} dpi raster.",
            )
    st.markdown(
        f"<div class='figure'><img src='data:image/svg+xml;base64,{b64}' "
        f"alt='Grid of {len(smiles)} chemical structures'/></div>",
        unsafe_allow_html=True,
    )
