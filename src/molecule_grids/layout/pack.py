"""Squeezed layouts: tight grid, text-like rows and free (outline-based) packing.

All packers work in drawing units and are deterministic (seeded). Molecules may be
turned by quarter turns; members of one ``group`` (analogues sharing a scaffold) always
share their turn so they stay aligned.

Inputs shared by the packers:

- ``dims[i, r]``: ink ``(w, h)`` of molecule ``i`` turned ``r`` quarter turns.
- ``capw[i]``: caption width (0 without a caption).
- ``group[i]``: group id of molecule ``i``.
"""

import math
import random

import numpy as np

SEED = 0
TOLERANCE = 0.10  # fewer grid columns only when they save more than this share of the height


def _anneal(state, cost, neighbour, iterations, rng):
    """Simulated annealing; returns the best state seen."""
    cur, cur_cost = state, cost(state)
    best, best_cost = cur, cur_cost
    t0, t1 = 0.02 * cur_cost, 1e-4 * cur_cost
    for k in range(iterations):
        t = t0 * (t1 / t0) ** (k / max(iterations - 1, 1))
        new = neighbour(cur, rng)
        new_cost = cost(new)
        if new_cost <= cur_cost or rng.random() < math.exp((cur_cost - new_cost) / t):
            cur, cur_cost = new, new_cost
            if cur_cost < best_cost:
                best, best_cost = cur, cur_cost
    return best, best_cost


def _iterations(n):
    return min(4000, 150 * n)


def _turn_group(rot, groups, rng):
    """Flip the quarter turn (0 <-> 1) of one random group."""
    g = rng.choice(groups)
    rot = dict(rot)
    rot[g] = 1 - rot[g]
    return rot


def _blocks(idx, group):
    """Molecules as blocks: each family (shared ``group``) together, in input order."""
    blocks = {}
    for i in idx:
        blocks.setdefault(group[i], []).append(i)
    return [tuple(b) for b in blocks.values()]


def _move(seq, rng):
    """Swap two items of ``seq`` or move one elsewhere."""
    seq = list(seq)
    if len(seq) > 1:
        a, b = rng.sample(range(len(seq)), 2)
        if rng.random() < 0.5:
            seq[a], seq[b] = seq[b], seq[a]
        else:
            seq.insert(b, seq.pop(a))
    return seq


def _box(i, r, dims, capw):
    w, h = dims[i, r]
    return max(w, capw[i]), h


def grid_pack(idx, columns, max_width, dims, capw, group, gap, strip, together=False):
    """Tight aligned grid: column widths and row heights fit their largest molecule.

    Tries column counts from ``columns`` down, keeping the grid within ``max_width``, and
    keeps the one with the smallest height; fewer columns only when they save more than
    ``TOLERANCE`` of it. With ``together``, each family fills consecutive cells in reading
    order (the search moves whole families).

    Returns
    -------
    list of list of tuple
        Rows of ``(i, r)``.
    """
    rng = random.Random(SEED)
    groups = sorted({group[i] for i in idx})
    n = len(idx)

    def size(state):
        C, cells, rot = state
        R = len(cells) // C
        colw, rowh = [0.0] * C, [0.0] * R
        for k, i in enumerate(cells):
            if i is not None:
                w, h = _box(i, rot[group[i]], dims, capw)
                r, c = divmod(k, C)
                colw[c] = max(colw[c], w)
                rowh[r] = max(rowh[r], h)
        width = gap + sum(w + gap for w in colw if w > 0)
        height = gap + sum(h + strip + gap for h in rowh if h > 0)
        return height + 10 * max(0.0, width - max_width)  # too wide: heavily penalised

    def from_blocks(C, blocks):
        cells = [i for b in blocks for i in b]
        return cells + [None] * (-len(cells) % C)

    def neighbour(state, rng):
        C, cells, rot = state
        if together:  # reorder whole families, then refill the cells in reading order
            if rng.random() < 0.2:
                return C, cells, _turn_group(rot, groups, rng)
            order = _move(_blocks([i for i in cells if i is not None], group), rng)
            return C, from_blocks(C, order), rot
        if rng.random() < 0.2 or len(cells) < 2:
            return C, cells, _turn_group(rot, groups, rng)
        a, b = rng.sample(range(len(cells)), 2)
        if cells[a] is None and cells[b] is None:
            return state
        cells = list(cells)
        cells[a], cells[b] = cells[b], cells[a]
        return C, cells, rot

    def start(C):
        """Best of two sort-based layouts: similar widths per column, or similar heights per row."""
        R = math.ceil(n / C)
        rot = dict.fromkeys(groups, 0)
        if together:  # families in input order, largest first
            blocks = sorted(_blocks(idx, group), key=len, reverse=True)
            return C, from_blocks(C, blocks), rot
        bw = sorted(idx, key=lambda i: -_box(i, 0, dims, capw)[0])
        bh = sorted(idx, key=lambda i: -dims[i, 0][1])
        a = [None] * (R * C)
        for c in range(C):
            for r, i in enumerate(sorted(bw[c * R : (c + 1) * R], key=lambda i: -dims[i, 0][1])):
                a[r * C + c] = i
        b = [None] * (R * C)
        for r in range(R):
            for c, i in enumerate(sorted(bh[r * C : (r + 1) * C], key=lambda i: -_box(i, 0, dims, capw)[0])):
                b[r * C + c] = i
        return min(((C, a, rot), (C, b, rot)), key=size)

    best, best_cost = None, None
    for C in range(min(columns, n), 0, -1):
        state, c = _anneal(start(C), size, neighbour, _iterations(n), rng)
        if best is None or c < (1 - TOLERANCE) * best_cost:
            best, best_cost = state, c
        if c > 2 * best_cost:  # fewer columns only get taller from here
            break
    C, cells, rot = best
    rows = [[(i, None if i is None else rot[group[i]]) for i in cells[r : r + C]] for r in range(0, len(cells), C)]
    rows = [row for row in rows if any(i is not None for i, _ in row)]
    if together:  # keep reading order, so families stay together
        return rows
    # Drop empty columns; put fuller columns and rows first so gaps end up at the bottom right
    # (whole columns and rows move, so the area is unchanged).
    keep = [c for c in range(C) if any(row[c][0] is not None for row in rows)]
    keep.sort(key=lambda c: sum(row[c][0] is None for row in rows))
    rows = [[row[c] for c in keep] for row in rows]
    return sorted(rows, key=lambda row: sum(i is None for i, _ in row))


def rows_pack(idx, max_width, dims, capw, group, gap, strip, together=False):
    """Rows like text: molecules flow left to right at their own widths, up to ``max_width``.

    With ``together``, each family stays consecutive (the search moves whole families).

    Returns
    -------
    list of list of tuple
        Rows of ``(i, r)``.
    """
    rng = random.Random(SEED)
    groups = sorted({group[i] for i in idx})

    def breaks(state):
        order, rot = state
        rows, row, x = [], [], gap
        for i in (i for block in order for i in block):
            w = _box(i, rot[group[i]], dims, capw)[0]
            if row and x + w + gap > max_width:
                rows.append(row)
                row, x = [], gap
            row.append(i)
            x += w + gap
        return rows + [row]

    def size(state):  # the width is fixed, so only the height counts
        _, rot = state
        return gap + sum(max(_box(i, rot[group[i]], dims, capw)[1] for i in row) + strip + gap for row in breaks(state))

    def neighbour(state, rng):
        order, rot = state
        u = rng.random()
        if u < 0.2:
            return order, _turn_group(rot, groups, rng)
        return _move(order, rng), rot

    # Shelf start: landscape turns, tallest first (next-fit decreasing height).
    rot = {g: 0 for g in groups}
    for i in idx:
        w, h = dims[i, 0]
        if h > w and all(group[j] != group[i] for j in idx if j != i):
            rot[group[i]] = 1
    blocks = _blocks(idx, group) if together else [(i,) for i in idx]
    start = (sorted(blocks, key=lambda b: -max(dims[i, rot[group[i]]][1] for i in b)), rot)
    (order, rot), _ = _anneal(start, size, neighbour, _iterations(len(idx)), rng)
    return [[(i, rot[group[i]]) for i in row] for row in breaks((order, rot))]


def free_pack(idx, max_width, masks, res, group, together=False):
    """Interlocking packing on occupancy masks, largest first, lowest fit first.

    Parameters
    ----------
    idx : list of int
        Molecules to place.
    max_width : float
        Available width, in drawing units.
    masks : dict
        ``masks[i, r]`` is a boolean array (rows, cols) of the dilated outline, molecule
        plus caption, at ``res`` drawing units per pixel.
    res : float
        Drawing units per mask pixel.
    group : list
        Group id of each molecule; a group's first placement fixes its turn.
    together : bool
        Place each family's members one after another, so they land side by side.

    Returns
    -------
    list of tuple
        ``(i, r, x, y)``: mask top-left corner in drawing units, relative to the packing area.
    """
    width = max(int(max_width // res), max(m.shape[1] for m in masks.values()))
    occ = np.zeros((64, width), dtype=bool)
    area = {i: int(masks[i, 0].sum()) for i in idx}
    if together:  # largest families first, members back to back
        blocks = sorted(_blocks(idx, group), key=lambda b: -sum(area[i] for i in b))
        order = [i for b in blocks for i in sorted(b, key=lambda i: -area[i])]
    else:
        order = sorted(idx, key=lambda i: -area[i])
    turned, placed = {}, []
    for i in order:
        options = [turned[group[i]]] if group[i] in turned else [0, 1, 2, 3]
        best = None
        for r in options:
            m = masks[i, r]
            h, w = m.shape
            if w > width:
                continue
            filled = np.flatnonzero(occ.any(axis=1))
            top = int(filled[-1]) + 1 if filled.size else 0
            need = top + h + 1
            if occ.shape[0] < need:
                occ = np.vstack([occ, np.zeros((need - occ.shape[0] + 64, width), dtype=bool)])
            # Skyline: the band where a fit can first appear, from the lowest gap to the top.
            sky = np.where(occ.any(axis=0), occ.shape[0] - np.argmax(occ[::-1], axis=0), 0)
            y0 = max(0, int(sky.min()) - h)
            band = occ[y0 : top + h, :]
            win = np.lib.stride_tricks.sliding_window_view(band, m.shape)
            hits = np.tensordot(win, m, axes=([2, 3], [0, 1]))
            free = np.argwhere(hits == 0)
            if not free.size:
                continue
            ys, xs = free[:, 0] + y0, free[:, 1]
            k = np.lexsort((xs, ys + h))  # lowest bottom edge, then leftmost
            cand = (int(ys[k[0]] + h), int(xs[k[0]]), int(ys[k[0]]), r)
            if best is None or cand < best:
                best = cand
        _, x, y, r = best
        m = masks[i, r]
        occ[y : y + m.shape[0], x : x + m.shape[1]] |= m
        turned.setdefault(group[i], r)
        placed.append((i, r, x * res, y * res))
    return placed
