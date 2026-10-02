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
TOLERANCE = 0.10  # fewer grid columns only when they save more than this share of the area


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


def _box(i, r, dims, capw):
    w, h = dims[i, r]
    return max(w, capw[i]), h


def grid_pack(idx, columns, dims, capw, group, gap, strip):
    """Tight aligned grid: column widths and row heights fit their largest molecule.

    Tries column counts from ``columns`` down to half of it and keeps fewer columns only
    when they save more than ``TOLERANCE`` of the area.

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
        used = [w for w in colw if w > 0]
        return (gap + sum(w + gap for w in used)) * (gap + sum(h + strip + gap for h in rowh if h > 0))

    def neighbour(state, rng):
        C, cells, rot = state
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
    for C in range(min(columns, n), max(1, columns // 2) - 1, -1):
        state, c = _anneal(start(C), size, neighbour, _iterations(n), rng)
        if best is None or c < (1 - TOLERANCE) * best_cost:
            best, best_cost = state, c
    C, cells, rot = best
    rows = [[(i, None if i is None else rot[group[i]]) for i in cells[r : r + C]] for r in range(0, len(cells), C)]
    rows = [row for row in rows if any(i is not None for i, _ in row)]
    # Drop empty columns; put fuller columns and rows first so gaps end up at the bottom right
    # (whole columns and rows move, so the area is unchanged).
    keep = [c for c in range(C) if any(row[c][0] is not None for row in rows)]
    keep.sort(key=lambda c: sum(row[c][0] is None for row in rows))
    rows = [[row[c] for c in keep] for row in rows]
    return sorted(rows, key=lambda row: sum(i is None for i, _ in row))


def rows_pack(idx, max_width, dims, capw, group, gap, strip):
    """Rows like text: molecules flow left to right at their own widths, up to ``max_width``.

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
        for i in order:
            w = _box(i, rot[group[i]], dims, capw)[0]
            if row and x + w + gap > max_width:
                rows.append(row)
                row, x = [], gap
            row.append(i)
            x += w + gap
        return rows + [row]

    def size(state):
        _, rot = state
        width = height = 0.0
        for row in breaks(state):
            boxes = [_box(i, rot[group[i]], dims, capw) for i in row]
            width = max(width, gap + sum(w + gap for w, _ in boxes))
            height += max(h for _, h in boxes) + strip + gap
        return width * (gap + height)

    def neighbour(state, rng):
        order, rot = state
        u = rng.random()
        if u < 0.2:
            return order, _turn_group(rot, groups, rng)
        order = list(order)
        a, b = rng.sample(range(len(order)), 2) if len(order) > 1 else (0, 0)
        if u < 0.6:
            order[a], order[b] = order[b], order[a]
        else:
            order.insert(b, order.pop(a))
        return order, rot

    # Shelf start: landscape turns, tallest first (next-fit decreasing height).
    rot = {g: 0 for g in groups}
    for i in idx:
        w, h = dims[i, 0]
        if h > w and all(group[j] != group[i] for j in idx if j != i):
            rot[group[i]] = 1
    start = (sorted(idx, key=lambda i: -dims[i, rot[group[i]]][1]), rot)
    (order, rot), _ = _anneal(start, size, neighbour, _iterations(len(idx)), rng)
    return [[(i, rot[group[i]]) for i in row] for row in breaks((order, rot))]


def free_pack(idx, max_width, masks, res, group):
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

    Returns
    -------
    list of tuple
        ``(i, r, x, y)``: mask top-left corner in drawing units, relative to the packing area.
    """
    width = max(int(max_width // res), max(m.shape[1] for m in masks.values()))
    occ = np.zeros((64, width), dtype=bool)
    order = sorted(idx, key=lambda i: -int(masks[i, 0].sum()))
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
