"""Assert that no two *rendered* text artists overlap in a matplotlib figure.

Panel layouts here were rearranged repeatedly and each rearrangement introduced a
collision only caught by eye. This makes the check mechanical.

Two classes of false positive are excluded deliberately:
  * tick labels for ticks outside the axes view limits - the artists exist and report
    visible=True, but are clipped and never drawn;
  * identical text drawn at the same place by twinned axes (ax.twinx() duplicates the
    shared x tick labels exactly on top of the originals).
"""
import matplotlib.text


def _live_ticklabels(ax):
    out = []
    for axis, lim, getter in ((ax.xaxis, ax.get_xlim(), ax.get_xticklabels),
                              (ax.yaxis, ax.get_ylim(), ax.get_yticklabels)):
        lo, hi = min(lim), max(lim)
        locs = axis.get_ticklocs()
        labs = getter()
        for loc, lab in zip(locs, labs):
            if lo - 1e-9 <= loc <= hi + 1e-9:
                out.append(lab)
    return out


def find_text_overlaps(fig, tol=2.0):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    items = []
    for ai, ax in enumerate(fig.axes):
        if not ax.axison:                       # axes used only to host a legend
            continue
        cand = [ax.title, ax.xaxis.label, ax.yaxis.label]
        cand += _live_ticklabels(ax)
        cand += [c for c in ax.get_children()
                 if isinstance(c, matplotlib.text.Text) and c not in (ax.title,)]
        for t in cand:
            if t is None or not t.get_visible() or not t.get_text().strip():
                continue
            if any(t is prev for prev, _, _ in items):
                continue
            try:
                bb = t.get_window_extent(renderer=rend)
            except Exception:
                continue
            if bb.width <= 0 or bb.height <= 0:
                continue
            items.append((t, f"ax{ai}:" + t.get_text().replace("\n", " / ")[:40], bb))

    out = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            ta, tb = items[i][0], items[j][0]
            a, b = items[i][2], items[j][2]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox <= tol or oy <= tol:
                continue
            # twinned axes redraw the shared tick labels exactly on top of each other
            if ta.get_text() == tb.get_text():
                inter = ox * oy
                if inter > 0.85 * min(a.width * a.height, b.width * b.height):
                    continue
            out.append((items[i][1], items[j][1], round(ox), round(oy)))
    return out, len(items)


def assert_no_text_overlap(fig, name=""):
    clashes, n = find_text_overlaps(fig)
    if clashes:
        print(f"  !! {name}: {len(clashes)} TEXT OVERLAP(S)")
        for a, b, ox, oy in clashes[:12]:
            print(f"     '{a}'  <->  '{b}'   overlap {ox}x{oy}px")
    else:
        print(f"  {name}: text-overlap check clean ({n} text items)")
    return clashes
