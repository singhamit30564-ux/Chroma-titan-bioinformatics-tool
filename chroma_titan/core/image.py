"""Bio-image analysis on numpy arrays (skimage-free).

Images are read from netpbm files (``.pgm``/``.ppm``, the only lossless plain
formats shipped with the examples) or built synthetically, and every operation is
implemented with numpy so it works headless.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np

from chroma_titan.core.io import resolve


def read_image(src) -> np.ndarray:
    """Read a PGM/PPM (P2/P3/P5/P6) file, or accept an ndarray directly."""
    if isinstance(src, np.ndarray):
        return src.astype(float)
    if hasattr(src, "shape"):
        return np.asarray(src, dtype=float)
    p = Path(resolve(src))
    data = p.read_bytes()
    if data[:2] in (b"P2", b"P3", b"P5", b"P6"):
        kind = data[:2].decode()
        tokens: list[bytes] = []
        pos = 0
        n = len(data)
        while len(tokens) < 4 and pos < n:
            while pos < n and data[pos:pos + 1].isspace():
                pos += 1
            if pos < n and data[pos:pos + 1] == b"#":
                while pos < n and data[pos:pos + 1] != b"\n":
                    pos += 1
                continue
            end = pos
            while end < n and not data[end:end + 1].isspace():
                end += 1
            if end == pos:
                break
            tokens.append(data[pos:end])
            pos = end
        if len(tokens) < 4:
            raise ValueError("truncated PGM/PPM header")
        w, h, maxv = int(tokens[1]), int(tokens[2]), int(tokens[3])
        if kind in ("P2", "P3"):
            nums = np.array([float(x) for x in data[pos:].split()], dtype=float)
            nch = 3 if len(nums) >= w * h * 3 else 1
            return nums[: w * h * nch].reshape(h, w, nch).astype(float)
        raw = data[pos + 1:]
        if kind == "P5":
            arr = np.frombuffer(raw[: w * h], dtype=np.uint8).astype(float).reshape(h, w)
            return arr
        return np.frombuffer(raw[: w * h * 3], dtype=np.uint8).astype(float).reshape(h, w, 3)
    if data[:4] == b"\x89PNG" or data[:2] == b"\xff\xd8":
        try:
            from PIL import Image
            import io as _io

            with Image.open(_io.BytesIO(data)) as im:
                return np.asarray(im, dtype=float)
        except Exception as exc:  # noqa: BLE001
            raise ValueError(f"cannot decode image: {exc}") from None
    # fall back: parse a whitespace matrix (useful for tests)
    rows = [ln for ln in data.decode("utf-8", "replace").splitlines() if ln.strip()]
    return np.array([[float(x) for x in ln.replace(",", " ").split()] for ln in rows], dtype=float)


def _count_ws_tokens(header: bytes) -> int:
    return len([t for t in re.split(rb"\s+", header.strip()) if t])


def write_pgm(arr: np.ndarray, path: str | None = None, maxval: int = 255) -> str:
    a = np.clip(arr, 0, maxval)
    h, w = a.shape[:2]
    if a.ndim == 3:
        text = f"P3\n{w} {h}\n{maxval}\n"
        flat = a.astype(int).flatten().tolist()
        lines = [" ".join(str(v) for v in flat[i:i + w * 3]) for i in range(0, len(flat), w * 3)]
        text += "\n".join(lines) + "\n"
    else:
        text = f"P2\n{w} {h}\n{maxval}\n"
        text += "\n".join(" ".join(str(int(v)) for v in row) for row in a) + "\n"
    if path:
        Path(path).write_text(text)
    return text


def to_gray(img: np.ndarray) -> np.ndarray:
    a = np.asarray(img, dtype=float)
    if a.ndim == 3:
        if a.shape[-1] == 1:
            return a[..., 0]
        if a.shape[-1] == 2:
            return a[..., :2].mean(axis=-1)
        if a.shape[-1] >= 3:
            return 0.2989 * a[..., 0] + 0.5870 * a[..., 1] + 0.1140 * a[..., 2]
        return a.reshape(a.shape[0], -1)
    return a


def normalize(img: np.ndarray, lo: float | None = None, hi: float | None = None,
              mode: str = "minmax") -> np.ndarray:
    a = img.astype(float)
    if mode == "zscore":
        m, s = float(a.mean()), float(a.std()) or 1.0
        return (a - m) / s
    if mode == "percentile":
        lo = lo if lo is not None else 1.0
        hi = hi if hi is not None else 99.0
        a_lo, a_hi = np.percentile(a, lo), np.percentile(a, hi)
    else:
        a_lo = float(np.min(a)) if lo is None else lo
        a_hi = float(np.max(a)) if hi is None else hi
    rng = (a_hi - a_lo) or 1.0
    return np.clip((a - a_lo) / rng, 0, 1)


def threshold_otsu(img: np.ndarray, nbins: int = 256) -> dict:
    g = to_gray(img).flatten()
    lo, hi = float(g.min()), float(g.max())
    hist, edges = np.histogram(g, bins=nbins, range=(lo, hi))
    hist = hist.astype(float)
    total = hist.sum() or 1.0
    p = hist / total
    omega = np.cumsum(p)
    mid = (edges[:-1] + edges[1:]) / 2
    mu = np.cumsum(p * mid)
    mu_t = mu[-1]
    denom = omega * (1 - omega)
    denom[denom == 0] = 1e-12
    sigma_b = (mu_t * omega - mu) ** 2 / denom
    k = int(np.argmax(sigma_b))
    thr = float(mid[k])
    return {"threshold": round(thr, 4), "nbins": nbins,
            "between_class_variance": round(float(sigma_b[k]), 6),
            "foreground_fraction": round(float((g > thr).mean()), 5)}


def binarize(img: np.ndarray, value: float | None = None, above: bool = True) -> np.ndarray:
    g = to_gray(img)
    t = float(np.mean(g)) if value is None else value
    return ((g > t) if above else (g <= t)).astype(float)


def gaussian_blur(img: np.ndarray, sigma: float = 1.0, radius: int = 0) -> np.ndarray:
    r = radius or max(1, int(round(3 * sigma)))
    x = np.arange(-r, r + 1)
    k = np.exp(-(x ** 2) / (2 * sigma * sigma))
    k /= k.sum()
    out = _convolve1d(_convolve1d(img, k, axis=0), k, axis=1)
    return out


def _convolve1d(a: np.ndarray, kernel: np.ndarray, axis: int = 1) -> np.ndarray:
    pad = len(kernel) // 2
    a = np.moveaxis(a, axis, -1)
    padded = np.pad(a, ((0, 0),) * (a.ndim - 1) + ((pad, pad),), mode="edge")
    out = np.zeros_like(a, dtype=float)
    for i, w in enumerate(kernel):
        out += w * padded[..., i:i + a.shape[-1]]
    return np.moveaxis(out, -1, axis)


def median_filter(img: np.ndarray, size: int = 3) -> np.ndarray:
    r = size // 2
    a = to_gray(img).astype(float)
    padded = np.pad(a, r, mode="edge")
    h, w = a.shape
    stack = np.stack([padded[i:i + h, j:j + w] for i in range(size) for j in range(size)], -1)
    return np.median(stack, axis=-1)


def sobel(img: np.ndarray) -> dict:
    g = to_gray(img).astype(float)
    kx = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], dtype=float)
    ky = kx.T
    gx = _convolve2d(g, kx)
    gy = _convolve2d(g, ky)
    return {"magnitude": np.hypot(gx, gy), "angle": np.degrees(np.arctan2(gy, gx)),
            "gradient_x": gx, "gradient_y": gy}


def _convolve2d(a: np.ndarray, k: np.ndarray) -> np.ndarray:
    r = k.shape[0] // 2
    p = np.pad(a, r, mode="edge")
    h, w = a.shape
    out = np.zeros_like(a, dtype=float)
    for i in range(k.shape[0]):
        for j in range(k.shape[1]):
            out += k[i, j] * p[i:i + h, j:j + w]
    return out


def laplacian(img: np.ndarray) -> np.ndarray:
    k = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]], dtype=float)
    return _convolve2d(to_gray(img), k)


def histogram(img: np.ndarray, bins: int = 256) -> list[dict]:
    g = to_gray(img).flatten()
    hist, edges = np.histogram(g, bins=bins)
    return [{"bin_min": round(float(edges[i]), 3), "bin_max": round(float(edges[i + 1]), 3),
             "count": int(hist[i]), "fraction": round(float(hist[i]) / max(1, g.size), 6)}
            for i in range(len(hist))]


def equalize_histogram(img: np.ndarray) -> np.ndarray:
    g = to_gray(img)
    lo, hi = float(g.min()), float(g.max())
    rng = (hi - lo) or 1.0
    q = np.clip((g - lo) / rng * 255, 0, 255).astype(np.uint64)
    hist = np.bincount(q.flatten(), minlength=256).astype(float)
    cdf = hist.cumsum()
    cdf = cdf / cdf[-1]
    return (cdf[q] * (hi - lo) + lo).astype(float)


def label(img: np.ndarray, connectivity: int = 2, background: float | None = None) -> dict:
    """Two-pass connected-component labelling (4- or 8-connectivity)."""
    a = to_gray(img)
    binary = (a > np.mean(a)) if background is None else (a != background)
    binary = binary.astype(np.int8)
    h, w = binary.shape
    lab = np.zeros((h, w), dtype=np.int64)
    parent = [0]
    back = ((-1, -1), (-1, 0), (-1, 1), (0, -1)) if connectivity == 2 else ((-1, 0), (0, -1))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    next_id = 1
    for i in range(h):
        for j in range(w):
            if not binary[i, j]:
                continue
            roots = []
            for di, dj in back:
                ni, nj = i + di, j + dj
                if 0 <= ni < h and 0 <= nj < w and lab[ni, nj]:
                    roots.append(find(int(lab[ni, nj])))
            if roots:
                m = min(roots)
                lab[i, j] = m
                for r in roots:
                    if r != m:
                        parent[max(r, m)] = min(r, m)
            else:
                lab[i, j] = next_id
                parent.append(next_id)
                next_id += 1
    remap: dict[int, int] = {}
    out = np.zeros_like(lab)
    for i in range(h):
        for j in range(w):
            if lab[i, j]:
                r = find(int(lab[i, j]))
                if r not in remap:
                    remap[r] = len(remap) + 1
                out[i, j] = remap[r]
    return {"labels": out, "n_objects": len(remap)}


def regionprops(labeled: np.ndarray, img: np.ndarray | None = None) -> list[dict]:
    rows = []
    h, w = labeled.shape
    for idx in range(1, int(labeled.max()) + 1):
        mask = labeled == idx
        ys, xs = np.nonzero(mask)
        area = int(mask.sum())
        if not area:
            continue
        cy, cx = float(ys.mean()), float(xs.mean())
        if len(ys) > 1:
            cov = np.cov(np.stack([xs - cx, ys - cy]))
            eig = np.linalg.eigvalsh(cov)
            major = 2 * math.sqrt(max(0.0, float(eig[-1])))
            minor = 2 * math.sqrt(max(1e-9, float(eig[0])))
            ecc = math.sqrt(max(0.0, 1 - (minor / major) ** 2)) if major else 0.0
        else:
            major = minor = 0.0
            ecc = 0.0
        intensity = float(np.mean(to_gray(img)[mask])) if img is not None else 1.0
        rows.append({
            "label": idx, "area": area, "centroid_x": round(cx, 3), "centroid_y": round(cy, 3),
            "bbox_xmin": int(xs.min()), "bbox_xmax": int(xs.max()),
            "bbox_ymin": int(ys.min()), "bbox_ymax": int(ys.max()),
            "eccentricity": round(ecc, 4), "major_axis_length": round(major, 4),
            "minor_axis_length": round(minor, 4), "mean_intensity": round(intensity, 4),
            "equivalent_diameter": round(2 * math.sqrt(area / math.pi), 4),
            "solidity": round(float(area / max(1, (xs.max() - xs.min() + 1) * (ys.max() - ys.min() + 1))), 4),
            "perimeter_estimate": round(float(np.sum(mask ^ np.roll(mask, 1, 0) | mask ^ np.roll(mask, 1, 1))), 1),
            "extent": round(area / max(1, w * h), 6),
        })
    return rows


def remove_small_objects(labeled: np.ndarray, min_area: int = 10) -> np.ndarray:
    out = labeled.copy()
    for idx in range(1, int(labeled.max()) + 1):
        if int((labeled == idx).sum()) < min_area:
            out[labeled == idx] = 0
    return out


def fill_holes(img: np.ndarray) -> np.ndarray:
    a = (to_gray(img) > np.mean(img)).astype(float)
    lab = label(1 - a)["labels"]
    holes = np.isin(lab, [i for i in range(1, int(lab.max()) + 1)
                          if not (np.nonzero(lab == i)[0].min() == 0
                                  or np.nonzero(lab == i)[0].max() == a.shape[0] - 1
                                  or np.nonzero(lab == i)[1].min() == 0
                                  or np.nonzero(lab == i)[1].max() == a.shape[1] - 1)])
    return np.clip(a + holes, 0, 1)


def binary_erosion(img: np.ndarray, iterations: int = 1) -> np.ndarray:
    a = (to_gray(img) > np.mean(img)).astype(float)
    for _ in range(iterations):
        m = np.ones_like(a)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                m = np.minimum(m, np.roll(np.roll(a, di, 0), dj, 1))
        a = m
    return a


def binary_dilation(img: np.ndarray, iterations: int = 1) -> np.ndarray:
    a = (to_gray(img) > np.mean(img)).astype(float)
    for _ in range(iterations):
        m = np.zeros_like(a)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                m = np.maximum(m, np.roll(np.roll(a, di, 0), dj, 1))
        a = m
    return a


def distance_transform(img: np.ndarray) -> np.ndarray:
    a = (to_gray(img) > np.mean(img)).astype(float)
    h, w = a.shape
    big = float(h + w)
    d = np.where(a > 0, 0.0, big)
    for i in range(h):
        for j in range(w):
            if a[i, j]:
                continue
            best = d[i, j]
            for di, dj in ((-1, 0), (0, -1), (-1, -1), (-1, 1)):
                ni, nj = i + di, j + dj
                if 0 <= ni < h and 0 <= nj < w:
                    best = min(best, d[ni, nj] + math.hypot(di, dj))
            d[i, j] = best
    for i in range(h - 1, -1, -1):
        for j in range(w - 1, -1, -1):
            if a[i, j]:
                continue
            best = d[i, j]
            for di, dj in ((1, 0), (0, 1), (1, 1), (1, -1)):
                ni, nj = i + di, j + dj
                if 0 <= ni < h and 0 <= nj < w:
                    best = min(best, d[ni, nj] + math.hypot(di, dj))
            d[i, j] = best
    return d


def crop(img: np.ndarray, top: int = 0, left: int = 0, height: int = 0, width: int = 0) -> np.ndarray:
    a = img
    h, w = a.shape[:2]
    height = height or h
    width = width or w
    return a[top:top + height, left:left + width]


def resize_nearest(img: np.ndarray, scale: float = 0.5, width: int = 0, height: int = 0) -> np.ndarray:
    a = img
    h, w = a.shape[:2]
    nh = height or max(1, int(round(h * scale)))
    nw = width or max(1, int(round(w * scale)))
    yi = (np.arange(nh) * h / nh).astype(int)
    xi = (np.arange(nw) * w / nw).astype(int)
    return a[np.ix_(yi, xi)]


def add_labels_overlay(labeled: np.ndarray, colors: int = 8) -> np.ndarray:
    """Coloured segmentation overlay (RGB float image)."""
    h, w = labeled.shape
    out = np.zeros((h, w, 3), dtype=float)
    palette = [(0.16, 0.62, 0.56), (0.91, 0.44, 0.32), (0.15, 0.27, 0.33), (0.96, 0.64, 0.38),
               (0.54, 0.69, 0.49), (0.91, 0.77, 0.42), (0.71, 0.40, 0.46), (0.43, 0.35, 0.48)]
    for idx in range(1, int(labeled.max()) + 1):
        c = palette[(idx - 1) % min(colors, len(palette))]
        out[labeled == idx] = c
    return out


def invert(img: np.ndarray, maxval: float | None = None) -> np.ndarray:
    a = img.astype(float)
    return (maxval if maxval is not None else float(a.max())) - a


def intensity_projection(img: np.ndarray, axis: int = 0, mode: str = "max") -> np.ndarray:
    img = np.asarray(img, dtype=float)
    a = to_gray(img) if (img.ndim == 3 and img.shape[-1] <= 4) else img
    if mode == "mean":
        return a.mean(axis=axis)
    if mode == "sum":
        return a.sum(axis=axis)
    if mode == "min":
        return a.min(axis=axis)
    return a.max(axis=axis)


def radial_profile(img: np.ndarray, bins: int = 32) -> list[dict]:
    g = to_gray(img)
    h, w = g.shape
    cy, cx = h / 2, w / 2
    yy, xx = np.mgrid[0:h, 0:w]
    r = np.hypot(yy - cy, xx - cx)
    edges = np.linspace(0, r.max(), bins + 1)
    idx = np.clip(np.digitize(r, edges) - 1, 0, bins - 1)
    out = []
    for b in range(bins):
        m = idx == b
        if m.any():
            out.append({"radius_min": round(float(edges[b]), 3), "radius_max": round(float(edges[b + 1]), 3),
                        "mean_intensity": round(float(g[m].mean()), 4),
                        "std": round(float(g[m].std()), 4), "pixels": int(m.sum())})
    return out


def texture_stats(img: np.ndarray, levels: int = 4) -> list[dict]:
    """Multi-scale Haar-like pyramid energy/entropy descriptor."""
    g = to_gray(img).astype(float)
    out: list[dict] = []
    cur = g
    for lvl in range(levels):
        h = cur.shape[0] - cur.shape[0] % 2
        w = cur.shape[1] - cur.shape[1] % 2
        if h < 2 or w < 2:
            break
        cur = cur[:h, :w]
        out.append({"level": lvl, "size": f"{h}x{w}",
                    "mean": round(float(cur.mean()), 4), "std": round(float(cur.std()), 4),
                    "energy": round(float((cur ** 2).mean()), 4),
                    "entropy": round(_entropy(cur), 5)})
        cur = (cur[0::2, 0::2] + cur[1::2, 0::2] + cur[0::2, 1::2] + cur[1::2, 1::2]) / 4
    return out


def _entropy(a: np.ndarray) -> float:
    hist = np.bincount(np.clip(a.astype(np.int64) - int(a.min()), 0, 255).flatten(), minlength=256).astype(float)
    p = hist / (hist.sum() or 1.0)
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def glcm(img: np.ndarray, distance: int = 1, angle: float = 0.0, levels: int = 8) -> list[dict]:
    """Grey-level co-occurrence matrix statistics (contrast, homogeneity, ...)."""
    g = to_gray(img).astype(float)
    lo, hi = float(g.min()), float(g.max()) or 1.0
    q = np.clip(((g - lo) / ((hi - lo) or 1.0) * (levels - 1)).round(), 0, levels - 1).astype(int)
    dy, dx = int(round(math.sin(angle) * distance)), int(round(math.cos(angle) * distance))
    n = q.shape[0] - abs(dy)
    m = q.shape[1] - abs(dx)
    if n <= 0 or m <= 0:
        return []
    a = q[:n, :m]
    b = q[abs(dy):abs(dy) + n, abs(dx):abs(dx) + m]
    M = np.zeros((levels, levels), dtype=float)
    np.add.at(M, (a, b), 1.0)
    M = M / (M.sum() or 1.0)
    ii, jj = np.meshgrid(np.arange(levels), np.arange(levels), indexing="ij")
    asm = float((M ** 2).sum())
    contrast = float((M * (ii - jj) ** 2).sum())
    ident = float(M.sum())
    corr_num = float((M * (ii - (M * ii).sum()) * (jj - (M * jj).sum())).sum())
    var_i = float((M * (ii - (M * ii).sum()) ** 2).sum()) or 1e-12
    var_j = float((M * (jj - (M * jj).sum()) ** 2).sum()) or 1e-12
    return [{"angle": angle, "distance": distance, "asm": round(asm, 6),
             "contrast": round(contrast, 6), "correlation": round(corr_num / math.sqrt(var_i * var_j), 6),
             "homogeneity": round(float((M / (1 + np.abs(ii - jj))).sum()), 6),
             "energy": round(math.sqrt(asm), 6), "idm": round(float((M / (1 + (ii - jj) ** 2)).sum()), 6),
             "probability": round(ident, 6)}]


def watershed_simple(img: np.ndarray, markers: int = 4) -> dict:
    """Marker-controlled region growing (watershed-lite) on a gradient image."""
    g = to_gray(img)
    grad = np.abs(laplacian(g))
    h, w = grad.shape
    flat = grad.flatten()
    q = np.clip((flat / (flat.max() or 1.0) * (markers - 1)).astype(int), 0, markers - 1)
    seeds = []
    for k in range(markers):
        idxs = np.nonzero(q == k)[0]
        if len(idxs):
            seeds.append(int(idxs[np.argsort(flat[idxs])[0]]))
    lab = np.full(h * w, -1, dtype=int)
    for i, s in enumerate(seeds):
        lab[s] = i
    frontier = list(seeds)
    while frontier:
        nxt = []
        for pos in frontier:
            y, x = divmod(pos, w)
            for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ny, nx = y + dy, x + dx
                if 0 <= ny < h and 0 <= nx < w and lab[ny * w + nx] < 0:
                    lab[ny * w + nx] = lab[pos]
                    nxt.append(ny * w + nx)
        frontier = nxt
    lab_img = lab.reshape(h, w).astype(float) + 1
    return {"labels": lab_img, "n_regions": len(seeds)}


def image_stats(img: np.ndarray) -> dict:
    g = to_gray(img)
    return {"width": int(g.shape[1]), "height": int(g.shape[0]),
            "channels": int(img.shape[2]) if img.ndim == 3 else 1,
            "min": round(float(g.min()), 4), "max": round(float(g.max()), 4),
            "mean": round(float(g.mean()), 4), "std": round(float(g.std()), 4),
            "median": round(float(np.median(g)), 4), "mse_brightness": round(float((g ** 2).mean()), 4),
            "dynamic_range": round(float(g.max() - g.min()), 4),
            "entropy_bits": round(_entropy(g), 4)}
