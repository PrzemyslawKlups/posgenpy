"""Synthetic test datasets, built from posgen's own generator scripts.

`examples/generateTestPos_easy.xml` (and `_hard.xml`) are the posgen scripts the
original posgenpy used to make its test cube: a 40 nm Fe matrix with 100
spherical Cu-Ni clusters cut into it. posgen is a C++ executable this package
no longer needs, so `generate_synthetic_pos` runs those scripts itself.

Only the handful of posgen elements the generator scripts use is understood,
and anything else raises rather than being skipped - a script this cannot run
faithfully should fail, not produce a plausible-looking wrong dataset. The
random numbers are numpy's, not posgen's, so the atoms differ from a posgen
run of the same script; the clusters (positions, radii, compositions) and the
matrix are the same.
"""
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

from .io import write_pos

__all__ = ["generate_synthetic_pos"]


def _xyz(element):
    """The x/y/z child values of a posgen element, as a (3,) array."""
    return np.array([float(element.find(axis).get("value")) for axis in "xyz"])


def _is_true(value):
    return str(value).strip().lower() in ("true", "1", "yes")


def _spatrand(element, rng):
    """Uniform random atoms in [0, bound]^3 at `density` atoms per nm^3."""

    bound = element.find("bound")
    size = np.array([float(bound.get(axis)) for axis in "xyz"])
    density = float(element.find("density").get("value"))

    # <atom mass=.../> and <count value=.../> come in pairs, in that order
    masses, weights = [], []
    for child in element:
        if child.tag == "atom":
            masses.append(float(child.get("mass")))
        elif child.tag == "count":
            weights.append(float(child.get("value")))
        elif child.tag not in ("bound", "density"):
            raise ValueError(f"<spatrand> child <{child.tag}> is not supported")

    if len(masses) != len(weights) or not masses:
        raise ValueError("<spatrand> needs one <count> after every <atom>")

    n = int(round(density * np.prod(size)))
    positions = rng.random((n, 3)) * size

    weights = np.asarray(weights)
    mc = rng.choice(np.asarray(masses), size=n, p=weights / weights.sum())

    return positions, mc


def _translation(element):
    """The total shift of a <geotransform>, which may only translate."""

    shift = np.zeros(3)
    for transform in element:
        if transform.tag != "translate":
            raise ValueError(
                f"<geotransform><{transform.tag}> is not supported, only <translate>"
            )
        shift = shift + _xyz(transform)
    return shift


def _clip_test(element):
    """A function telling, for an (n, 3) array, which atoms a <clip> keeps."""

    shapes = list(element)
    if len(shapes) != 1:
        raise ValueError("<clip> must hold exactly one shape")
    shape = shapes[0]

    if shape.tag == "sphere":
        radius = float(shape.get("radius"))
        centre = _xyz(shape)

        def inside(positions):
            # one cheap pass over x first: a 1 nm sphere in a 40 nm cube only
            # needs the full distance for the few atoms in its slab
            near = np.flatnonzero(np.abs(positions[:, 0] - centre[0]) <= radius)
            result = np.zeros(len(positions), dtype=bool)
            result[near] = ((positions[near] - centre) ** 2).sum(axis=1) <= radius**2
            return result

    elif shape.tag == "box":
        corners = list(shape)
        if [c.tag for c in corners] != ["x", "y", "z", "x", "y", "z"]:
            raise ValueError("<clip><box> needs two corners as x, y, z, x, y, z")
        values = np.array([float(c.get("value")) for c in corners]).reshape(2, 3)
        low, high = values.min(axis=0), values.max(axis=0)

        def inside(positions):
            return np.all((positions >= low) & (positions <= high), axis=1)

    else:
        raise ValueError(f"<clip><{shape.tag}> is not supported")

    # invert keeps what is outside the shape - a hole
    if _is_true(element.get("invert", "false")):
        return lambda positions: ~inside(positions)
    return inside


def generate_synthetic_pos(xml_file, pos_file=None, seed: int = 0) -> pd.DataFrame:
    """Run a posgen generator script and write the .pos file it describes.

    The script is run in document order on one working set of atoms:
    `<spatrand>` adds random atoms, `<geotransform><translate>` moves them,
    `<clip>` keeps the atoms inside a sphere or box (outside, with
    `invert="true"`), `<store>` keeps a copy under a name (and empties the
    working set with `clear="true"`), `<recall>` adds a stored copy back, and
    `<possave>` writes the .pos file.

    Parameters
    ----------
    xml_file : str or Path
        The posgen script, e.g. ``examples/generateTestPos_easy.xml``.
    pos_file : str or Path, optional
        Where to write the result. Defaults to the script's own
        `<possave file=...>`, next to the script.
    seed : int
        Seed for the random atom positions and labels, so the same seed
        always gives the same file.

    Returns
    -------
    pandas.DataFrame
        The atoms written, in the shape `read_pos_to_df` returns
        (ionIdx, x, y, z, mc).
    """

    xml_file = Path(xml_file)
    root = ET.parse(xml_file).getroot()
    if root.tag != "posscript":
        raise ValueError(f"{xml_file} is not a posgen script (<{root.tag}>)")

    rng = np.random.default_rng(seed)

    # The working set is a list of (positions, mc) chunks, joined only when it
    # is saved. The generator scripts recall the whole running dataset once
    # per cluster, and a list makes that an append instead of a copy of every
    # atom. Every step builds new arrays rather than editing these, so a store
    # can hold the same chunks the working set does.
    working = []
    stores = {}
    saved = None

    for element in root:
        tag = element.tag

        if tag == "version":
            continue

        elif tag == "spatrand":
            working.append(_spatrand(element, rng))

        elif tag == "geotransform":
            shift = _translation(element)
            working = [(positions + shift, mc) for positions, mc in working]

        elif tag == "clip":
            keeps = _clip_test(element)
            clipped = []
            for positions, mc in working:
                keep = keeps(positions)
                clipped.append(
                    (positions, mc) if keep.all() else (positions[keep], mc[keep])
                )
            working = clipped

        elif tag == "store":
            stores[element.get("name")] = list(working)
            if _is_true(element.get("clear", "false")):
                working = []

        elif tag == "recall":
            name = element.get("name")
            if name not in stores:
                raise ValueError(f"<recall name={name!r}> before any <store> of it")
            working = working + stores[name]

        elif tag == "possave":
            positions = np.concatenate([p for p, _ in working] or [np.zeros((0, 3))])
            mc = np.concatenate([m for _, m in working] or [np.zeros(0)])
            saved = pd.DataFrame(
                {
                    "ionIdx": np.arange(1, len(mc) + 1),
                    "x": positions[:, 0].astype(np.float32),
                    "y": positions[:, 1].astype(np.float32),
                    "z": positions[:, 2].astype(np.float32),
                    "mc": mc.astype(np.float32),
                }
            )
            target = Path(pos_file) if pos_file else xml_file.parent / element.get("file")
            write_pos(saved, target)

        else:
            raise ValueError(f"posgen element <{tag}> is not supported")

    if saved is None:
        raise ValueError(f"{xml_file} has no <possave>, so it writes nothing")

    return saved
