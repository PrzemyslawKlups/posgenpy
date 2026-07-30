"""Reading .pos / .rrng files and labelling atoms with ion species.

Pure-Python replacement for the .pos/.rrng handling previously split between
``preposgenpy.py`` and ``legacy/clustering_corrections.py``'s ``read_pos``/``read_rrng``.
"""
import re

import numpy as np
import pandas as pd

# ===========================================================================
# Opening POS file


def read_pos_to_df(filename):
    """Load a binary .pos file into a DataFrame with columns ionIdx, x, y, z, mc."""

    data = np.fromfile(filename, dtype=">f4")
    data = data.reshape((-1, 4))

    # convert to native-endian float32
    data = data.astype(np.float32)

    df = pd.DataFrame(
        {
            "ionIdx": np.arange(1, len(data) + 1),
            "x": data[:, 0],
            "y": data[:, 1],
            "z": data[:, 2],
            "mc": data[:, 3],
        }
    )

    return df


# ===========================================================================
# Opening RRNG file and merging data


def range_info_from_rrng(rrng_line):
    """
    Parse a single RRNG range line.

    Returns
    -------
    mc_begin : float
    mc_end : float
    ion_name : str
    ion_volume : float
    ion_color : str
        Hex colour string from the RRNG file.
    formula : dict[str, int]
        Element symbol -> atom count for this range, e.g. {"Fe": 2, "O": 3}.
    """

    rrng_line = rrng_line.strip()

    if not rrng_line.startswith("Range"):
        raise ValueError("Parsed line is not an RRNG range.")

    # Mass range
    mass_match = re.search(
        r'=\s*([\d.]+)\s+([\d.]+)',
        rrng_line
    )

    mc_begin = float(mass_match.group(1))
    mc_end = float(mass_match.group(2))

    # Volume
    vol_match = re.search(r'Vol:([\d.]+)', rrng_line)
    ion_volume = float(vol_match.group(1))

    # Keep original RRNG hex colour
    color_match = re.search(r'Color:([0-9A-Fa-f]{6})', rrng_line)
    ion_color = color_match.group(1).upper()

    # Element composition, e.g. "Fe:1 O:1" -> {"Fe": 1, "O": 1}.
    # Always parsed (regardless of an explicit Name: alias) so a formula is
    # available even for ranges given a custom Name.
    species_matches = re.findall(
        r'([A-Z][a-z]?):(\d+)',
        rrng_line
    )
    formula = {element: int(count) for element, count in species_matches}

    # Handle Name:Fm style entries
    name_match = re.search(r'Name:(\w+)', rrng_line)

    if name_match:
        ion_name = name_match.group(1)

    else:
        parts = []
        for element, count in formula.items():
            parts.append(element if count == 1 else f"{element}{count}")
        ion_name = "".join(parts)

    if not formula:
        # No parseable element composition (e.g. a bare custom Name with no
        # species list) - fall back to treating the name itself as a
        # single-atom "element" so downstream composition counting still works.
        formula = {ion_name: 1}

    return mc_begin, mc_end, ion_name, ion_volume, ion_color, formula


def read_rrng(filename):
    """Read a .rrng (IVAS/AP Suite "hash-style") range file."""

    rows = []

    with open(filename) as f:
        for line in f:
            if line.startswith("Range"):
                rows.append(range_info_from_rrng(line))

    return pd.DataFrame(
        rows,
        columns=["mc_begin", "mc_end", "ion", "volume", "color", "formula"]
    )


def read_apt_range_file(filename):
    """
    Read an AP Suite-style .RRNG export (the '[Ions]'/'[Ranges]' CSV format),
    as a fallback for files where read_rrng() finds no "Range..." lines.

    Returns the same schema as read_rrng(): mc_begin, mc_end, ion, volume,
    color, formula. Note this format doesn't carry a per-ion volume, so
    `volume` is NaN.
    """

    this_range_file = pd.read_csv(filename)

    not_range_part = True
    rows = []

    for row in this_range_file['[Ions]']:

        if "[Ranges]" in row:
            not_range_part = False
            continue
        if not_range_part:
            continue
        elif 'Number=' in row:
            continue

        split_row = row.split(sep=' ')
        split_row = [x for x in split_row if x != '']

        lower = split_row.pop(0)
        lower = lower.split('=')[-1]
        upper = split_row.pop(0)
        name = split_row.pop(0).replace("Name:", "")
        colour = split_row.pop(-1).replace("Color:", "")
        comp = ' '.join(split_row)

        species_matches = re.findall(r'([A-Z][a-z]?):(\d+)', comp)
        formula = {element: int(count) for element, count in species_matches}
        if not formula:
            formula = {name: 1}

        rows.append(
            {
                "mc_begin": float(lower),
                "mc_end": float(upper),
                "ion": name,
                "volume": np.nan,
                "color": colour.upper(),
                "formula": formula,
            }
        )

    return pd.DataFrame(
        rows,
        columns=["mc_begin", "mc_end", "ion", "volume", "color", "formula"]
    )


def build_ion_formula_map(rrng_df):
    """Map each unique ion label to its {element: count} formula.

    If the same ion label appears in more than one range (e.g. isotopes),
    the first range's formula wins.
    """

    formula_map = {}
    for ion, formula in zip(rrng_df["ion"], rrng_df["formula"]):
        formula_map.setdefault(ion, formula)
    return formula_map


def get_ion_colours_from_rrng(rrng_df):
    """Map each unique ion label to its IVAS hex colour, e.g. {"Fe": "#FF0000"}."""

    colours = rrng_df.loc[:, ["ion", "color"]].drop_duplicates()
    return {row.ion: f"#{row.color}" for row in colours.itertuples()}


def assign_rrng_to_pos(pos_df, rrng_df):
    """Label each atom in pos_df with its ion name + colour, by mc range.

    Ranges are tested in file order and the FIRST match wins, matching
    libatomprobe's `RangeFile::getIonID` (`src/io/ranges.cpp:2871`), which
    returns on its first hit. Bounds are inclusive at both ends. This only
    matters for range files with overlapping ranges - posgen does not reject
    those, it just takes the lowest-indexed match.
    """

    pos_df = pos_df.copy()

    pos_df["ion"] = None
    pos_df["color"] = None

    mc = pos_df["mc"].to_numpy()
    unassigned = np.ones(len(pos_df), dtype=bool)

    for low, high, ion, color in zip(
        rrng_df["mc_begin"],
        rrng_df["mc_end"],
        rrng_df["ion"],
        rrng_df["color"],
    ):

        mask = (mc >= low) & (mc <= high) & unassigned

        pos_df.loc[mask, "ion"] = ion
        pos_df.loc[mask, "color"] = color
        unassigned &= ~mask

    return pos_df
