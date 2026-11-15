"""Reading .pos / .apt / .rrng files and labelling atoms with ion species.

Pure-Python replacement for the .pos/.rrng handling previously split between
``preposgenpy.py`` and ``legacy/clustering_corrections.py``'s ``read_pos``/``read_rrng``.
"""
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ===========================================================================
# Opening POS file


def _pos_frame(xyz, mc):
    """The DataFrame every reader returns: ionIdx (1-based), x, y, z, mc."""

    # native-endian float32, whatever byte order the file used
    xyz = np.asarray(xyz).astype(np.float32)
    mc = np.asarray(mc).astype(np.float32)

    return pd.DataFrame(
        {
            "ionIdx": np.arange(1, len(mc) + 1),
            "x": xyz[:, 0],
            "y": xyz[:, 1],
            "z": xyz[:, 2],
            "mc": mc,
        }
    )


def read_pos_to_df(filename):
    """Load a .pos file into a DataFrame with columns ionIdx, x, y, z, mc.

    A file ending in `.apt` (any case) is read as an AP Suite .apt file
    instead - see `read_apt_to_df` - and comes back as the same frame, so
    either can be passed here. The format is picked by the extension, not by
    the file's first bytes: a .pos has no header, and "APT\\0" read as a .pos
    float is a valid x coordinate (~13.02 nm).
    """

    if Path(filename).suffix.lower() == ".apt":
        return read_apt_to_df(filename)

    data = np.fromfile(filename, dtype=">f4")
    data = data.reshape((-1, 4))

    return _pos_frame(data[:, :3], data[:, 3])


def write_pos(pos_df, filename):
    """Write the x, y, z, mc columns of `pos_df` as a binary .pos file.

    The inverse of `read_pos_to_df`: one big-endian float32 record of
    (x, y, z, mc) per atom, which is all the .pos format holds - `ionIdx` and
    any other column are not written.
    """

    records = np.column_stack(
        [pos_df[column].to_numpy() for column in ("x", "y", "z", "mc")]
    )
    records.astype(">f4").tofile(filename)


# ===========================================================================
# Reading and writing APT (AP Suite) files

# AP Suite 6 ".apt" layout, as read by NOMAD's ifes_apt_tc_data_modeling
# (apt/apt6_headers.py, apt6_sections.py, apt6_reader.py): a 540-byte file
# header, then sections, each a 148-byte header followed by its data.
# Everything is little-endian; strings are fixed-length, zero-padded UTF-16.
_APT_FILE_HEADER = np.dtype([
    ("cSignature", "S4"),
    ("iHeaderSize", "<i4"),
    ("iHeaderVersion", "<i4"),
    ("wcFileName", "<u2", (256,)),
    ("ftCreationTime", "<u8"),       # Windows FILETIME: 100 ns since 1601
    ("llIonCount", "<u8"),
])
_APT_SECTION_HEADER = np.dtype([
    ("cSignature", "S4"),
    ("iHeaderSize", "<i4"),
    ("iHeaderVersion", "<i4"),
    ("wcSectionType", "<u2", (32,)),
    ("iSectionVersion", "<i4"),
    ("eRelationshipType", "<u4"),    # 1 = one record per ion
    ("eRecordType", "<u4"),          # 1 = fixed-size records
    ("eRecordDataType", "<u4"),      # 3 = IEEE float
    ("iDataTypeSize", "<i4"),        # bits
    ("iRecordSize", "<i4"),          # bytes
    ("wcDataUnit", "<u2", (16,)),
    ("llRecordCount", "<u8"),
    ("llByteCount", "<u8"),
])
assert _APT_FILE_HEADER.itemsize == 540 and _APT_SECTION_HEADER.itemsize == 148

# seconds between the FILETIME epoch (1601-01-01) and the Unix one
_FILETIME_EPOCH_OFFSET = 11_644_473_600


def _utf16(text, length):
    units = np.frombuffer(text.encode("utf-16-le"), dtype="<u2")[: length - 1]
    out = np.zeros(length, dtype="<u2")
    out[: len(units)] = units
    return out


def _apt_section(section_type, unit, n_ions, values_per_ion, header_size):
    header = np.zeros(1, dtype=_APT_SECTION_HEADER)
    header["cSignature"] = b"SEC\x00"
    header["iHeaderSize"] = header_size
    header["iHeaderVersion"] = 2
    header["wcSectionType"] = _utf16(section_type, 32)
    header["iSectionVersion"] = 1
    header["eRelationshipType"] = 1
    header["eRecordType"] = 1
    header["eRecordDataType"] = 3
    header["iDataTypeSize"] = 32
    header["iRecordSize"] = 4 * values_per_ion
    header["wcDataUnit"] = _utf16(unit, 16)
    header["llRecordCount"] = n_ions
    header["llByteCount"] = n_ions * 4 * values_per_ion
    return header


def write_apt(pos_df, filename):
    """Write the x, y, z, mc columns of `pos_df` as an AP Suite .apt file.

    The file holds the two sections a .pos file carries - "Position" (nm) and
    "Mass" (Da), as float32 - so it is the same data as `write_pos` in the
    format AP Suite 6 saves its own experiments in.

    The layout follows NOMAD's open reader of the format
    (`ifes_apt_tc_data_modeling`), not AMETEK documentation, and has not been
    checked against AP Suite itself. In particular, the "Position" section
    carries six extra floats before its data, which that reader documents as
    the point cloud's min/max bounds "in x, y, z direction respectively" -
    written here as xmin, xmax, ymin, ymax, zmin, zmax. If AP Suite will not
    open a file, `write_pos` writes the same atoms in the older format it
    also imports.
    """

    positions = pos_df[["x", "y", "z"]].to_numpy(dtype="<f4")
    mass = pos_df["mc"].to_numpy(dtype="<f4")
    n_ions = len(mass)
    if n_ions == 0:
        # the format has no empty file: readers reject llIonCount = 0
        raise ValueError("an .apt file must hold at least one ion")

    header = np.zeros(1, dtype=_APT_FILE_HEADER)
    header["cSignature"] = b"APT\x00"
    header["iHeaderSize"] = _APT_FILE_HEADER.itemsize
    header["iHeaderVersion"] = 2
    header["wcFileName"] = _utf16(Path(filename).name, 256)
    header["ftCreationTime"] = int((time.time() + _FILETIME_EPOCH_OFFSET) * 10**7)
    header["llIonCount"] = n_ions

    bounds = np.column_stack(
        [positions.min(axis=0), positions.max(axis=0)]
    ).ravel().astype("<f4")

    with open(filename, "wb") as f:
        f.write(header.tobytes())
        # the six bounds count as part of the Position section's header
        f.write(_apt_section("Position", "nm", n_ions, 3, 148 + 6 * 4).tobytes())
        f.write(bounds.tobytes())
        f.write(positions.tobytes())
        f.write(_apt_section("Mass", "Da", n_ions, 1, 148).tobytes())
        f.write(mass.tobytes())


def _from_utf16(units):
    """Decode a zero-padded UTF-16 header field - the inverse of `_utf16`."""

    units = np.asarray(units, dtype="<u2")
    zeros = np.flatnonzero(units == 0)
    if len(zeros):
        units = units[: zeros[0]]
    return units.tobytes().decode("utf-16-le", errors="replace")


def _apt_sections(f, filename, start, file_size):
    """Walk an .apt file's sections: name -> (section header, data offset).

    Each section is a header of `iHeaderSize` bytes (148, or more when it
    carries extra fields - "Position" has its six bounds there) followed by
    `llByteCount` bytes of data, and the next section starts right after.
    """

    sections = {}
    offset = start

    while offset < file_size:
        if offset + _APT_SECTION_HEADER.itemsize > file_size:
            raise ValueError(f"{filename} is truncated inside a section header")

        f.seek(offset)
        section = np.fromfile(f, dtype=_APT_SECTION_HEADER, count=1)[0]
        if bytes(section["cSignature"]) != b"SEC":
            raise ValueError(f"{filename} has no section header at byte {offset}")
        header_size = int(section["iHeaderSize"])
        if header_size < _APT_SECTION_HEADER.itemsize:
            raise ValueError(
                f"{filename} has a section header of {header_size} bytes at "
                f"byte {offset}"
            )

        data_offset = offset + header_size
        offset = data_offset + int(section["llByteCount"])
        if offset > file_size:
            raise ValueError(f"{filename} is truncated inside a section's data")

        sections.setdefault(_from_utf16(section["wcSectionType"]), (section, data_offset))

    return sections


def _apt_float_section(f, filename, sections, name, values_per_ion, n_ions):
    """One of the file's per-ion float sections, as an (n_ions, values_per_ion) array."""

    if name not in sections:
        raise ValueError(
            f"{filename} has no {name!r} section (it has: {', '.join(sections)})"
        )

    section, data_offset = sections[name]
    bits = int(section["iDataTypeSize"])
    if (
        section["eRecordDataType"] != 3
        or bits not in (32, 64)
        or section["iRecordSize"] != values_per_ion * bits // 8
    ):
        raise ValueError(
            f"the {name!r} section of {filename} is not {values_per_ion} "
            f"float(s) per ion"
        )
    if section["llRecordCount"] != n_ions:
        raise ValueError(
            f"the {name!r} section of {filename} has {section['llRecordCount']} "
            f"records for {n_ions} ions"
        )
    if section["llByteCount"] < n_ions * section["iRecordSize"]:
        raise ValueError(f"the {name!r} section of {filename} is short of data")

    f.seek(data_offset)
    data = np.fromfile(f, dtype=f"<f{bits // 8}", count=n_ions * values_per_ion)
    return data.reshape(n_ions, values_per_ion)


def read_apt_to_df(filename):
    """Load an AP Suite .apt file into the frame `read_pos_to_df` returns.

    An .apt holds one section per measured quantity - in an AP Suite 6
    export: Voltage, Delta Pulse, pulse, erate, Temp, Mass, Noise, tofc,
    Multiplicity, Vap, Detector Coordinates, Position and Epos ToF. The
    columns are taken from "Position" (x, y, z in nm) and "Mass" (mc in Da),
    the two a .pos file holds. Every other section is skipped without being
    read, so a large file costs no more memory than its .pos would.

    The layout is the one `write_apt` writes (see the notes there). It
    matches an AP Suite 6 export, the six floats of "Position"'s header
    included (xmin, xmax, ymin, ymax, zmin, zmax). `read_pos_to_df` calls
    this for any file ending in `.apt`.

    Raises
    ------
    ValueError
        If the file is not an .apt file, is truncated, or has no
        Position/Mass section with one float record per ion.
    """

    file_size = Path(filename).stat().st_size

    with open(filename, "rb") as f:
        header = np.fromfile(f, dtype=_APT_FILE_HEADER, count=1)
        if len(header) == 0 or bytes(header[0]["cSignature"]) != b"APT":
            raise ValueError(f"{filename} is not an AP Suite .apt file")
        header = header[0]

        header_size = int(header["iHeaderSize"])
        if header_size < _APT_FILE_HEADER.itemsize:
            raise ValueError(f"{filename} has a file header of {header_size} bytes")

        n_ions = int(header["llIonCount"])
        sections = _apt_sections(f, filename, header_size, file_size)
        xyz = _apt_float_section(f, filename, sections, "Position", 3, n_ions)
        mc = _apt_float_section(f, filename, sections, "Mass", 1, n_ions)

    return _pos_frame(xyz, mc[:, 0])


# ===========================================================================
# The final clusters as posgen's per-atom files


def export_cluster_pos_files(
    clustered_atoms, pos_df, output_dir, prefix, apt=False, id_offset=1
):
    """Write a cluster search out as the three files posgen wrote for it.

    So the clusters can be opened in AP Suite, IVAS or anything else that
    reads .pos files:

    | file | holds |
    |---|---|
    | `{prefix}_cluster.pos` | every clustered atom - core, border and the enveloped bulk - with its own mass |
    | `{prefix}_clusterID.pos` | the same atoms in the same order, with `cluster_id + id_offset` in the mass column |
    | `{prefix}_unclustered.pos` | every other atom of `pos_df`, unranged ones included |

    `id_offset=1` is posgen's own `<clusterid offset="1">`, so cluster 0 of
    `cluster_stats` is "1 Da" in the clusterID file and no cluster sits at
    0 Da. The cluster and unclustered files together are exactly the atoms of
    `pos_df`. Clustered atoms are grouped by cluster, in pos-file order within
    each.

    Parameters
    ----------
    clustered_atoms : pandas.DataFrame
        `run_clustering` output: ionIdx, x, y, z, mc and cluster_id (-1 =
        unclustered).
    pos_df : pandas.DataFrame
        The whole tip as `read_pos_to_df` returns it - the source of the
        unclustered atoms, since `run_clustering` drops unranged ones.
    apt : bool
        Also write each file as `.apt` (see `write_apt`). An empty file is
        skipped, since the format cannot hold zero ions.

    Returns
    -------
    dict
        File kind ("cluster", "clusterID", "unclustered", plus "*_apt") ->
        path written.
    """

    for name, frame in (("clustered_atoms", clustered_atoms), ("pos_df", pos_df)):
        if "ionIdx" not in frame.columns:
            raise ValueError(f"{name} needs an ionIdx column to match atoms by")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    clustered = clustered_atoms.loc[clustered_atoms["cluster_id"] != -1]
    clustered = clustered.sort_values(["cluster_id", "ionIdx"], kind="stable")

    # ionIdx is only a row number, which any other tip has too - so the
    # clustered atoms must also sit where pos_df puts those rows
    source = pos_df.set_index("ionIdx").reindex(clustered["ionIdx"])
    if source["x"].isna().any() or not np.allclose(
        source[["x", "y", "z"]].to_numpy(), clustered[["x", "y", "z"]].to_numpy()
    ):
        raise ValueError(
            "clustered_atoms holds atoms that are not in pos_df - they are "
            "not from the same tip"
        )

    unclustered = pos_df.loc[~pos_df["ionIdx"].isin(clustered["ionIdx"])]

    ids = clustered[["x", "y", "z"]].copy()
    ids["mc"] = clustered["cluster_id"].to_numpy() + id_offset

    written = {}
    for kind, frame in (
        ("cluster", clustered), ("clusterID", ids), ("unclustered", unclustered)
    ):
        path = output_dir / f"{prefix}_{kind}.pos"
        write_pos(frame, path)
        written[kind] = path

        if apt and len(frame):
            apt_path = path.with_suffix(".apt")
            write_apt(frame, apt_path)
            written[f"{kind}_apt"] = apt_path

    return written


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
