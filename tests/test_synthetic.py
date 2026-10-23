"""Tests for the .pos writer and the posgen-script generator of the test cube.

The example notebooks run on a cube this generator makes, so a mistake here
would show up as a plausible clustering of the wrong data. These pin what the
generator scripts rely on: holes are really empty, clusters really fill them,
and the counts and compositions follow the script.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import posgenpy as pgp

EXAMPLES = Path(__file__).parent.parent / "examples"

# a 10 nm matrix of masses 1:3 = 1:9 with one 1.5 nm hole, refilled by a
# pure-mass-2 cluster - the same pattern generateTestPos_*.xml repeats 100 times
TINY_SCRIPT = """<!DOCTYPE posscript SYSTEM "posscript.dtd">
<posscript>
  <version value="0.0.1"/>
  <spatrand>
    <bound x="10" y="10" z="10"/>
    <density value="20"/>
    <atom index="(000)" mass="1.0"/>
    <count value="1"/>
    <atom index="(000)" mass="3.0"/>
    <count value="9"/>
  </spatrand>
  <geotransform><translate><x value="-5"/><y value="-5"/><z value="-5"/></translate></geotransform>
  <clip invert="true"><sphere radius="1.5"><x value="1"/><y value="2"/><z value="-1"/></sphere></clip>
  <store name="running" forceram="true" clear="true"/>

  <spatrand>
    <bound x="3" y="3" z="3"/>
    <density value="20"/>
    <atom index="(000)" mass="2.0"/>
    <count value="1"/>
  </spatrand>
  <geotransform><translate><x value="-1.5"/><y value="-1.5"/><z value="-1.5"/></translate></geotransform>
  <clip><sphere radius="1.5"><x value="0"/><y value="0"/><z value="0"/></sphere></clip>
  <geotransform><translate><x value="1"/><y value="2"/><z value="-1"/></translate></geotransform>
  <recall name="running"/>
  <store name="running" clear="true" forceram="true"/>

  <recall name="running"/>
  <clip><box><x value="-5"/><y value="-5"/><z value="-5"/><x value="5"/><y value="5"/><z value="5"/></box></clip>
  <possave file="tiny.pos"/>
</posscript>
"""


@pytest.fixture
def tiny(tmp_path):
    script = tmp_path / "tiny.xml"
    script.write_text(TINY_SCRIPT)
    return script


def test_write_pos_round_trips_through_read_pos(tmp_path):
    rng = np.random.default_rng(1)
    original = pd.DataFrame(rng.normal(size=(50, 4)), columns=["x", "y", "z", "mc"])

    pgp.write_pos(original, tmp_path / "round.pos")
    read = pgp.read_pos_to_df(tmp_path / "round.pos")

    assert list(read.columns) == ["ionIdx", "x", "y", "z", "mc"]
    np.testing.assert_array_equal(
        read[["x", "y", "z", "mc"]].to_numpy(),
        original.to_numpy().astype(np.float32),
    )


def test_hole_is_empty_and_the_cluster_fills_it(tiny):
    atoms = pgp.generate_synthetic_pos(tiny, seed=0)
    positions = atoms[["x", "y", "z"]].to_numpy(float)
    in_sphere = ((positions - [1, 2, -1]) ** 2).sum(axis=1) <= 1.5**2

    # no matrix atom (mass 1 or 3) inside, no cluster atom (mass 2) outside
    assert set(atoms.loc[in_sphere, "mc"]) == {2.0}
    assert not (atoms.loc[~in_sphere, "mc"] == 2.0).any()

    # density x volume, to within the counting noise of one random draw
    sphere_volume = 4 / 3 * np.pi * 1.5**3
    assert in_sphere.sum() == pytest.approx(20 * sphere_volume, rel=0.15)
    assert len(atoms) == pytest.approx(20 * 1000, rel=0.02)


def test_composition_and_reproducibility(tiny, tmp_path):
    atoms = pgp.generate_synthetic_pos(tiny, seed=0)
    matrix = atoms.loc[atoms["mc"] != 2.0, "mc"]
    assert (matrix == 1.0).mean() == pytest.approx(0.1, abs=0.01)

    # the default output is the script's own <possave>, next to the script
    assert (tiny.parent / "tiny.pos").exists()

    again = pgp.generate_synthetic_pos(tiny, tmp_path / "again.pos", seed=0)
    pd.testing.assert_frame_equal(atoms, again)
    other = pgp.generate_synthetic_pos(tiny, tmp_path / "other.pos", seed=1)
    assert not atoms[["x", "y", "z"]].equals(other[["x", "y", "z"]])


def test_unsupported_element_raises(tmp_path):
    script = tmp_path / "rotate.xml"
    script.write_text(TINY_SCRIPT.replace(
        "<translate><x value=\"-5\"/><y value=\"-5\"/><z value=\"-5\"/></translate>",
        "<rotate angle=\"10\"/>",
    ))
    with pytest.raises(ValueError, match="rotate"):
        pgp.generate_synthetic_pos(script)


def test_easy_example_cube(tmp_path):
    """The script the example notebooks run on: a 40 nm cube at 30 atoms/nm^3."""
    atoms = pgp.generate_synthetic_pos(
        EXAMPLES / "generateTestPos_easy.xml", tmp_path / "easy.pos"
    )

    # posgen's own run of this script wrote 1,919,955 atoms, 11,919 Cu and
    # 12,010 Ni (the old repo's posgenpyTest_easy.pos)
    assert len(atoms) == pytest.approx(1_919_955, rel=0.001)
    counts = atoms["mc"].value_counts()
    assert counts[1.0] == pytest.approx(11_919, rel=0.05)
    assert counts[2.0] == pytest.approx(12_010, rel=0.05)

    positions = atoms[["x", "y", "z"]].to_numpy()
    assert positions.min() >= -20 and positions.max() <= 20
