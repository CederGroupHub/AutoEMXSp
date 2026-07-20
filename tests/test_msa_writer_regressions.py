#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regression tests for MSA writing and acquisition-time spectrum persistence."""

from pathlib import Path

from autoemx.utils.legacy.spectrum_pointer_writer import write_spectrum_pointer_file


def test_write_spectrum_pointer_file_preserves_template_tail(tmp_path: Path):
    template_path = tmp_path / "EM_metadata.msa"
    template_path.write_text(
        "".join(
            [
                "#FORMAT      : EMSA/MAS Spectral Data File\n",
                "#VERSION     : 1.0\n",
                "#NPOINTS     : 2\n",
                "#LIVETIME    : 1.00000000\n",
                "#REALTIME    : 2.00000000\n",
                "#SPECTRUM\n",
                "1.0\n",
                "2.0\n",
                "#ENDOFDATA   : \n",
                "TRAILER_LINE\n",
            ]
        ),
        encoding="utf-8",
    )
    output_path = tmp_path / "written.msa"

    write_spectrum_pointer_file(
        str(output_path),
        spectrum_vals=[10.0, 20.0, 30.0],
        xperchan=0.01,
        offset=-0.03,
        sample_result_dir=str(tmp_path),
        template_filename="EM_metadata.msa",
        live_time=7.5,
        real_time=8.5,
    )

    written_text = output_path.read_text(encoding="utf-8")

    assert "#ENDOFDATA   :" in written_text
    assert "TRAILER_LINE" in written_text
    assert "10.0\n" in written_text
    assert "20.0\n" in written_text
    assert "30.0\n" in written_text
    assert "#NPOINTS     : 3\n" in written_text
    assert "#LIVETIME    : 7.50000000\n" in written_text
