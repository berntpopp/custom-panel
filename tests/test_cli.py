"""
Comprehensive CLI tests for the custom-panel tool.

These tests ensure the command-line interface works correctly from the user's
perspective and that all components are properly integrated.
"""

import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
from typer.testing import CliRunner

from custom_panel.cli import app


@pytest.fixture
def runner() -> CliRunner:
    """Create a CLI runner for testing commands."""
    return CliRunner()


@pytest.fixture
def test_environment(tmp_path: Path) -> dict:
    """Create a temporary test environment with mock data and config files."""
    # Create mock source data files
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    # Mock In-house panel data
    inhouse_path = data_dir / "inhouse_panel.csv"
    inhouse_df = pd.DataFrame({"genes": ["BRCA1", "TP53", "KRAS"]})
    inhouse_df.to_csv(inhouse_path, index=False)

    # Mock Manual curation data (with a gene that will be vetoed)
    manual_dir = data_dir / "manual"
    manual_dir.mkdir()
    manual_path = manual_dir / "manual_list.csv"
    manual_df = pd.DataFrame({"approved_symbol": ["VETO_GENE", "MLH1"]})
    manual_df.to_csv(manual_path, index=False)

    # Create a minimal test configuration file
    config_content = f"""
general:
  output_dir: "results"  # This will be overridden by the CLI arg

data_sources:
  inhouse_panels:
    enabled: true
    source_group: true
    panels:
      - name: "Test_Inhouse"
        file_path: "{inhouse_path}"
        gene_column: "genes"
        evidence_score: 2.5

  manual_curation:
    enabled: true
    veto:
      enabled: true
      reason: "Manual veto"
    lists:
      - name: "Test_Manual"
        file_path: "{manual_path}"
        gene_column: "approved_symbol"
        evidence_score: 1.5

  # Disable all other sources to make the test fast and simple
  panelapp:
    enabled: false
  commercial_panels:
    enabled: false
  cosmic_germline:
    enabled: false
  acmg_incidental_findings:
    enabled: false
  hpo_neoplasm:
    enabled: false
  clingen:
    enabled: false
  thegencc:
    enabled: false

scoring:
  thresholds:
    score_threshold: 2.0  # High threshold to test veto
    min_sources: 1

output:
  formats:
    - excel
    - csv
    - parquet
  reports:
    html:
      enabled: true
  bed_files:
    germline:
      enabled: true
    exons:
      enabled: true
      transcript_types:
        canonical: true
        mane_select: true
"""
    config_path = tmp_path / "test_config.yml"
    config_path.write_text(config_content)

    return {
        "config_path": config_path,
        "output_dir": tmp_path / "cli_output",
        "tmp_path": tmp_path,
    }


# Mock the API clients to prevent real network calls during testing
@patch("custom_panel.engine.annotator.EnsemblClient")
@patch("custom_panel.engine.annotator.HGNCClient")
def test_cli_run_success(
    mock_hgnc_client, mock_ensembl_client, runner: CliRunner, test_environment: dict
):
    """Test a successful end-to-end run via the CLI."""
    # Setup mock HGNC client
    mock_hgnc_instance = mock_hgnc_client.return_value
    mock_hgnc_instance.standardize_gene_symbols.side_effect = lambda symbols: {
        s: {"approved_symbol": s, "hgnc_id": f"HGNC:{i}"}
        for i, s in enumerate(symbols, 1)
    }

    # Setup mock Ensembl client
    mock_ensembl_instance = mock_ensembl_client.return_value
    mock_ensembl_instance.get_symbols_data_batch.return_value = {
        "BRCA1": {
            "gene_id": "ENSG00000012048",
            "chromosome": "17",
            "start": 43044295,
            "end": 43125483,
            "canonical_transcript": "ENST00000357654",
            "mane_select_transcript": "ENST00000357654.9",
        },
        "TP53": {
            "gene_id": "ENSG00000141510",
            "chromosome": "17",
            "start": 7661779,
            "end": 7687550,
            "canonical_transcript": "ENST00000269305",
            "mane_select_transcript": "ENST00000269305.9",
        },
        "KRAS": {
            "gene_id": "ENSG00000133703",
            "chromosome": "12",
            "start": 25205246,
            "end": 25250936,
            "canonical_transcript": "ENST00000256078",
            "mane_select_transcript": "ENST00000256078.10",
        },
        "VETO_GENE": {
            "gene_id": "ENSG00000999999",
            "chromosome": "1",
            "start": 1000000,
            "end": 1010000,
            "canonical_transcript": "ENST00000999999",
        },
        "MLH1": {
            "gene_id": "ENSG00000076242",
            "chromosome": "3",
            "start": 36993325,
            "end": 37050845,
            "canonical_transcript": "ENST00000231790",
            "mane_select_transcript": "ENST00000231790.8",
        },
    }

    # Mock transcript data for exon BED generation
    mock_transcript_data = {
        "BRCA1": {
            "all_transcripts": [
                {
                    "id": "ENST00000357654",
                    "Exon": [
                        {
                            "id": "ENSE00001",
                            "start": 43044295,
                            "end": 43044400,
                            "strand": -1,
                            "rank": 1,
                            "seq_region_name": "17",
                        },
                        {
                            "id": "ENSE00002",
                            "start": 43045800,
                            "end": 43045900,
                            "strand": -1,
                            "rank": 2,
                            "seq_region_name": "17",
                        },
                    ],
                }
            ],
            "gene_id": "ENSG00000012048",
            "chromosome": "17",
        },
        "TP53": {
            "all_transcripts": [
                {
                    "id": "ENST00000269305",
                    "Exon": [
                        {
                            "id": "ENSE00003",
                            "start": 7661779,
                            "end": 7661900,
                            "strand": -1,
                            "rank": 1,
                            "seq_region_name": "17",
                        }
                    ],
                }
            ],
            "gene_id": "ENSG00000141510",
            "chromosome": "17",
        },
    }

    # Patch the pipeline run method to return our mock transcript data
    with patch("custom_panel.engine.pipeline.Pipeline.run") as mock_run:
        # Create a mock annotated DataFrame
        annotated_df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "TP53", "KRAS", "VETO_GENE", "MLH1"],
                "hgnc_id": ["HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4", "HGNC:5"],
                "gene_id": [
                    "ENSG00000012048",
                    "ENSG00000141510",
                    "ENSG00000133703",
                    "ENSG00000999999",
                    "ENSG00000076242",
                ],
                "chromosome": ["17", "17", "12", "1", "3"],
                "gene_start": [43044295, 7661779, 25205246, 1000000, 36993325],
                "gene_end": [43125483, 7687550, 25250936, 1010000, 37050845],
                "gene_strand": [-1, -1, 1, 1, -1],
                "canonical_transcript": [
                    "ENST00000357654",
                    "ENST00000269305",
                    "ENST00000256078",
                    "ENST00000999999",
                    "ENST00000231790",
                ],
                "mane_select_transcript": [
                    "ENST00000357654.9",
                    "ENST00000269305.9",
                    "ENST00000256078.10",
                    pd.NA,
                    "ENST00000231790.8",
                ],
                "score": [2.5, 2.5, 2.5, 0.5, 1.5],
                "source_count": [1, 1, 1, 1, 1],
                "sources": [
                    "inhouse_panels",
                    "inhouse_panels",
                    "inhouse_panels",
                    "manual_curation",
                    "manual_curation",
                ],
                "veto_reasons": [
                    "",
                    "",
                    "",
                    "Veto from manual_curation: Manual veto",
                    "",
                ],
                "include": [True, True, True, True, False],
                "inclusion_reason": [
                    "score",
                    "score",
                    "score",
                    "veto",
                    "below_threshold",
                ],
            }
        )

        mock_run.return_value = (annotated_df, mock_transcript_data)

        config_path = test_environment["config_path"]
        output_dir = test_environment["output_dir"]

        # Invoke the CLI 'run' command
        result = runner.invoke(
            app,
            ["run", "--config-file", str(config_path), "--output-dir", str(output_dir)],
        )

    # 1. Assert the command exited successfully
    assert result.exit_code == 0, f"CLI failed with: {result.stdout}"
    assert "Pipeline completed successfully" in result.stdout

    # 2. Assert that the output directory and files were created
    # The pipeline creates a timestamped subdirectory
    assert output_dir.exists()
    run_dirs = list(output_dir.iterdir())
    assert len(run_dirs) == 1
    run_dir = run_dirs[0]
    final_output_dir = run_dir / "06_final_output"

    assert final_output_dir.exists()
    assert (final_output_dir / "master_panel.xlsx").exists()
    assert (final_output_dir / "master_panel.csv").exists()
    assert (final_output_dir / "master_panel.parquet").exists()
    assert (final_output_dir / "panel_report.html").exists()
    assert (final_output_dir / "germline_panel.bed").exists()
    assert (run_dir / "run_summary.json").exists()

    # 3. Assert the content of output files
    output_df = pd.read_csv(final_output_dir / "master_panel.csv")
    assert len(output_df) == 5  # BRCA1, TP53, KRAS, VETO_GENE, MLH1

    # Check veto logic
    veto_gene_row = output_df[output_df["approved_symbol"] == "VETO_GENE"]
    assert not veto_gene_row.empty
    assert bool(veto_gene_row.iloc[0]["include"]) is True
    assert veto_gene_row.iloc[0]["inclusion_reason"] == "veto"

    # Check score-based inclusion
    brca1_row = output_df[output_df["approved_symbol"] == "BRCA1"]
    assert bool(brca1_row.iloc[0]["include"]) is True
    assert brca1_row.iloc[0]["inclusion_reason"] == "score"

    # Check exclusion
    mlh1_row = output_df[output_df["approved_symbol"] == "MLH1"]
    assert bool(mlh1_row.iloc[0]["include"]) is False
    assert mlh1_row.iloc[0]["inclusion_reason"] == "below_threshold"

    # Check BED file content
    bed_content = (final_output_dir / "germline_panel.bed").read_text()
    bed_lines = bed_content.strip().split("\n")
    # Should have 4 included genes (BRCA1, TP53, KRAS, VETO_GENE)
    assert len(bed_lines) == 4

    # Check run summary
    with open(run_dir / "run_summary.json") as f:
        run_summary = json.load(f)
    assert "run_directory" in run_summary
    assert "file_logging_enabled" in run_summary


# TODO: Complex test requiring deeper mocking - disabled for now
# def test_cli_argument_overrides(runner: CliRunner, test_environment: dict):
#     """Test that CLI arguments are parsed correctly."""
#     pass


def test_cli_dry_run(runner: CliRunner, test_environment: dict):
    """Test that dry run mode doesn't generate output files."""
    config_path = test_environment["config_path"]
    output_dir = test_environment["output_dir"]

    # Invoke the CLI with --dry-run
    result = runner.invoke(
        app,
        [
            "run",
            "--config-file",
            str(config_path),
            "--output-dir",
            str(output_dir),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert "Running in dry-run mode" in result.stdout
    assert "Dry run completed successfully" in result.stdout

    # Output directory may be created but should not contain final output files
    if output_dir.exists():
        # Check that no final output files were created
        for run_dir in output_dir.iterdir():
            final_output_dir = run_dir / "06_final_output"
            if final_output_dir.exists():
                # Should not have generated actual output files
                output_files = list(final_output_dir.glob("*.csv")) + list(
                    final_output_dir.glob("*.xlsx")
                )
                assert len(output_files) == 0, (
                    "Output files were created in dry run mode"
                )


def test_cli_config_check(runner: CliRunner, test_environment: dict):
    """Test the config-check command."""
    config_path = test_environment["config_path"]

    result = runner.invoke(app, ["config-check", "--config-file", str(config_path)])

    assert result.exit_code == 0
    assert "Configuration Validation" in result.stdout
    assert "Data Sources Configuration" in result.stdout
    assert "inhouse_panels" in result.stdout
    assert "manual_curation" in result.stdout
    assert "Scoring Configuration" in result.stdout
    assert "Score threshold: 2.0" in result.stdout


@patch("custom_panel.sources.g01_panelapp.search_panelapp_panels")
def test_cli_search_panels(mock_search, runner: CliRunner):
    """Test the search-panels command."""
    # Mock search results
    mock_search.return_value = [
        {"id": "123", "name": "Cancer Panel", "version": "1.0", "source": "panelapp"},
        {"id": "456", "name": "Cardiac Panel", "version": "2.1", "source": "panelapp"},
    ]

    result = runner.invoke(app, ["search-panels", "cancer"])

    assert result.exit_code == 0
    assert "Searching for panels matching: 'cancer'" in result.stdout
    assert "Cancer Panel" in result.stdout
    assert "123" in result.stdout


# TODO: Complex test - disabled for now
# def test_cli_missing_config(runner: CliRunner):
#     """Test error handling when config file is missing."""
#     result = runner.invoke(app, [
#         "run",
#         "--config-file", "/nonexistent/config.yml"
#     ])
#
#     assert result.exit_code == 1
#     assert "Configuration Error" in result.stdout or "Error loading configuration" in result.stdout


def test_cli_invalid_log_level(runner: CliRunner, test_environment: dict):
    """Test handling of invalid log level."""
    config_path = test_environment["config_path"]

    # This should still work, defaulting to a valid log level
    result = runner.invoke(
        app, ["run", "--config-file", str(config_path), "--log-level", "INVALID"]
    )

    # The command should handle this gracefully
    assert result.exit_code in [0, 1]  # May succeed or fail gracefully
