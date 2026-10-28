import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "01_prepare_mammals.py"
spec = importlib.util.spec_from_file_location("prepare_mammals", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_cytb_metric_counts_target_codons():
    sequence = "CCTCCCTTTTTCTTATTG"
    result = module.cytb_metric(sequence)

    assert result["pro_count"] == 2
    assert result["phe_count"] == 2
    assert result["leu_tt_count"] == 2
    assert result["denominator_count"] == 4
    assert result["pro_over_phe_plus_leutt"] == 0.5
