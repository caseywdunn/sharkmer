import sharkmer_viewer
import os
import tempfile


def test_test():
    # Make sure testing works
    assert 1 == 1


def test_sharkmer_viewer():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    histo_path = os.path.join(current_dir, 'data', 'Cordagalma.histo')
    stats_path = os.path.join(current_dir, 'data', 'Cordagalma.stats')
    with tempfile.NamedTemporaryFile(suffix='.html', delete=False) as f:
        out_path = f.name
    try:
        return_val = sharkmer_viewer.create_report(
            histo_path, stats_path, out_path, 'Cordagalma', 0
        )
        assert return_val == 0
        assert os.path.exists(out_path)
    finally:
        if os.path.exists(out_path):
            os.unlink(out_path)