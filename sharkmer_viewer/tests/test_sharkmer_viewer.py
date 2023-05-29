import sharkmer_viewer
import os


def test_test():
    # Make sure testing works
    assert 1 == 1


def test_sharkmer_viewer():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(current_dir, 'data', 'Cordagalma.histo')
    return_val = sharkmer_viewer.create_report(data_path)
    assert return_val == 0