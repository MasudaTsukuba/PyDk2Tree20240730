"""
test_basic.py
2024/8/6, T. Masuda
Amagasa Laboratory, University of Tsukuba
"""
import pytest
from src.DKTree import DKTree


def test_basic():
    tree = DKTree()
    tree.insert_entries(10)
    tree.add_edge(1, 2)
    ans1 = tree.report_edge(1, 2)
    tree.remove_edge(1, 2)
    ans2 = tree.report_edge(1, 2)
    assert ans1
    assert not ans2
    pass

