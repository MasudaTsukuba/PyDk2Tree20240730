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


def test_series():
    nr_of_nodes = 1000
    tree = DKTree()
    tree.insert_entries(1002)
    all_nodes = [i+1 for i in range(nr_of_nodes + 1)]
    for i in range(nr_of_nodes):
        tree.add_edge(i + 1, i + 2)
    ans1 = tree.report_all_edges(all_nodes, all_nodes)
    assert len(ans1) == nr_of_nodes
    ans2 = tree.report_edge(1000, 1001)
    assert ans2
    ans3 = tree.report_edge(1001, 1000)
    assert not ans3


def test_series2():
    nr_of_nodes = 1000
    tree = DKTree()
    tree.insert_entries(1002)
    all_nodes = [i+1 for i in range(nr_of_nodes + 1)]
    for i in range(nr_of_nodes):
        tree.add_edge(i + 1, i + 1)
    ans1 = tree.report_all_edges(all_nodes, all_nodes)
    assert len(ans1) == nr_of_nodes
    ans2 = tree.report_edge(1000, 1000)
    assert ans2
    ans3 = tree.report_edge(1001, 1000)
    assert not ans3
