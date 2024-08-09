"""
test_file.py
2024/8/6, T. Masuda
Amagasa Laboratory, University of Tsukuba
"""
from src.DKTree import DKTree


def test_file():
    number_of_nodes = 11
    all_nodes = [i for i in range(number_of_nodes)]
    tree = DKTree()
    tree.insert_entries(number_of_nodes + 1)
    with open('sample_input_20240730.txt', 'r') as file:
        for line in file:
            values = line.split(' ')
            a = int(values[0])
            b = int(values[1])
            tree.add_edge(a, b)
    ans1 = tree.report_edge(1, 2)
    assert ans1

    ans2 = tree.report_edge(2, 1)
    assert not ans2

    findings = tree.report_all_edges(all_nodes, all_nodes)
    assert len(findings) == 9

    tree.remove_edge(1, 2)
    ans3 = tree.report_edge(1, 2)
    assert not ans3

    findings = tree.report_all_edges(all_nodes, all_nodes)
    assert len(findings) == 8
    pass

