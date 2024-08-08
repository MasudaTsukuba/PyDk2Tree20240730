"""
MyTest20240806.py
2024/8/6, T.Masuda
Amagasa Laboratory, University of Tsukuba
"""

from src.DKTree import DKTree


def my_test():
    tree = DKTree()
    tree.insert_entries(10)
    tree.add_edge(1, 2)
    ans1 = tree.report_edge(1, 2)
    tree.remove_edge(1, 2)
    ans2 = tree.report_edge(1, 2)
    pass


if __name__ == '__main__':
    my_test()
