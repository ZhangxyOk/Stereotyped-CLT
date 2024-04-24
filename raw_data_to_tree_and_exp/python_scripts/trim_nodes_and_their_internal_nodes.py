import sys 
import os
from ete3 import Tree
import argparse

def Parsers():
    parser = argparse.ArgumentParser(description="Keep nodes, remove other nodes and their internal nodes")
    parser.add_argument('-t', '--tree', type=str, nargs='?', required=True, help="tree file path")
    parser.add_argument('-n', '--nodes', type=str, nargs='?',required=True, help="nodes. seperate by , ")
    parser.add_argument('-o', '--out', type=str, nargs='?', required=True, help="out tree file")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    options = Parsers()
    tree = Tree(options.tree)
    node_list = list(options.nodes.strip().split(","))
    tree.prune(node_list)
    tree.write(outfile=options.out)