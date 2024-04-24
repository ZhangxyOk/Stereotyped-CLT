import os
import sys
import argparse
from ete3 import Tree

def create_balanced_tree_from_scratch(nodes):
    t = Tree()
    all_leaves = 0
    while all_leaves < nodes:
        for each_leaf in t.get_leaves():
            each_leaf.add_child()
            each_leaf.add_child()
            all_leaves = len(t.get_leaves())
            #print(all_leaves)
            if all_leaves > nodes:
                break
      ## add node name
    i = 1
    for l in t.get_leaves():
        l.name = "node_" + str(i)
        i += 1
    return t


def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nodes', type=int, nargs='?', required=True, help="Input the number of leaves need")
    parser.add_argument('-o', '--out', type=str, nargs='?', required=True, help="PATH of the output tree file")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    options = Parsers()
    out_tree = create_balanced_tree_from_scratch(options.nodes)
    out_tree.write(outfile=options.out)