from Bio import Phylo
from gd_utils import get_tree_from_str


def process_lengths_tree(tree):
    out_l = []
    for clade in tree.find_clades():
        clade.branch_length = None
        clade.confidence = None
    Phylo.write(tree, "tmp.txt", "newick")
    with open("tmp.txt", "r", encoding="utf-8") as f_in:
        for el in f_in:
            out = el.replace(":0.00000", "")
            out_l.append(out)
    return out_l


def get_sim_st(in_path):
    with open(in_path) as f:
        return f.readlines()


def trim_sim_st(in_path):
    st_tree = Phylo.read(in_path + "s_tree.trees", "newick")
    # species = set([])
    # for leaf in st_tree.get_terminals():
    #     species.add(leaf.name)
    st = process_lengths_tree(st_tree)[0]
    return st


def prepare_gt_trees(in_path, n):
    val = ""
    gt = []
    for i in range(1, n):
        if i <= 9:
            val = "00"+str(i)
        elif i <= 99:
            val = "0"+str(i)
        else:
            val = str(i)
        gt_tree = Phylo.read(in_path + "g_trees" + val + ".trees", "newick")
        for leaf in gt_tree.get_terminals():
            leaf.name = leaf.name.split("_")[0]
        gt.append(process_lengths_tree(gt_tree)[0])
    return gt


if __name__ == "__main__":

    # path = "wgd-simulated/n20-wgd-1/1/"
    # path = "wgd-simulated/n20-wgd-2/1/"
    # path = "wgd-simulated/n20-wgd-3/1/"
    # path = "wgd-simulated/n20-wgd-4/1/"
    path = "wgd-simulated/n20-wgd-5/1/"
    # st_path = "wgd-simulated/s_tree.trees"
    # st_tree = Phylo.read(st_path, "newick")
    # st = process_lengths_tree(st_tree)[0]
    # # sp tree for simphy in n20-wgd-1/1/ consist of a duplicated part ((1_1,10_1),3_1),((1_2,10_2),3_2)) => ((1,10),3)
    st = "((11,((4,(((7,8),9),((1,10),3))),(((14,12),2),((6,13),(15,5))))),(((18,17),19),(16,20)))"

    gt = prepare_gt_trees(path, 100)

    with open("tests/sim5.txt", "w", encoding="utf-8") as f_out:
        f_out.write("100\n")
        for g in gt:
            f_out.write(g)
        f_out.write("#st\n")
        f_out.write(st)



