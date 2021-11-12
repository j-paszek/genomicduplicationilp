#!/usr/bin/env python
# coding: utf-8

from Bio import Phylo
from io import StringIO


def get_tree_from_str(tree):  # input tree in string format
    return Phylo.read(StringIO(tree), "newick")


# get number of a node; to address the matrix in the model
def get_clade_number(node, st):
    num_clade_list = [(i, j) for i, j in enumerate(st.find_clades(order='level'))]
    ll = [el for el in num_clade_list if el[1] == node]
    return ll[0][0]


# input:  node (from gene tree), species tree
# output: least-common ancestor mapping of a (gt) node in a species tree (node in st)
def lca(node, st):
    tmp_list = []  # FIXME: to properly use biopython common_ancestor
    for el in node.get_terminals():  # very wild conversion is done
        dd = {"name": el.name}  # instead of inputting a list of leaves
        tmp_list.append(dd)  # we input a list of corresponding dict with leave names
    return st.common_ancestor(*tmp_list)


# input : node, species tree
# output: True if node is a duplication
def is_duplication(node, st):
    for child in node:
        if lca(child, st) == lca(node, st):
            return True
    return False


# input : node from gene tree, species tree
# output: True if node is a speciation that is an ancestor of a duplication
def is_an_ancestor_of_a_duplication(node, st):
    for ino in node.get_nonterminals():
        if is_duplication(ino, st):
            return True
    return False


# input:  gene tree node, and a species tree
#         interval_limit - to limit movement of nodes
# output: list that corresponds to an interval of possible mappings in FHS model (nodes)
def get_interval(gt_node, st, interval_limit=0):
    bottom = lca(gt_node, st)  # lca is a bottom interval end
    if is_duplication(gt_node, st) or is_an_ancestor_of_a_duplication(gt_node, st):  # upper end is root
        interval = st.get_path(bottom)[::-1]  # FIXME : get_path excludes root
        interval.append(st.root)  # see above
        excl = 0
    else:  # NOTE: We can limit mappings by not considering variants for speciations under duplications
        interval = [bottom]
        excl = 1  # exclude
    # test_int(interval)
    interval_numbers = [get_clade_number(el, st) for el in interval]
    if interval_limit > 0:  # if parameter is set we bound the interval
        if interval_limit < len(interval_numbers):  # if interval is longer than the bound
            # node numbers ascend from 0-the root
            return sorted(interval_numbers)[-interval_limit:], excl
    return interval_numbers, excl


# input:  gene tree, and a species tree (trees)
# output: list of lists; every list corresponds to an interval
def intervals_for_gt(gt, st, interval_limit=0):
    ino = gt.get_nonterminals()
    odp = []
    for gt_node in ino:
        interval, excl = get_interval(gt_node, st, interval_limit)
        odp.append((interval, gt_node, excl))
    return odp


# input:  two duplication nodes,
# output: true if those nodes are comparable duplications
def are_comparable_fhs(n1, n2):
    a = n2.is_parent_of(n1)
    b = n1.is_parent_of(n2)
    if not (a | b):  # check if nodes are comparable
        return False
    return True


# input:  species tree node index i; number of duplications; and two lists[j], that for a duplication j consist of:
#         dcomp - list of numbers of comparable duplications
#         dint  - interval - list of allowed nodes in a species tree
# output: list of chains; chain is a list of comparable duplications
def get_chains_for_sp_node(i, dupno, dcomp, dint, exclude):
    ret = []
    for j in range(dupno):
        if not exclude[j]:
            # we ignore chains that start in a speciation not above any duplication
            # we start that chain from its ancestor duplication
            # simplifies model by limiting the number of chains; with all chains also works
            if i in dint[j]:  # if not, duplication j cannot be mapped to node i
                if not dcomp[j]:  # if there are no comparable duplications we return only duplication j
                    ret.append([j])
                else:
                    # exclude speciations that can not be converted into duplications
                    dup_l = [m for m in dcomp[j] if not exclude[m]]
                    # (+) exclude nodes that are not mapped into i - due to limited interval model extension
                    # in model extension leaf duplication mapping can be limited so it is possible
                    # that i is not in dint[k]
                    ll = [k for k in dup_l if i in dint[k]]
                    lll = [k for k in ll if k > j]  # to check if is a leaf duplication (see below)
                    if not lll:  # add current leaf duplication j to the chain of comparable duplications
                        ret.append(ll + [j])  # l instead of dcomp[j] simplifies model
    return ret


def print_clade(clade):
    lll = [str(x) for x in clade.get_terminals()]
    return ';'.join(lll)


def process_gt(gtree, st, dno, dlist, dcomp, dint, issp, lcam, childr, exclude, interval_limit=0):
    gt = get_tree_from_str(gtree)
    int_l = intervals_for_gt(gt, st, interval_limit)
    start = dno
    for (subl, node, excl) in int_l:
        dint.append(subl)  # elements are added on position dno
        dlist.append(node)
        exclude.append(excl)
        dcomp.append([])  # we put empty list - further filled with comparable duplication numbers
        issp.append(not is_duplication(node, st))
        lcam.append(get_clade_number(lca(node, st), st))
        dno = dno + 1
    # comparable intervals are only in the same gene tree
    for j in range(start, dno - 1):
        for k in range(j + 1, dno):
            if are_comparable_fhs(dlist[j], dlist[k]):
                dcomp[j].append(k)
                dcomp[k].append(j)
    # set childr - list of non-terminal children - element = id of a node
    for j in range(start, dno):
        childclades = [clade for clade in int_l[j - start][1].clades if not clade.is_terminal()]
        nlist = [x[1] for x in int_l]  # nodes in gt - dup or sp/dup
        childr.append([start + nlist.index(c) for c in childclades])
    # note: childr,... - numbers of nodes in gene tree - pre-order
    return int_l, dno, dlist, dcomp, dint, issp, lcam, childr, exclude
