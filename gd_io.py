from gd_utils import get_tree_from_str

def read_interval_file(filename):
    gt_list, st_list = [], []
    try:
        t = [line.strip() for line in open(filename, "r").readlines() if len(line.strip()) and line.strip()[0] != '#']
        gt_list, st_list = t[1:-1], t[-1]
    except IOError as e:
        print(filename, " I/O error({0}): {1}".format(e.errno, e.strerror))
        quit()
    return gt_list, st_list


def check(species, genetreelist):
    counter = 0
    for g in genetreelist:
        ok = False
        gt = get_tree_from_str(g)
        leaves = [x.name for x in gt.get_terminals()]
        if species[0] in leaves:
            ok = True
        if species[1] in leaves:
            ok = True
        if ok:
            counter = counter + 1
    return counter


if __name__ == "__main__":
    gtl, st = read_interval_file("tests/treefam.txt")
    stree = get_tree_from_str(st)
    species = [x.name for x in stree.get_terminals()]

    # for s in species:
    #     print(s, check([s], gtl))
    print(check(["ORYSA","ARATH"], gtl))