import gurobipy as gp
from gurobipy import *
import pandas as pd
from gd_utils import get_tree_from_str, get_chains_for_sp_node, print_clade, process_gt


def run_optimizer(stree, gf, interval_limit=0, max_converted_spec=-1, time_limit=-1, outfile="out", debug=0):
    solution, support = 0, {}
    try:
        dno = 0  # enumerate duplications; current number of duplication <=> interval
        dlist = []  # list of duplication nodes (order as above)
        dcomp = []  # list - for a corresponding node - list of comparable duplication numbers
        dint = []  # list - for a corresponding node - interval - list of allowed nodes in species tree
        issp = []  # list - for a corresponding node - True if speciation
        lcam = []  # list - for a corresponding node - index of a st node; its lca-mapping
        childr = []  # list - for a corresponding node - list of dno indexes of children
        exclude = []  # list - for a corresponding node - 1 if exclude; speciation not above duplication
        # to manage the support of a genomic duplication by gene trees
        mapfgtid = []  # list - for a corresponding node - id of a gene tree
        gtid = 0  # id of a gene tree; sets starting id of a gene tree

        # Species tree from input into constraints
        st = get_tree_from_str(stree)
        n = len(st.get_terminals())
        m = sum(len(get_tree_from_str(gt).get_terminals()) for gt in gf)
        maxdup = m - len(gf)
        print("Size of input data: n =", n, ", m =", m, ", maxdup =", maxdup)

        # Create a new model
        m = gp.Model("FHS model ME clustering")

        # ************* Create variables ****************
        # dup          - every variable corresponds to a species tree node
        # dup[i] = k   - will be equal to ME_score = k for that node
        dup = m.addVars(range(2 * n - 1), vtype=GRB.INTEGER, name="MEscore for node")

        # mapf         - every variable corresponds to an internal gene tree node
        # mapf[i] = k  - node i from gene tree is mapped to a node k in the species tree
        #    The position is represented by species tree node number, which is assigned in level-order.
        #    The root node number is 0, and any ancestor number is lesser than the number of its descendant.
        mapf = m.addVars(maxdup, vtype=GRB.INTEGER, name="mapping")

        # inttab[i,j]     - intervals are stored here
        # i - potdup      - for every duplication there is a corresponding interval - in total potdup intervals
        # j - 2*n-1       - interval spans over a subset of species tree nodes; total = 2*n-1;
        # inttab[i,j] = 1 - if node j in species tree S is in an interval of a duplication i
        inttab = m.addVars(maxdup, 2 * n - 1, vtype=GRB.BINARY, name="interval")

        # chain                  - every variable corresponds to a chain that corresponds to a species tree node
        # chain[i,...            - for a species tree node i
        #                        - variable value = sum( (inttab[d, i]) for d in chain);
        #                          duplications d from chain are comparable;
        #                          inttab[d, i] == 1 iff d is mapped to i
        #                        - there can be at most dupno different chains
        # get_chains_for_sp_node - returns n chains of duplications
        # chain[i,j] = k         - if j < n; k is a number of duplications from a chain j mapped to i
        #                          otherwise k==0
        # max_ function can be applied to only variables; thus chain variables have to be created to compute MEscore;
        chain = m.addVars(range(2 * n - 1), maxdup, vtype=GRB.INTEGER, name="chain")

        # variable to remove speciations that are not convereted into duplications from duplication count
        specrem = m.addVars(maxdup, vtype=GRB.BINARY, name="not converted speciation")
        # ifspecmaptolca[x]=0 if speciation x is mapped to its lca-mapping, ifspecmaptolca[x]<>0 if not
        ifspecmaptolca = m.addVars(maxdup, vtype=GRB.INTEGER, name="if mapped to lca")
        # ifspecmaptolcadec[x]=1 if speciation x is mapped to its lca-mapping, ifspecmaptolcadec[x]=0 if not
        ifspecmaptolcadec = m.addVars(maxdup, vtype=GRB.BINARY, name="if mapped to lca dec")

        # compare lca mapping # NOTE: error for <0 values (!)
        iffchtolca = m.addVars(maxdup, vtype=GRB.INTEGER, name="if first child mapped to my lca")
        # iffchtolcadec[x]=1 if 1st child of speciation x is mapped to its lca-mapping;iffchtolcadec[x]=0 if not
        iffchtolcadec = m.addVars(maxdup, vtype=GRB.BINARY, name="child1 if mapped to lca dec")
        iffchtolcaandp = m.addVars(maxdup, vtype=GRB.BINARY, name="child1 of x and x mapped to lca of x")

        # and for the second child (we assume binary trees)
        ifschtolca = m.addVars(maxdup, vtype=GRB.INTEGER, name="if second child mapped to my lca")
        ifschtolcadec = m.addVars(maxdup, vtype=GRB.BINARY, name="child2 if mapped to lca dec")
        ifschtolcaandp = m.addVars(maxdup, vtype=GRB.BINARY, name="child2 of x and x mapped to lca of x")

        # variables to limit number of converted speciations
        ifconverted = m.addVars(maxdup, vtype=GRB.BINARY, name="if spec converted to dup")

        # ************* Objective        ****************
        # Set objective
        m.setObjective(dup.sum(), GRB.MINIMIZE)

        # ************* Constraints      ****************

        if debug == 3:
            print("Data processing")

        # constraints for monotonicity -
        for g in gf:
            (l, dno, dlist, dcomp, dint, issp, lcam, childr, exclude) = \
                process_gt(g, st, dno, dlist, dcomp, dint, issp, lcam, childr, exclude, interval_limit)
            for lelem in l:  #
                mapfgtid.append(gtid)
            gtid = gtid + 1

        if debug == 3:
            print("Interval constraints")

        # constraints for intervals    - duplication is placed at exactly one node from all interval nodes
        for i in range(dno):
            # constraint to control that duplication is placed at exactly one node from all interval nodes
            m.addConstr(sum((inttab[i, k] for k in dint[i])) == 1, "interval" + str(i))
            # and it cannot be placed outside an interval
            m.addConstr(sum((inttab[i, k] for k in range(2 * n - 1) if k not in dint[i])) == 0,
                        "outside of an interval" + str(i))

        if debug == 3:
            print("Monotonicity constraints")

        # constraints for monotonicty
        for i in range(dno):
            m.addConstr(mapf[i] == sum((k * inttab[i, k] for k in dint[i])), "position" + str(i))
            for j in dcomp[i]:  # monotonicty apply only to all comparable duplications
                if dlist[j].is_parent_of(dlist[i]):
                    m.addConstr(mapf[j] <= mapf[i], "monotonicity" + str(j) + " " + str(i))

        #         # FIXME: for t3 test to force optimal solution
        #         m.addConstr(inttab[5,0] == 1, name="no penalty")

        if debug == 3:
            print("Not-converted speciation constraints")
        # set constraints for speciations that are not converted
        for i in range(dno):
            if not issp[i]:  # a duplication - no penalty
                m.addConstr(specrem[i] == 0, name="no penalty, duplication")
                m.addConstr(ifconverted[i] == 0, name="no conversion, duplication")
            elif exclude[i]:  # speciation that cannot be converted into duplication
                # below constraints are enough; this is to simplify the model
                m.addConstr(specrem[i] == 1, name="penalty, speciation cannot be converted")
                m.addConstr(ifconverted[i] == 0, name="no conversion")
            else:
                # for a speciation v
                # if v is not mapped into lca-mapping of v, then v is a duplication;
                #     (mapf[i] - lcam[i]) = 0 => mapped to lca-mapping
                #     note (!) : lcam[i] >= mapf[i] from the definition of mapping
                m.addConstr(ifspecmaptolca[i] == (lcam[i] - mapf[i]), name="lca map check")
                m.addConstr((ifspecmaptolcadec[i] == 0) >> (ifspecmaptolca[i] >= 1), name="lca map no")
                #  >> (ifspecmaptolca[i] * 2 >= 1)
                m.addConstr((ifspecmaptolcadec[i] == 1) >> (ifspecmaptolca[i] <= 0), name="lca map yes")
                #  >> (ifspecmaptolca[i] * 2 <= 1)

                # v is converted into duplication; (mapf[i] - lcam[i]) <> 0
                m.addConstr((ifspecmaptolcadec[i] == 0) >> (specrem[i] == 0), name="no penalty, lca of v")
                m.addConstr((ifspecmaptolcadec[i] == 0) >> (ifconverted[i] == 1), name="conversion")

                # check if child of v is mapped into lca-mapping of v--
                if len(childr[i]) == 2:
                    # first
                    m.addConstr(iffchtolca[i] == (mapf[childr[i][0]] - mapf[i]), name="child1 lca map check")
                    m.addConstr((iffchtolcadec[i] == 0) >> (iffchtolca[i] * 2 >= 1), name="child1 lca map no")
                    m.addConstr((iffchtolcadec[i] == 1) >> (iffchtolca[i] * 2 <= 1), name="child1 lca map yes")

                    m.addConstr(iffchtolcaandp[i] == and_(iffchtolcadec[i], ifspecmaptolcadec[i]),
                                name="no conversion, child1 enforce duplication")
                    m.addConstr((iffchtolcaandp[i] == 1) >> (specrem[i] == 0), name="no penalty, child1")
                    m.addConstr((iffchtolcaandp[i] == 1) >> (ifconverted[i] == 1), name="conversion")
                    # second
                    m.addConstr(ifschtolca[i] == (mapf[childr[i][1]] - mapf[i]), name="child2 lca map check")
                    m.addConstr((ifschtolcadec[i] == 0) >> (ifschtolca[i] * 2 >= 1), name="child2 lca map no")
                    m.addConstr((ifschtolcadec[i] == 1) >> (ifschtolca[i] * 2 <= 1), name="child2 lca map yes")

                    m.addConstr(ifschtolcaandp[i] == and_(ifschtolcadec[i], ifspecmaptolcadec[i]),
                                name="no conversion, child2 enforce duplication")
                    m.addConstr((ifschtolcaandp[i] == 1) >> (specrem[i] == 0), name="no penalty, child2")
                    m.addConstr((ifschtolcaandp[i] == 1) >> (ifconverted[i] == 1), name="conversion")

                elif len(childr[i]) == 1:  # only one non-terminal child
                    # NOTE: cannot use abs_(mapf[childr[i][0]] - lcam[i]) in ILP model
                    # lcam[i] replace with mapf[i] and (ifspecmaptolcadec[i] == 1)
                    m.addConstr(iffchtolca[i] == (mapf[childr[i][0]] - mapf[i]), name="child1 lca map check")
                    m.addConstr((iffchtolcadec[i] == 0) >> (iffchtolca[i] * 2 >= 1), name="child1 lca map no")
                    m.addConstr((iffchtolcadec[i] == 1) >> (iffchtolca[i] * 2 <= 1), name="child1 lca map yes")

                    m.addConstr(iffchtolcaandp[i] == and_(iffchtolcadec[i], ifspecmaptolcadec[i]),
                                name="no conversion, child1 enforce duplication")
                    m.addConstr((iffchtolcaandp[i] == 1) >> (specrem[i] == 0), name="no penalty, child1")
                    m.addConstr((iffchtolcaandp[i] == 1) >> (ifconverted[i] == 1), name="conversion")

                # -- if yes, then v is a duplication

        if debug == 3:
            print("Chain constraints")

        # constrains to set the values of duplication chain variables according to intervals (inttab)
        if debug == 2:
            print("dcomp ", dcomp)
        for i in range(2 * n - 1):
            chains = get_chains_for_sp_node(i, maxdup, dcomp, dint, exclude)
            if debug == 2:
                print("chains ", i, chains)
            c = 0
            for el in chains:
                # NOTE: if the speciation is not converted into duplication then it is at the top of the chain
                m.addConstr(chain[i, c] == sum((inttab[k, i] - specrem[k] * inttab[k, i]) for k in el),
                            name="chain" + str(i) + "," + str(c))
                c = c + 1
            for j in range(c, maxdup):
                m.addConstr(chain[i, j] == 0, name="chain" + str(i) + "," + str(j))

        # constrains to compute k - the MEscore of single species tree node i
        #                         - k is equal to the length of maximal chain of comparable duplications mapped to i
        for i in range(2 * n - 1):
            m.addConstr(dup[i] == max_(chain[i, j] for j in range(maxdup)), name="k for node(" + str(i) + ")")

        if max_converted_spec > -1:
            m.addConstr(sum(ifconverted[i] for i in range(dno)) <= max_converted_spec, name="conversion Limit")

        if debug == 3:
            print("Begin optimize")
        if debug == 3:
            m.write("out.lp")

        # Optimize model
        if not time_limit == -1:
            m.Params.TimeLimit = time_limit
        m.optimize()

        if m.objVal == float('+inf'):   # no feasible solution (due to termination)
            return solution, support

        if debug == 3:
            print(dup)

        # ******************* creating the output ***********************
        # OLD
        # # calculates the support of duplication episodes by gene trees (number of gene trees)
        # d = {}  # key - species tree node, value - set of supporting gene trees
        # for jj in range(dno):
        #     if not specrem[jj].x:  # if specrem = 1 then we have a speciation
        #         if mapf[jj].x not in d:
        #             d[mapf[jj].x] = {mapfgtid[jj]}
        #         else:
        #             d[mapf[jj].x].add(mapfgtid[jj])
        # if counter in d:
        #     data_list.append((print_clade(clade), dup[counter].x, len(d[counter])))
        # else:
        #     data_list.append((print_clade(clade), dup[counter].x, math.nan))

        max_chain_bound = int(m.objVal)
        dd = {}  # key - species tree node, value - list of sets of supporting gene trees along the chain
        for i in range(2 * n - 1):
            chains = get_chains_for_sp_node(i, maxdup, dcomp, dint, exclude)
            if i not in dd:
                dd[i] = []
            for el in chains:                    # el is a duplication chain from a gene tree
                c = 0
                for k in el:                     # k is a duplication in a gene tree
                    if mapf[k].x == i:           # if duplication k is mapped to a species tree node i
                        if not specrem[k].x:     # we have to exclude top nodes if they are not converted speciations
                            if len(dd[i]) > c:
                                dd[i][c].add(mapfgtid[k])    # add support (gt of k) for chain height (c)
                                c = c + 1
                            elif len(dd[i]) == c:
                                dd[i].append({mapfgtid[k]})  # insert support (gt of k) for new height (c)
                                c = c + 1
        print("Gene tree support:")
        print(dd)

        solution, support = m.objVal, dd

        # prints the duplication episodes assigned to particular nodes
        counter = 0
        data_list = []
        for clade in st.find_clades(order='level'):
            tmp_list = [print_clade(clade), dup[counter].x]
            for n in range(max_chain_bound):
                if n < len(dd[counter]):
                    tmp_list.append(len(dd[counter][n]))
                    tmp_list.append(dd[counter][n])
                else:
                    tmp_list.append(math.nan)
                    tmp_list.append(math.nan)
            data_list.append(tmp_list)
            counter = counter + 1

        column_list = ["Node Cluster", "ME score"]
        for n in range(max_chain_bound):
            column_list.append("Level " + str(n) + " support")
            column_list.append("Level " + str(n) + " trees")

        df = pd.DataFrame(data_list, columns=column_list)
        df.to_csv(outfile + ".csv")

        # ****************************************************************

        print('Total MEscore: %g' % m.objVal)

        if debug > 1:
            for v in m.getVars():
                if "MEscore" in v.varName:
                    print('%s %g' % (v.varName, v.x))
                if debug == 3:
                    print('%s %g' % (v.varName, v.x))

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError:
        print('Encountered an attribute error')

    return solution, support

