import ete3
import numpy

import logging
logger = logging.getLogger(__name__)


def num_desc_leaves(parent):
    if parent.is_leaf():
        parent.add_feature("desc_leaves", 1)
        return 1
    else:
        leaves = 0
        for child in parent.children:
            leaves += num_desc_leaves(child)
        parent.add_feature("desc_leaves", leaves)
        return leaves


def accumulate_ed_scores(parent, running_sum, all_scores, root=True):
    if not root:
        running_sum += parent.dist/parent.desc_leaves

    parent.add_feature("ed_score", running_sum)

    if parent.is_leaf():
        all_scores[(parent.name, parent.domain, parent.kingdom, parent.phylum, parent.clas, parent.order, parent.family)] = running_sum
    else:
        for child in parent.children:
            accumulate_ed_scores(child, running_sum, all_scores, False)


def get_ed_scores(dated_tre):
    # first, compute branch lengths
    for node in dated_tre.traverse():
        if node.up:
            if node.is_leaf():
                this_date = 0
            else:
                this_date = node.date

            parent_date = node.up.date

            node.dist = parent_date - this_date

    # then, label each node with the number of leaves below it
    num_desc_leaves(dated_tre)

    # finally, recurse through tree computing ED scores
    all_scores = {}
    accumulate_ed_scores(dated_tre, 0, all_scores)

    return all_scores


def add_ed_scores(dated_tre, existing_scores):
    logger.info("Computing ED scores.")
    new_scores = get_ed_scores(dated_tre)
    for key in new_scores:
        if key not in existing_scores:
            existing_scores[key] = []

        existing_scores[key].append(new_scores[key])


def write_ed_scores(scores, filename):
    logger.info("Writing ED score distributions.")
    fout = open(filename, "w")
    fout.write("leaf_label,domain,kingdom,phylum,class,order,family,mean,min,25pct,median,75pct,max,N\n")
    for key in scores:
        fout.write("%s,%s,%s,%s,%s,%s,%s,%f,%f,%f,%f,%f,%f,%d\n" % (key[0] if "," not in key[0] else ('"' + key[0] + '"'),
                                                        key[1],
                                                        key[2],
                                                        key[3],
                                                        key[4],
                                                        key[5],
                                                        key[6],
                                                        numpy.mean(scores[key]),
                                                        numpy.min(scores[key]),
                                                        numpy.percentile(scores[key],25),
                                                        numpy.percentile(scores[key],50),
                                                        numpy.percentile(scores[key],75),
                                                        numpy.max(scores[key]),
                                                        len(scores[key])))


def compute_rf_distances(output_folder, output_tree_filename, n):
    logger.info("Computing Robinson-Foulds distances.")

    fout = open("%s/rf_distances.csv" % (output_folder),"w")

    results = []
    for i in range(n):
        tre1 = ete3.Tree(open("%s/%s_%d.tre" % (output_folder, output_tree_filename, i+1),"r").read(),format=1)
        results.append([])
        for j in range(n):
            if j > i:
                logger.info("Computing RF distance between tree %d and tree %d." % (i+1,j+1))
                tre2 = ete3.Tree(open("%s/%s_%d.tre" % (output_folder, output_tree_filename, j+1),"r").read(),format=1)
                x = tre1.robinson_foulds(tre2)
                results[-1].append(x[0])
            else:
                results[-1].append(None)


    for i in range(n):
        for j in range(n):
            if results[i][j]:
                fout.write("%f" % (results[i][j]))
            else:
                fout.write("0")

            if j != n-1:
                fout.write(",")
            else:
                fout.write("\n")
