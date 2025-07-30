import gc
import ete3
import numpy as np

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


def accumulate_ed_scores(parent, running_sum, all_scores, get_dict=True, root=True):
    if not root:
        running_sum += parent.dist/parent.desc_leaves

    parent.add_feature("ed_score", running_sum)

    if parent.is_leaf() and get_dict:
        all_scores[(parent.name, parent.domain, parent.kingdom, parent.phylum, parent.clas, parent.order, parent.family)] = running_sum
    else:
        for child in parent.children:
            accumulate_ed_scores(child, running_sum, all_scores, get_dict, False)


def get_ed_scores(dated_tre, get_dict=True):
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
    if get_dict:
        all_scores = {}
        accumulate_ed_scores(dated_tre, 0, all_scores, True)

        return all_scores
    else:
        accumulate_ed_scores(dated_tre, 0, {}, False)


def add_ed_scores(dated_tre, existing_scores):
    logger.info("Computing ED scores.")
    new_scores = get_ed_scores(dated_tre)
    for key in new_scores:
        if key not in existing_scores:
            existing_scores[key] = []

        existing_scores[key].append(new_scores[key])


def write_ed_scores(filename, scores):
    logger.info("Writing ED score distributions.")
    fout = open(filename, "w")
    fout.write("leaf_name\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tmean\tmin\t25pct\tmedian\t75pct\tmax\tN\n")
    for key in scores:
        fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n" % (key[0],
                                                                                 key[1],
                                                                                 key[2],
                                                                                 key[3],
                                                                                 key[4],
                                                                                 key[5],
                                                                                 key[6],
                                                                                 np.mean(scores[key]),
                                                                                 np.min(scores[key]),
                                                                                 np.percentile(scores[key],25),
                                                                                 np.percentile(scores[key],50),
                                                                                 np.percentile(scores[key],75),
                                                                                 np.max(scores[key]),
                                                                                 len(scores[key])))
    fout.close()


def compute_pd(parent):
    """Compute total PD below each node. Assume the tree has branch lengths."""
    pd = 0
    for child in parent.children:
        pd += compute_pd(child)

    parent.add_feature("pd", pd)

    return pd + parent.dist


def save_pd_for_clades(tre, pd_clades, pd_dict, dates_dict, spp_dict):
    for node in tre.traverse():
        if node.name in pd_clades:
            pd_dict[node.name].append(node.pd)
            dates_dict[node.name].append(node.date)
            spp_dict[node.name].append(len(node.get_leaves()))


def write_pd_dists(filename, pd_dict, dates_dict, spp_dict):
    def list_to_tab_str(ls):
        ls_str = ""
        for item in ls[:-1]:
            ls_str += str(item)
            ls_str += '\t'
        ls_str += str(ls[-1])

        return ls_str

    fout = open("%s_phyla.txt" % filename, "w")
    fout.write("Clade\tspp.\tMean\tMin\t2.5pct\t16pct\t25pct\t50pct\t75pct\t84pct\t97.5pct\tMax\n")
    for clade in pd_dict:
        fout.write("%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n" % (clade,
                                                                             spp_dict[clade][0],
                                                                             np.mean(pd_dict[clade]),
                                                                             np.percentile(pd_dict[clade], 0),
                                                                             np.percentile(pd_dict[clade], 2.5),
                                                                             np.percentile(pd_dict[clade], 16),
                                                                             np.percentile(pd_dict[clade], 25),
                                                                             np.percentile(pd_dict[clade], 50),
                                                                             np.percentile(pd_dict[clade], 75),
                                                                             np.percentile(pd_dict[clade], 84),
                                                                             np.percentile(pd_dict[clade], 97.5),
                                                                             np.percentile(pd_dict[clade], 100),
                                                                             list_to_tab_str(pd_dict[clade])))

    fout.close()

    fout = open("%s_dates.txt" % filename, "w")
    fout.write("Clade\tMean\tMin\t2.5pct\t16pct\t25pct\t50pct\t75pct\t84pct\t97.5pct\tMax\n")
    for clade in dates_dict:
        fout.write("%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n" % (clade,
                                                                             spp_dict[clade][0],
                                                                             np.mean(dates_dict[clade]),
                                                                             np.percentile(dates_dict[clade], 0),
                                                                             np.percentile(dates_dict[clade], 2.5),
                                                                             np.percentile(dates_dict[clade], 16),
                                                                             np.percentile(dates_dict[clade], 25),
                                                                             np.percentile(dates_dict[clade], 50),
                                                                             np.percentile(dates_dict[clade], 75),
                                                                             np.percentile(dates_dict[clade], 84),
                                                                             np.percentile(dates_dict[clade], 97.5),
                                                                             np.percentile(dates_dict[clade], 100),
                                                                             list_to_tab_str(dates_dict[clade])))

    fout.close()


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
                fout.write("%f" % (x[0]))
                del tre2
                gc.collect()
            else:
                results[-1].append(None)
                if j == i:
                    fout.write("0")
                else:
                    fout.write("%f" % (results[j][i]))

            if j != n-1:
                fout.write(",")
            else:
                fout.write("\n")

        del tre1
        gc.collect()

    fout.close()


def date_labelling_guo(parent):
    """Recurse in postorder through the tree, labelling each node with the oldest date below it
    and the path length (number of nodes) to that date.
    Return value is: [oldest date found so far below this node, path length to it]
    """
    if parent.is_leaf():
        if parent.date != 0:
            print("Leaf node with non-zero date")
        oldest_paths = [[parent.date, 1]]
        parent.oldest_paths = oldest_paths
    else:
        oldest_paths = []
        for child in parent.children:
            child_paths = date_labelling_guo(child)
            for path in child_paths:
                oldest_paths.append(path)

        if parent.date is None:
            parent.oldest_paths = oldest_paths
            for path in oldest_paths:
                path[1] += 1
        else:
            oldest_paths = [[parent.date, 1]]
            parent.oldest_paths = oldest_paths

    return oldest_paths


def compute_ed_scores_guo(parent, ancestral_date, root=True):
    ed_scores = []

    if not root:
        for path in parent.oldest_paths:
            if path[1] != 0:
                ed_scores.append(((ancestral_date-path[0])/path[1]) / parent.desc_leaves)

    parent.add_feature("ed_scores_guo", ed_scores)

    for child in parent.children:
        if parent.date:
            compute_ed_scores_guo(child, parent.date, False)
        else:
            compute_ed_scores_guo(child, ancestral_date, False)


def sum_ed_scores_guo(parent, root=True):
    i = 0
    if not root:
        for child in parent.children:
            for j in range(len(child.ed_scores_guo)):
                if len(parent.ed_scores_guo) == 1:
                    child.ed_scores_guo[j] += parent.ed_scores_guo[0]
                else:
                    child.ed_scores_guo[j] += parent.ed_scores_guo[i]
                i += 1

    for child in parent.children:
        sum_ed_scores_guo(child, False)


def lineages_through_time(tre):
    """Assumes all nodes are dated and have 0, 1 or 2 children. If leaf nodes have a date of 0, they
    are extant species; if they have a date of more than zero, they are extinct."""

    dates_list = []
    nonimputed_dates = []
    for node in tre.traverse():
        if node.date > 0:
            dates_list.append((node.date, len(node.children)-1))

            if not node.imputed_date:
                nonimputed_dates.append(-node.date)

    dates_list.sort(reverse=True)

    lineages = [2]
    dates = [-dates_list[0][0]]
    for i in range(1,len(dates_list)):
        lineages.append(lineages[-1] + dates_list[i][1])
        dates.append(-dates_list[i][0])

    return dates, lineages, nonimputed_dates
