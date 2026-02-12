# BSD 3-Clause License

# Copyright (c) 2025, Jonathan David Duke

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import gc
import ete3
import numpy as np

import logging
logger = logging.getLogger(__name__)


def label_desc_leaves(parent):
    """Label each node with the number of leaf nodes below it."""

    if parent.is_leaf:
        parent.add_prop("num_leaves", 1)
        return 1
    else:
        leaves = 0
        for child in parent.children:
            leaves += label_desc_leaves(child)
        parent.add_prop("num_leaves", leaves)
        return leaves


def compute_ed_scores(parent, scores_dict=None, running_sum=None, root=True):
    """Recursive function to accumulate ED scores across the tree, i.e. starting at the root,
    sum up branch lengths weighted by 1/the number of descendant leaves. Assumes that
    all nodes have a "num_leaves" attribute; if not then call label_desc_leaves() first.
    Assumes tree has branch lengths, not just dates."
    """

    if root:
        running_sum = 0
    else:
        running_sum += parent.dist/parent.props["num_leaves"]

    parent.add_prop("ed_score", running_sum)

    if parent.is_leaf and scores_dict is not None:
        scores_dict[(parent.name,
                    parent.props["domain"],
                    parent.props["kingdom"],
                    parent.props["phylum"],
                    parent.props["clas"],
                    parent.props["order"],
                    parent.props["family"])] = running_sum
    else:
        for child in parent.children:
            compute_ed_scores(child, scores_dict, running_sum, False)


def append_ed_scores(dated_tre, existing_scores):
    """Compute ED scores for the dated_tre, and add the score for each species to
    the dictionary existing_scores. The dictionary is keyed on leaf names, e.g. for the
    Open Tree the keys will be of the form Genus_name_ott1234, and the values
    are a list of ED scores for the species. Assumes tree has branch lengths.
    """

    logger.info("Computing ED scores.")
    if "num_leaves" not in dated_tre.props:
        label_desc_leaves(dated_tre)

    new_scores = {}
    compute_ed_scores(dated_tre, 0, new_scores, True)

    for key in new_scores:
        if key not in existing_scores:
            existing_scores[key] = []

        existing_scores[key].append(new_scores[key])


def write_ed_scores(filename, scores):
    logger.info("Writing ED score distributions.")
    fout = open(filename, "w")
    fout.write("leaf_name\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tmean\tmin\t2.5pct\t25pct\tmedian\t75pct\t97.5pct\tmax\tN\n")
    for key in scores:
        fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n" % (key[0],
                                                                                         key[1],
                                                                                         key[2],
                                                                                         key[3],
                                                                                         key[4],
                                                                                         key[5],
                                                                                         key[6],
                                                                                         np.mean(scores[key]),
                                                                                         np.min(scores[key]),
                                                                                         np.percentile(scores[key],2.5),
                                                                                         np.percentile(scores[key],25),
                                                                                         np.percentile(scores[key],50),
                                                                                         np.percentile(scores[key],75),
                                                                                         np.percentile(scores[key],97.5),
                                                                                         np.max(scores[key]),
                                                                                         len(scores[key])))
    fout.close()


def compute_pd(parent, root_call=True):
    """Compute total PD below each node. Assume the tree has branch lengths."""
    pd = 0
    for child in parent.children:
        pd += compute_pd(child, False)

    parent.add_prop("pd", pd)

    if root_call:
        return pd
    else:
        return pd + parent.dist


def save_pd_for_clades(tre, pd_clades, pd_dict, dates_dict, spp_dict):
    for node in tre.traverse():
        if node.name in pd_clades:
            pd_dict[node.name].append(node.props["pd"])
            dates_dict[node.name].append(node.props["date"])
            spp_dict[node.name].append(len(node))


def write_pd_dists(filename, pd_dict, dates_dict, spp_dict):
    def list_to_tab_str(ls):
        ls_str = ""
        for item in ls[:-1]:
            ls_str += str(item)
            ls_str += '\t'
        ls_str += str(ls[-1])

        return ls_str

    fout = open("%s_pd_for_clades.txt" % filename, "w")
    fout.write("Clade\tspp.\tMean\tMin\t2.5pct\t16pct\t25pct\tmedian\t75pct\t84pct\t97.5pct\tMax\n")
    for clade in pd_dict:
        if len(spp_dict[clade]) > 0:
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
        if len(spp_dict[clade]) > 0:
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


def assign_iucn_status(tre, iucn_filename="config/latest_iucn_2025.csv"):
    import csv

    iucn_lookup = {}
    with open(iucn_filename, newline='') as csvfile:
        rdr = csv.reader(csvfile)
        for idx, line in enumerate(rdr):
            if idx == 0:
                # first line has column headings
                continue
            iucn_lookup[int(line[0])] = line[1]

    for node in tre.leaves():
        name_parts = node.name.split('_')
        ottid = int(name_parts[-1][3:])

        if ottid in iucn_lookup:
            status = iucn_lookup[ottid]
            if status == "LR/nt" or status == "LR/cd":
                status = "NT"
            if status == "LR/lc":
                status = "LC"

            node.add_prop("iucn_status", iucn_lookup[ottid])


def assign_extinction_risks(tre, rng=None, lookup_table="config/p_extinction.csv", randomise_risk=False, missing_value=None):
    """Add extinction probabilities to leaves in the tree. Assumes iucn_status is labelled
    on the leaves in the tree.
    If randomise_risk=True, draw an extinction risk from the distribution rather than using the median.
    If missing_value is not None, use the value given instead of drawing a random value.
    """

    import csv

    pext_lookup = {"ALL": [[], None]}

    with open(lookup_table, newline='') as csvfile:
        rdr = csv.reader(csvfile)
        for idx, line in enumerate(rdr):
            if idx == 0:
                # first line has column headings
                continue
            if line[0] not in pext_lookup:
                pext_lookup[line[0]] = [[], None]

            pext_lookup[line[0]][0].append(float(line[1]))
            pext_lookup['ALL'][0].append(float(line[1]))

    for key in pext_lookup:
        pext_lookup[key][1] = np.median(pext_lookup[key][0])

    for leaf in tre.leaves():
        if "iucn_status" in leaf.props and leaf.props["iucn_status"] in pext_lookup:
            status = leaf.props["iucn_status"]
            if randomise_risk:
                # pick value from distribution
                random_idx = rng.integers(len(pext_lookup[status][0]))
                risk = pext_lookup[status][0][random_idx]
            else:
                # use median
                risk = pext_lookup[status][1]
        else:
            if missing_value is not None:
                risk = missing_value
            else:
                random_idx = rng.integers(len(pext_lookup["ALL"][0]))
                risk = pext_lookup["ALL"][0][random_idx]

        leaf.add_prop("pext", risk)


def compute_edge2_scores(tre, scores_dict=None):
    """Compute ED2 and EDGE2 scores for the input tre. Assumes each leaf node has a property "pext" which
    contains its extinction probability. Leaf nodes will be labelled with properties "ed2_score" and "edge2_score"
    containing their scores.
    """

    # recursively compute and label the products of extinction risks below each internal node
    def label_pext_products(parent):
        if parent.is_leaf:
            return parent.props["pext"]
        else:
            product = 1

            for child in parent.children:
                product *= label_pext_products(child)

            parent.add_prop("pext_product", product)

            return product

    label_pext_products(tre)

    # then push those products down the tree, multiplying them by branch lengths. For a leaf node we compute
    # ED2 by taking the terminal branch length and adding (sum of products above / own extinction risk)
    for node in tre.traverse(strategy="preorder"):
        if node.is_leaf:
            node.add_prop("ed2_score", node.dist + node.up.props["ed2_intermediate"]/node.props["pext"])
            node.add_prop("edge2_score", node.props["pext"] * node.props["ed2_score"])

            if scores_dict is not None:
                key = (node.props["domain"] if node.props["domain"] is not None else "no_domain",
                       node.props["kingdom"] if node.props["kingdom"] is not None else "no_kingdom",
                       node.props["phylum"] if node.props["phylum"] is not None else "no_phylum",
                       node.name)

                scores_dict.setdefault(key,[[], []])[0].append(node.props["ed2_score"])
                scores_dict[key][1].append(node.props["edge2_score"])
        else:
            if node is tre:
                if node.dist is None:
                    node.add_prop("ed2_intermediate", 0)
                else:
                    node.add_prop("ed2_intermediate", node.dist * node.props["pext_product"])
            else:
                node.add_prop("ed2_intermediate", node.up.props["ed2_intermediate"] + node.dist * node.props["pext_product"])


def write_edge2_scores(filename, scores_dict, per_phylum=None):
    if per_phylum:
        keys = [(key[0], key[1], key[2], -np.mean(scores_dict[key][1]), key[3]) for key in scores_dict]
        keys.sort()

        dict_to_write = {}

        current_key = None
        count = 0
        for i, sorted_key in enumerate(keys):
            key = (sorted_key[0], sorted_key[1], sorted_key[2], sorted_key[4])
            if current_key is None:
                count = 1
                dict_to_write[key] = scores_dict[key]
                current_key = key
            elif key[0] == current_key[0] and key[1] == current_key[1] and key[2] == current_key[2]:
                if count < per_phylum:
                    dict_to_write[key] = scores_dict[key]
                    count += 1
            else:
                count = 1
                dict_to_write[key] = scores_dict[key]
                current_key = key

    else:
        dict_to_write = scores_dict

    with open(filename, "w") as fout:
        # column headings
        col_headings =  "leaf_name,domain,kingdom,phylum,N,"
        col_headings += "ED2_mean,min,2.5pct,median,97.5pct,max,"
        col_headings += "EDGE2_mean,min,2.5pct,median,97.5pct,max\n"
        fout.write(col_headings)

        for key in dict_to_write:
            string_to_write = "{:s},{:s},{:s},{:s},{:d},"
            string_to_write += "{:f},{:f},{:f},{:f},{:f},{:f},"
            string_to_write += "{:f},{:f},{:f},{:f},{:f},{:f}\n"

            fout.write(string_to_write.format(key[3] if key[3] is not None else "",
                                              key[0] if key[0] is not None else "",
                                              key[1] if key[1] is not None else "",
                                              key[2] if key[2] is not None else "",
                                              len(scores_dict[key][0]),
                                              np.mean(dict_to_write[key][0]),
                                              np.min(dict_to_write[key][0]),
                                              np.percentile(dict_to_write[key][0],2.5),
                                              np.percentile(dict_to_write[key][0],50),
                                              np.percentile(dict_to_write[key][0],97.5),
                                              np.max(dict_to_write[key][0]),
                                              np.mean(dict_to_write[key][1]),
                                              np.min(dict_to_write[key][1]),
                                              np.percentile(dict_to_write[key][1],2.5),
                                              np.percentile(dict_to_write[key][1],50),
                                              np.percentile(dict_to_write[key][1],97.5),
                                              np.max(dict_to_write[key][1])
                                              ))


def compute_evoh(tre, rho):
    def evoh_beta(rho, edge_length):
        return np.exp(-rho * edge_length)

    def label_evoh_p(parent, rho):
        if parent.is_leaf:
            parent.add_prop("evoh_p", 1)
            return 1
        else:
            product = 1
            for child in parent.children:
                product *= (1 - label_evoh_p(child, rho) * evoh_beta(rho, child.dist))

            p = 1 - product

            parent.add_prop("evoh_p", p)

            return p

    label_evoh_p(tre, rho)

    gamma_rho = 0
    for node in tre.traverse(strategy="preorder"):
        if node is tre:
            # don't include stem - assume we are complete tree and the root node is already the origin of life
            continue

        gamma_rho += (1 - evoh_beta(rho, node.dist)) * node.props["evoh_p"]

    gamma_rho /= (1 - np.exp(-rho))

    return gamma_rho


def compute_gamma(tre):
    """Gamma statistic from Pybus & Harvey (2000). Will be normally distributed with mean 0 for
    a tree with branch lengths that fit a constant-rate pure birth model.
    Here, assume a bifurcating tree with all nodes are dated, and dated consistently,
    with leaf nodes that have date 0. Variable names and subscripts are taken from the Pybus paper.
    """

    dates_list = []
    for node in tre.traverse():
        if node.props["date"] > 0:
            dates_list.append(node.props["date"])

    dates_list.sort(reverse=True)

    g = []
    for i in range(1, len(dates_list)):
        g.append(dates_list[i-1] - dates_list[i])
    g.append(dates_list[-1])

    n = len(tre)

    T = 0
    for j in range(2, n+1):
        # j = 2, ..., n
        T += j * g[j-2]

    tmpsum = 0
    for i in range(2, n):
        # i = 2, ..., n-1
        for k in range(2, i+1):
            # k = 2, ..., i
            tmpsum += k * g[k-2]

    gamma = ( (1/(n-2) * tmpsum) - T/2 ) / ( T * np.sqrt(1 / (12*(n-2))) )

    return gamma


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
    For use in interpolation algorithm from Guo et al. (2025).
    """
    if parent.is_leaf:
        if parent.props["date"] != 0:
            print("Leaf node with non-zero date")
        oldest_paths = [[parent.props["date"], 1]]
        parent.oldest_paths = oldest_paths
    else:
        oldest_paths = []
        for child in parent.children:
            child_paths = date_labelling_guo(child)
            for path in child_paths:
                oldest_paths.append(path)

        if parent.props["date"] is None:
            parent.oldest_paths = oldest_paths
            for path in oldest_paths:
                path[1] += 1
        else:
            oldest_paths = [[parent.props["date"], 1]]
            parent.oldest_paths = oldest_paths

    return oldest_paths


def compute_ed_scores_guo(parent, ancestral_date, root=True):
    """Implementation of ED interpolation algorithm from Guo et al. (2025)."""

    ed_scores = []

    if not root:
        for path in parent.oldest_paths:
            if path[1] != 0:
                ed_scores.append(((ancestral_date-path[0])/path[1]) / parent.props["desc_leaves"])

    parent.add_prop("ed_scores_guo", ed_scores)

    for child in parent.children:
        if parent.props["date"]:
            compute_ed_scores_guo(child, parent.props["date"], False)
        else:
            compute_ed_scores_guo(child, ancestral_date, False)


def sum_ed_scores_guo(parent, root=True):
    """Compute PD by summing ED computed my method from Guo et al. (2025)."""

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
        if node.props["date"] > 0:
            dates_list.append((node.props["date"], len(node.children)-1))

            if not node.imputed_date:
                nonimputed_dates.append(-node.props["date"])

    dates_list.sort(reverse=True)

    lineages = [2]
    dates = [-dates_list[0][0]]
    for i in range(1,len(dates_list)):
        lineages.append(lineages[-1] + dates_list[i][1])
        dates.append(-dates_list[i][0])

    return dates, lineages, nonimputed_dates
