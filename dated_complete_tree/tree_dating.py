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


import copy
import numpy as np

import logging
logger = logging.getLogger(__name__)


def assign_dates(tre, dates, sample_dates=False, rng=None, num_sources=False, use_youngest_dates=False, use_oldest_dates=False):
    """Assign chronosynth dates to nodes in the tree.
    If sample_dates=True, first create a dictionary of all date sources whoe keys are the sources and values their (random) priority ordering.
    For each dated node, assign the date of the highest-priority source."""

    # remove any dates of 0.0
    zero_dates = []
    for ott_name in dates["node_ages"]:
        for i, source in enumerate(dates["node_ages"][ott_name]):
            if float(source["age"]) < 0.000001:
                zero_dates.append((ott_name, i))

    for ott_name, i in zero_dates:
        del dates["node_ages"][ott_name][i]
        if len(dates["node_ages"][ott_name]) == 0:
            del dates["node_ages"][ott_name]

    # randomly prioritise all date sources
    if sample_dates:
        all_sources = set()
        for ott_name in dates["node_ages"]:
            for source in dates["node_ages"][ott_name]:
                all_sources.add(source["source_id"])

        source_priorities = list(all_sources)
        rng.shuffle(source_priorities)

        priority_dict = {}
        for i, sourceid in enumerate(source_priorities):
            priority_dict[sourceid] = i

    # date_sources = {}

    # assign dates to any node with at least one source
    for node in tre.traverse(strategy="preorder"):
        ott_name = node.name.split('_')[-1]
        date = None
        sources = None

        if node.is_leaf:
            # leaf nodes are assumed to be extant; assign date of zero
            date = 0
        elif ott_name in dates["node_ages"]:
            ages = [float(source["age"]) for source in dates["node_ages"][ott_name]]

            if num_sources:
                node.add_prop("date_sourceids", set([source["source_id"] for source in dates["node_ages"][ott_name]]))

            # date_sources[node] = copy.copy(ages)

            if sample_dates:
                # sample dates based on source prioritisation
                sources = [source["source_id"] for source in dates["node_ages"][ott_name]]
                source_choice = None
                for i, sourceid in enumerate(sources):
                    if source_choice is None:
                        source_choice = (i, sourceid)
                        logger.debug("Initial source choice for %s is %s, priority %d." % (ott_name, sourceid, priority_dict[sourceid]))
                    else:
                        if priority_dict[sourceid] < priority_dict[source_choice[1]]:
                            logger.debug("Superseding source choice for %s with %s, priority %d vs %d." % (ott_name, sourceid, priority_dict[sourceid], priority_dict[source_choice[1]]))
                            source_choice = (i, sourceid)

                date = round(ages[source_choice[0]], 5)

            elif use_youngest_dates:
                date = round(min(ages),5)
            elif use_oldest_dates:
                date = round(max(ages),5)
            else:
                date = round(np.median(ages), 5)

            if date < 0.000001:
                logger.warning("Warning: computed age of zero on interior node %s; instead, setting date for this node to None" % node.name)

        node.add_prop("date", date)
        node.add_prop("imputed_date", False)
        node.add_prop("imputation_type", 0)

    # return date_sources


def remove_inconsistent_dates(parent, mrad=None):
    """Recurse in preorder through the tree, removing the date info from any nodes that have dates older
    than a date found on an ancestor. (mrad = most recent ancestor date)
    Implemented to ape BLADJ; prefer use of label_older_descendants followed by dq_date_removal.
    """

    if mrad is None:
        # must be root node
        if not parent.props["date"]:
            raise Exception("Require root date.")
        next_mrad = parent.props["date"]
    else:
        next_mrad = mrad
        # if we have a date, check it against the most recent ancestor
        if parent.props["date"]:
            if round(parent.props["date"],6) >= round(mrad,6):
                # if it is older than an ancestor, throw away the date information at this node
                parent.props["date"] = None
                # logger.info("Removing inconsistent date at node %s." % parent.name)
            else:
                next_mrad = parent.props["date"]

    for child in parent.children:
        remove_inconsistent_dates(child, next_mrad)


def strip_undated_nodes(tre):
    """Return a copy of the input tree cut down to only dated nodes."""

    stripped_tre = tre.copy()

    nodes_to_remove = []
    for node in stripped_tre.traverse(strategy="preorder"):
        if node == stripped_tre:
            continue

        if not node.props["date"]:
            nodes_to_remove.append(node)

    for node in nodes_to_remove:
        parent = node.up

        for child in node.get_children():
            parent.add_child(child.detach())

        node.detach()

    return stripped_tre


def label_older_descendants(parent):
    """Labels each node with a list of descendant nodes that have dates older than its own."""

    descendants = []

    for child in parent.children:
        child_descendants = label_older_descendants(child)
        descendants.extend(child_descendants)
        # ignore node if is has no date (None) or a date of 0
        if child.props["date"]:
            descendants.append(child)

    # store list of descendant nodes with an older date than this node
    if parent.props["date"]:
        parent.props["older_descendants"] = [
            child for child in descendants if round(child.props["date"],6) >= round(parent.props["date"],6)
        ]

    return descendants


def build_dq_dict(tre):
    """Take a tree labelled by the label_older_descendants function and builds a dictionary. The keys
    are the inconsistent nodes. The values are a list: [set of older descendants,
                                                        set of younger ancestors,
                                                        pointer to equivalent node in whole tree]
    """

    dq_dict = {}

    for node in tre.traverse(strategy="preorder"):
        if node == tre:
            # ignore root node; we assume this is correct
            continue

        if node.props["date"] and len(node.props["older_descendants"]) > 0:
            for desc_node in node.props["older_descendants"]:
                if node.name not in dq_dict:
                    dq_dict[node.name] = [set([]), set([]), node]
                if desc_node.name not in dq_dict:
                    dq_dict[desc_node.name] = [set([]), set([]), desc_node]

                dq_dict[node.name][0].add(desc_node.name)
                dq_dict[desc_node.name][1].add(node.name)

    return dq_dict


def dq_date_removal(tre):
    """Remove dates from nodes such that as few dates are removed as possible to create a consistent set of dates
    on the tree.
    """

    dq_dict = build_dq_dict(tre)

    dq_counts = {}
    for key_node_name in dq_dict:
        for node_name in dq_dict[key_node_name][0]:
            if node_name not in dq_counts:
                dq_counts[node_name] = 0
            dq_counts[node_name] += 1

        for node_name in dq_dict[key_node_name][1]:
            if node_name not in dq_counts:
                dq_counts[node_name] = 0
            dq_counts[node_name] += 1

    dq_counts_info = []
    for node_name in dq_counts:
        dq_counts_info.append([dq_counts[node_name], dq_dict[node_name][2]])

    def tuple_value(x):
        # used for picking a max value: i.e. highest count, with tie breakers being:
        #  1. favour removing dates on mrca nodes rather than named nodes;
        #  2. favour removing the tied node whose paired node is younger
        if x[1].name[:4] == "mrca":
            mrca = 1
        else:
            mrca = 0

        return (x[0], mrca, -round(x[1].props["date"],6))


    if len(dq_counts_info) == 0:
        logger.info("In dq_date_removal: No inconsistent dates to remove.")
        return

    max_count = max(dq_counts_info, key=tuple_value)

    while max_count[0] > 0:
        node_to_blank = max_count[1]

        node_to_blank.props["date"] = None

        dq_counts_info.remove(max_count)

        for j in range(len(dq_counts_info)):
            key_node_name = dq_counts_info[j][1].name

            if node_to_blank.name in dq_dict[key_node_name][0]:
                dq_dict[key_node_name][0].remove(node_to_blank.name)
                dq_counts_info[j][0] -= 1
                continue

            if node_to_blank.name in dq_dict[key_node_name][1]:
                dq_dict[key_node_name][1].remove(node_to_blank.name)
                dq_counts_info[j][0] -= 1

        max_count = max(dq_counts_info, key=tuple_value)


def use_shortest_path(x):
    """Key function for use in max function in date_labelling, so that max uses the shortest path as a tiebreaker."""
    return [x[0], -x[1]]


def date_labelling(parent):
    """Recurse in postorder through the tree, labelling each node with the oldest date below it
    and the path length (number of nodes) to that date. If the oldest date is a tie (usually 0),
    we must choose a path to it. This computes both the longest path (most nodes) and the shortest.

    Return value is: [oldest date found so far below this node, longest path length to it],
                     [oldest date found so far below this node, shortest path length to it]
    """
    if parent.is_leaf:
        if parent.props["date"] != 0:
            logger.warning("Leaf node %s has non-zero date." % parent.name)
    else:
        oldest_path_long = [0, -1e8]
        oldest_path_short = [0, 1e8]
        for child in parent.children:
            new_path_long, new_path_short = date_labelling(child)
            oldest_path_long = max(oldest_path_long, new_path_long)
            oldest_path_short = max(oldest_path_short, new_path_short, key=use_shortest_path)

    if parent.props["date"] is None:
        parent.props["oldest_path_long"] = copy.copy(oldest_path_long)
        parent.props["oldest_path_short"] = copy.copy(oldest_path_short)
        oldest_path_long[1] += 1
        oldest_path_short[1] += 1
    else:
        # if there is a date, doesn't matter what was received; just return this date
        oldest_path_long = [parent.props["date"], 0]
        oldest_path_short = [parent.props["date"], 0]
        parent.props["oldest_path_long"] = copy.copy(oldest_path_long)
        parent.props["oldest_path_short"] = copy.copy(oldest_path_short)
        oldest_path_long[1] += 1
        oldest_path_short[1] += 1

    return oldest_path_long, oldest_path_short


def impute_clade_birth_model(choices, dates, rng):
    # assumes a bifurcating tree
    i = 1
    while len(choices) > 0:
        probs = np.array([c.props["num_leaves"]-1 for c in choices])
        probs = probs / np.sum(probs)

        next_node_idx = rng.choice(range(len(choices)), p=probs)
        next_node = choices[next_node_idx]

        next_node.props["date"] = dates[i]
        next_node.props["imputation_type"] = 4
        next_node.props["imputed_date"] = True

        choices.remove(next_node)

        if not next_node.children[0].is_leaf:
            choices.append(next_node.children[0])
        if not next_node.children[1].is_leaf:
            choices.append(next_node.children[1])

        i += 1


def impute_missing_dates(tre, l=1, m=0, use_logN_model=False, use_birth_model=False, rng=None, counts=[0,0]):
    """Traverse the tree in preorder, giving each undated node a date spaced along the path between between
    its parent (which always has a date, since this is preorder traversal) and the oldest date found below
    it (as labelled by the date_labelling function). Assumes root node is dated.

    When the oldest date beneath is a tie (usually because the oldest date is 0), the tie can be broken
    by using the longest path or the shortest path. This function computes both versions, then interpolates
    between the two solutions based on the parameter l:
        date = l * date_from_longest_path + (1-l) * date_from_shortest_path

    In addition, interpolation along a path can use equal spacing (m=0), or spacing that biases dates older
    (m > 0) or younger (m < 0). Uses spacing along an exponential function, i.e. y = exp(m*x).
    Values of m between -2 and 2 are pretty sensible.
    """
    if use_logN_model or use_birth_model:
        def label_pct_dates(parent):
            if parent.is_leaf:
                results = [0,0,1]
            elif not parent.props["date"]:
                results = [1,0,0]
            else:
                results = [1,1,0]

            for child in parent.children:
                new_results = label_pct_dates(child)
                results[0] += new_results[0]
                results[1] += new_results[1]
                results[2] += new_results[2]

            parent.add_prop("child_tree_size", results[0])
            parent.add_prop("num_dates", results[1])
            parent.add_prop("num_leaves", results[2])

            return results

        label_pct_dates(tre)

    for node in tre.traverse(strategy="preorder"):
        if node is tre:
            node.props["imputed_date"] = False
            continue

        if node.props["date"] is None:
            if use_logN_model and node.props["oldest_path_long"][0] == 0:
                counts[0] += 1
                if node.props["num_leaves"] > 1 and node.up.props["num_leaves"] > 1 and node.up.props["num_leaves"] > node.props["num_leaves"]:
                    node.props["date"] = node.up.props["date"] * np.log(node.props["num_leaves"])/np.log(node.up.props["num_leaves"])
                else:
                    # need backup option of standard BLADJ in case of o---o---o situation or where num_leaves is the same for
                    # both parent and child
                    node.props["date"] = node.up.props["date"] - (node.up.props["date"] - node.props["oldest_path_long"][0]) / (node.props["oldest_path_long"][1]+1)
                node.props["imputation_type"] = 3
            elif use_birth_model and node.props["oldest_path_long"][0] == 0:
                # birth model, assuming that the leaves are infinitesimally close to the next speciation event

                # one lineage to model or two?
                if not node.up.children[0].is_leaf and node.up.children[0].props["oldest_path_long"][0] == 0 and not node.up.children[1].is_leaf and node.up.children[1].props["oldest_path_long"][0] == 0:
                    lineages = 2
                    leaves = node.up.props["num_leaves"]
                    choices = [node.up.children[0], node.up.children[1]]
                else:
                    lineages = 1
                    leaves = node.props["num_leaves"]
                    choices = [node]

                crown_date = node.up.props["date"]
                birth_rate = (np.log(lineages) - np.log(leaves+1)) / crown_date

                # dates[0] will be crown date; we do not use this
                dates = np.log(np.arange(lineages, leaves+1) / (leaves+1)) / birth_rate

                impute_clade_birth_model(choices, dates, rng)
            else:
                if node.up.props["imputed_date"]:
                    if node.up.props["oldest_path_long"][0] == node.props["oldest_path_long"][0] and node.up.props["oldest_path_long"][1] == node.props["oldest_path_long"][1]+1:
                        node.props["date_above_long"] = node.up.props["date_above_long"]
                        node.props["mu_spacing_long"] = node.up.props["mu_spacing_long"]
                    else:
                        node.props["date_above_long"] = node.up.props["date_long"]
                        mu_spacing_long  = np.exp(m * np.linspace(0, 1, node.props["oldest_path_long"][1]+1))
                        node.props["mu_spacing_long"]  = np.cumsum( mu_spacing_long  / np.sum(mu_spacing_long) )

                    if node.up.props["oldest_path_short"][0] == node.props["oldest_path_short"][0] and node.up.props["oldest_path_short"][1] == node.props["oldest_path_short"][1]+1:
                        node.props["date_above_short"] = node.up.props["date_above_short"]
                        node.props["mu_spacing_short"] = node.up.props["mu_spacing_short"]
                    else:
                        node.props["date_above_short"] = node.up.props["date_short"]
                        mu_spacing_short = np.exp(m * np.linspace(0, 1, node.props["oldest_path_short"][1]+1))
                        node.props["mu_spacing_short"] = np.cumsum( mu_spacing_short / np.sum(mu_spacing_short) )

                else:
                    node.props["date_above_long"] = node.up.props["date"]
                    node.props["date_above_short"] = node.up.props["date"]

                    mu_spacing_long  = np.exp(m * np.linspace(0, 1, node.props["oldest_path_long"][1]+1))
                    mu_spacing_short = np.exp(m * np.linspace(0, 1, node.props["oldest_path_short"][1]+1))

                    node.props["mu_spacing_long"]  = np.cumsum( mu_spacing_long  / np.sum(mu_spacing_long) )
                    node.props["mu_spacing_short"] = np.cumsum( mu_spacing_short / np.sum(mu_spacing_short) )

                node.props["date_long"]  = node.props["date_above_long"] - (node.props["date_above_long"] - node.props["oldest_path_long"][0]) * node.props["mu_spacing_long"][-(node.props["oldest_path_long"][1]+1)]
                node.props["date_short"] = node.props["date_above_short"] - (node.props["date_above_short"] - node.props["oldest_path_short"][0]) * node.props["mu_spacing_short"][-(node.props["oldest_path_short"][1]+1)]

                counts[1] += 1
                node.props["date"] = l*node.props["date_long"] + (1-l)*node.props["date_short"]

                if l == 1 and m == 0 and node.props["oldest_path_long"][0] == 0:
                    node.props["imputation_type"] = 1
                elif l == 0 and m == 0 and node.props["oldest_path_long"][0] == 0:
                    node.props["imputation_type"] = 2
                elif node.props["oldest_path_long"][0] == 0:
                    node.props["imputation_type"] = 5
                else:
                    node.props["imputation_type"] = 6

            node.props["imputed_date"] = True


def round_to_4sf(x):
    """Round to four sig figs; intended for ages in Mya, so assumes input is less than 10000 (would be
    older than the solar system).
    """
    if x > 1000:
        return round(x)
    elif x > 100:
        return round(x,1)
    elif x > 10:
        return round(x,2)
    elif x > 1:
        return round(x,3)
    else:
        return round(x,4)


def compute_branch_lengths(tre, round_numbers=False):
    """Fill in 'dist' field with branch lengths. Intended for a fully dated tree."""
    for node in tre.traverse(strategy="preorder"):
        if node.up:
            dist = node.up.props["date"] - node.props["date"]
            if dist < 0:
                logger.error("Error: Negative branch length above %s" % node.name)

            if round_numbers:
                node.dist = round_to_4sf(dist)
            else:
                node.dist = dist


def write_tree_with_branch_lengths(tre, filename):
    """Write out dated tree in Newick format (suitable for OneZoom). Branch lengths rounded
    to 4 sig figs to save space in the text file."""
    # compute_branch_lengths(tre, round_numbers=False)

    from ete4.parser.newick import DIST, NAME
    MY_DIST = DIST.copy()
    MY_DIST['write'] = lambda x: "%.7f" % float(x)
    custom_parser = {
        'leaf':     [NAME, MY_DIST],
        'internal': [NAME, MY_DIST]
    }

    tre.write(outfile=filename,
                 parser=custom_parser,
                 format_root_node=True)


def compute_dates(parent):
    """Fill in 'date' field with dates, given a tree with branch lengths in units of time. Assume leaf nodes have a date of 0.
    """
    if parent.is_leaf:
        parent.props["date"] = 0
    else:
        for child in parent.children:
            compute_dates(child)
            parent.props["date"] = child.props["date"] + child.dist
