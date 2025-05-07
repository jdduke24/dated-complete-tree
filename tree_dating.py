import copy

import logging
logger = logging.getLogger(__name__)


def remove_inconsistent_dates(parent, mrad):
    """Recurse in preorder through the tree, removing the date info from any nodes that have dates older
    than a date found on an ancestor. (mrad = most recent ancestor date)
    """
    next_mrad = mrad
    # if we have a date, check it against the most recent ancestor
    if parent.date:
        if parent.date > mrad:
            # if it is older than an ancestor, throw away the date information at this node
            parent.date = None
            logger.info("Removing inconsistent date at node %s." % parent.name)
        else:
            next_mrad = parent.date

    for child in parent.children:
        remove_inconsistent_dates(child, next_mrad)


def date_labelling(parent):
    """Recurse in postorder through the tree, labelling each node with the oldest date below it
    and the path length (number of nodes) to that date.
    Return value is: [oldest date found so far below this node, path length to it]
    """
    if parent.is_leaf():
        if parent.date != 0:
            logger.warning("Leaf node %s has non-zero date." % parent.name)
        oldest_path = [parent.date, 0]
    else:
        oldest_path = [0, 0]
        for child in parent.children:
            oldest_path = max(oldest_path, date_labelling(child))

    parent.add_feature("oldest_path", copy.deepcopy(oldest_path))

    if parent.date is None:
        oldest_path[1] += 1
    else:
        oldest_path = [parent.date, 1]

    return oldest_path


def impute_missing_dates(tre):
    """Traverse the tree in preorder, giving each undated node a date spaced along the path between between
    its parent (which always has a date, since this is preorder traversal) and the oldest date found below
    it (as labelled by the date_labelling function). Assumes root node is dated.
    """
    for node in tre.traverse(strategy="preorder"):
        if node.is_root():
            node.add_feature("imputed_date",False)
            continue

        if node.date is None:
            node.date = node.up.date - (node.up.date - node.oldest_path[0]) / (node.oldest_path[1]+1)
            node.imputed_date = True


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
    for node in tre.traverse():
        if node.up:
            dist = node.up.date - node.date
            if dist < 0:
                logger.warning("Warning: Negative branch length above %s" % node.name)
            if round_numbers:
                node.dist = round_to_4sf(node.up.date - node.date)
            else:
                node.dist = node.up.date - node.date


def write_tree_with_branch_lengths(tre, filename):
    """Write out dated tree in Newick format (suitable for OneZoom). Branch lengths rounded
    to 4 sig figs to save space in the text file."""
    compute_branch_lengths(tre, round_numbers=True)

    tre.write(outfile=filename,
                 format=1,
                 format_root_node=True)


def compute_dates(tre, root_date):
    """Fill in 'date' field with dates, given a tree with branch lengths in units of time. Intended for a tree given
    branch lengths by phylocom bladj. Requires the date used by bladj on the root node.
    """
    tre.add_feature("date", root_date)
    tre.add_feature("imputed_date", False)
    for node in tre.traverse(strategy="preorder"):
        if node.up:
            node.add_feature("date", node.up.date - node.dist)
            node.add_feature("imputed_date", False)
