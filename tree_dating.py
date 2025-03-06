import copy

import logging
logger = logging.getLogger(__name__)

# mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
def remove_inconsistent_dates(parent, mrad):
    next_mrad = mrad
    # if we have a date, check it against the most recent ancestor
    if parent.date:
        if parent.date > mrad:
            # if it is older than an ancestor, throw away all the date information at this node
            parent.date = None
            logger.info("Removing inconsistent date at node %s." % parent.name)
        else:
            next_mrad = parent.date

    for child in parent.children:
        remove_inconsistent_dates(child, next_mrad)


# label dates with nearest dates below them and path lengths (number of nodes) to those dates
# return value is [oldest date found so far below this node, path length to it]
def date_labelling(parent):
    if parent.is_leaf():
        if parent.date != 0:
            logger.warning("Leaf node %s has non-zero date." % parent.name)
        oldest_path = [parent.date,0]
    else:
        oldest_path = [0,0]
        for child in parent.children:
            oldest_path = max(oldest_path,date_labelling(child))

    parent.add_feature("oldest_path", copy.deepcopy(oldest_path))

    if parent.date is None:
        oldest_path[1] += 1
    else:
        oldest_path = [parent.date,1]

    return oldest_path


# impute dates by dating each undated node between its parent and the oldest date found below it
def impute_missing_dates(tre):
    for node in tre.traverse(strategy="preorder"):
        if node.is_root():
            node.add_feature("imputed_date",False)
            continue

        if node.date is None:
            node.date = node.up.date - (node.up.date - node.oldest_path[0]) / (node.oldest_path[1]+1)
            node.imputed_date = True


def compute_branch_lengths(tre, round_numbers=False):
    for node in tre.traverse():
        if node.up:
            dist = node.up.date - node.date
            if dist < 0:
                logger.warning("Warning: Negative branch length above %s" % node.name)
            if round_numbers:
                node.dist = round_to_4sf(node.up.date - node.date)
            else:
                node.dist = node.up.date - node.date


def round_to_4sf(x):
    """Round to four sig figs: assumes input is less than 10000."""
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


def write_tree_with_branch_lengths(tre, filename):
    compute_branch_lengths(tre, round_numbers=True)

    tre.write(outfile=filename,
                 format=1,
                 format_root_node=True)
