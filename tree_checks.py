from taxonomy_utils import tx_levels
import logging
logger = logging.getLogger(__name__)


def get_phylogeny_only(tre):
    if tre.ph_tx != "PH":
        raise Exception("Input tree must begin with phylogeny.")

    new_tree = tre.copy()

    to_delete = []
    for node in new_tree.traverse():
        if node.ph_tx != "PH":
            to_delete.append(node)

    for node in to_delete:
        node.detach()
        del node

    import gc
    gc.collect()

    return new_tree


def check_bifurcating(tre, filename=None):
    if filename:
        fout = open(filename, "wt")

    one_child_nodes = []
    polytomy_nodes = []
    for node in tre.traverse():
        if len(node.children) > 2:
            polytomy_nodes.append((node.name, node.tx_level, len(node.children)))

        if len(node.children) == 1:
            one_child_nodes.append((node.name, node.tx_level, node.children[0].tx_level))

    for name, rank, desc_rank in one_child_nodes:
        if filename:
            fout.write("%s,%s,%s\n" % (name, rank, desc_rank))
        else:
            print(name, rank)

    for name, rank, num_children in polytomy_nodes:
        if filename:
            fout.write("%s,%s,%d\n" % (name, rank, num_children))
        else:
            print(name, rank)

    print("There are %d one-child nodes." % len(one_child_nodes))

    print("There are %d polytomies." % len(polytomy_nodes))

    if filename:
        fout.close()


def count_subspecies(tre, return_subsp=False):
    count = 0
    subsp = set()
    for node in tre.traverse():
        if tx_levels[node.tx_level] > 0 and tx_levels[node.tx_level] < 3:
            count += 1
            if return_subsp:
                subsp.add(node)

    if return_subsp:
        return subsp
    else:
        return count


def check_taxonomy_order(parent, file=None, current_rank=None):
    if current_rank is None:
        current_rank = parent.tx_level

    if tx_levels[parent.tx_level] > 0 or tx_levels[current_rank] < 0:
        current_rank = parent.tx_level

    errors = 0
    for child in parent.children:
        if tx_levels[current_rank] > 0 and tx_levels[child.tx_level] > 0:
            if tx_levels[child.tx_level] > tx_levels[current_rank]:
                ancestor = child.up
                while(tx_levels[ancestor.tx_level] < 0 or tx_levels[ancestor.tx_level] > tx_levels[child.tx_level]):
                    ancestor = ancestor.up
                if file:
                    file.write("Error: ancestor %s has rank %s but descendant %s has rank %s\n" % (ancestor.name, current_rank, child.name, child.tx_level))
                else:
                    print("Error: ancestor %s has rank %s but descendant %s has rank %s" % (ancestor.name, current_rank, child.name, child.tx_level))
                errors += 1
        check_taxonomy_order(child, file, current_rank)


def find_phy_under_tax(tre, file):
    for node in tre.traverse(strategy="preorder"):
        if node.up:
            if node.up.ph_tx == "TX" and node.ph_tx == "PH":
                file.write("PH node %s, %s sits under TX node %s, %s\n" % (node.name, node.tx_level, node.up.name, node.up.tx_level))


def get_tx_levels(tre):
    tx_levels = set()
    for node in tre.traverse(strategy="preorder"):
        tx_levels.add(node.tx_level)
    return tx_levels


def check_zero_dates(tre):
    for node in tre.traverse():
        if node.is_leaf() and node.date is None:
            print(node.name, "is leaf with date None")
        if node.is_leaf() and node.date is not None and node.date > 0:
            print(node.name, "is leaf with date > 0")
        if not node.is_leaf() and node.date == 0:
            print(node.name, "is interior node with date 0")


# mrad = most recent ancestor date of the current node - if the current node is longer ago than this, that's bad
def check_inconsistent_dates(parent, fout, mrad=None, rounded=True):
    if mrad is None:
        if not parent.date:
            raise Exception("Require root date.")
        next_mrad = parent.date

    else:
        next_mrad = mrad
        # if we have a date, check it against the most recent ancestor
        if parent.date:
            if rounded:
                if round(parent.date,6) >= round(mrad,6):
                    fout.write("Found inconsistent date at %s: date of %f is later than %f\n" % (parent.name, parent.date, mrad))
                else:
                    next_mrad = parent.date
            else:
                if parent.date >= mrad:
                    fout.write("Found inconsistent date at %s: date of %f is later than %f\n" % (parent.name, parent.date, mrad))
                else:
                    next_mrad = parent.date

    for child in parent.children:
        check_inconsistent_dates(child, fout, next_mrad, rounded)
