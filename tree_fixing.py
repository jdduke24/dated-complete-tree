from taxonomy_utils import tx_levels
import ete3
import random

import logging
logger = logging.getLogger(__name__)

# Helper functions:

# create a new node
def create_node(name):
    new_node = ete3.TreeNode(name=name)
    new_node.add_feature("tx_level", "mrca")
    new_node.add_feature("desc_year", None)
    new_node.add_feature("ancestral_rank", None)
    new_node.add_feature("desc_rank", None)
    new_node.add_feature("ph_tx", "IN")
    new_node.add_feature("date", None)
    new_node.add_feature("imputed_date", False)
    new_node.add_feature("info", None)

    return new_node


# insert a new node between the parent and the child
def insert_node(parent, child, name):
    new_node = create_node(name)
    new_node.ancestral_rank = parent.ancestral_rank
    new_node.desc_rank = child.desc_rank

    new_node.add_child(child.detach())
    parent.add_child(new_node)

    return new_node


def remove_tree_below(node, root=True, to_delete=None):
    """Delete the whole tree beneath this node (not including this node)."""
    if node is None:
        return

    if root:
        to_delete = []

    for child in node.children:
        remove_tree_below(child, False, to_delete)
        to_delete.append(child)

    if root:
        for nd in to_delete:
            logger.debug("Removing children of %s: Deleted %s." % (node.name, nd.name))
            nd.detach()
            del nd


def remove_node_and_parents(node, subspecies_only=True):
    """Delete this node, and clean up the tree above it, such any mrca, no-rank, or below-species parents that have no children other than the chain
    being deleted are also deleted.
    """
    if node is None:
        return

    parent = node.up
    if len(node.children) == 0:
        logger.debug("Deleted %s. Parent is %s" % (node.name, parent.name if parent else "None"))
        node.detach()
        del node
        if subspecies_only:
            if parent is not None and tx_levels[parent.tx_level] < tx_levels["species"]:
                remove_node_and_parents(parent, subspecies_only)
            elif parent is not None and len(parent.children) == 0 and parent.date is None:
                logger.info("While fixing, set %s to date 0 from None" % (parent.name))
                parent.date = 0
        else:
            remove_node_and_parents(parent, subspecies_only)


def remove_subspecies(tre):
    """Remove any subspecies, varieties, formae etc from the tree. If there are subspecies with no species node, keep one of the subspecies
    and promote it to species rank.
    """

    # Prune away anything below species
    species_nodes = []
    species_names = set()
    for node in tre.traverse():
        if node.tx_level == "species":
            species_nodes.append(node)
            species_names.add(node.species_name)

    for node in species_nodes:
        logger.debug("Removing all nodes below species node %s." % (node.name))
        remove_tree_below(node)
        node.date = 0

    # Deal with remaining nodes of rank below species.
    # First get lists of such nodes, in dictionaries where keys are the species names.
    subsp_dict = {}
    for node in tre.traverse():
        # if node.tx_level == "no rank - terminal":
        #     if node.species_name not in norank_terminals:
        #         norank_terminals[node.species_name] = []
        #     norank_terminals[node.species_name].append(node)

        if (node.tx_level == "no rank - terminal" or
                node.tx_level == "subspecies" or
                node.tx_level == "varietas" or
                node.tx_level == "variety" or
                node.tx_level == "forma" or
                node.tx_level == "infraspecificname"):
            if node.species_name not in subsp_dict:
                subsp_dict[node.species_name] = []
            subsp_dict[node.species_name].append(node)


    # For other sub-species ranks:
    # Where there are multiple subspecies (or varieties etc) of the same species, but no node representing the species itself,
    # remove all but one of the subspecies and "promote" the remaining subspecies into a species.
    for sp in subsp_dict:
        logger.debug("List for %s has length %d." % (sp, len(subsp_dict[sp])))
        if sp in species_names:
            # if there is a species node with this name elsewhere in the tree, delete all subspecies
            logger.info("Species %s already exists; deleting entire list of subspecies." % (sp))
            keep = None
        else:
            # otherwise, choose one of the subspecies nodes to represent this species

            # 1. only one node has this species name, so make that a species
            if len(subsp_dict[sp]) == 1:
                logger.info("Node %s kept and promoted to species from %s because it is the only example of species %s." % (subsp_dict[sp][0].name, subsp_dict[sp][0].tx_level, sp))
                subsp_dict[sp][0].tx_level = "species (promoted)"
                subsp_dict[sp][0].ancestral_rank = "species (promoted)"
                subsp_dict[sp][0].desc_rank = "species (promoted)"
                if subsp_dict[sp][0].date is None:
                    subsp_dict[sp][0].date = 0
                continue

            sp_found = None

            # 2. otherwise, look to see if one of the nodes has a name with repeated species name - Genus_speciesname_speciesname - if so, keep that
            for i in range(len(subsp_dict[sp])):
                nm = subsp_dict[sp][i].name.replace("subsp.","")
                nm = nm.replace("var.","")
                nm = nm.replace("f.","")
                nm = nm.replace("  "," ")
                nm = nm.replace("__","_")

                if ' ' in nm:
                    nm_parts = nm.split(' ')[:-1]
                else:
                    nm_parts = nm.split('_')[:-1]

                if len(nm_parts) > 2 and nm_parts[1] == nm_parts[2]:
                    logger.info("Node %s kept and promoted to species from %s based on repetition of species name." % (subsp_dict[sp][i].name, subsp_dict[sp][i].tx_level))
                    subsp_dict[sp][i].tx_level = "species (promoted)"
                    subsp_dict[sp][i].ancestral_rank = "species (promoted)"
                    subsp_dict[sp][i].desc_rank = "species (promoted)"
                    if subsp_dict[sp][i].date is None:
                        subsp_dict[sp][i].date = 0
                    sp_found = i
                    break

            # 3. finally, if there wasn't a good reason to keep a particular node from the duplicates, pick one at random to keep and promote it to species
            if sp_found is None:
                keep = random.choice(subsp_dict[sp])
                logger.info("Node %s kept and promoted to species from %s based on random choice." % (keep.name, keep.tx_level))
                keep.tx_level = "species (promoted)"
                keep.ancestral_rank = "species (promoted)"
                keep.desc_rank = "species (promoted)"
                if keep.date is None:
                    keep.date = 0
            else:
                keep = subsp_dict[sp][i]


        # delete all the duplicates that didn't get chosen above, cleaning up below and above as well
        for node in subsp_dict[sp]:
            if node is None:
                continue
            if node is not keep:
                remove_tree_below(node)
                remove_node_and_parents(node)


def fix_polyphyly(tofix_dict, expand_parent_backbones=False):
    """Go through dictionary of nodes to be "fixed" and move each into the backbone associated with it."""

    keys_to_remove = []
    for key in tofix_dict:
        # if there is no backbone, we don't need to (can't) fix this
        backbone_size = len(tofix_dict[key][1])
        if backbone_size == 0:
            keys_to_remove.append(key)
            continue

        # detach everything about to be moved from its parent
        inserts_size = len(tofix_dict[key][0])
        for i in range(inserts_size):
            if len(tofix_dict[key][0][i].up.children) == 1:
                if not tofix_dict[key][0][i].up.date:
                    logger.info("While fixing, set %s to date 0 from None" % (tofix_dict[key][0][i].up.name))
                    tofix_dict[key][0][i].up.date = 0

            tofix_dict[key][0][i].detach()

    for key in keys_to_remove:
        del tofix_dict[key]

    # Go through the candidate nodes in random order and insert them at random places into the backbone
    keys = list(tofix_dict.keys())
    random.shuffle(keys) # NB: with a big list here, the number of random permutations is so great that most can never be generated!
    for key in keys:
        inserts = list(range(len(tofix_dict[key][0])))
        random.shuffle(inserts)

        # dictionary of genera already seen so we can keep species of the same genus together
        genera_found = {}

        count = 0
        for i in inserts:
            node_to_move = tofix_dict[key][0][i]

            if (node_to_move.info != "OTH FIX" or                    # - only OTH FIX nodes might have multiple species from the same genus
                    tx_levels[node_to_move.tx_level] == tx_levels["no rank"] or        # - nodes of no rank don't have a genus
                    tx_levels[node_to_move.tx_level] >  tx_levels["species group"] or  # - if it has rank higher than species group, it can't be a species from a genus we have seen before
                    node_to_move.genus_name not in genera_found):    # - if it is not in the list then we haven't see it, duh

                # we haven't seen this genus before - so this node can be put in a random place in the backbone
                backbone_size = len(tofix_dict[key][1])
                child = tofix_dict[key][1][random.randint(0,backbone_size-1)]
                if node_to_move.tx_level == "species":
                    genera_found[node_to_move.genus_name] = node_to_move
            else:
                # we have seen this genus before - so put this node with the others of the same genus to ensure monophyly
                child = genera_found[node_to_move.genus_name]

            parent = child.up

            if len(parent.children) >= 2:
                new_internal_node = insert_node(parent, child, "mrcat")
                new_internal_node.add_child(node_to_move)
                if expand_parent_backbones:
                    for root in tofix_dict[key][2]:
                        for level in tx_levels:
                            if tx_levels[level] < tx_levels[child.ancestral_rank]:
                                # if this is a backbone for a rank *below* that of the ancestral rank of the new node,
                                # the new node can be added to this backbone.
                                # for example, if the new node has ancestral rank family, that could be a potential
                                # backbone node for an OTH FIX node that is a subspecies, species, or genus, but not any higher.
                                if (root, level) in tofix_dict:
                                    tofix_dict[(root,level)][1].append(new_internal_node)
                else:
                    tofix_dict[key][1].append(new_internal_node)
            else:
                parent.add_child(node_to_move)

            if expand_parent_backbones:
                for root in tofix_dict[key][2]:
                    for level in tx_levels:
                        if tx_levels[level] < tx_levels[child.ancestral_rank]:
                            # as above
                            if (root, level) in tofix_dict:
                                tofix_dict[(root,level)][1].append(node_to_move)
            else:
                tofix_dict[key][1].append(node_to_move)

            count += 1


# fix polytomy directly beneath given parent
def fix_polytomy(parent):
    if len(parent.children) <= 2:
        raise("Error: Trying to fix a non-polytomy")

    # pick which children to move
    nodes_to_move = random.sample(parent.children, len(parent.children)-2)

    # detach them
    for node_to_move in nodes_to_move:
        node_to_move.detach()

    # a sibling is a node above which we might want to insert the node we are moving
    possible_siblings = list(parent.children)
    possible_siblings.append(parent)

    for i in range(len(nodes_to_move)):
        node_to_move = nodes_to_move[i]

        new_sibling = random.choice(possible_siblings)

        if new_sibling is parent:
            new_node = create_node("mrcap")
            new_node.ancestral_rank = parent.ancestral_rank
            new_node.desc_rank = parent.desc_rank
            current_children = list(parent.children)
            for child in current_children:
                new_node.add_child(child.detach())

            parent.add_child(new_node)
            parent.add_child(node_to_move)

        else:

            new_node = insert_node(new_sibling.up, new_sibling, "mrcap")
            new_node.add_child(node_to_move)

        possible_siblings.append(new_node)
        possible_siblings.append(node_to_move)


def fix_all_polytomies(tre):
    polytomies = []
    for node in tre.traverse(strategy='preorder'):
        if len(node.children) > 2:
            polytomies.append(node)

    for node in polytomies:
        fix_polytomy(node)


def delete_one_child_nodes(tre):
    one_child_nodes = []
    for node in tre.traverse(strategy="preorder"):
        if len(node.children) == 1:
            one_child_nodes.append(node)

    for node in one_child_nodes:
        node.up.add_child(node.children[0].detach())
        node.detach()
        del node
