from taxonomy_utils import tx_levels
from taxonomy_utils import get_genus_and_species
import ete4

import logging
logger = logging.getLogger(__name__)

# Helper functions:

def create_node(name, forced_insert=False):
    """Create a new node with metadata."""
    new_node = ete4.Tree()
    new_node.name = name

    new_node.props["tx_level"] = "mrca"
    new_node.props["ancestral_rank"] = None
    new_node.props["desc_rank"] = None
    if forced_insert:
        new_node.props["ph_tx"] = "FI"
    else:
        new_node.props["ph_tx"] = "IN"
    new_node.props["info"] = None

    return new_node


def insert_node(parent, child, name, forced_insert=False):
    """Insert a new node, with metadata, between the parent and the child."""
    new_node = create_node(name, forced_insert)
    new_node.props["ancestral_rank"] = parent.props["ancestral_rank"]
    new_node.props["desc_rank"] = child.props["desc_rank"]

    new_node.add_child(child.detach())
    parent.add_child(new_node)

    return new_node


def insert_node_below(parent, name):
    """Insert a new node below the parent, and move the current children of the parent to the newly-inserted node"""
    new_node = create_node(name)

    new_node.props["ancestral_rank"] = parent.props["ancestral_rank"]

    highest_desc_rank = "species"
    for child in parent.children:
        if tx_levels[child.props["tx_level"]] > tx_levels[highest_desc_rank]:
            highest_desc_rank = child.props["tx_level"]
    new_node.props["desc_rank"] = highest_desc_rank

    orig_children = parent.children.copy()
    for child in orig_children:
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
    being deleted are also deleted. If subspecies_only=False, will clean up all parents with no children, regardless of rank.
    """
    if node is None:
        return

    parent = node.up
    if len(node.children) == 0:
        logger.debug("Deleted %s. Parent is %s" % (node.name, parent.name if parent else "None"))
        node.detach()
        del node
        if subspecies_only:
            if parent is not None and tx_levels[parent.props["tx_level"]] < tx_levels["species"]:
                remove_node_and_parents(parent, subspecies_only)
        else:
            remove_node_and_parents(parent, subspecies_only)


# Tree topology adjustment functions:

def remove_subspecies(tre, rng):
    """Remove any subspecies, varieties, formae etc from the tree. If there are subspecies with no species node, keep one of the subspecies
    and promote it to species rank.
    """

    # Prune away anything below species
    species_nodes = []
    species_names = set()
    for node in tre.traverse():
        if tx_levels[node.props["tx_level"]] == tx_levels["species"]:
            species_nodes.append(node)
            species_names.add(node.props["species_name"])

    for node in species_nodes:
        if len(node.children) > 0:
            logger.debug("Removing all nodes below species node %s." % (node.name))
            remove_tree_below(node)

    # Deal with remaining nodes of rank below species.
    # First get lists of such nodes, in dictionaries where keys are the species names.
    subsp_dict = {}
    for node in tre.traverse():
        if (node.props["tx_level"] == "no rank - terminal" or
                node.props["tx_level"] == "subspecies" or
                node.props["tx_level"] == "varietas" or
                node.props["tx_level"] == "variety" or
                node.props["tx_level"] == "forma" or
                node.props["tx_level"] == "infraspecificname"):
            if node.props["species_name"] not in subsp_dict:
                subsp_dict[node.props["species_name"]] = []
            subsp_dict[node.props["species_name"]].append(node)


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
                logger.info("Node %s kept and promoted to species from %s because it is the only example of species %s." % (subsp_dict[sp][0].name, subsp_dict[sp][0].props["tx_level"], sp))
                subsp_dict[sp][0].props["tx_level"] = "species (promoted)"
                subsp_dict[sp][0].props["ancestral_rank"] = "species (promoted)"
                subsp_dict[sp][0].props["desc_rank"] = "species (promoted)"

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
                    logger.info("Node %s kept and promoted to species from %s based on repetition of species name." % (subsp_dict[sp][i].name, subsp_dict[sp][i].props["tx_level"]))
                    subsp_dict[sp][i].props["tx_level"] = "species (promoted)"
                    subsp_dict[sp][i].props["ancestral_rank"] = "species (promoted)"
                    subsp_dict[sp][i].props["desc_rank"] = "species (promoted)"

                    sp_found = i
                    break

            # 3. finally, if there wasn't a good reason to keep a particular node from the duplicates, pick one at random to keep and promote it to species
            if sp_found is None:
                keep = subsp_dict[sp][rng.integers(len(subsp_dict[sp]))]
                logger.info("Node %s kept and promoted to species from %s based on random choice." % (keep.name, keep.props["tx_level"]))
                keep.props["tx_level"] = "species (promoted)"
                keep.props["ancestral_rank"] = "species (promoted)"
                keep.props["desc_rank"] = "species (promoted)"
            else:
                keep = subsp_dict[sp][i]


        # delete all the duplicates that didn't get chosen above, cleaning up below and above as well
        for node in subsp_dict[sp]:
            if node is None:
                continue
            if node is not keep:
                remove_tree_below(node)
                remove_node_and_parents(node)


def impute_species_into_empty_taxa(tre):
    """Finds empty higher-than-species taxa, and imputes a representative random species into them."""
    new_parents = []
    for node in tre.traverse():
        if node.is_leaf and tx_levels[node.props["tx_level"]] != tx_levels["species"]:
            new_parents.append(node)

    for node in new_parents:
        nm_parts = node.name.split('_')
        sp_name = ""
        for part in nm_parts[:-1]:
            sp_name += (part + "_")
        sp_name += "sp._ott0000"

        new_node = create_node(sp_name)
        genus, species = get_genus_and_species(sp_name)
        new_node.add_prop("genus_name", genus)
        new_node.add_prop("species_name", species)
        new_node.props["tx_level"] = "species (imputed)"
        new_node.props["ph_tx"] = "IM"

        node.add_child(new_node)

        logger.info("Added representative species %s as a child of %s node %s. %s %s" % (sp_name, node.props["tx_level"], node.name, genus, species))


def remove_nonspecies_leaves(tre):
    """Remove any leaf nodes that are of rank above species, e.g. empty genera."""
    to_remove = []
    for node in tre.traverse():
        if node.is_leaf and tx_levels[node.props["tx_level"]] != tx_levels["species"]:
            logger.info("Non-species leaf %s removed, rank %s, %s." % (node.name, node.props["tx_level"], node.props["ph_tx"]))
            to_remove.append(node)

    for node in to_remove:
        remove_node_and_parents(node, subspecies_only=False)


def key_to_node(key):
    if type(key) is tuple:
        return key[0]
    elif type(key) is ete4.core.tree.Tree:
        return key
    else:
        return None


def fix_polyphyly(tofix_dict, rng, expand_parent_backbones=False):
    """Go through dictionary of nodes to be "fixed" and move each into the backbone associated with it."""

    keys_to_remove = []
    for key in tofix_dict:
        # if there is no backbone, or only the root of the clade is in the backbone, we don't need to fix this
        backbone_size = len(tofix_dict[key][1])
        if backbone_size == 0 or (backbone_size == 1 and tofix_dict[key][1][0] is key_to_node(key)):
            keys_to_remove.append(key)
            continue

        # detach everything about to be moved from its parent
        inserts_size = len(tofix_dict[key][0])
        for i in range(inserts_size):
            tofix_dict[key][0][i].detach()

    for key in keys_to_remove:
        del tofix_dict[key]

    # Go through the candidate nodes in random order and insert them at random places into the backbone
    keys = list(tofix_dict.keys())
    rng.shuffle(keys) # NB: with a big list here, the number of random permutations is so great that most can never be generated!
    for key in keys:
        inserts = list(range(len(tofix_dict[key][0])))
        rng.shuffle(inserts)

        # dictionary of genera already seen so we can keep species of the same genus together
        genera_found = {}

        count = 0
        for i in inserts:
            node_to_move = tofix_dict[key][0][i]

            if (node_to_move.props["info"] != "OTH FIX" or                    # - only OTH FIX nodes might have multiple species from the same genus
                    tx_levels[node_to_move.props["tx_level"]] == tx_levels["no rank"] or        # - nodes of no rank don't have a genus
                    tx_levels[node_to_move.props["tx_level"]] >  tx_levels["species group"] or  # - if it has rank higher than species group, it can't be a species from a genus we have seen before
                    node_to_move.props["genus_name"] not in genera_found):    # - if it is not in the list then we haven't see it, duh

                # we haven't seen this genus before - so this node can be put in a random place in the backbone
                backbone_size = len(tofix_dict[key][1])
                child = tofix_dict[key][1][rng.integers(backbone_size)]
                if tx_levels[node_to_move.props["tx_level"]] == tx_levels["species"]:
                    genera_found[node_to_move.props["genus_name"]] = node_to_move
            else:
                # this is a species and we have seen its genus before - therefore this is a species in a non-monophyletic genus,
                # so put this node with the others of the same genus to ensure monophyly
                child = genera_found[node_to_move.props["genus_name"]]

            if type(key_to_node(key)) == ete4.core.tree.Tree and child is key_to_node(key):
                # the "child" chosen is the root node of the group - we want to insert a new node *below* this
                if len(child.children) < 2:
                    child.add_child(node_to_move)
                    new_internal_node = None
                else:
                    new_internal_node = insert_node_below(child, "mrcaimp")
                    child.add_child(node_to_move)
            else:
                # not the root node; insert above, unless the parent has only 1 child in which case just attach directly
                parent = child.up
                if len(parent.children) < 2:
                    parent.add_child(node_to_move)
                    new_internal_node = None
                else:
                    new_internal_node = insert_node(parent, child, "mrcaimp")
                    new_internal_node.add_child(node_to_move)

            # if we created a new node, add it to the backbone
            if new_internal_node:
                if expand_parent_backbones:
                    for root in tofix_dict[key][2]:
                        for level in tx_levels:
                            if tx_levels[level] < tx_levels[child.props["ancestral_rank"]]:
                                # if this is a backbone for a rank *below* that of the ancestral rank of the new node,
                                # the new node can be added to this backbone.
                                # for example, if the new node has ancestral rank family, that could be a potential
                                # backbone node for an OTH FIX node that is a subspecies, species, or genus, but not any higher.
                                if (root, level) in tofix_dict:
                                    tofix_dict[(root,level)][1].append(new_internal_node)
                else:
                    tofix_dict[key][1].append(new_internal_node)

            # add the node we moved to the backbone
            if expand_parent_backbones:
                for root in tofix_dict[key][2]:
                    for level in tx_levels:
                        if tx_levels[level] < tx_levels[child.props["ancestral_rank"]]:
                            # as above
                            if (root, level) in tofix_dict:
                                tofix_dict[(root,level)][1].append(node_to_move)
            else:
                tofix_dict[key][1].append(node_to_move)

            count += 1



def fix_polytomy(parent, rng):
    """Fix polytomy directly beneath given parent: choose uniformly at random from possible topologies."""
    if len(parent.children) <= 2:
        raise Exception("Error: Trying to fix a non-polytomy")

    # pick which children to move - all but 2
    nodes_to_move = parent.get_children()
    del nodes_to_move[rng.integers(len(nodes_to_move))]
    del nodes_to_move[rng.integers(len(nodes_to_move))]

    # detach them
    for node_to_move in nodes_to_move:
        node_to_move.detach()

    # a sibling is a node above which we might want to insert the node we are moving
    possible_siblings = list(parent.children)
    possible_siblings.append(parent)

    for i in range(len(nodes_to_move)):
        node_to_move = nodes_to_move[i]

        new_sibling = possible_siblings[rng.integers(len(possible_siblings))]

        if new_sibling is parent:
            new_node = create_node("mrcapoly")
            new_node.props["ancestral_rank"] = parent.props["ancestral_rank"]
            new_node.props["desc_rank"] = parent.props["desc_rank"]
            current_children = list(parent.children)
            for child in current_children:
                new_node.add_child(child.detach())

            parent.add_child(new_node)
            parent.add_child(node_to_move)

        else:

            new_node = insert_node(new_sibling.up, new_sibling, "mrcapoly")
            new_node.add_child(node_to_move)

            if tx_levels[node_to_move.props["desc_rank"]] > tx_levels[new_node.props["desc_rank"]]:
                new_node.props["desc_rank"] = node_to_move.props["desc_rank"]
            if tx_levels[node_to_move.props["desc_rank"]] > tx_levels[new_node.up.props["desc_rank"]]:
                new_node.up.props["desc_rank"] = node_to_move.props["desc_rank"]

        possible_siblings.append(new_node)
        possible_siblings.append(node_to_move)


def fix_all_polytomies(tre, rng):
    """Fix at random all polytomies in the tree."""
    polytomies = []
    for node in tre.traverse(strategy='preorder'):
        if len(node.children) > 2:
            polytomies.append(node)

    for node in polytomies:
        fix_polytomy(node, rng)


def delete_one_child_nodes(tre, maintain_branch_lengths=False):
    """Strip out nodes with only one child. If maintain_branch_lengths=True, add the branch length above
    the deleted node to the branch below it.
    """
    one_child_nodes = []
    for node in tre.traverse(strategy="preorder"):
        if len(node.children) == 1:
            one_child_nodes.append(node)

    for node in one_child_nodes:
        if not node.up:
            # if root no has one one child, delete it and move the root to the original root's child
            tre = tre.children[0]
            tre.detach()
        else:
            if maintain_branch_lengths:
                node.children[0].dist += node.dist
            node.up.add_child(node.children[0].detach())
            node.detach()
        del node

    return tre


def strip_birds(tre, ejm_birds_filename="config/OTT_crosswalk_2024.csv"):
    import csv

    desired_ottids = set()

    with open(ejm_birds_filename, newline='') as csvfile:
        rdr = csv.reader(csvfile)
        for idx, line in enumerate(rdr):
            if idx == 0:
                # first line has column headings
                continue
            desired_ottids.add(int(line[8]))

    aves_root = None
    for node in tre.search_nodes(name="Aves_ott81461"):
        aves_root = node
        break

    to_remove = []
    if aves_root is not None:
        for node in aves_root.traverse():
            if tx_levels[node.props["tx_level"]] == tx_levels["species"] or tx_levels[node.props["tx_level"]] == tx_levels["subspecies"]:
                nm_parts = node.name.split('_')
                ottid = int(nm_parts[-1][3:])
                if ottid not in desired_ottids:
                    to_remove.append(node)
                else:
                    desired_ottids.remove(ottid)

    for node in to_remove:
        if tx_levels[node.props["tx_level"]] == tx_levels["subspecies"]:
            remove_node_and_parents(node, subspecies_only=True)

    for node in to_remove:
        if tx_levels[node.props["tx_level"]] == tx_levels["species"] and len(node.children) == 0:
            remove_node_and_parents(node, subspecies_only=False)


def strip_turtles(tre, turtles_filename="config/turtle_checklist_ott.csv"):
    import csv

    desired_ottids = set()

    with open(turtles_filename, newline='') as csvfile:
        rdr = csv.reader(csvfile)
        for idx, line in enumerate(rdr):
            if idx == 0:
                # first line has column headings
                continue
            desired_ottids.add(int(line[1]))

    turtles_root = None
    for node in tre.search_nodes(name="Testudines_ott639666"):
        turtles_root = node
        break

    to_remove = []
    if turtles_root is not None:
        for node in turtles_root.traverse():
            if tx_levels[node.props["tx_level"]] == tx_levels["species"] or tx_levels[node.props["tx_level"]] == tx_levels["subspecies"]:
                nm_parts = node.name.split('_')
                ottid = int(nm_parts[-1][3:])
                if ottid not in desired_ottids:
                    to_remove.append(node)
                else:
                    desired_ottids.remove(ottid)

    for node in to_remove:
        if tx_levels[node.props["tx_level"]] == tx_levels["subspecies"]:
            remove_node_and_parents(node, subspecies_only=True)

    for node in to_remove:
        if tx_levels[node.props["tx_level"]] == tx_levels["species"] and len(node.children) == 0:
            remove_node_and_parents(node, subspecies_only=False)


def fix_taxonomy_ordering(tre, filename="config/taxonomy_fixes.csv"):
    """Adjust the tree to ensure correct taxonomic ordering, using hand-curated config file.
    Type 0: remove the rank from the ancestor node.
    Type 1: remove the rank from the descendant node.
    Type 2: the inconsistent descendant node has only one child, so just delete the node. Assume no branch lengths are set yet.
    """

    import csv

    with open(filename, newline='') as csvfile:
        rdr = csv.reader(csvfile)
        for idx, line in enumerate(rdr):
            if idx == 0:
                # first line has column headings
                continue
            if line[2] == "0":
                node_to_find = None
                for node in tre.search_nodes(name=line[0]):
                    node_to_find = node
                    break
                if node_to_find is None:
                    raise Exception("Node in taxonomy ordering config is not in tree.")

                node.props["tx_level"] = "no rank"

            elif line[2] == "1":
                node_to_find = None
                for node in tre.search_nodes(name=line[1]):
                    node_to_find = node
                    break
                if node_to_find is None:
                    raise Exception("Node in taxonomy ordering config is not in tree.")

                node.props["tx_level"] = "no rank"

            if line[2] == "2":
                node_to_find = None
                for node in tre.search_nodes(name=line[1]):
                    node_to_find = node
                    break
                if node_to_find is None:
                    raise Exception("Node in taxonomy ordering config is not in tree.")

                # assume this isn't the root node
                node_to_find.up.add_child(node.children[0].detach())
                node_to_find.detach()
                del node_to_find


def forced_taxa_moves(tre, filename="config/forced_taxa_moves.csv"):
    """Force some taxa to be sisters of other taxa.
    Initially created to force Eukaryota to be a sister of Archaea.
    """

    import csv

    with open(filename, newline='') as csvfile:
        rdr = csv.reader(csvfile)
        for idx, line in enumerate(rdr):
            if idx == 0:
                # first line has column headings
                continue

            node_to_move = None
            for node in tre.search_nodes(name=line[0]):
                node_to_move = node
                break
            if node_to_move is None:
                raise Exception("Node in forced taxa moves config is not in tree.")

            sister_node = None
            for node in tre.search_nodes(name=line[1]):
                sister_node = node
                break
            if sister_node is None:
                raise Exception("Node in forced taxa moves config is not in tree.")

            new_internal_node = insert_node(sister_node.up, sister_node, "mrcaimp", forced_insert=True)
            new_internal_node.add_child(node_to_move.detach())

            # pretend this was phylogeny so it doesn't get moved again
            node_to_move.props["ph_tx"] = "PH"
