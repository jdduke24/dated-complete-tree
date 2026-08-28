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


import ete4

from .taxonomy_utils import tx_levels

from .tree_fixing_utils import create_node
from .tree_fixing_utils import insert_node
from .tree_fixing_utils import insert_node_below
from .tree_fixing_utils import remove_node_and_parents
from .tree_fixing_utils import open_config_csv

import logging
logger = logging.getLogger(__name__)


def remove_nonspecies_leaves(tre):
    """Remove any leaf nodes that are of rank above species, e.g. empty genera."""
    to_remove = []
    for node in tre.traverse(strategy="preorder"):
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


def fix_remaining_polytomies(tre, rng):
    """Fix at random all polytomies in the tree."""
    polytomies = []
    for node in tre.traverse(strategy="preorder"):
        if len(node.children) > 2:
            polytomies.append(node)

    for node in polytomies:
        fix_polytomy(node, rng)


def delete_one_child_nodes(tre, maintain_branch_lengths=True):
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
            if maintain_branch_lengths and node.children[0].dist and node.dist:
                node.children[0].dist += node.dist
            node.up.add_child(node.children[0].detach())
            node.detach()
        del node

    return tre


def forced_taxa_moves(tre, filename=None):
    """Force some taxa to be sisters of other taxa.
    Initially created to force Eukaryota to be a sister of Archaea.
    """

    import csv

    with open_config_csv(filename, "forced_taxa_moves.csv") as csvfile:
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
