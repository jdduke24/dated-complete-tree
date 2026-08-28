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


from .taxonomy_utils import tx_levels

import ete4
from importlib.resources import files

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


def key_to_node(key):
    if type(key) is tuple:
        return key[0]
    elif type(key) is ete4.core.tree.Tree:
        return key
    else:
        return None


def open_config_csv(filename, default_basename):
    """Open a config CSV in text mode with ``newline=''``.

    If ``filename`` is ``None``, open the file ``default_basename`` from the
    package's bundled ``config/`` directory. Otherwise treat ``filename`` as a
    user-supplied filesystem path.
    """
    if filename is None:
        return (files(__package__) / "config" / default_basename).open(newline="")
    return open(filename, newline="")
