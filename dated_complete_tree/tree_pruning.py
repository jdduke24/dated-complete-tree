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
from .taxonomy_utils import get_genus_and_species

from .tree_fixing_utils import create_node
from .tree_fixing_utils import insert_node
from .tree_fixing_utils import remove_tree_below
from .tree_fixing_utils import remove_node_and_parents
from .tree_fixing_utils import open_config_csv

import logging
logger = logging.getLogger(__name__)


# Tree pruning functions

def remove_based_on_props(tre, props_list):
    nodes_to_remove = set()
    for node in tre.traverse(strategy='preorder'):
        for prop in props_list:
            if prop in node.props:
                nodes_to_remove.add(node)

    for node in nodes_to_remove:
        logger.info("Deleting extinct or unwanted node %s and tree below it." % node.name)
        remove_tree_below(node)
        remove_node_and_parents(node, False)


def remove_subspecies(tre, rng):
    """Remove any subspecies, varieties, formae etc from the tree. If there are subspecies with no species node, keep one of the subspecies
    and promote it to species rank.
    """

    # Prune away anything below species
    species_nodes = []
    species_names = set()
    for node in tre.traverse(strategy="preorder"):
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
    for node in tre.traverse(strategy="preorder"):
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
    for node in tre.traverse(strategy="preorder"):
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


def strip_birds(tre, ejm_birds_filename=None):
    import csv

    desired_ottids = set()

    with open_config_csv(ejm_birds_filename, "birds.csv") as csvfile:
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
        for node in aves_root.traverse(strategy="preorder"):
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


def strip_turtles(tre, turtles_filename=None):
    import csv

    desired_ottids = set()

    with open_config_csv(turtles_filename, "turtles.csv") as csvfile:
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
        for node in turtles_root.traverse(strategy="preorder"):
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
