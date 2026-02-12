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


from chronosynth import chronogram as cg
import ete4
import tree_fixing
import tree_dating
from taxonomy_utils import tx_levels
from taxonomy_utils import get_genus_and_species
import json
import re

import logging
logger = logging.getLogger(__name__)


def load_metadata(date_cache="chronosynth_date_info/node_ages.json",
                  annotations="opentree16.1_tree/annotations.json",
                  taxonomy="ott3.7.3/taxonomy.tsv",
                  force_dates_refresh=False,
                  force_no_dates_refresh=False):
    """Load metadata.
    Takes:
    date_cache -- string specifying path for a JSON file which will be used by Chronosynth to cache date information from phylogenies.
    annotations -- string specifying path for the annotations.json file from a Open Tree tree release.
    taxonomy -- string specifying path for the taxonomy.tsv file from an Open Tree Taxonomy release.

    Returns:
    dates -- JSON dictionary containing dates from phylogenies.
    phylogeny_nodes -- set of OTT ids and ancestor node ids for nodes in the Open Tree of Life that are in an underlying phylogeny,
                       rather than only in taxonomy.
    taxa -- dictionary; keys are integer ids from the Open Tree Taxonomy;
                        values are a tuple: (taxonomic rank string e.g. 'species', True/False flag specifying whether the taxon is extinct).
    """

    # get all dates from phylesystem studies
    if date_cache is not None:
        dates = cg.build_synth_node_source_ages(cache_file_path=date_cache, force_reload=force_dates_refresh, force_no_reload=force_no_dates_refresh)
        sources_to_delete = set(["ot_1250@tree2"])
        deletions = []
        for ott_name in dates["node_ages"]:
            for i, source in enumerate(dates["node_ages"][ott_name]):
                if source["source_id"] in sources_to_delete:
                    deletions.append((ott_name, i))

        deletions.sort(reverse=True)

        for ott_name, i in deletions:
            del dates["node_ages"][ott_name][i]
            if len(dates["node_ages"][ott_name]) == 0:
                del dates["node_ages"][ott_name]


    # get phylogeny annotations - if a node is not here, it is only from taxonomy
    if annotations is not None:
        annotations_dict = json.load(open(annotations,"r"))
        phylogeny_nodes = set(annotations_dict["nodes"].keys())

    # get dictionary of taxa labels (e.g. "family", "genus", "species") by OTT uid
    taxa = {}
    taxonomy_file = open(taxonomy,"r")
    for line in taxonomy_file.readlines()[1:]:
        components = line.split("|")

        flags = components[6].strip()
        flags = flags.split(",")
        extinct = False
        if "extinct" in flags or "extinct_inherited" in flags:
            extinct = True

        taxa[int(components[0].strip())] = (components[3].strip(), extinct)

    if date_cache is not None and annotations is not None:
        return dates, phylogeny_nodes, taxa
    elif annotations is not None:
        return phylogeny_nodes, taxa
    else:
        return taxa

def build_and_annotate_tree(phylogeny_nodes,
                            taxa,
                            tree_filename="opentree16.1_tree/labelled_supertree/labelled_supertree_ottnames.tre",
                            keep_all_dates=False,
                            ignore_extinct=True,
                            has_branch_lengths=False,
                            suppress_logging=False):

    """Build an ETE3 tree containing the Open Tree of Life and annotate it with extra information to be used in topology resolution and dating.
    Removes nodes marked as "extinct" or "extinct_inherited" in the Open Tree Taxonomy."

    Takes:
    phylogeny_nodes, taxa -- the outputs of load_metadata().
    tree_filename -- string specifying path to a newick tree to be loaded.

    Returns:
    tre -- the root node of the tree.
    """

    tree_string = open(tree_filename,"r").read()
    tree_string = tree_string.replace("''", "")  # remove '' in names (the tree loader will remove the quotes which wrap single-quoted names)
    tree_string = tree_string.replace('"', "")   # remove all double-quote characters
    tree_string = tree_string.replace('“', "")   # remove all double-quote characters
    tree_string = tree_string.replace('”', "")   # remove all double-quote characters
    tree_string = tree_string.replace("\u2009", "")   # weird unicode space - remove
    tree_string = tree_string.replace("\xa0", " ")    # weird unicode space - replace with normal space
    tree_string = tree_string.replace("о", "o")    # weird unicode o - replace with normal o
    tree_string = tree_string.replace("с", "c")    # weird unicode c - replace with normal c

    if not has_branch_lengths:
        tree_string = tree_string.replace(":", "")  # remove : in names

    tree_string = tree_string.replace(" ", "_")    # replace all spaces with underscores - so Newick format works without quoted names

    tre = ete4.Tree(tree_string, parser=1)

    if has_branch_lengths and not tre.name:
        tre.name = "mrca_root"

    if not suppress_logging:
        logger.info("ETE4 tree loaded. Beginning annotation")

    count = 0

    extinct_nodes = set()
    # add features to tree nodes: taxonomy level, degree, date (Mya)
    for node in tre.traverse(strategy='preorder'):
        # get OTT identifier (if the node is not only a most recent common ancestor) and look it up in the OT taxonomy
        node.name = node.name.strip("'")
        node.name = node.name.strip("_")
        node.name = node.name.replace(',', "")   # remove commas in node names
        node.name = node.name.replace('(', "")   # remove brackets in node names
        node.name = node.name.replace(')', "")   # remove brackets in node names
        node.name = node.name.replace('[', "")   # remove brackets in node names
        node.name = node.name.replace(']', "")   # remove brackets in node names

        if node.name[:4] == "mrca":
            tx_level = "mrca"
            ott_name = node.name
            if not ignore_extinct:
                node.add_prop("extinct", False)
        elif node.name[:11] == "uncultured_" or node.name[:11] == "Uncultured_" or node.name[:13] == "unidentified_" or node.name[:13] == "Unidentified_" or "intergeneric_hybrids" in node.name:
            extinct_nodes.add(node)
            ott_name = node.name
            continue
        else:
            if ' ' in node.name:
                ott_name = node.name.split(' ')[-1]
            else:
                ott_name = node.name.split('_')[-1]
            ott_uid = int(ott_name[3:])

            if ott_uid == 0:
                # must be reading in a tree with my imputed ott0000 species in it
                tx_level = "species"
            else:

                if ignore_extinct:
                    if taxa[ott_uid][1]:
                        # extinct or extinct_inherited
                        extinct_nodes.add(node)
                        continue
                else:
                    if taxa[ott_uid][1] or node.name[:1] == 'x_':
                        node.add_prop("extinct", True)
                    else:
                        node.add_prop("extinct", False)

                tx_level = taxa[ott_uid][0]

                if re.search(r".*_f\._.*", node.name) and tx_level == "no rank - terminal":
                    tx_level = "forma"
                    if not suppress_logging:
                        logger.info("Based on name, %s given rank 'forma' rather than 'no rank - terminal'." % node.name)
                elif re.search(r".*_var\._.*", node.name) and tx_level == "no rank - terminal":
                    tx_level = "variety"
                    if not suppress_logging:
                        logger.info("Based on name, %s given rank 'variety' rather than 'no rank - terminal'." % node.name)

        node.add_prop("tx_level", tx_level)

        if phylogeny_nodes is not None:
            if ott_name in phylogeny_nodes or has_branch_lengths:
                node.add_prop("ph_tx", "PH")
            else:
                node.add_prop("ph_tx", "TX")
                if not ignore_extinct and node.extinct:
                    # ignore extinct nodes from taxonomy
                    extinct_nodes.add(node)
                    continue

        if has_branch_lengths:
            node.add_prop("date", None)
            node.add_prop("imputed_date", False)
            node.add_prop("imputation_type", 0)

        if not ignore_extinct and node.extinct and node.date is None:
            # ignore extinct nodes with no date
            extinct_nodes.add(node)
            continue

        if tx_levels[tx_level] > 0 and tx_levels[tx_level] <= tx_levels["species group"]:
            genus, species = get_genus_and_species(node.name, ignore_extinct)
            if species == "extinct":
                extinct_nodes.add(node)
                continue
            node.add_prop("genus_name", genus)
            node.add_prop("species_name", species)

        node.add_prop("info", None)

        if not suppress_logging:
            count += 1
            if count % 10000 == 0:
                logger.info("%d nodes annotated" % count)

    if not suppress_logging:
        logger.info("Annotation complete. %d nodes annotated." % count)

    # remove things we noted as extinct (may include non-extinct nodes that we want to delete for other reasons)
    for node in extinct_nodes:
        if not suppress_logging:
            logger.info("Deleting extinct or unwanted node %s and tree below it." % node.name)
        tree_fixing.remove_tree_below(node)
        tree_fixing.remove_node_and_parents(node, False)

    if has_branch_lengths:
        tree_dating.compute_dates(tre)
        tree_dating.compute_branch_lengths(tre) # may modify existing branch lengths! - ensures tree is ultrametric

    return tre
