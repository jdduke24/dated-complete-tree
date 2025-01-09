from chronosynth import chronogram as cg
import ete3
import tree_fixing
from taxonomy_utils import tx_levels
from taxonomy_utils import get_genus_and_species
import json
import re

import logging
logger = logging.getLogger(__name__)

##############################################
# Get metadata for tree -
#  - dates from Chronosynth/DateLife
#  - phylogeny annotations from Open Tree of Life
#  - taxa labels from Open Tree Taxonomy

def load_metadata(date_cache="chronosynth_date_info/node_ages.json",
                      phylogeny="opentree14.9_tree/annotations.json",
                      taxonomy="ott3.6/taxonomy.tsv",
                      descr_dates="oz_data/descr_dates.csv"):
    """Load metadata.
    Takes:
    date_cache -- string specifying path for a JSON file which will be used by Chronosynth to cache date information from phylogenies.
    phylogeny -- string specifying path for the annotations.json file from a Open Tree tree release.
    taxonomy -- string specifying path for the taxonomy.tsv file from an Open Tree Taxonomy release.
    desc_dates - string specifying a CSV file containing the dates of scientific descriptions, by OTT id.

    Returns:
    dates -- JSON dictionary containing dates from phylogenies.
    phylogeny_nodes -- set of OTT ids and ancestor node ids for nodes in the Open Tree of Life that are in an underlying phylogeny,
                       rather than only in taxonomy.
    taxa -- dictionary; keys are ids from the Open Tree Taxonomy; values are a tuple with the taxonomic level string (e.g. 'species') and a True/False
            flag specifying whether the taxon is extinct.
    descr_years -- dictionary: keys are OTT ids as integers; values are integer years of scientific description of the taxon, if available.
    """

    # get all dates from phylesystem studies
    dates = cg.build_synth_node_source_ages(cache_file_path=date_cache)

    # get phylogeny annotations - if a node is not here, it is only from taxonomy
    annotations = json.load(open(phylogeny,"r"))
    phylogeny_nodes = set(annotations["nodes"].keys())

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

    descr_years = {}
    descr_file = open(descr_dates,"r")
    for line in descr_file.readlines():
        components = line.split(",")
        descr_years[int(components[0].strip())] = int(components[2].strip())

    return dates, phylogeny_nodes, taxa, descr_years

##############################################


##############################################
# load in whole tree and write out subset of it for analysis
# (would be better if we can get a ete tree directly rather than going via dendropy)

def write_subtree(supertree_filename="opentree14.9_tree/labelled_supertree/labelled_supertree_ottnames.tre",
                 subtree_root_node="Passeriformes ott1041547",
                 output_filename="my_unresolved_passeriformes.tre"):

    import dendropy

    # get subtree from downloaded full OTL, including taxa labels and OTT uids, and write out in newick form
    supertree_str = open(supertree_filename,"r").read()
    supertree = dendropy.Tree.get_from_string(supertree_str, schema='newick')

    supertree_copy = dendropy.Tree(supertree)
    new_root_node = supertree_copy.find_node_with_label(subtree_root_node)

    new_root_node.parent_node = None
    subtree = dendropy.Tree(seed_node=new_root_node)

    subtree.write(path=output_filename,schema="newick")

##############################################


##############################################
# create ete3 tree from my subtree

def build_and_annotate_tree(dates,
                            phylogeny_nodes,
                            taxa,
                            descr_years,
                            tree_filename="opentree14.9_tree/labelled_supertree/labelled_supertree_ottnames.tre"):
    """Build an ETE3 tree containing the Open Tree of Life and annotate it with extra information to be used in topology resolution and dating.
    Removes nodes marked as "extinct" or "extinct_inherited" in the Open Tree Taxonomy."

    Takes:
    dates, phylogeny_nodes, taxa, descr_years -- the outputs of load_metadata().
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

    tree_string = tree_string.replace(" ", "_")    # replace all spaces with underscores - so Newick format works without quoted names

    tre = ete3.Tree(tree_string,format=1,quoted_node_names=True)

    logger.info("ETE3 tree loaded. Beginning annotation")

    count = 0

    fout = open("names_all.csv",'w')

    extinct_nodes = set()
    # add features to tree nodes: taxonomy level, degree, date (Mya, *not* branch length)
    for node in tre.traverse(strategy='preorder'):
        # get OTT identifier (if the node is not only a most recent common ancestor) and look it up in the OT taxonomy
        node.name = node.name.strip("'")
        node.name = node.name.strip("_")
        descr_year = None
        if node.name[:4] == "mrca":
            tx_level = "mrca"
            ott_name = node.name
        elif node.name[:14] == "uncultured_ott" or node.name[:16] == "unidentified_ott":
            extinct_nodes.add(node)
            ott_name = node.name
            continue
        else:
            if ' ' in node.name:
                ott_name = node.name.split(' ')[-1]
            else:
                ott_name = node.name.split('_')[-1]
            ott_uid = int(ott_name[3:])

            if taxa[ott_uid][1]:
                # extinct, ignore
                extinct_nodes.add(node)
                continue

            tx_level = taxa[ott_uid][0]

            if re.search(r".*_f\._.*", node.name) and tx_level == "no rank - terminal":
                tx_level = "forma"
                logger.info("Based on name, %s given rank 'forma' rather than 'no rank - terminal'." % node.name)
            elif re.search(r".*_var\._.*", node.name) and tx_level == "no rank - terminal":
                tx_level = "variety"
                logger.info("Based on name, %s given rank 'variety' rather than 'no rank - terminal'." % node.name)

            if ott_uid in descr_years:
                descr_year = descr_years[ott_uid]

        node.add_feature("tx_level", tx_level)
        node.add_feature("descr_year", descr_year)

        if ott_name in phylogeny_nodes:
            node.add_feature("ph_tx", "PH")
            #del phylogeny_nodes[ott_name]
        else:
            node.add_feature("ph_tx", "TX")

        date = None
        if ott_name in dates['node_ages']:
            ages = [float(source['age']) for source in dates['node_ages'][ott_name]]
            ages.sort()

            num_ages = len(ages)
            midpoint = int((num_ages-1)/2)
            if num_ages % 2 == 0:
                # assume num_ages > 0 or we wouldn't have got here
                median_age = (ages[midpoint] + ages[midpoint+1]) / 2
            else:
                median_age = ages[midpoint]

            if node.is_leaf():
                logger.warning("Warning: %s is dated leaf of age %f" % (node.name, median_age))

            if not node.is_leaf() and median_age < 0.000001:
                logger.warning("Warning: computed date of zero on interior node %s; instead, setting date for this node to None" % node.name)
            else:
                date = median_age
        elif node.is_leaf():
            date = 0
        node.add_feature("date", date)
        node.add_feature("imputed_date", False)

        if tx_levels[tx_level] > 0 and tx_levels[tx_level] <= tx_levels["species group"]:
            genus, species = get_genus_and_species(node.name,fout)
            if species == "extinct":
                extinct_nodes.add(node)
                continue
            node.add_feature("genus_name", genus)
            node.add_feature("species_name", species)

        node.add_feature("info", None)

        count += 1
        if count % 10000 == 0:
            logger.info("%d nodes annotated" % count)

    logger.info("Annotation complete. %d nodes annotated." % count)
    fout.close()

    # remove things we noted as extinct (may include non-extinct nodes that we want to delete for other reasons)
    for node in extinct_nodes:
        tree_fixing.remove_tree_below(node)
        tree_fixing.remove_node_and_parents(node, False)


    def add_anc_ranks(parent, current_rank=None, plants=False, animals=False, domain=None, kingdom=None, phylum=None, clas=None, order=None, family=None):
        """Add ancestral ranks. If this node has no rank, this is the rank of the first node above it which
        does have a rank. If this node has a rank, that rank is both its ancestral and descendant rank.
        Also annotate every node with the domain, kingdom, phylum, class, order and family it is under (if any).
        """
        if current_rank is None:
            # only occurs at root (first call of function)
            current_rank = parent.tx_level
            parent.add_feature("ancestral_rank", parent.tx_level)

        if tx_levels[parent.tx_level] > 0 or tx_levels[current_rank] < 0:
            current_rank = parent.tx_level

        if parent.tx_level == "domain":
            domain = parent.name
        elif parent.tx_level == "kingdom":
            kingdom = parent.name
        elif parent.tx_level == "phylum":
            phylum = parent.name
        elif parent.tx_level == "class":
            clas = parent.name
        elif parent.tx_level == "order":
            order = parent.name
        elif parent.tx_level == "family":
            family = parent.name

        parent.add_feature("domain", domain)
        parent.add_feature("kingdom", kingdom)
        parent.add_feature("phylum", phylum)
        parent.add_feature("clas", clas)
        parent.add_feature("order", order)
        parent.add_feature("family", family)

        for child in parent.children:
            if tx_levels[child.tx_level] < 0:
                if plants and current_rank == "section":
                    child.add_feature("ancestral_rank", "plantsection")
                elif plants and current_rank == "subsection":
                    child.add_feature("ancestral_rank", "plantsubsection")
                elif animals and current_rank == "section":
                    child.add_feature("ancestral_rank", "animalsection")
                elif animals and current_rank == "subsection":
                    child.add_feature("ancestral_rank", "animalsubsection")
                else:
                    child.add_feature("ancestral_rank", current_rank)
            else:
                if plants:
                    if child.tx_level == "section":
                        child.tx_level = "plantsection"
                    elif child.tx_level == "subsection":
                        child.tx_level = "plantsubsection"
                elif animals:
                    if child.tx_level == "section":
                        child.tx_level = "animalsection"
                    elif child.tx_level == "subsection":
                        child.tx_level = "animalsubsection"

                child.add_feature("ancestral_rank", child.tx_level)

            if child.name == "Metazoa_ott691846":
                # animals
                add_anc_ranks(child, current_rank, False, True, domain, kingdom, phylum, clas, order, family)
            elif child.name == "mrcaott2ott148" or parent.name == "Fungi_ott352914":
                # plants and fungi
                add_anc_ranks(child, current_rank, True, False, domain, kingdom, phylum, clas, order, family)
            else:
                add_anc_ranks(child, current_rank, plants, animals, domain, kingdom, phylum, clas, order, family)

    logger.info("Adding ancestral ranks.")
    add_anc_ranks(tre)


    def add_desc_ranks(parent, plants=False, animals=False):
        """Add descendant ranks. If this node has no rank, this is the rank of the first node below it which
        does have a rank. If this node has a rank, that rank is both its ancestral and descendant rank.
        """
        if parent.is_leaf():
            parent.add_feature("desc_rank",parent.tx_level)
            return parent.tx_level
        else:
            if parent.name == "Metazoa_ott691846":
                # animals
                animals = True
                plants = False
            elif parent.name == "mrcaott2ott148" or parent.name == "Fungi_ott352914":
                # plants and fungi
                animals = False
                plants = True

            max_desc_rank = "subspecies"
            for child in parent.children:
                child_rank = add_desc_ranks(child, plants, animals)

                if plants:
                    if child_rank == "section":
                        child_rank = "plantsection"
                    elif child_rank == "subsection":
                        child_rank = "plantsubsection"
                elif animals:
                    if child_rank == "section":
                        child_rank = "animalsection"
                    elif child_rank == "subsection":
                        child_rank = "animalsubsection"

                if tx_levels[child_rank] > tx_levels[max_desc_rank]:
                    max_desc_rank = child_rank

            if tx_levels[parent.tx_level] > 0:
                parent.add_feature("desc_rank",parent.tx_level)
                return parent.tx_level

            parent.add_feature("desc_rank",max_desc_rank)
            return max_desc_rank

    logger.info("Adding descendant ranks.")
    add_desc_ranks(tre)

    return tre
