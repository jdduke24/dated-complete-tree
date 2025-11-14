from taxonomy_utils import tx_levels

import logging
logger = logging.getLogger(__name__)


def add_anc_ranks(parent, current_rank=None, plants=False, animals=False, domain=None, kingdom=None, phylum=None, clas=None, order=None, family=None):
    """Add ancestral ranks. If this node has no rank, this is the rank of the first node above it which
    does have a rank. If this node has a rank, that rank is both its ancestral and descendant rank.
    Also annotate every node with the domain, kingdom, phylum, class, order and family it is under (if any).
    """
    if current_rank is None:
        # only occurs at root (first call of function)
        current_rank = parent.props["tx_level"]
        parent.add_prop("ancestral_rank", parent.props["tx_level"])

    if tx_levels[parent.props["tx_level"]] > 0 or tx_levels[current_rank] < 0:
        current_rank = parent.props["tx_level"]

    if parent.props["tx_level"] == "domain":
        domain = parent.name
    elif parent.props["tx_level"] == "kingdom":
        kingdom = parent.name
    elif parent.props["tx_level"] == "phylum":
        phylum = parent.name
    elif parent.props["tx_level"] == "class":
        clas = parent.name
    elif parent.props["tx_level"] == "order":
        order = parent.name
    elif parent.props["tx_level"] == "family":
        family = parent.name

    parent.add_prop("domain", domain)
    parent.add_prop("kingdom", kingdom)
    parent.add_prop("phylum", phylum)
    parent.add_prop("clas", clas)
    parent.add_prop("order", order)
    parent.add_prop("family", family)

    for child in parent.children:
        if tx_levels[child.props["tx_level"]] < 0:
            if plants and current_rank == "section":
                child.add_prop("ancestral_rank", "plantsection")
            elif plants and current_rank == "subsection":
                child.add_prop("ancestral_rank", "plantsubsection")
            elif animals and current_rank == "section":
                child.add_prop("ancestral_rank", "animalsection")
            elif animals and current_rank == "subsection":
                child.add_prop("ancestral_rank", "animalsubsection")
            else:
                child.add_prop("ancestral_rank", current_rank)
        else:
            if plants:
                if child.props["tx_level"] == "section":
                    child.props["tx_level"] = "plantsection"
                elif child.props["tx_level"] == "subsection":
                    child.props["tx_level"] = "plantsubsection"
            elif animals:
                if child.props["tx_level"] == "section":
                    child.props["tx_level"] = "animalsection"
                elif child.props["tx_level"] == "subsection":
                    child.props["tx_level"] = "animalsubsection"

            child.add_prop("ancestral_rank", child.props["tx_level"])

        if child.name == "Metazoa_ott691846":
            # animals
            add_anc_ranks(child, current_rank, False, True, domain, kingdom, phylum, clas, order, family)
        elif child.name == "mrcaott2ott148" or parent.name == "Fungi_ott352914":
            # plants and fungi
            add_anc_ranks(child, current_rank, True, False, domain, kingdom, phylum, clas, order, family)
        else:
            add_anc_ranks(child, current_rank, plants, animals, domain, kingdom, phylum, clas, order, family)


def add_desc_ranks(parent, plants=False, animals=False):
    """Add descendant ranks. If this node has no rank, this is the rank of the first node below it which
    does have a rank. If this node has a rank, that rank is both its ancestral and descendant rank.
    """
    if parent.is_leaf:
        parent.add_prop("desc_rank",parent.props["tx_level"])
        return parent.props["tx_level"]
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

        if tx_levels[parent.props["tx_level"]] > 0:
            parent.add_prop("desc_rank",parent.props["tx_level"])
            return parent.props["tx_level"]

        parent.add_prop("desc_rank",max_desc_rank)
        return max_desc_rank


def populate_genus_dict(parent, genus_dict, nmp_genus_dict, genus_root, kingdom="Other", root=True):
    """Recurse through the tree, populating:
    genus_dict: step 1, identify taxonomic species placed directly under phylogenetic genus nodes, along with their associated backbone.
    nmp_genus_dict: step 2, identify non-monophyletic genera, ie taxonomic species which could be moved into more than one monophyletic group.
    """

    # these variables keep track of the names of genera found, and the number of species or genus nodes found that are not beneath a genus node.
    # they contain info only for the current parent and (at the end of the for loop) everything beneath it in the tree.

    if parent.props["tx_level"] == "kingdom":
        kingdom = parent.name

    if root and parent.props["ph_tx"] == "PH" and tx_levels[parent.props["tx_level"]] == tx_levels["genus"]:
        # if the root node of the whole tree is itself a genus, set this as the current genus root
        genus_root = parent

        genus_dict[genus_root] = [[],[genus_root]]

    species_found = set()
    genera_found = set()

    for child in parent.children:
        # keep track of all genera and species found
        if tx_levels[child.props["tx_level"]] == tx_levels["species"]:
            if child.props["genus_name"] not in genera_found:
                # not seen this genus before
                genera_found.add("%s/%s" % (kingdom, child.props["genus_name"]))

            if child.props["species_name"] not in species_found:
                # not seen this species before
                species_found.add("%s/%s" % (kingdom, child.props["species_name"]))

        # set up lists for collecting info about this genus/species
        if tx_levels[child.props["tx_level"]] == tx_levels["genus"] and child.props["ph_tx"] == "PH":
            # dict for step 1
            genus_dict[child] = [[],[child]]

        elif tx_levels[child.props["tx_level"]] == tx_levels["species"] and "%s/%s" % (kingdom, child.props["genus_name"]) not in nmp_genus_dict:
            # dict for step 2
            nmp_genus_dict["%s/%s" % (kingdom, child.props["genus_name"])] = [[],[]]


        if child.props["ph_tx"] == "TX":
            # nodes from a taxonomy may need fixing
            if tx_levels[child.props["tx_level"]] == tx_levels["species"]:
                if tx_levels[parent.props["tx_level"]] == tx_levels["genus"] and len(parent.children) > 2:
                    # we found a species from a taxonomy. If it is directly attached to its genus node and
                    # the genus node has more than 2 children, it needs fixing
                    genus_dict[genus_root][0].append(child)
                    child.props["info"] = "GR FIX %s" % genus_root.name
                elif genus_root:
                    # if we are either deeper in the genus, or directly under a genus node with fewer than 3 children,
                    # we don't need to do anything with this node
                    child.props["info"] = "GR NONE %s" % genus_root.name # mark so they aren't picked up in step 4
                else:
                    # otherwise, we have a species from a taxonomy that is not below any genus node
                    #  -- need to collect up to be moved into a non-monophyletic group,
                    #     but if it is attached directly to some subgenus or species group, leave it where it is
                    if tx_levels[parent.props["tx_level"]] <= tx_levels["species"] or tx_levels[parent.props["tx_level"]] >= tx_levels["genus"]:
                        nmp_genus_dict["%s/%s" % (kingdom, child.props["genus_name"])][0].append(child)
                        child.props["info"] = "NMP FIX %s/%s" % (kingdom, child.props["genus_name"])
            elif tx_levels[child.props["tx_level"]] != tx_levels["genus"]:
                # don't recurse from species or genus nodes from taxonomy - but anything of higher rank will be recursed through
                new_genera_found, new_species_found = populate_genus_dict(child, genus_dict, nmp_genus_dict, genus_root, kingdom, False)
                genera_found.update(new_genera_found)
                species_found.update(new_species_found)

        else: # child.props["ph_tx"] == "PH":
            if tx_levels[child.props["tx_level"]] == tx_levels["species"]:
                if genus_root:
                    # if we are below a genus node, add to step 1 backbone
                    genus_dict[genus_root][1].append(child)
                    child.props["info"] = "GR BKB %s" % genus_root.name
                else:
                    # or if not, add to step 2 backbone
                    nmp_genus_dict["%s/%s" % (kingdom, child.props["genus_name"])][1].append(child)
                    child.props["info"] = "NMP BKB %s/%s" % (kingdom, child.props["genus_name"])
            elif tx_levels[child.props["tx_level"]] == tx_levels["genus"]:
                # genus from a phylogeny - recurse down from here
                new_genera_found, new_species_found = populate_genus_dict(child, genus_dict, nmp_genus_dict, child, kingdom, False)
                genera_found.update(new_genera_found)
                species_found.update(new_species_found)
            else:
                # not a species or genus.
                # if genus_root exists then we are somewhere under a genus node - add this node to step 1 backbone
                if genus_root:
                    genus_dict[genus_root][1].append(child)
                    child.props["info"] = "GR BKB %s" % genus_root.name
                new_genera_found, new_species_found = populate_genus_dict(child, genus_dict, nmp_genus_dict, genus_root, kingdom, False)
                genera_found.update(new_genera_found)
                species_found.update(new_species_found)


    if len(genera_found) == 1:
        # if we found species from only one genus below this node and some nodes that are in the step 2 backbone, add
        # this node also to the backbone
        genus = genera_found.copy()
        genus_name = genus.pop()
        if genus_name in nmp_genus_dict and len(nmp_genus_dict[genus_name][1]) > 0:
            nmp_genus_dict[genus_name][1].append(parent)
            if parent.props["info"] is None:
                parent.props["info"] = "NMP EX BKB %s" % genus_name
            else:
                parent.props["info"] += "NMP EX BKB %s" % genus_name

    return (genera_found, species_found)



####################################################################################
# find remaining nodes to fix
def populate_tofix_dict(parent, tofix_dict, nmp_genus_dict, kingdom="Other"):
    """Recurse through tree labelling taxonomic nodes that need moving into a phylogenetic backbone. Ignores
    nodes under genus nodes (these are taken care of previously with GR labels), and species already collected
    up as part of a non-monophyletic genus.
    """
    if parent.props["tx_level"] == "kingdom":
        kingdom = parent.name

    # can't fix original root node, so it's ok to go straight into children
    for child in parent.children:
        if child.props["info"] is not None and child.props["info"][:2] == "GR":
            # this is already in a genus node backbone - ignore and don't recurse further
            continue

        if child.props["ph_tx"] == "TX":
            # could be a node we need to fix
            if tx_levels[child.props["tx_level"]] == tx_levels["species"] and parent.props["info"] is not None and parent.props["info"][:3] == "NMP":
                # this is a species below a NMP node - ignore and don't recurse further
                continue
            if tx_levels[child.props["tx_level"]] == tx_levels["species"] and child.props["info"] is not None and child.props["info"][:3] == "NMP":
                genus_key = "%s/%s" % (kingdom, child.props["genus_name"])
                if genus_key in nmp_genus_dict:
                    if len(nmp_genus_dict[genus_key][1]) > 0:
                        # if we already found this species in a non-monophyletic genus and we also
                        # found a backbone for that genus, we will fix it there rather than here
                        continue
                    else:
                        # but if there isn't a backbone for this genus, it gets moved into the "other" fixing category
                        del nmp_genus_dict[genus_key]

            # otherwise, note parent node below which we will construct a backbone
            if parent not in tofix_dict:
                tofix_dict[parent] = {}
                parent.props["info"] = "OTH PARENT"

            if child.props["tx_level"] == "no rank":
                used_rank = child.props["desc_rank"]
            else:
                used_rank = child.props["tx_level"]

            if used_rank not in tofix_dict[parent]:
                tofix_dict[parent][used_rank] = [[], [parent], None]


            tofix_dict[parent][used_rank][0].append(child)
            child.props["info"] = "OTH FIX"

            continue

        populate_tofix_dict(child, tofix_dict, nmp_genus_dict, kingdom)


def oth_backbone_type(child, parent, backbone_rank):
    """Determine if the node "child" should be part of a backbone for a taxonomic node.
    A backbone node could be:
        1. an MRCA node whose ancestral rank is higher than we are constructing a backbone for.
        2. A non-MRCA node of equal or higher rank than we are constructing a backbone for.
        3. A non-MRCA node whose rank is less than the backbone rank, but with a parent whose ancestral rank
          is higher than both this node and the backbone rank.
    """
    if (tx_levels[child.props["tx_level"]] < 0 and tx_levels[child.props["ancestral_rank"]] > tx_levels[backbone_rank]):
        return 1
    elif tx_levels[child.props["tx_level"]] >= tx_levels[backbone_rank]:
        return 2
    elif parent and (tx_levels[child.props["tx_level"]] > 0 and
                     tx_levels[child.props["tx_level"]] < tx_levels[backbone_rank] and
                     tx_levels[parent.props["ancestral_rank"]] > tx_levels[child.props["tx_level"]] and
                     tx_levels[parent.props["ancestral_rank"]] > tx_levels[backbone_rank]):
        return 3
    else:
        return None


def populate_tofix_bkb(parent, tofix_dict, current_roots, root_call=True):
    """Populate the possible backbone for each OTH FIX node. The backbone consists of anywhere under the same parent
    that is not already labelled and GR or NMP, and is of an appropriate rank.

    The keys of the tofix_dict are nodes labelled as OTH PARENT.
    Each entry in tofix_dict then has a list of ranks, each of which will contain the backbone for an OTH FIX
    node of that rank.

    The list current_roots maintains a list of OTH PARENT keys for the tofix_dict, into which the current node may
    be inserted.
    """

    if root_call:
        # first call only
        if parent.props["info"] and parent.props["info"][:10] == "OTH PARENT" and (tx_levels[parent.props["tx_level"]] < 0  or tx_levels[parent.props["tx_level"]] > tx_levels["genus"]):
            # current_roots maintains a list of this node and the OTH PARENT nodes above this one. These nodes
            # will have their backbones populated as we go along. If this node is suitable for a backbone of the
            # current parent, it will be suitable for parents above this one too.
            # Assume nothing is going to go above the root.
            current_roots.append(parent)
            for rank in tofix_dict[parent]:
                tofix_dict[parent][rank][2] = list(current_roots)


    for child in parent.children:
        if ((child.props["info"] is None or child.props["info"] == "OTH PARENT") and (child.props["ph_tx"] == "PH" or child.props["ph_tx"] == "IN" or child.props["ph_tx"] == "FI")):
            for root in current_roots:
                for rank in tofix_dict[root]:
                    bkb_type = oth_backbone_type(child, parent, rank)
                    if bkb_type:
                        tofix_dict[root][rank][1].append(child)
                        if child.props["info"] is None:
                            child.props["info"] = "OTH BKB %d %s %s" % (bkb_type, root.name, rank)
                        else:
                            child.props["info"] += "\nOTH BKB %d %s %s" % (bkb_type, root.name, rank)

            if tx_levels[child.props["tx_level"]] < 0  or tx_levels[child.props["tx_level"]] > tx_levels["genus"]:
                # never recurse through a genus or below, or through an inserted node (which at this stage must
                # have come from NMP fixing, so going any further would be interrupting a monophyletic group - unless
                # it was added by a forced taxa move ie ph_tx == "FI").
                if child.props["info"] and child.props["info"][:10] == "OTH PARENT":
                    current_roots.append(child)
                    for rank in tofix_dict[child]:
                        tofix_dict[child][rank][2] = list(current_roots)

                if child.props["ph_tx"] == "PH" or child.props["ph_tx"] == "FI":
                    populate_tofix_bkb(child, tofix_dict, current_roots, False)

                if child.props["info"] and child.props["info"][:10] == "OTH PARENT":
                    current_roots.pop()

        elif ((child.props["info"] is None or child.props["info"][:3] == "NMP") and
              (parent.props["info"] is not None and parent.props["info"][:3] == "OTH") and
              (child.props["ph_tx"] == "PH" or child.props["ph_tx"] == "IN")):
            # root node of a NMP group can be in backbone, because we insert above it
            for root in current_roots:
                for rank in tofix_dict[root]:
                    if rank == "no rank":
                        # for the purpose of backbone-building only, treat an OTH FIX node of "no rank" like "genus",
                        # i.e. an OTH FIX node of no rank cannot be put below a genus
                        used_rank = "genus"
                    else:
                        used_rank = rank
                    bkb_type = oth_backbone_type(child, parent, used_rank)
                    if bkb_type:
                        tofix_dict[root][rank][1].append(child)
                        if child.props["info"] is None:
                            child.props["info"] = "OTH EX BKB %d %s %s" % (bkb_type, root.name, rank)
                        else:
                            child.props["info"] += "\nOTH EX BKB %d %s %s" % (bkb_type, root.name, rank)


def process_tofix_bkb(tofix_dict):
    """Process the output of populate_tofix_bkb into the format needed for tree_fixing.fix_polyphyly."""
    new_dict = {}
    for root in tofix_dict:
        for rank in tofix_dict[root]:
            new_dict[(root,rank)] = tofix_dict[root][rank]

    return new_dict
