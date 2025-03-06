from taxonomy_utils import tx_levels

import logging
logger = logging.getLogger(__name__)

def populate_genus_dict(parent, genus_dict, nmp_genus_dict, genus_root):
    """Recurse through the tree, populating:
    genus_dict: step 1, identify taxonomy-only species placed directly under genus nodes, along with their associated backbone.
    nmp_genus_dict: step 2, identify non-monophyletic genera, ie taxonomy-only species which could be moved into more than one monophyletic group.
    """

    # these variables keep track of the names of genera found, and the number of species or genus nodes found that are not beneath a genus node.
    # they contain info only for the current parent and (at the end of the for loop) everything beneath it in the tree.

    species_found = set()
    genera_found = set()

    for child in parent.children:
        # keep track of all genera and species found
        if child.tx_level == "species":
            if child.genus_name not in genera_found:
                # not seen this genus before
                genera_found.add(child.genus_name)

            if child.species_name not in species_found:
                # not seen this species before
                species_found.add(child.species_name)

        # set up lists for collecting info about this genus/species
        if child.tx_level == "genus" and child.ph_tx == "PH":
            # dict for step 1
            genus_dict[child] = [[],[]]
        elif child.tx_level == "species" and child.genus_name not in nmp_genus_dict:
            # dict for step 2
            nmp_genus_dict[child.genus_name] = [[],[]]


        if child.ph_tx == "TX":
            # nodes from a taxonomy may need fixing
            if child.tx_level == "species":
                if parent.tx_level == "genus" and len(parent.children) > 2:
                    # we found a species from a taxonomy. If it is directly attached to its genus node and
                    # the genus node has more than 2 children, it needs fixing
                    genus_dict[genus_root][0].append(child)
                    child.info = "GR FIX %s" % genus_root.name
                elif genus_root:
                    # if we are either deeper in the genus, or directly under a genus node with fewer than 3 children,
                    # we don't need to do anything with this node
                    child.info = "GR NONE %s" % genus_root.name # mark so they aren't picked up in step 4
                else:
                    # otherwise, we have a species from a taxonomy that is not below any genus node
                    #  -- need to collect up to be moved into a non-monophyletic group,
                    #     but if it is attached directly to some subgenus or species group, leave it where it is
                    if tx_levels[parent.tx_level] <= tx_levels["species"] or tx_levels[parent.tx_level] >= tx_levels["genus"]:
                        nmp_genus_dict[child.genus_name][0].append(child)
                        child.info = "NMP FIX %s" % child.genus_name
            elif child.tx_level != "genus":
                # don't recurse from species or genus nodes from taxonomy - but anything of higher rank will be recursed through
                new_genera_found, new_species_found = populate_genus_dict(child, genus_dict, nmp_genus_dict, genus_root)
                genera_found.update(new_genera_found)
                species_found.update(new_species_found)

        else: # child.ph_tx == "PH":
            if child.tx_level == "species":
                if genus_root:
                    # if we are below a genus node, add to step 1 backbone
                    genus_dict[genus_root][1].append(child)
                    child.info = "GR BKB %s" % genus_root.name
                else:
                    # or if not, add to step 2 backbone
                    nmp_genus_dict[child.genus_name][1].append(child)
                    child.info = "NMP BKB %s" % child.genus_name
            elif child.tx_level == "genus":
                # genus from a phylogeny - recurse down from here
                new_genera_found, new_species_found = populate_genus_dict(child, genus_dict, nmp_genus_dict, child)
                genera_found.update(new_genera_found)
                species_found.update(new_species_found)
            else:
                # not a species or genus.
                # if genus_root exists then we are somewhere under a genus node - add this node to step 1 backbone
                if genus_root:
                    genus_dict[genus_root][1].append(child)
                    child.info = "GR BKB %s" % genus_root.name
                new_genera_found, new_species_found = populate_genus_dict(child, genus_dict, nmp_genus_dict, genus_root)
                genera_found.update(new_genera_found)
                species_found.update(new_species_found)


    if len(genera_found) == 1:
        # if we found species from only one genus below this node and some nodes that are in the step 2 backbone, add
        # this node also to the backbone
        genus = genera_found.copy()
        genus_name = genus.pop()
        if genus_name in nmp_genus_dict and len(nmp_genus_dict[genus_name][1]) > 0:
            nmp_genus_dict[genus_name][1].append(parent)
            if parent.info is None:
                parent.info = "NMP PAR BKB %s" % genus_name
            else:
                parent.info += "NMP PAR BKB %s" % genus_name

    return (genera_found, species_found)



####################################################################################
# find remaining nodes to fix
def populate_tofix_dict(parent, tofix_dict, nmp_genus_dict):
    # can't fix original root node, so fine for recursion to go straight into children
    for child in parent.children:
        if child.info is not None and child.info[:2] == "GR":
            # this is already in a genus node backbone - ignore and don't recurse further
            continue

        if child.ph_tx == "TX":
            # could be a node we need to fix
            if child.tx_level == "species" and parent.info is not None and parent.info[:3] == "NMP":
                # this is a species below a NMP node - ignore and don't recurse further
                continue
            if child.tx_level == "species" and child.info is not None and child.info[:3] == "NMP":
                if child.genus_name in nmp_genus_dict:
                    if len(nmp_genus_dict[child.genus_name][1]) > 0:
                        # if we already found this species in a non-monophyletic genus and we also
                        # found a backbone for that genus, we will fix it there rather than here
                        continue
                    else:
                        # but if there isn't a backbone for this genus, it gets moved into the "other" fixing category
                        del nmp_genus_dict[child.genus_name]

            # otherwise, note parent node below which we will construct a backbone
            if parent not in tofix_dict:
                tofix_dict[parent] = {}
                parent.info = "OTH PARENT"

            if child.tx_level == "no rank":
                used_rank = child.desc_rank
            else:
                used_rank = child.tx_level

            if used_rank not in tofix_dict[parent]:
                tofix_dict[parent][used_rank] = [[],[],None]

            tofix_dict[parent][used_rank][0].append(child)
            child.info = "OTH FIX"

            continue

        populate_tofix_dict(child, tofix_dict, nmp_genus_dict)

# and here, populate the possible backbone for each node to fix, ie anywhere under the same parent that
# is not already labelled, since labelled nodes belong to some monophyletic group already
def populate_tofix_bkb(parent, tofix_dict, current_roots, root_call=True):
    if root_call:
        # first call only; repeat of the contents of the loop below
        if parent.info == "OTH PARENT" and (tx_levels[parent.tx_level] < 0  or tx_levels[parent.tx_level] >= tx_levels["subtribe"]):
            # never recurse through a genus or below
            current_roots.append(parent)
            populate_tofix_bkb(parent, tofix_dict, current_roots, False)
            for rank in tofix_dict[parent]:
                tofix_dict[parent][rank][2] = list(current_roots)
            current_roots.pop()

        if (parent.info is None or parent.info == "OTH PARENT") and parent.ph_tx == "PH":
            for root in current_roots:
                for rank in tofix_dict[root]:
                    if (tx_levels[parent.tx_level] < 0 and tx_levels[parent.ancestral_rank] > tx_levels[rank]) or tx_levels[parent.tx_level] >= tx_levels[rank]:
                        tofix_dict[root][rank][1].append(parent)
                        if parent.info is None:
                            parent.info = "OTH BKB %s %s" % (root.name, rank)
                        else:
                            parent.info += "OTH BKB %s %s" % (root.name, rank)
            if tx_levels[parent.tx_level] < 0  or tx_levels[parent.tx_level] >= tx_levels["subtribe"]:
                # never recurse through a genus or below
                populate_tofix_bkb(parent, tofix_dict, current_roots, False)


    for child in parent.children:
        if child.info == "OTH PARENT" and (tx_levels[child.tx_level] < 0  or tx_levels[child.tx_level] >= tx_levels["subtribe"]):
            # never recurse through a genus or below
            current_roots.append(child)
            populate_tofix_bkb(child, tofix_dict, current_roots, False)
            for rank in tofix_dict[child]:
                tofix_dict[child][rank][2] = list(current_roots)
            current_roots.pop()

        if (child.info is None or child.info == "OTH PARENT") and child.ph_tx == "PH":
            for root in current_roots:
                for rank in tofix_dict[root]:
                    if ((tx_levels[child.tx_level] < 0 and tx_levels[child.ancestral_rank] > tx_levels[rank])
                            or tx_levels[child.tx_level] >= tx_levels[rank]
                            or (tx_levels[child.tx_level] > 0 and tx_levels[child.tx_level] < tx_levels[rank] and tx_levels[parent.ancestral_rank] > tx_levels[child.tx_level] and tx_levels[parent.ancestral_rank] > tx_levels[rank])):
                        tofix_dict[root][rank][1].append(child)
                        if child.info is None:
                            child.info = "OTH BKB %s %s" % (root.name, rank)
                        else:
                            child.info += "OTH BKB %s %s" % (root.name, rank)
            if tx_levels[child.tx_level] < 0  or tx_levels[child.tx_level] >= tx_levels["subtribe"]:
                # never recurse through a genus or below
                populate_tofix_bkb(child, tofix_dict, current_roots, False)

        elif ((child.info is None or child.info[:3] == "NMP")
                  and (parent.info is not None and parent.info[:3] == "OTH")
                  and (child.ph_tx == "PH" or child.ph_tx == "IN")):
            # root node of a NMP group can be in backbone, because we insert above it
            for root in current_roots:
                for rank in tofix_dict[root]:
                    if rank == "no rank":
                        # for the purpose of backbone-building only, treat an OTH FIX node of "no rank" like "genus",
                        # i.e. an OTH FIX node of no rank cannot be put below a genus
                        used_rank = "genus"
                    else:
                        used_rank = rank
                    if ((tx_levels[child.tx_level] < 0 and tx_levels[child.ancestral_rank] > tx_levels[used_rank])
                            or tx_levels[child.tx_level] >= tx_levels[used_rank]
                            or (tx_levels[child.tx_level] > 0 and tx_levels[child.tx_level] < tx_levels[rank] and tx_levels[parent.ancestral_rank] > tx_levels[child.tx_level] and tx_levels[parent.ancestral_rank] > tx_levels[rank])):
                        tofix_dict[root][rank][1].append(child)
                        if child.info is None:
                            child.info = "OTH EX BKB %s %s" % (root.name, rank)
                        else:
                            child.info += "OTH EX BKB %s %s" % (root.name, rank)


def process_tofix_bkb(tofix_dict):
    new_dict = {}
    for root in tofix_dict:
        for rank in tofix_dict[root]:
            new_dict[(root,rank)] = tofix_dict[root][rank]

    return new_dict
