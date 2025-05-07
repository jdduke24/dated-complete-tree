import os
import sys
import gc

import tree_loading
import tree_labelling
import tree_fixing
import tree_dating

import argparse
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename="main.log", filemode="w", force=True, level=logging.DEBUG)

#####################################################################################################################


def generate_trees(args):
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    sys.setrecursionlimit(10000)

    # Load metadata for tree from Open Tree, Chronosynth and OneZoom
    dates, phylogeny_nodes, taxa, descr_years = tree_loading.load_metadata(date_cache=args.date_cache,
                                                                           phylogeny=args.phylogeny,
                                                                           taxonomy=args.taxonomy,
                                                                           descr_dates=args.descr_dates)

    # Create ETE3 tree structure for entire Open Tree of Life, with my annotations
    whole_tre_unmodified = tree_loading.build_and_annotate_tree(dates,
                                                                phylogeny_nodes,
                                                                taxa,
                                                                descr_years,
                                                                tree_filename=args.supertree)

    if args.compute_ed:
        import tree_metrics
        ed_scores = {}

    if args.compute_rf:
        args.maintain_species_set = True

    if args.maintain_species_set:
        import random
        import datetime
        f_seeds = open("%s/seeds.txt" % args.output_folder, "wt")

    for n in range(args.num_trees):
        # Copy tree - we will change the copy, and keep the original unchanged so we can restore it next iteration without
        # reloading everything
        whole_tre = whole_tre_unmodified.copy()

        # Remove anything below species level (this includes promoting some subspecies to species if they are the only
        # example of their species)
        if args.maintain_species_set:
            random.seed(100)
        tree_fixing.remove_subspecies(whole_tre)
        tree_fixing.remove_nonspecies_leaves(whole_tre)

        #####################################################################################################################

        # Now, fixing the topology:
        if args.maintain_species_set:
            sd = datetime.datetime.now().timestamp()
            f_seeds.write("%d,%f\n" % (n+1,sd))
            random.seed(sd)
        # First, do labelling for steps 1-3:
        #  - 1-2 are independent of each other; step 3 collects up nodes not labelled in 1-2.
        #  - tree is only labelled at this stage; modifications are made in tree_fixing functions.
        genus_dict = {}       # step 1, nodes below genus nodes
        nmp_genus_dict = {}   # step 2, non-monophyletic genera
        tree_labelling.populate_genus_dict(whole_tre, genus_dict, nmp_genus_dict, None)

        tofix_dict = {}       # step 3, all other nodes from taxonomy (not phylogenies) to
                              # be moved to a suitable place in the tree, such that we
                              # generated a plausible hypothetical tree
        tree_labelling.populate_tofix_dict(whole_tre, tofix_dict, nmp_genus_dict)

        # Second, fix the topology based on the labels.
        # Fix steps 1 and 2.
        tree_fixing.fix_polyphyly(genus_dict)
        tree_fixing.fix_polyphyly(nmp_genus_dict)

        tree_fixing.remove_nonspecies_leaves(whole_tre)

        # Find and label backbone for step 3, after steps 1 an 2 already fixed.
        tree_labelling.populate_tofix_bkb(whole_tre, tofix_dict, [])
        fix_dict = tree_labelling.process_tofix_bkb(tofix_dict)

        # Finally, fix step 3.
        tree_fixing.fix_polyphyly(fix_dict, expand_parent_backbones=True)

        tree_fixing.remove_nonspecies_leaves(whole_tre)

        # Last of all, polytomy resolution.
        tree_fixing.fix_all_polytomies(whole_tre)

        # Optional - mainly here because it's good for OneZoom: remove one-child nodes. Gives a fully bifurcating tree.
        tree_fixing.delete_one_child_nodes(whole_tre)

        #####################################################################################################################

        # Date imputation
        tree_dating.remove_inconsistent_dates(whole_tre, whole_tre.date)
        tree_dating.date_labelling(whole_tre)
        tree_dating.impute_missing_dates(whole_tre)

        # All nodes now dated - set dists in ete and write out tree.
        tree_dating.write_tree_with_branch_lengths(whole_tre, filename="%s/%s_%d.tre" % (args.output_folder, args.output_tree_filename, n+1))

        if args.compute_ed:
            tree_metrics.add_ed_scores(whole_tre, ed_scores)

        del whole_tre
        gc.collect()


    if args.compute_ed:
        tree_metrics.write_ed_scores(ed_scores, filename="%s/ed_scores.csv" % (args.output_folder))

    if args.compute_rf:
        tree_metrics.compute_rf_distances(args.output_folder, args.output_tree_filename, args.num_trees)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a set of dated trees of all life, based on the Open Tree of Life and Chronosynth"
            "Optionally, also generate evolutionary distinctiveness scores for the trees"
        )
    )

    parser.add_argument("--num_trees",
                        help="How many trees to generate",
                        type=int,
                        default=1)

    parser.add_argument("--output_folder",
                        help="Path of folder where output trees will be written in Newick format",
                        default="output")

    parser.add_argument("--output_tree_filename",
                        help="Filename for output trees, e.g. 'dated_tree' would result in trees named 'dated_tree_1', 'dated_tree_2' etc.",
                        default="dated_tree")

    parser.add_argument("--supertree",
                        help="Path of the labelled_supertree_ottnames.tre file from the Open Tree of Life",
                        default="opentree14.9_tree/labelled_supertree/labelled_supertree_ottnames.tre")

    parser.add_argument("--date_cache",
                        help="Path of the date cache generated by Chronosynth",
                        default="chronosynth_date_info/node_ages.json")

    parser.add_argument("--phylogeny",
                        help="Path of the annotations.json file from the Open Tree of Life",
                        default="opentree14.9_tree/annotations.json")

    parser.add_argument("--taxonomy",
                        help="Path of the taxonomy.tsv file from the Open Tree Taxonomy",
                        default="ott3.6/taxonomy.tsv")

    parser.add_argument("--descr_dates",
                        help="Path of a csv file containing description dates for species",
                        default="oz_data/descr_dates.csv")

    parser.add_argument("--maintain_species_set",
                        help="Flag: whether to keep the set of species the same across all trees, enabling computation of distance metrics between trees. Default: False.",
                        action="store_true")

    parser.add_argument("--compute_ed",
                        help="Flag: whether to compute a distribution of ED scores. A csv file summarising the scores will be placed in the output folder. Default: False.",
                        action="store_true")

    parser.add_argument("--compute_rf",
                        help="Flag: whether to compute Robinson-Foulds distances between all pairs of trees. A csv file of distances be placed in the output folder. Default: False.",
                        action="store_true")

    args = parser.parse_args()
    generate_trees(args)


if __name__ == "__main__":
    main()
