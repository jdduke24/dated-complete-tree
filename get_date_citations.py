import os
os.chdir('/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/python/dated-complete-tree')

import json

dates = json.load(open('chronosynth_date_info/node_ages.json'))

my_dates = {}
my_sources = {}
for node in dates['node_ages']:
    my_dates[node] = []
    my_sources[node] = []
    for dat in dates['node_ages'][node]:
        my_dates[node].append(float(dat['age']))
        my_sources[node].append(dat['source_id'])

all_source_trees = set()
for node in my_sources:
    for source_tree in my_sources[node]:
        all_source_trees.add(source_tree)


seen_trees = set()
refs_text = []
for source_tree in all_source_trees:
    file_name = source_tree.split('@')[0]

    if file_name in seen_trees:
        continue
    else:
        seen_trees.add(file_name)

    phyle_folder = file_name[:3] + file_name[-2:]

    tree_name = source_tree.split('@')[1]

    print(source_tree, phyle_folder, file_name, tree_name)

    tree_dict = json.load(open("/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/phylesystem-1/study/%s/%s/%s.json" % (phyle_folder,
                                                                                                                                          file_name,
                                                                                                                                          file_name),
                                                                                                                                          'r'))

    ref_text = tree_dict["nexml"]["^ot:studyPublicationReference"].strip()
    ref_text = ref_text.replace('\n', '')

    refs_text.append(ref_text + '\n')


refs_text.sort()

fout = open('date_source_refs.txt', 'wt')
fout.writelines(refs_text)
fout.close()
