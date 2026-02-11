import json

file = 'ot_61/ot_3061/ot_3061'
tree = 'tree3'

tree_dict = json.load(open("/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/phylesystem-1/study/%s.json" % file,'r'))

for node in tree_dict['nexml']['treesById']['trees1']['treeById'][tree]['edgeBySourceId']:
    for edge in tree_dict['nexml']['treesById']['trees1']['treeById'][tree]['edgeBySourceId'][node]:
        print(edge, tree_dict['nexml']['treesById']['trees1']['treeById'][tree]['edgeBySourceId'][node][edge]['@length'])
        dist = float(tree_dict['nexml']['treesById']['trees1']['treeById'][tree]['edgeBySourceId'][node][edge]['@length'])
        tree_dict['nexml']['treesById']['trees1']['treeById'][tree]['edgeBySourceId'][node][edge]['@length'] = round(dist*100,6)

json.dump(tree_dict, open("/Users/Jonathan/Library/CloudStorage/Dropbox/Imperial/Tree_of_Life/Open_Tree/phylesystem-1/study/%s_new.json" % file,'w'),
          indent=0,
          sort_keys=True,
          separators=(',', ': '))
