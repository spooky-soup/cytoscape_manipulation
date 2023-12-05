from Bio import Phylo
import networkx
import py4cytoscape as p4c

if __name__ == "__main__":
    tree = Phylo.read('./OG_trees/OG0000955.nwk', 'newick')
    print(tree)
    net = Phylo.to_networkx(tree)
    print(net)
    cyjs = networkx.cytoscape_data(net, name='test tree')
    print(cyjs)
    cnt_id = 0
    for node in cyjs['elements']['nodes']:
        node['data']['value'] = str(node['data']['value'])
        node['data']['id'] = node['data']['id'] + str(cnt_id)
        cnt_id += 1
    for edge in cyjs['elements']['edges']:
        edge['data']['source'] = str(edge['data']['source'])
        edge['data']['target'] = str(edge['data']['target'])
    p4c.create_network_from_cytoscapejs(cyjs, title='tree')
