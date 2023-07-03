import os

import pandas as pd
import py4cytoscape as p4c
import numpy as np

CYS_PATH = "/home/spooky/РНФ_сети/test_network_merged_byname.cys"
if __name__ == '__main__':
    path = os.path.join(CYS_PATH)
    p4c.open_session(path)
    network = "Merged Network"
    p4c.cytoscape_ping()
    edge_table = p4c.get_table_columns(table="edge", network=network)
    node_table = p4c.get_table_columns(table="node", network=network)
    og_list = node_table['OG'].unique().tolist()
    og_list = list(filter(lambda n: 'OG' in n, og_list))
    og_relation = pd.DataFrame(index=og_list)
    for og in og_list:
        print(og)
        neighbours = p4c.get_first_neighbors(node_names=(node_table[node_table['OG'] == og]).SUID.tolist(),
                                             as_nested_list=True,
                                             network=network)
        neighbours = list(filter(lambda n: 'OG' not in n, np.hstack([p[1] for p in neighbours])))
        neighbours_og = []
        for n in neighbours:
            neighbours_og.append((node_table[node_table['name'] == n])['OG'].values[0])

        unique_og = np.unique(neighbours_og)
        col = {}
        for n_og in unique_og:
            col[n_og] = neighbours_og.count(n_og)
        og_relation[og] = pd.Series(col)
    print(og_relation)
    og_relation.to_csv(path_or_buf='./orthogroups_connectivity.csv')
