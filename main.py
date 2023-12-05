import logging

import py4cytoscape
import requests.exceptions
from py4cytoscape import CyError
from tables_instruments import *

logging.basicConfig(level=logging.INFO, filename="py_log.log", filemode="w")


if __name__ == "__main__":
    # Запускаем цитоскейп:
    try:
        p4c.cytoscape_ping()
        p4c.cytoscape_version_info()
    except requests.exceptions.RequestException:
        logging.critical("In cytoscape_ping: Cannot find local or remote Cytoscape. Start Cytoscape and then proceed.")
        exit()
    except CyError:
        logging.critical("In cytoscape_ping: Error connecting to CyREST or version is unsupported.")
        exit()

    # Открываем тестовую сессию (находится в одной папке с этим кодом)
    path = os.path.join("test_network_AD.cys")
    p4c.open_session(path)
    # Загружаем табличку с ортогруппами (находится в одной папке с этим кодом)
    og_df = read_and_reformat_ortho_table("table_orthogroups.csv")
    # Получаем из цитоскейпа таблицы узлов и ребер для каждой сети
    network_list = p4c.get_network_list()
    network_number = p4c.get_network_count()
    node_tables = {}
    for n in range(network_number):
        node_tables[network_list[n]] = p4c.get_table_columns(table="node", network=network_list[n])

    edge_tables = {}
    for n in range(network_number):
        edge_tables[network_list[n]] = p4c.get_table_columns(table="edge", network=network_list[n],
                                                             columns=['name', 'shared name', 'node1_string_id', 'node2_string_id'])
    # костыль для тестовой сессии
    network_list.remove('string_interactions_short-18.tsv')
    all_og_nodes = []
    all_og_edges = []
    # Производим над каждой сетью манипуляции для подгрузки ортогрупп в node tables в цитоскейп
    # Описание функций лежит в tables_instruments.py
    for network in network_list:
        # TOD_: заменить маппинг имен на встроенную функцию map из py4cytoscape?

        # Добавляем ID генов другого формата в таблицу узлов из таблицы ребер
        names = get_names_relation(edge_tables[network])
        set_names_in_node_table(node_tables[network], names)
        # Присваиваем генам соответствующие ортогруппы
        # TODO: проверить правильность порядка ортогрупп после set
        specie_og_list = list(set(assign_og_to_gene(node_tables[network], og_df)))
        # Обновляем в видовой сети таблицу узлов (добавился столбец name)
        # TODO: посмотреть, можно ли переписать на загрузку только одного столбца, а не всей таблицы
        load_altered_node_table_to_cyto(node_tables[network], network)
        added_og_nodes = p4c.add_cy_nodes(specie_og_list, network=network)
        all_og_nodes += added_og_nodes
        # Подгружаем столбец "is_OG", значения в котором равны True если это узел ортогруппы, False - если это ген
        load_is_og_col(node_tables[network], specie_og_list, network)
        og_edges_list = node_tables[network][['name', 'OG']].dropna().values.tolist()
        # Дезайн нодов
        p4c.set_node_shape_bypass(specie_og_list, 'ELLIPSE', network=network)
        p4c.set_node_color_bypass(specie_og_list, '#FFA753', network=network)
        # Добавляем связи между генами и ортогруппами
        added_og_edges = p4c.add_cy_edges(og_edges_list, edge_type='is in orthogroup', network=network)
        # Дезайн ребер
        added_og_edges = [dict(i)['SUID'] for i in added_og_edges]
        all_og_edges += added_og_edges
        p4c.set_edge_color_bypass(added_og_edges, '#FFA753', network=network)
        # p4c.select_nodes(specie_og_list, preserve_current_selection=False, by_col='shared name', network=network)
        # Применить раскладку
        # py4cytoscape.layout_network(layout_name='attribute-circle', network=network)
    # Merge tables
    p4c.merge_networks(sources=network_list, title='all_species_merged')
    # Дезайн общей сети
    p4c.set_current_network(network='all_species_merged')
    p4c.set_node_shape_bypass(all_og_nodes, 'ELLIPSE', network='all_species_merged')
    p4c.set_node_color_bypass(all_og_nodes, '#FFA753', network='all_species_merged')
    p4c.set_edge_color_bypass(all_og_edges, '#FFA753', network='all_species_merged')
    p4c.layout_network()
