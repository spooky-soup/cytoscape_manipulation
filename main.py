import logging
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
        edge_tables[network_list[n]] = p4c.get_table_columns(table="edge", network=network_list[n], columns=['name',
                                                                                                             'node1_string_id',
                                                                                                             'node2_string_id'])

    # Производим над каждой сетью манипуляции для подгрузки ортогрупп в node tables в цитоскейп
    # Описание функций лежит в tables_instruments.py
    for network in network_list:
        if network == 'string_interactions_short-18.tsv':
            continue
        names = get_names_relation(edge_tables[network])
        set_names_in_node_table(node_tables[network], names)
        assign_og_to_gene(node_tables[network], og_df)
        load_altered_node_table_to_cyto(node_tables[network], network)
