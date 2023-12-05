import numpy as np
import pandas as pd
import os
import py4cytoscape as p4c
# Функции для работы с разныии таблицами


# "Плоский" список из вложенных списков
def flatten(xs):
    result = []
    if isinstance(xs, (list, tuple)):
        for x in xs:
            result.extend(flatten(x))
    elif pd.isna(xs):
        pass
    else:
        result.append(xs)
    return result


# Изменение формата таблицы с ортогруппами
# Строчки - ортогруппы, столбцы - виды (по номерам)
def read_and_reformat_ortho_table(path: str) -> pd.DataFrame:
    # Файл с ортогруппами
    og_path = os.path.join(path)
    og_df = pd.read_csv(og_path)
    og_df = pd.DataFrame(np.vstack([og_df.columns, og_df]))
    og_df.set_index(og_df.iloc[:, 0], inplace=True, drop=True)
    og_df.drop(columns=0, inplace=True)
    og_dict = {}
    for r in og_df.iterrows():
        og_dict.update({r[0]: sorted(flatten([list(s.str.split(', ')) for s in r[1:]]))})

    res_df = pd.DataFrame(index=list(og_dict.keys()))
    for o in og_dict.keys():
        for p in og_dict[o]:
            n = int(p.split('.')[0])
            if n not in res_df.columns:
                res_df[n] = list([] for _ in range(len(og_dict)))
            res_df[n][o].append(p)

    return res_df


# Из таблиц цитоскейпа получаем соотношение тривиального названия и ID
def get_names_relation(edge_table: pd.DataFrame) -> dict:
    cnt = 0
    edge_table['name1'] = ['']*len(edge_table.index)
    edge_table['name2'] = ['']*len(edge_table.index)
    for s in edge_table.loc[:, 'name']:
        gene1, gene2 = s.split(' (interacts with) ')
        edge_table.loc[edge_table.index[cnt], 'name1'] = gene1
        edge_table.loc[edge_table.index[cnt], 'name2'] = gene2
        cnt += 1
    names = dict(zip(edge_table.loc[:, 'name1'], edge_table.loc[:, 'node1_string_id']))
    names.update(dict(zip(edge_table.loc[:, 'name2'], edge_table.loc[:, 'node2_string_id'])))
    return names


# Добавляем в таблицу узлов ID
def set_names_in_node_table(node_table: pd.DataFrame, relation: dict):
    node_table['string_id'] = ['']*len(node_table.index)
    for i in node_table.index:
        name = node_table.loc[i, 'name']
        node_table.loc[i, 'string_id'] = relation[name]


# Подгружаем в node_table столбец "is_OG" и ноды ортогрупп
def load_is_og_col(node_table: pd.DataFrame, og_list: list, network: str):
    is_og_dict = {'name': list(node_table['name']) + og_list, 'is_og': [False]*len(node_table.index) + [True]*len(og_list)}
    og_df = pd.DataFrame(data=is_og_dict, index=list(node_table.index) + og_list)
    p4c.load_table_data(og_df, table='node', data_key_column='name', table_key_column='name', network=network)


# Добавляем в таблицу узлов ортогруппы
# Возвращает список ортогрупп для текущего вида
def assign_og_to_gene(node_table: pd.DataFrame, og_table: pd.DataFrame):
    node_table['OG'] = ['']*len(node_table.index)
    node_table['is_OG'] = [False]*len(node_table.index)
    specie = int(node_table.loc[node_table.index[0], 'string_id'].split('.')[0])
    node_table['specie'] = [specie]*len(node_table.index)
    node_table['name'] = str(specie) + ". " + node_table['shared name']
    node_table.set_index('string_id', inplace=True)
    og_list = []
    for i in og_table.index:
        for gene in og_table.loc[i, specie]:
            if gene not in node_table.index:
                continue
            node_table.loc[gene, 'OG'] = i
            og_list.append(i)
    node_table['OG'].replace('', np.nan, inplace=True)
    node_table.reset_index(inplace=True)
    node_table.set_index('SUID', inplace=True)
    return og_list


# Загружаем видоизмененную таблицу узлов обратно в цитоскейп
def load_altered_node_table_to_cyto(node_table: pd.DataFrame, network_name: str):
    p4c.load_table_data(node_table, table='node', network=network_name,
                        data_key_column='shared name', table_key_column='shared name')
