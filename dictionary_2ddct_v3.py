import sys
import pandas as pd
import numpy as np
import pymysql
from pymongo import MongoClient


def main():
    database_type = sys.argv[1]
    dict_gene_ddct1 = Dictionary2ddctVersion3(database_type)
    dict_gene_ddct = dict_gene_ddct1.get_dict_2ddct_groups()

    print(dict_gene_ddct)


class Dictionary2ddctVersion3:
    def __init__(self, database_type):
        self._database_type = database_type
        self.df_complete_2_ddct = self.df_2_delta_delta()

    def df_2_delta_delta(self):
        df = self.collect_data_from_database()
        df_processed = self.process_df_for_operations(df)
        array_housekeeping_column = self.get_array_housekeeping_column(df_processed)
        df_dct = self.get_delta_ct(df_processed, array_housekeeping_column)
        df_2_ddct = self.apply_2_delta_delta(df_dct)
        df_complete_2_ddct = self.add_categorical_columns_to_df(df, df_2_ddct)

        return df_complete_2_ddct

    @staticmethod
    def create_df_from_mongodb():
        host = 'localhost'
        port = 27017
        client = MongoClient(host, port)
        mongo_database = client.database_Teresa
        collection = mongo_database.liver_expression_data
        records_of_collection = list(collection.find())
        df = pd.DataFrame(records_of_collection)

        return df

    @staticmethod
    def create_df_from_mysql():
        host = 'localhost'
        user = 'root'
        password = 'password'
        database = 'Gene_Expression'
        mysql_db = pymysql.connect(host, user, password, database)
        data_query = 'SELECT * FROM gene_ct'
        df = pd.read_sql(data_query, con=mysql_db)

        return df

    def collect_data_from_database(self):
        if self._database_type == 'mongo':
            df = self.create_df_from_mongodb()

            return df

        else:
            df = self.create_df_from_mysql()

            return df

    @staticmethod
    def process_df_for_operations(df):
        df_rounded = df.iloc[:, 3:].round(2)
        df_processed = df_rounded.replace({0.00: np.nan})

        return df_processed

    @staticmethod
    def get_array_housekeeping_column(df_processed):
        array_housekeeping_column = df_processed.iloc[:, 0].to_list()

        return array_housekeeping_column

    @staticmethod
    def get_delta_ct(df_processed, array_housekeeping):
        df_dct = df_processed.iloc[:, 1:].subtract(array_housekeeping, axis=0)

        return df_dct

    @staticmethod
    def apply_2_delta_delta(df_dct):
        df_2_ddct = df_dct.apply(lambda x: 2 ** (-x))

        return df_2_ddct

    def add_categorical_columns_to_df(self, df, df_2_ddct):
        if self._database_type == 'mongo':
            df_no_id_column = df.drop(columns=['_id'])
        else:
            df_no_id_column = df.drop(columns=['id'])
        df_renamed_exp_group_column = df_no_id_column.rename(columns={'exp_group': 'Group'})
        categorical_columns = df_renamed_exp_group_column.iloc[:, :2]
        parts_for_aggrupation = [categorical_columns, df_2_ddct]
        df_complete_2_ddct = pd.concat(parts_for_aggrupation, axis=1)

        return df_complete_2_ddct

    @staticmethod
    def mask_method(df_only_numeric_values, var):
        q1 = df_only_numeric_values[var].quantile(0.25)
        q3 = df_only_numeric_values[var].quantile(0.75)
        iqr = q3 - q1
        mask = ((df_only_numeric_values[var] > q1 - (1.5 * iqr)) &
                (df_only_numeric_values[var] < q3 + (1.5 * iqr)))

        return mask

    def get_dict_2ddct(self, df_group):
        df_only_numeric_values = df_group.iloc[:, 2:]
        dict_2ddct = {}
        for var in df_only_numeric_values:
            mask_var = self.mask_method(df_only_numeric_values, var)
            df_mask = df_only_numeric_values[mask_var]
            column_var = df_mask[var]
            dict_2ddct[var] = list(column_var)

        return dict_2ddct

    def get_dict_2ddct_groups(self):
        dict_2ddct_groups = {}
        array_experimental_groups = list(set(self.df_complete_2_ddct['Group']))
        for group in array_experimental_groups:
            mask_group = self.df_complete_2_ddct['Group'] == group
            df_group = self.df_complete_2_ddct[mask_group]
            dict_2ddct_groups[group] = self.get_dict_2ddct(df_group)

        return dict_2ddct_groups


if __name__ == '__main__':
    main()
