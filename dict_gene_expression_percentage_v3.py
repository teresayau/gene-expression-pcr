import sys
import numpy as np
from dictionary_2ddct_v3 import Dictionary2ddctVersion3


def main():
    database_type = sys.argv[1]
    dict_gene_expression1 = DictGeneExpressionPercentageVersion3(database_type)
    dict_expression_percent_gene_group = dict_gene_expression1.get_dict_expression_percent_group_gene()
    dict_expression_percent_group_gene = dict_gene_expression1.get_dict_expression_percent_gene_group()

    print(dict_expression_percent_gene_group, dict_expression_percent_group_gene)


class DictGeneExpressionPercentageVersion3:
    def __init__(self, database_type):
        self.dict_2ddct = Dictionary2ddctVersion3(database_type).get_dict_2ddct_groups()
        self.dict_expression_percent_group_gene = self.get_dict_expression_percent_group_gene()

    def get_dict_2ddct_mean_control(self):
        control_group = self.dict_2ddct.get('Control')
        dict_2ddct_mean = {}
        for gene, array_ct in control_group.items():
            ct_mean = np.mean(array_ct)
            dict_2ddct_mean[gene] = ct_mean

        return dict_2ddct_mean

    def get_dict_expression_percent(self, group_data):
        dict_mean_control = self.get_dict_2ddct_mean_control()
        dict_data_percentage = {}
        for gene, array_ct in group_data.items():
            if gene in dict_mean_control:
                division_of_ct_values = array_ct / dict_mean_control.get(gene)
                percentage_calculation = division_of_ct_values * 100
                dict_data_percentage[gene] = percentage_calculation

        return dict_data_percentage

    def get_dict_expression_percent_group_gene(self):
        dict_data_percentage_groups = {}
        for group_name, group_data in self.dict_2ddct.items():
            dict_data_percentage_groups[group_name] = self.get_dict_expression_percent(group_data)

        return dict_data_percentage_groups

    def get_array_genes(self):
        group_1 = list(self.dict_expression_percent_group_gene.keys())[0]
        array_genes = self.dict_expression_percent_group_gene.get(group_1).keys()

        return array_genes

    def get_dict_expression_percent_per_gene(self, var):
        dict_expression_percent_per_gene = {}
        for group_name, group_data in self.dict_expression_percent_group_gene.items():
            var_on_group = group_data[var]
            dict_expression_percent_per_gene[group_name] = var_on_group

        return dict_expression_percent_per_gene

    def get_dict_expression_percent_gene_group(self):
        array_genes = self.get_array_genes()
        dict_expression_percent_gene_group = {}
        for var in array_genes:
            dict_expression_percent_gene_group[var] = self.get_dict_expression_percent_per_gene(var)

        return dict_expression_percent_gene_group


if __name__ == '__main__':
    main()
