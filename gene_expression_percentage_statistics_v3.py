import sys
import numpy as np
from scipy import stats
import datetime
import pandas as pd
from dict_gene_expression_percentage_v3 import DictGeneExpressionPercentageVersion3


def main():
    database_type = sys.argv[1]
    gene_expression_statistics1 = GeneExpressionPercentageStatisticsVersion3(database_type)
    gene_expression_statistics1.file_dataframes_statistics()


class GeneExpressionPercentageStatisticsVersion3:
    def __init__(self, database_type):
        self.dict_percentage_group_gene = \
            DictGeneExpressionPercentageVersion3(
                database_type).get_dict_expression_percent_group_gene()
        self.dict_percentage_gene_group = \
            DictGeneExpressionPercentageVersion3(
                database_type).get_dict_expression_percent_gene_group()
        self.array_genes = DictGeneExpressionPercentageVersion3(
            database_type).get_array_genes()

    def file_dataframes_statistics(self):
        df_mean = self.convert_dict_mean_to_dataframe()
        df_univariate_test = self.convert_dict_univariate_to_dataframe()
        self.save_excel_file_dataframes(df_mean, df_univariate_test)

    @staticmethod
    def get_dict_percentage_mean(group_data):
        dict_percentage_mean = {}
        for gene, array_percentage in group_data.items():
            mean_percent = np.mean(array_percentage)
            mean_percent_rounded = mean_percent.round(2)
            dict_percentage_mean[gene] = mean_percent_rounded

        return dict_percentage_mean

    def get_dict_percentage_mean_groups(self):
        dict_percentage_mean_groups = {}
        for group_name, group_data in self.dict_percentage_group_gene.items():
            dict_percentage_mean_groups[group_name] = self.get_dict_percentage_mean(group_data)

        return dict_percentage_mean_groups

    @staticmethod
    def normality_test_method(array_percentage):
        array_length = len(array_percentage)

        if array_length >= 8:
            k2, normaltest_p = stats.normaltest(array_percentage)

            if normaltest_p > 0.05:
                return True

        return False

    def get_dict_normal_test(self, group_data):
        dict_normal_test = {}
        for gene, array_percentage in group_data.items():
            dict_normal_test[gene] = self.normality_test_method(array_percentage)

        return dict_normal_test

    def get_dict_normal_test_groups(self):

        dict_normal_test_groups = {}
        for group_name, group_data in self.dict_percentage_group_gene.items():
            dict_normal_test_groups[group_name] = self.get_dict_normal_test(group_data)

        return dict_normal_test_groups

    @staticmethod
    def levene_test_method(gene_data):
        array1, array2, array3, array4, array5 = gene_data.values()
        w_value, levene_p = stats.levene(array1,
                                         array2,
                                         array3,
                                         array4,
                                         array5)
        if levene_p > 0.05:
            return True
        else:
            return False

    def get_dict_levene_test(self):
        dict_levene_test = {}
        for gene, gene_data in self.dict_percentage_gene_group.items():
            dict_levene_test[gene] = self.levene_test_method(gene_data)
        return dict_levene_test

    def combine_normaltest_dict_by_group(self, gene):
        dict_normaltest = self.get_dict_normal_test_groups()
        normality = True
        for group_name, group_data in dict_normaltest.items():
            normality_in_group = group_data[gene]
            normality = normality and normality_in_group

        return normality

    def get_dict_normaltest_combined(self):
        dict_normaltest = self.get_dict_normal_test_groups()
        dict_normaltest_combined = {}
        group1 = list(dict_normaltest.keys())[0]
        array_genes_normaltest = dict_normaltest.get(group1).keys()

        for gene in array_genes_normaltest:
            combined_normaltest = self.combine_normaltest_dict_by_group(gene)
            dict_normaltest_combined[gene] = combined_normaltest

        return dict_normaltest_combined

    @staticmethod
    def oneway_anova_method(gene_data):
        array1, array2, array3, array4, array5 = gene_data.values()
        f_value, oneway_p = stats.f_oneway(array1,
                                           array2,
                                           array3,
                                           array4,
                                           array5)
        if oneway_p < 0.05:
            return 'oneway_p: {}, Reject Ho'.format(oneway_p)
        else:
            return 'oneway_p: {}, Accept Ho'.format(oneway_p)

    @staticmethod
    def kruskal_wallis_method(gene_data):
        array1, array2, array3, array4, array5 = gene_data.values()
        h_value, kruskal_p = stats.kruskal(array1,
                                           array2,
                                           array3,
                                           array4,
                                           array5)
        if kruskal_p < 0.05:
            return 'kruskal_p: {}, Reject Ho'.format(kruskal_p)
        else:
            return 'kruskal_p: {}, Accept Ho'.format(kruskal_p)

    def dict_stats_test_application(self):
        dict_normal = self.get_dict_normaltest_combined()
        dict_levene = self.get_dict_levene_test()
        dict_stats_test = {}
        for gene in self.array_genes:
            gene_data = self.dict_percentage_gene_group.get(gene)
            if dict_normal.get(gene) is True and dict_levene.get(gene) is True:
                stats_test = self.oneway_anova_method(gene_data)
            else:
                stats_test = self.kruskal_wallis_method(gene_data)

            dict_stats_test[gene] = stats_test

        return dict_stats_test

    def convert_dict_mean_to_dataframe(self):
        dict_mean = self.get_dict_percentage_mean_groups()
        return pd.DataFrame.from_dict(dict_mean, orient='index')

    def convert_dict_univariate_to_dataframe(self):
        dict_univariate_test = self.dict_stats_test_application()
        return pd.DataFrame.from_dict(dict_univariate_test, orient='index')

    @staticmethod
    def save_excel_file_dataframes(df_mean, df_univariate):
        excel_writer = pd.ExcelWriter(('gene_expression_stats' +
                                       datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') +
                                       '.xlsx'), engine='xlsxwriter')
        df_mean.to_excel(excel_writer, sheet_name='mean')
        df_univariate.to_excel(excel_writer, sheet_name='univariate_test', header=False)
        excel_writer.save()

    # @staticmethod
    # def save_file_dataframes(df_mean, df_univariate_test):
    #     new_file_name = 'gene_expression_stats'
    #     file_df_mean_univariate = Path(new_file_name)
    #     file_df_mean_univariate.touch()
    #     file_df_mean_univariate.write_text('Gene expression mean \n' + str(df_mean) + '\n\n' +
    #                                        'Gene expression univariate test \n' +
    #                                        str(df_univariate_test))
    #     file_df_mean_univariate.rename(new_file_name + '_' +
    #                                    datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') +
    #                                    '.txt')
    #
    #     return file_df_mean_univariate


if __name__ == '__main__':
    main()
