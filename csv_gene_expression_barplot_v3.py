import sys
import matplotlib.pyplot as plt
from dict_gene_expression_percentage_v3 import DictGeneExpressionPercentageVersion3
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import datetime


def main():
    database_type = sys.argv[1]
    gene_expression_plot1 = CsvGeneExpressionBarplot3(database_type)
    gene_expression_plot1.plot_gene_expression_data()


class CsvGeneExpressionBarplot3:
    def __init__(self, database_type):
        self.dict_percentage_gene_group = \
            DictGeneExpressionPercentageVersion3(
                database_type).get_dict_expression_percent_gene_group()

    def plot_gene_expression_data(self):
        dict_mean = self.get_dict_mean_gene_group()
        dict_sem = self.get_dict_sem_gene_group()
        self.barplot_gene_expression(dict_mean, dict_sem)

    @staticmethod
    def get_dict_mean(gene_data):
        dict_mean = {}
        for group_name, group_data in gene_data.items():
            mean_group = group_data.mean()
            dict_mean[group_name] = mean_group

        return dict_mean

    def get_dict_mean_gene_group(self):
        dict_mean_gene_group = {}
        for gene_name, gene_data in self.dict_percentage_gene_group.items():
            mean_group_gene = self.get_dict_mean(gene_data)
            dict_mean_gene_group[gene_name] = mean_group_gene

        return dict_mean_gene_group

    @staticmethod
    def get_dict_sem_group(gene_data):
        dict_sem = {}
        for group_name, group_data in gene_data.items():
            sem_group = stats.sem(group_data)
            dict_sem[group_name] = sem_group

        return dict_sem

    def get_dict_sem_gene_group(self):
        dict_sem_gene_group = {}
        for gene_name, gene_data in self.dict_percentage_gene_group.items():
            sem_group_gene = self.get_dict_sem_group(gene_data)
            dict_sem_gene_group[gene_name] = sem_group_gene

        return dict_sem_gene_group

    @staticmethod
    def barplot_gene_expression(dict_mean, dict_sem):
        with PdfPages('gene_expression_' +
                      datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S') +
                      '.pdf') as pdf_barplots:
            array_genes_barplot = dict_mean.keys()
            for gene in array_genes_barplot:
                mean_gene_groups = dict_mean.get(gene).values()
                sem_gene_groups = dict_sem.get(gene).values()
                x_axis = dict_mean.get(gene).keys()
                plt.clf()
                plt.bar(x_axis, mean_gene_groups,
                        yerr=sem_gene_groups,
                        ecolor='black',
                        capsize=5,
                        edgecolor='black')
                plt.title(gene)
                pdf_barplots.savefig()


if __name__ == '__main__':
    main()