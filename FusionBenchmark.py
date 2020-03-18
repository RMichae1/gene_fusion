import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from os.path import join


class FusionBenchmark:
    """
    Computes Benchmarking values and creates figures for the provided gene-fusion tsv
    """

    def __init__(self, fusion_tsv, out_path, truth_path):
        """
        :param fusion_df: tsv containing all called fusions (Caller, Sample, Name, LeftBreak, RightBreak)
        :param out_path: path to write output to
        :param truth_path: path to truth files
        """
        self.output = out_path
        self.stats_dir = join(self.output, "stats")
        self.fig_dir = join(self.output, "figures")
        if not os.path.isdir(self.stats_dir):
            os.mkdir(self.stats_dir)
        if not os.path.isdir(self.fig_dir):
            os.mkdir(self.fig_dir)
        self.fusion_df = pd.read_csv(fusion_tsv, sep="\t")
        self.truth_path = truth_path
        self.truth50_file = join(self.truth_path, "sim_50.truth_set.dat")
        self.truth101_file = join(self.truth_path, "sim_101.truth_set.dat")
        self.truth50 = pd.read_csv(self.truth50_file, sep="|", header=None, names=("Sample", "FusionName"))
        self.truth101 = pd.read_csv(self.truth101_file, sep="|", header=None, names=("Sample", "FusionName"))
        self.true_fusions50 = self.true_fusions_df(self.truth50)
        self.true_fusions101 = None
        self.benchmark50_df = self.compute_stats50()
        # self.true_fusion101 = self.true_fusions_df(self.truth101

    def true_fusions_df(self, truth_df):
        """
        returns aggregated df of correctly called fusions per caller in with True,
        False column
        """
        self.fusion_df["Fusion_Truth"] = self.fusion_df["#FusionName"].isin(truth_df["FusionName"])
        selected_df = self.fusion_df.groupby(["FusionCaller", "Sample", "Fusion_Truth"])[
            ["#FusionName"]].count().unstack(fill_value=0).stack().reset_index()
        return selected_df

    def write_true_fusions50(self):
        self.true_fusions50.to_csv(join(self.stats_dir, "true_fusions50.tsv"), sep="\t")

    def write_true_fusions101(self):
        self.true_fusions101.to_csv(join(self.stats_dir, "true_fusions101.tsv"), sep="\t")

    def compute_stats50(self):
        benchmark_df = pd.DataFrame(columns=["FusionCaller", "Sample", "TP", "FP",
                                             "FN", "Precision"])
        caller_sample_df = self.fusion_df.groupby(["FusionCaller", "Sample"]).first().reset_index()
        benchmark_df["FusionCaller"] = caller_sample_df["FusionCaller"]
        benchmark_df["Sample"] = caller_sample_df["Sample"]
        benchmark_df["TP"] = self.true_fusions50[self.true_fusions50["Fusion_Truth"]==True].reset_index()["#FusionName"]
        benchmark_df["FP"] = self.true_fusions50[self.true_fusions50["Fusion_Truth"]==False].reset_index()["#FusionName"]
        benchmark_df["Precision"] = self.compute_precision(tp=benchmark_df["TP"], fp=benchmark_df["FP"])

        false_negatives_df = self.false_negatives(truth_df=self.truth50)
        # this df is not sorted, therefore merging was necessary to get the right values in place
        benchmark_df = pd.merge(benchmark_df, false_negatives_df, how="left", on=["FusionCaller", "Sample"])
        benchmark_df = benchmark_df.dropna(axis=1, how='all')
        benchmark_df.rename(columns={'FN_y': 'FN'}, inplace=True)

        benchmark_df["Recall"] = self.compute_tpr(tp=benchmark_df["TP"], fn=benchmark_df["FN"])
        pos_bench = benchmark_df[(benchmark_df["Precision"] > 0) & (benchmark_df["Recall"] > 0)]
        benchmark_df["F1"] = self.f1(prec=pos_bench["Precision"], recall=pos_bench["Recall"])
        benchmark_df = benchmark_df.fillna(0)
        return benchmark_df

    def write_stats50(self):
        self.benchmark50_df.to_csv(join(self.stats_dir, "caller_stats.tsv", sep="\t"))

    def false_negatives(self, truth_df: pd.DataFrame) -> pd.Series:
        """
        Returns fusions per caller that were not detected
        :param truth_df: truth dataframe
        :return: Series of missed calls counted
        """
        false_neg_df = pd.DataFrame(columns=["FusionCaller", "Sample", "FusionName"])
        for caller in self.fusion_df["FusionCaller"].unique():
            caller_selection = self.fusion_df[self.fusion_df["FusionCaller"]==caller]
            false_negatives = truth_df[~truth_df["FusionName"].isin(
                caller_selection["#FusionName"])].groupby("Sample").count().reset_index()
            # caller is not set as column, as this is implicit from the call
            false_negatives["FusionCaller"] = caller
            false_neg_df = false_neg_df.append(false_negatives, ignore_index=True)
        false_neg_df.rename(columns={'FusionName': 'FN'}, inplace=True)
        return false_neg_df

    @staticmethod
    def compute_tpr(tp, fn):
        return tp / (tp + fn)

    @staticmethod
    def compute_precision(tp, fp):
        return tp / (tp + fp)

    @staticmethod
    def f1(prec, recall):
        f_score = 2*((prec * recall)/(prec + recall))
        return f_score

    def plot_true_fusions50(self):
        self.plot_true_fusions(true_fusions_df=self.true_fusions50, savefig=True)

    def plot_true_fusions101(self):
        self.plot_true_fusions(true_fusions_df=self.true_fusions101, savefig=True)

    def plot_true_fusions(self, true_fusions_df, savefig=False):
        """
        Plots true fusions barplots stacked for each sample of each caller in the fusion tsv
        :param truth_file:
        :param savefig:
        :return:
        """
        # we have two values per sample, therefore deselect every other
        labels = ["{}-{}".format(caller, sample) for i, (caller, sample)
                  in enumerate(zip(true_fusions_df["FusionCaller"], true_fusions_df["Sample"]))
                  if i % 2 == 0]
        ind = range(len(true_fusions_df["Sample"].unique()) * len(true_fusions_df["FusionCaller"].unique()))

        p1 = plt.bar(ind, np.array(true_fusions_df[true_fusions_df["Fusion_Truth"] == False]["#FusionName"]), width=0.35)
        p2 = plt.bar(ind, np.array(true_fusions_df[true_fusions_df["Fusion_Truth"] == True]["#FusionName"]),
                     bottom=np.array(true_fusions_df[true_fusions_df["Fusion_Truth"] == False]["#FusionName"]),
                     width=0.35)

        plt.ylabel("counts")
        plt.title("Fusions Detected by Caller")
        plt.xticks(ind, labels, rotation=20)
        plt.legend((p1[0], p2[0]), ('False Positives', 'True Positives'))
        plt.tight_layout()
        if savefig:
            plt.savefig(join(self.output, "figures", "true_fusion.png"))
        plt.show()
        plt.clf()

    def plot_stats(self, savefig=False):
        boxplot = self.benchmark50_df.boxplot(column=["Precision", "Recall", "F1"], by="FusionCaller")
        if savefig:
            plt.savefig(join(self.output, "figures", "stats_boxplot.png"))
        plt.show()


if __name__ == '__main__':
    benchmark = FusionBenchmark(fusion_tsv="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output/fusion_benchmark.tsv",
                                out_path="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output",
                                truth_path="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/")
    benchmark.plot_true_fusions50()
    benchmark.write_true_fusions50()
    print(benchmark.benchmark50_df)
    benchmark.plot_stats(savefig=True)