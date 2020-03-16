import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
        self.fusion_df = pd.read_csv(fusion_tsv, sep="\t")
        self.truth_path = truth_path
        self.truth50 = np.loadtxt(join(self.truth_path, "sim_50.truth_set.dat"), dtype=str, delimiter="|")
        self.truth101 = np.loadtxt(join(self.truth_path, "sim_101.truth_set.dat"), dtype=str, delimiter="|")
        self.true_fusions50 = self.true_fusions_df(self.truth50)

    def true_fusions_df(self, truth_df):
        """
        returns aggregated df of correctly called fusions per caller in with True,
        False column
        """
        print(self.fusion_df)
        self.fusion_df["Fusion_Truth"] = self.fusion_df["#FusionName"].isin(truth_df["FusionName"])
        selected_df = self.fusion_df.groupby(["FusionCaller", "Sample", "Fusion_Truth"])[
            ["#FusionName"]].count().unstack(fill_value=0).stack().reset_index()
        return selected_df

    def false_negatives_df(self):
        """
        Returns the missed fusions per Caller
        """
        selected_df = None
        return selected_df

    def compute_sensitivity(self):
        tp = 0
        fn = 1
        S = tp / (tp + fn)
        return S

    def compute_loss(self):
        loss = 0
        return loss

    def compute_accuracy(self):
        acc = 0
        return acc

    def compute_f1(self):
        f1 = (2*prec)

    def plot_true_fusions50(self):
        self.plot_true_fusions(truth_file=self.truth50, savefig=True)

    def plot_true_fusions101(self):
        self.plot_true_fusions(truth_file=self.truth101, savefig=True)

    def plot_true_fusions(self, truth_file, savefig=False):
        """
        Plots true fusions barplots stacked for each sample of each caller in the fusion tsv
        :param truth_file:
        :param savefig:
        :return:
        """
        # we have two values per sample, therefore deselect every other
        labels = ["{}-{}".format(caller, sample) for i, (caller, sample)
                  in enumerate(zip(truth_file["FusionCaller"], truth_file["Sample"]))
                  if i % 2 == 0]
        ind = range(len(truth_file["Sample"].unique()) * len(truth_file["FusionCaller"].unique()))

        p1 = plt.bar(ind, np.array(truth_file[truth_file["Fusion_Truth"] == False]["#FusionName"]), width=0.35)
        p2 = plt.bar(ind, np.array(truth_file[truth_file["Fusion_Truth"] == True]["#FusionName"]),
                     bottom=np.array(truth_file[truth_file["Fusion_Truth"] == False]["#FusionName"]),
                     width=0.35)

        plt.ylabel("counts")
        plt.title("Fusions Detected by Caller")
        plt.xticks(ind, labels, rotation=20)
        plt.legend((p1[0], p2[0]), ('False Positives', 'True Positives'))
        if savefig:
            plt.savefig(join(self.out_path, "figures", "true_fusion.png"))
        plt.show()


if __name__ == '__main__':
    benchmark = FusionBenchmark(fusion_tsv="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output/fusion_benchmark.tsv",
                                out_path="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output",
                                truth_path="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/")
    benchmark.plot_true_fusions()