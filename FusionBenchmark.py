import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os.path import join


class FusionBenchmark:

    def __init__(self, fusion_df, out_path, truth_path):
        self.output = out_path
        self.fusion_df = fusion_df
        self.truth_path = truth_path
        self.truth50 = np.loadtxt(self.truth_path, dtype=str, delimiter="|")
        self.truth101 = np.loadtxt(self.truth_path, dtype=str, delimiter="|")

    def true_fusions_df(self):
        """
        returns aggregated df of correctly called fusions per caller in with True,
        False column
        """
        selected_df = self.fusion_df.groupby(["FusionCaller",
                                              "Sample",
                                              "Fusion_Truth"])[["#FusionName"]].count().unstack(
            fill_value=0).stack().reset_index()
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
        pass

    def compute_accuracy(self):
        pass

    def compute_f1(self):
        pass

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

        p1 = plt.bar(ind, np.array(truth_file[truth_file["Fusion_Truth"] == False]["#FusionName"]),
                     width=0.35)
        p2 = plt.bar(ind, np.array(truth_file[truth_file["Fusion_Truth"] == True]["#FusionName"]),
                     bottom=np.array(truth_file[truth_file["Fusion_Truth"] == False]["#FusionName"]),
                     width=0.35)

        plt.ylabel("counts")
        plt.title("Fusions Detected by Caller")
        plt.xticks(ind, labels, rotation=20)
        plt.legend((p1[0], p2[0]), ('False Positives', 'True Positives'))
        if savefig:
            plt.savefig(join(self.out_path, "true_fusion.png"))
        plt.show()