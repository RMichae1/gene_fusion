import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
from os.path import join


class FusionBenchmark:
    """
    Computes Benchmarking values and creates figures for the provided gene-fusion tsv
    E.g. This can be called with the 50nt or 101nt benchmark tsv file
    """

    def __init__(self, fusion_tsv, out_path, truth_file, analysis_name="sim50"):
        """
        :param fusion_df: tsv containing all called fusions (Caller, Sample, Name, LeftBreak, RightBreak)
        :param out_path: path to write output to
        :param truth_path: path to truth files
        """
        self.output = out_path
        self.analysis_name = analysis_name
        self.stats_dir = join(self.output, "stats_{}".format(analysis_name))
        self.fig_dir = join(self.output, "figures_{}".format(analysis_name))
        if not os.path.isdir(self.stats_dir):
            os.mkdir(self.stats_dir)
        if not os.path.isdir(self.fig_dir):
            os.mkdir(self.fig_dir)
        self.chr_ix = self.set_chr_ix()
        self.fusion_df = pd.read_csv(fusion_tsv, sep="\t")
        if self.analysis_name == "sim50":
            self.samples = ["sim_adipose", "sim_brain", "sim_colon", "sim_heart", "sim_testis"]
        elif self.analysis_name == "sim101":
            self.samples = ["sim1_reads", "sim2_reads", "sim3_reads", "sim4_reads", "sim5_reads"]
        else:
            raise ValueError("USAGE: Specify 50nt or 101nt analysis.")
        self.subselect_analysis()
        self.truth = pd.read_csv(truth_file, sep="|", header=None, names=("Sample", "FusionName"))
        self.true_fusions_df = self.find_true_fusions(self.truth)
        self.benchmark_df = self.compute_stats()
        self.breakpoints_df = self.get_breakpoint_dist(self.fusion_df)

    @staticmethod
    def set_chr_ix():
        chr_ix = ["chr{}".format(i + 1) for i in range(22)]
        chr_ix.append("chrX")
        chr_ix.append("chrY")
        return chr_ix

    def subselect_analysis(self):
        self.fusion_df = self.fusion_df[self.fusion_df["Sample"].isin(self.samples)]

    def find_true_fusions(self, truth_df) -> pd.DataFrame:
        """
        returns aggregated df of correctly called fusions per caller in with True,
        False column
        """
        self.fusion_df["Fusion_Truth"] = self.fusion_df["#FusionName"].isin(truth_df["FusionName"])
        selected_df = self.fusion_df.groupby(["FusionCaller", "Sample", "Fusion_Truth"])[
            ["#FusionName"]].count().unstack(fill_value=0).stack().reset_index()
        return selected_df

    def write_true_fusions(self) -> pd.DataFrame:
        self.true_fusions_df.to_csv(join(self.stats_dir, "true_fusions.tsv"), sep="\t")

    def compute_stats(self):
        """
        Compute statistics of individual fusioncallers, this is
        #True Positives, #False Positives, #False Negatives, Precision, Recall, F1 Score
        :return: df with calculated stats
        """
        benchmark_df = pd.DataFrame(columns=["FusionCaller", "Sample", "TP", "FP",
                                             "FN", "Precision"])
        caller_sample_df = self.fusion_df.groupby(["FusionCaller", "Sample"]).first().reset_index()
        benchmark_df["FusionCaller"] = caller_sample_df["FusionCaller"]
        benchmark_df["Sample"] = caller_sample_df["Sample"]
        benchmark_df["TP"] = self.true_fusions_df[self.true_fusions_df["Fusion_Truth"]==True].reset_index()["#FusionName"]
        benchmark_df["FP"] = self.true_fusions_df[self.true_fusions_df["Fusion_Truth"]==False].reset_index()["#FusionName"]
        benchmark_df["Precision"] = self.compute_precision(tp=benchmark_df["TP"], fp=benchmark_df["FP"])

        false_negatives_df = self.false_negatives(truth_df=self.truth)
        # this df is not sorted, therefore merging was necessary to get the right values in place
        benchmark_df = pd.merge(benchmark_df, false_negatives_df, how="left", on=["FusionCaller", "Sample"])
        benchmark_df = benchmark_df.dropna(axis=1, how='all')
        benchmark_df.rename(columns={'FN_y': 'FN'}, inplace=True)

        benchmark_df["Recall"] = self.compute_tpr(tp=benchmark_df["TP"], fn=benchmark_df["FN"])
        pos_bench = benchmark_df[(benchmark_df["Precision"] > 0) & (benchmark_df["Recall"] > 0)]
        benchmark_df["F1"] = self.f1(prec=pos_bench["Precision"], recall=pos_bench["Recall"])
        benchmark_df = benchmark_df.fillna(0)
        return benchmark_df

    def write_stats(self):
        stats_name = "caller_stats_{}.tsv".format(self.analysis_name)
        self.benchmark_df.to_csv(join(self.stats_dir, stats_name, sep="\t"))

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

    def get_breakpoint_dist(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        calculates breakpoints given chromosomal positions over all fusions
        :param df: dataframe with fusions to compute break-dist over
        :return: dataframe of counts of breakpoint locations
        """

        breakpoints_df = pd.DataFrame(np.zeros((len(self.chr_ix), len(self.chr_ix))),
                                      columns=self.chr_ix, index=self.chr_ix)
        # expand on dataframe to get chr locations of the Breaks
        df[["LChr", "LPos", "LType"]] = df.LeftBreakpoint.str.split(":", expand=True)
        df[["RChr", "RPos", "RType"]] = df.RightBreakpoint.str.split(":", expand=True)
        # fill df with breaks
        for left, right in zip(df["LChr"], df["RChr"]):
            # account for unclean indices
            if len(left) > 5:
                left = left[:5]
            if len(right) > 5:
                right = right[:5]
            if "_" in left:
                left = left.split("_")[0]
            if "_" in right:
                right = right.split("_")[0]
            breakpoints_df.loc[left, right] += 1
        return breakpoints_df

    @staticmethod
    def compute_tpr(tp, fn):
        """
        calculates true positive rate
        :param tp: True positives
        :param fn: False negatives
        :return: TPR
        """
        return tp / (tp + fn)

    @staticmethod
    def compute_precision(tp, fp):
        """
        calculates precision
        :param tp: True positives
        :param fp: false positives
        :return: precision
        """
        return tp / (tp + fp)

    @staticmethod
    def f1(prec, recall):
        """
        calculates unweighted F-Score, as a measure between precision and recall
        :param prec: precision
        :param recall: recall
        :return: f1-score
        """
        f_score = 2*((prec * recall)/(prec + recall))
        return f_score

    def plot_true_fusions(self, savefig=False):
        """
        Plots true fusions barplots stacked for each sample of each caller in the fusion tsv
        :param truth_file:
        :param savefig:
        :return:
        """
        # we have two values per sample, therefore deselect every other
        labels = ["{}-{}".format(call_sample[0], call_sample[1]) for i, call_sample
                  in enumerate(list(zip(self.true_fusions_df.FusionCaller, self.true_fusions_df.Sample)))
                  if i % 2 == 0]
        ind = range(len(self.true_fusions_df["Sample"].unique()) * len(self.true_fusions_df["FusionCaller"].unique()))

        p1 = plt.bar(ind, np.array(self.true_fusions_df[
                                       self.true_fusions_df["Fusion_Truth"]==False]["#FusionName"]), width=0.35)
        p2 = plt.bar(ind, np.array(self.true_fusions_df[
                                       self.true_fusions_df["Fusion_Truth"]==True]["#FusionName"]),
                     bottom=np.array(self.true_fusions_df[
                                         self.true_fusions_df["Fusion_Truth"]==False]["#FusionName"]), width=0.35)

        plt.ylabel("counts")
        plt.title("Fusions Detected by Caller")
        plt.xticks(ind, labels, rotation=90, fontsize=8)
        plt.legend((p1[0], p2[0]), ('False Positives', 'True Positives'))
        plt.tight_layout()
        if savefig:
            filename = "true_fusions_{}.png".format(self.analysis_name)
            plt.savefig(join(self.output, "figures_{}".format(self.analysis_name),  filename))
        plt.show()
        plt.clf()

    def plot_stats(self, savefig=False) -> None:
        boxplot = self.benchmark_df.boxplot(column=["Precision", "Recall", "F1"],
                                            by="FusionCaller", rot=15)
        plt.suptitle("")
        if savefig:
            plt.savefig(join(self.output, "figures_{}".format(self.analysis_name), "stats_boxplot.png"))
        plt.show()

    def plot_break_heatmap_comparison(self, caller1: str =None, caller2: str =None,
                                      savefig: bool=False) -> None:
        """
        Plots heatmap over breakpoint distribution for two callers in the given analysis
        :param caller1:
        :param caller2:
        :param savefig: boolean if plot should be written to output
        :return: None
        """
        caller1_df = self.fusion_df[self.fusion_df["FusionCaller"] == caller1]
        df1 = self.get_breakpoint_dist(caller1_df)
        caller2_df = self.fusion_df[self.fusion_df["FusionCaller"] == caller2]
        df2 = self.get_breakpoint_dist(caller2_df)

        chr_range = np.arange(len(self.chr_ix))

        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(5, 5))

        im1 = ax0.imshow(df1)
        im2 = ax1.imshow(df2)

        fig.colorbar(im1, ax=ax1, orientation="horizontal")

        ax0.set_xticks(chr_range)
        ax0.set_xticklabels(list(df1))
        ax0.set_yticks(chr_range)
        ax0.set_yticklabels(list(df2))

        ax1.set_xticks(chr_range)
        ax1.set_xticklabels(list(df2))
        ax1.set_yticks(chr_range)
        ax1.set_yticklabels([])

        plt.setp(ax0.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=8)
        plt.setp(ax0.get_yticklabels(), fontsize=8)
        plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=8)
        plt.setp(ax1.get_yticklabels(), fontsize=8)

        if savefig:
            filename = "break_heatmap_{}_{}.png".format(caller1, caller2)
            plt.savefig(join(self.output, "figures_{}".format(self.analysis_name), filename))

        fig.suptitle("\n\n\nFusion Breakpoints {} and {}".format(caller1, caller2))
        fig.tight_layout()
        fig.show()

    def plot_break_distribution(self, savefig: bool = False, caller: str = None) -> None:
        """
        Plots multiplot over distribution of all breakpoints, if caller is specified a subselection
        for only the caller data will be used
        :param savefig: boolean if output should be written to file
        :param caller: optional subselection
        :return: None
        """
        if caller:
            # save calculated breakpoint distribution and calculate the one for the subselect
            backup_df = self.breakpoints_df.copy()
            caller_selection = self.fusion_df[self.fusion_df["FusionCaller"] == caller]
            self.breakpoints_df = self.get_breakpoint_dist(caller_selection)

        left_breaks = self.breakpoints_df.sum(axis=1)
        right_breaks = self.breakpoints_df.sum()

        chr_range = np.arange(len(self.chr_ix))

        fig = plt.figure(figsize=(15, 15))
        gs = gridspec.GridSpec(2, 2, height_ratios=[4, 1], width_ratios=[3, 1])
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
        ax2 = plt.subplot(gs[2])
        ax3 = plt.subplot(gs[3])
        # fig.subplots_adjust(top=0.85)

        im = ax0.imshow(self.breakpoints_df, aspect="auto")
        fig.colorbar(im, ax=ax1)

        # plot heatmap of the distribution
        ax0.set_title("Fusion Distributions over breakpoints (l|r)")
        ax0.set_xticks(chr_range)
        ax0.set_xticklabels([])
        ax0.set_yticks(chr_range)
        ax0.set_yticklabels(list(self.breakpoints_df.index))

        # add text to fields
        for i in chr_range:
            for j in chr_range:
                text = ax0.text(j, i, self.breakpoints_df.iloc[i, j],
                                ha="center", va="center", color="w")

        # add barplot for right break distribution
        ax2.bar(chr_range, right_breaks, align="edge")
        ax2.set_xlim(0, len(chr_range))
        ax2.set_ylim(0, max(right_breaks) + 5)
        ax2.set_xticks(chr_range)
        ax2.set_xticklabels(list(self.breakpoints_df))

        # add barplot for left break distribution
        # reindex reverse sort of left breakage otherwise reversed
        ax1.barh(chr_range, left_breaks.iloc[::-1], align="edge")
        ax1.set_ylim(0, len(chr_range))
        ax1.set_xlim(0, max(left_breaks) + 5)
        ax1.set_xticks(np.arange(0, max(left_breaks) + 5, step=20, dtype=int))
        ax1.set_xticklabels(np.arange(0, max(left_breaks) + 5, step=20, dtype=int))
        ax1.set_yticks([])

        ax3.axis("off")

        plt.setp(ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=9)

        if savefig:
            filename = "break_dist_{}{}.png".format(self.analysis_name, caller)
            plt.savefig(join(self.output, "figures_{}".format(self.analysis_name), filename))

        if caller:
            self.breakpoints_df = backup_df.copy()

        fig.tight_layout()
        fig.show()


if __name__ == '__main__':
    # first analysis is for simulated data 50nt basis
    benchmark50 = FusionBenchmark(
        fusion_tsv="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output/fusion_benchmark.tsv",
        out_path="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output",
        truth_file="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/sim_50.truth_set.dat",
        analysis_name="sim50")
    benchmark50.plot_true_fusions(savefig=True)
    benchmark50.write_true_fusions()
    benchmark50.write_stats()
    print(benchmark50.benchmark_df)
    benchmark50.plot_stats(savefig=True)
    benchmark50.plot_break_distribution(savefig=True)
    benchmark50.plot_break_heatmap_comparison(caller1="FusionMap", caller2="STAR-Fusion", savefig=True)
    benchmark50.plot_break_distribution(savefig=True, caller="FusionMap")
    benchmark50.plot_break_distribution(savefig=True, caller="STAR-Fusion")

    # second analysis for longer 101nt samples
    benchmark101 = FusionBenchmark(
        fusion_tsv="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output/fusion_benchmark.tsv",
        out_path="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output",
        truth_file="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/sim_101.truth_set.dat",
        analysis_name="sim101")
    benchmark101.plot_true_fusions(savefig=True)
    benchmark101.write_true_fusions()
    benchmark101.write_stats()
    benchmark101.plot_stats(savefig=True)
    benchmark101.plot_break_heatmap_comparison(caller1="FusionMap", caller2="STAR-Fusion", savefig=True)
    benchmark101.plot_break_distribution(savefig=True)