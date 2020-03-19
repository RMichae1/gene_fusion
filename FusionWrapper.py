import pandas as pd
import os
from typing import List


class FusionWrapper:
    """
    Wraps benchmark tsv output into one format and writes to tsv
    """

    def __init__(self, output_path):
        """
        Takes path of benchmark outputs which are wrapped
        :param output_path:
        """
        self.output_path = output_path
        self.fusion_file = "/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output/{sample}/{caller}/{filename}"
        self.fusion_header = ["FusionCaller", "Sample", "#FusionName", "LeftBreakpoint",
                              "RightBreakpoint"]
        self.star_header = ['#FusionName', 'JunctionReadCount', 'SpanningFragCount',
                            'SpliceType', 'LeftGene', 'LeftBreakpoint', 'RightGene',
                            'RightBreakpoint', 'JunctionReads', 'SpanningFrags',
                            'LargeAnchorSupport', 'FFPM', 'LeftBreakDinuc',
                            'LeftBreakEntropy', 'RightBreakDinuc', 'RightBreakEntropy',
                            'annots', "Sample"]
        self.fmap_header = ['FusionID', 'UniqueCuttingPositionCount',
                            'SeedCount', 'RescuedCount', 'Strand', 'Chromosome1',
                            'Position1', 'Chromosome2', 'Position2', 'KnownGene1',
                            'KnownTranscript1', 'KnownExonNumber1', 'KnownTranscriptStrand1',
                            'KnownGene2', 'KnownTranscript2', 'KnownExonNumber2',
                            'KnownTranscriptStrand2', 'FusionJunctionSequence', 'FusionGene',
                            'SplicePattern', 'SplicePatternClass', 'FrameShift',
                            'FrameShiftClass', 'Distance', 'OnExonBoundary', 'Filter', "Sample"]
        self.fusion_df = self.combine_dfs()

    def combine_dfs(self) -> pd.DataFrame:
        """
        Combine all dataframe of available wrapped callers
        :return: final dataframe of all callers and all samples
        """
        combined_df = pd.DataFrame(columns=self.fusion_header)
        star_fusion_df = self.wrap_star()
        fmap_df = self.wrap_fmap()
        combined_df = combined_df.append(star_fusion_df).append(fmap_df)
        return combined_df

    def wrap_star(self) -> pd.DataFrame:
        star_df = self.read_star_fusion()
        wrapped_star_df = self.wrap_star_fusion(star_df)
        return wrapped_star_df

    def wrap_fmap(self) -> pd.DataFrame:
        fmap_df = self.read_fusion_map()
        unified_df = self.wrap_fmap_fusion(fmap_df)
        return unified_df

    def load_samples_to_df(self, output_df: pd.DataFrame, caller: str, fusion_tsv: str, header: List[str] = None) -> pd.DataFrame:
        """
        iterate over samples and load complete dataframe for each caller
        :param output_df: dataframe which is appended to
        :param caller: caller that generated the output
        :param fusion_tsv: file containing the fusion-data
        :param header: optional in case header needs overwriting
        :return: df of all samples
        """
        for sample in os.listdir(self.output_path):
            # replace when not working with simulated data
            if "sim" not in sample or "_reads" in sample:
                continue
            sample_file = self.fusion_file.format(sample=sample, caller=caller, filename=fusion_tsv)
            tmp_df = pd.read_csv(sample_file, sep="\t")
            # append sample name
            tmp_df["Sample"] = str(sample)
            # overwrite header in case sample names prevent merging
            if header:
                tmp_df.columns = header
            output_df = output_df.append(tmp_df)
        return output_df

    def read_star_fusion(self) -> pd.DataFrame:
        """
        read and merge Star-Fusion output
        :return: merged Star-Fusion dataframe
        """
        empty_star_df = pd.DataFrame(columns=self.star_header)
        merged_star_df = self.load_samples_to_df(output_df=empty_star_df, caller="star_fusion",
                                                 fusion_tsv="star-fusion.fusion_predictions.tsv")
        return merged_star_df

    def wrap_star_fusion(self, star_df: pd.DataFrame) -> pd.DataFrame:
        wrapped_df = pd.DataFrame(columns=self.fusion_header)
        wrapped_df["#FusionName"] = star_df["#FusionName"]
        wrapped_df["LeftBreakpoint"] = star_df["LeftBreakpoint"]
        wrapped_df["RightBreakpoint"] = star_df["RightBreakpoint"]
        wrapped_df["Sample"] = star_df["Sample"]
        wrapped_df["FusionCaller"] = "STAR-Fusion"
        return wrapped_df

    def read_fusion_map(self) -> pd.DataFrame:
        """
        read and merge FusionMap output
        :return: merged fusionmap dataframe
        """
        # FusionMap header elements are dependent on sample input name (Count columns)
        empty_fmap_df = pd.DataFrame(columns=self.fmap_header)
        merged_fmap_df = self.load_samples_to_df(output_df=empty_fmap_df, caller="fmap",
                                                 fusion_tsv="FusionDetection.FusionReport.Table.txt",
                                                 header=self.fmap_header)
        return merged_fmap_df

    def wrap_fmap_fusion(self, fmap_df: pd.DataFrame) -> pd.DataFrame:
        wrapped_df = pd.DataFrame(columns=self.fusion_header)
        wrapped_df["#FusionName"] = fmap_df["FusionGene"].str.replace(".->.", "--", regex=True)
        wrapped_df["LeftBreakpoint"] = "chr" + fmap_df["Chromosome1"].astype(str) + ":" + fmap_df[
            "Position1"].astype(str) + ":" + fmap_df["Strand"].astype(str).str[0]
        wrapped_df["RightBreakpoint"] = "chr" + fmap_df["Chromosome2"].astype(str) + ":" + fmap_df[
            "Position2"].astype(str) + ":" + fmap_df["Strand"].astype(str).str[1]
        wrapped_df["Sample"] = fmap_df["Sample"]
        wrapped_df["FusionCaller"] = "FusionMap"
        return wrapped_df

    def get_fusion_df(self) -> pd.DataFrame:
        return self.fusion_df

    def set_fusion_df(self, fusion_df) -> None:
        self.fusion_df = fusion_df

    def write_to_file(self) -> None:
        filepath = os.path.join(self.output_path, "fusion_benchmark.tsv")
        self.fusion_df.to_csv(filepath, sep="\t", index_label="CallIndex")


if __name__ == '__main__':
    wrapper = FusionWrapper(output_path="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output")
    print("Merged dataframe ... \n {}".format(wrapper.get_fusion_df()))
    wrapper.write_to_file()