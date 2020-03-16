import pandas as pd
import os


class FusionWrapper:

    def __init__(self, output_path):
        self.output_path = output_path
        self.fusion_file = "/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output/{sample}/{caller}/{filename}"
        self.fusion_df = self.combine_dfs()
        self.star_header = ['#FusionName', 'JunctionReadCount', 'SpanningFragCount',
                       'SpliceType', 'LeftGene', 'LeftBreakpoint', 'RightGene',
                       'RightBreakpoint', 'JunctionReads', 'SpanningFrags',
                       'LargeAnchorSupport', 'FFPM', 'LeftBreakDinuc',
                       'LeftBreakEntropy', 'RightBreakDinuc', 'RightBreakEntropy', 'annots']
        self.fmap_header = ['FusionID', 'UniqueCuttingPositionCount',
                       'SeedCount', 'RescuedCount', 'Strand', 'Chromosome1',
                       'Position1', 'Chromosome2', 'Position2', 'KnownGene1',
                       'KnownTranscript1', 'KnownExonNumber1', 'KnownTranscriptStrand1',
                       'KnownGene2', 'KnownTranscript2', 'KnownExonNumber2',
                       'KnownTranscriptStrand2', 'FusionJunctionSequence', 'FusionGene',
                       'SplicePattern', 'SplicePatternClass', 'FrameShift',
                       'FrameShiftClass', 'Distance', 'OnExonBoundary', 'Filter']
        self.fusion_header = ["FusionCaller", "Sample", "#FusionName", "LeftBreakpoint", "RightBreakpoint"]

    def combine_dfs(self) -> pd.DataFrame:
        combined_df = pd.DataFrame(columns=self.fusion_header)
        star_fusion_df = self.wrap_star()
        fmap_df = self.wrap_fmap()
        combined_df = combined_df.append(star_fusion_df).append(fmap_df)
        return combined_df

    def wrap_star(self) -> pd.DataFrame:
        star_df = self.read_star_fusion()
        wrapped_star_df = self.wrap_star(star_df)
        return wrapped_star_df

    def wrap_fmap(self) -> pd.DataFrame:
        fmap_df = self.read_fusion_map()
        unified_df = self.wrap_fmap(fmap_df)
        return unified_df

    def load_samples_to_df(self, output_df, fusion_tsv, header=None) -> pd.DataFrame:
        for sample in os.listdir(self.output_path):
            if "sim" not in sample:
                continue
            sample_file = self.fusion_file.format(sample=sample, caller="star_fusion",
                                                  filename="star-fusion.fusion_predictions.tsv")
            tmp_df = pd.read_csv(sample_file, sep="\t")
            # overwrite header in case sample names prevent merging
            if header:
                tmp_df.rename(columns=header)
            output_df.append(tmp_df)
        return output_df

    def read_star_fusion(self) -> pd.DataFrame:
        empty_star_df = pd.DataFrame(columns=self.star_header)
        merged_star_df = self.load_samples_to_df(output_df=empty_star_df,
                                                 fusion_tsv="star-fusion.fusion_predictions.tsv")
        return merged_star_df

    def wrap_star_fusion(self, star_df) -> pd.DataFrame:
        wrapped_df = pd.DataFrame(columns=self.read_fusion_map)
        wrapped_df["#FusionName"] = star_df["#FusionName"]
        wrapped_df["LeftBreakpoint"] = star_df["LeftBreakpoint"]
        wrapped_df["RightBreakpoint"] = star_df["RightBreakpoint"]
        wrapped_df["FusionCaller"] = "STAR-Fusion"
        wrapped_df["Sample"] = "sim_adipose" # TODO assign this dynamically when loading
        return wrapped_df

    def read_fusion_map(self):
        # FusionMap header elements are dependent on sample input name (Count columns)
        empty_fmap_df = pd.DataFrame(columns=self.fmap_header)
        merged_fmap_df = self.load_samples_to_df(output_df=empty_fmap_df,
                                                 fusion_tsv="FusionDetection.FusionReport.Table.txt",
                                                 header=self.fmap_header)
        return merged_fmap_df

    def wrap_fmap(self, fmap_df) -> pd.DataFrame:

        wrapped_df = pd.DataFrame(columns=self.fusion_header)
        wrapped_df["#FusionName"] = fmap_df["FusionGene"].str.replace(".->.", "--", regex=True)
        wrapped_df["LeftBreakpoint"] = "chr" + fmap_df["Chromosome1"].astype(str) + ":" + fmap_df[
            "Position1"].astype(str) + ":" + fmap_df["Strand"].astype(str).str[0]
        wrapped_df["RightBreakpoint"] = "chr" + fmap_df["Chromosome2"].astype(str) + ":" + fmap_df[
            "Position2"].astype(str) + ":" + fmap_df["Strand"].astype(str).str[1]
        wrapped_df["Sample"] = "sim_adipose"
        wrapped_df["FusionCaller"] = "FusionMap"
        return wrapped_df

    def get_fusion_df(self):
        return self.fusion_df

    def write_to_file(self):
        filepath = os.path.join(self.output_path, "fusion_benchmark.csv")
        self.fusion_header.to_tsv(filepath)


if __name__ == '__main__':
    wrapper = FusionWrapper(output_path="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/output")
    print(wrapper.get_fusion_df())