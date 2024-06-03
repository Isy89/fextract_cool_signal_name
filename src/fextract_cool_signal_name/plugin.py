import logging
import pathlib
from typing import Any
from typing import List
from typing import Optional

import click
import matplotlib
import pandas as pd
from lbfextract.utils import generate_time_stamp
from lbfextract.utils_classes import Signal

import lbfextract.fextract
from lbfextract.core import App
from lbfextract.fextract.schemas import Config, AppExtraConfig, ReadFetcherConfig, SingleSignalTransformerConfig, \
    SignalSummarizer

import pysam
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

dict_parameters = {
    100: {
        "figsize": (20, 10),
        "fontsize_title": 20,
        "fontsize_xlabel": 15,
        "fontsize_ylabel": 15
    },
    300: {
        "figsize": (50, 10),
        "fontsize_title": 40,
        "fontsize_xlabel": 30,
        "fontsize_ylabel": 30
    },
    500: {
        "figsize": (70, 10),
        "fontsize_title": 50,
        "fontsize_xlabel": 35,
        "fontsize_ylabel": 35
    },
    1000: {
        "figsize": (100, 10),
        "fontsize_title": 55,
        "fontsize_xlabel": 40,
        "fontsize_ylabel": 40
    }
}


def select_parameters(reads):
    read_length = len(reads)

    match read_length:
        case x if x <= 100:
            return dict_parameters[100]
        case x if x <= 300:
            return dict_parameters[300]
        case x if x <= 500:
            return dict_parameters[500]
        case x if x <= 1000:
            return dict_parameters[1000]
        case _:
            raise ValueError("Read length exceeds the available configurations")


class FextractHooks:
    @lbfextract.hookimpl
    def transform_single_intervals(self, transformed_reads: pd.DataFrame, config: Any,
                                   extra_config: Any) -> Signal:
        """
        :param transformed_reads: ReadsPerIntervalContainer containing a list of ReadsPerInterval which are
            basically lists with information about start and end of the interval
        :param config: config specific to the function
        :param extra_config: config containing context information plus extra parameters
        """

        def count_rads_in_interval(x: pysam.libcalignmentfile.IteratorRowRegion) -> int:
            return len(list(x))

        transformed_reads.index = transformed_reads.Chromosome.astype("str") + "-" + transformed_reads.Start.astype(
            "str") + "-" + transformed_reads.End.astype("str")

        reads_per_interval_pd_series = transformed_reads.reads_per_interval.apply(lambda x: count_rads_in_interval(x))

        return Signal(reads_per_interval_pd_series.values,
                      metadata=reads_per_interval_pd_series.index.copy(),
                      tags=tuple(["cool_signal", ]))

    @lbfextract.hookimpl
    def transform_all_intervals(self, single_intervals_transformed_reads: Signal, config: Any,
                                extra_config: Any) -> Signal:
        """
        :param single_intervals_transformed_reads: Signal object containing the signals per interval
        :param config: config specific to the function
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        return single_intervals_transformed_reads

    @lbfextract.hookimpl
    def plot_signal(self, signal: Signal, config: Any, extra_config: Any) -> matplotlib.figure.Figure:
        """
        :param signal: Signal object containing the signals per interval
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        reads_per_interval_pd_series = pd.Series(signal.array, index=signal.metadata)
        signal_type = "_".join(signal.tags) if signal.tags else ""

        try:
            parameters = select_parameters(reads_per_interval_pd_series)
        except ValueError:
            logger.warning("Plot not generated. Too many intervals to be plot")
            fig, ax = plt.subplots()
            return fig

        fig, ax = plt.subplots(figsize=parameters["figsize"])
        reads_per_interval_pd_series.plot(kind='bar', ax=ax)
        ax.set_title("Reads per interval", fontsize=parameters["fontsize_title"])
        ax.set_xlabel("Intervals", fontsize=parameters["fontsize_xlabel"])
        ax.set_ylabel("Reads count", fontsize=parameters["fontsize_ylabel"])

        fig.savefig(
            extra_config.ctx["output_path"] /
            f"{generate_time_stamp()}__{extra_config.ctx['id']}__{signal_type}_signal_plot.pdf",
            dpi=300)
        return fig


class CliHook:
    @lbfextract.hookimpl_cli
    def get_command(self) -> click.Command | List[click.Command]:
        @click.command()
        @click.option('--path_to_bam', type=click.Path(exists=False,
                                                       file_okay=True,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the bam file to be used')
        @click.option('--path_to_bed', type=click.Path(exists=False,
                                                       file_okay=True,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the bed file to be used')
        @click.option('--output_path', type=click.Path(exists=False,
                                                       file_okay=False,
                                                       dir_okay=True,
                                                       writable=True,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the output directory')
        @click.option("--skip_read_fetching", is_flag=True, show_default=True,
                      help='Boolean flag. When it is set, the fetching of the reads is skipped and the latest'
                           'timestamp of this run (identified by the id) is retrieved')
        @click.option("--exp_id", default=None, type=str, show_default=True,
                      help="run id")
        @click.option("--window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted around the middle point of an"
                           " interval present in the bedfile")
        @click.option("--flanking_window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted after the window")
        @click.option("--extra_bases", default=2000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted from the bamfile when removing the "
                           "unused bases to be sure to get all the proper paires, which may be mapping up to 2000 bs")
        @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                      help="number of intervals to be used to extract the signal, if it is higher then the provided"
                           "intervals, all the intervals will be used")
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min", "skip"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: (Undefined, Undefined, Undefined, Undefined)")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--flip_based_on_strand", is_flag=True,
                      show_default=True,
                      default=False,
                      help="flip the signal based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        # Here you can add the options you want to be available in the command line (they should match the ones in the
        # function)
        def extract_cool_signal_name(
                path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
                output_path: pathlib.Path,
                skip_read_fetching,
                window: int,
                flanking_window: int,
                extra_bases: int,
                n_binding_sites: int,
                summarization_method: str,
                cores: int,
                exp_id: Optional[str],
                flip_based_on_strand: bool,
                gc_correction_tag: Optional[str]
        ):
            """
            Cool signal plugin. Example of plugin for LBFextract.
            This plugin extracts the number of reads per interval provided in a BED file.
            """
            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites
            }
            reads_transformer_config = {}
            single_signal_transformer_config = {
                "signal_transformer": "coverage",
                "flip_based_on_strand": flip_based_on_strand,
                "gc_correction": True if gc_correction_tag else False,
                "tag": gc_correction_tag

            }
            transform_all_intervals_config = {
                "summarization_method": summarization_method
            }
            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            res = App(plugins_name=["cool_signal_name", ],
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=Config(reads_transformer_config),
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=SignalSummarizer(transform_all_intervals_config),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

        return extract_cool_signal_name


hook = FextractHooks()
hook_cli = CliHook()
