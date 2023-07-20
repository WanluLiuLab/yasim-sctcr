"""
TODO: docs

.. versionadded:: 0.1.0
"""
import glob
import os
from collections import defaultdict

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Mapping, Any, Type, Optional, Dict
from yasim.helper.depth_io import DepthType, read_depth, DepthParsingException
from yasim.helper.rna_seq import validate_adapter_args, run_rna_seq
from yasim.llrg_adapter import BaseLLRGAdapter

_lh = get_logger(__name__)


def sc_rna_seq_frontend(
        transcriptome_fasta_dir: str,
        output_fastq_dir_path: str,
        depth_dir_path: str,
        jobs: int,
        simulator_name: str,
        adapter_args: Mapping[str, Any],
        assembler_args: Mapping[str, Any],
        adapter_class: Type[BaseLLRGAdapter],
        is_pair_end: bool,
        llrg_executable_path: Optional[str] = None,
        not_perform_assemble: bool = False
) -> int:
    """
    :return: Exit Value

    .. versionadded:: 3.1.5
    """
    adapter_args = validate_adapter_args(
        adapter_args=adapter_args,
        adapter_class=adapter_class,
        llrg_executable_path=llrg_executable_path
    )

    barcode_depth_data_dict: Dict[str, DepthType] = {}
    for depth_file_path in tqdm(glob.glob(os.path.join(depth_dir_path, "*")), desc="Validating depth input"):
        try:
            depth_data = read_depth(depth_file_path)
        except DepthParsingException:
            _lh.error(f"RNA SEQ: Failed to parse depth file {depth_file_path}")
            return 1
        barcode = os.path.basename(depth_file_path).split(".")[0]
        if barcode in barcode_depth_data_dict:
            _lh.error(f"RNA SEQ: Duplicated barcode {barcode}")
            return 1
        barcode_depth_data_dict[barcode] = depth_data

    full_exception_dict = defaultdict(lambda: 0)
    for barcode, depth_data in barcode_depth_data_dict.items():
        exception_dict = run_rna_seq(
            transcriptome_fasta_dir=transcriptome_fasta_dir,
            output_fastq_prefix=os.path.join(output_fastq_dir_path, barcode),
            depth_data=depth_data,
            jobs=jobs,
            simulator_name=simulator_name,
            adapter_args=adapter_args,
            assembler_args=assembler_args,
            adapter_class=adapter_class,
            is_pair_end=is_pair_end,
            llrg_executable_path=llrg_executable_path,
            not_perform_assemble=not_perform_assemble,
            show_tqdm=False
        )
        for k, v in exception_dict.items():
            full_exception_dict[k] += v

    _lh.info(f"RNA SEQ: Status of errors: {dict(full_exception_dict)}")
    _lh.info("RNA SEQ: Simulation finished successfully")
    return 0
