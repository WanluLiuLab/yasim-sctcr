import copy
import functools
import itertools
import json
import operator
import random

import numpy as np

import yasim_sctcr
from labw_utils.bioutils.algorithm.alignment import SmithWatermanAligner
from labw_utils.bioutils.algorithm.sequence import (
    translate_cdna,
    TRANSL_TABLES,
    TRANSL_TABLES_NT,
)
from labw_utils.commonutils.serializer.json import AbstractJSONSerializable

from labw_utils.commonutils.lwio.safe_io import get_reader
from labw_utils.typing_importer import List, Tuple, Dict, Mapping, Any, Optional, Literal, Final

TCRTranslationTableType = List[Tuple[str, str, str, str]]
"""
[(REAL_NT, TRANS_AA, REAL_AA, STATUS)]
"""

Cdr3DeletionTableType = Dict[str, Dict[int, int]]
"""
{gene_type: {n_to_del: freq}}
"""


def align(ref_nt_seq: str, imgt_aa_seq: str) -> TCRTranslationTableType:
    """
    Align TCR on reference to AA sequence in IMGT Database,

    :param ref_nt_seq: Reference TCR NT in reference FASTA.
    :param imgt_aa_seq: IMGT predicted AA.
    :return: Alignment Result
    """
    swas: List[SmithWatermanAligner] = []
    translation_table = []
    for i in (0, 1, 2):
        seq = ref_nt_seq[i:]
        seq = seq[: len(seq) - len(seq) % 3]
        swa = SmithWatermanAligner(
            translate_cdna(seq),
            imgt_aa_seq,
            is_global=False,
            mismatch_score=-6,
            indel_score=-6,
        )
        swas.append(swa)
    mn = np.argmax(list(swa.score for swa in swas))
    backtrack = swas[mn].get_backtrack()[0].splitlines()[1:]
    """
    Should be like:
    
    MRLVARVTVFLTFGTIIDAKTTQPTSMDCAEGRAANLPCNHSTISGNEYVYWYRQIHSQGPQYIIHGLKNNETNEMASLIITEDRKSSTLILPHATLRDTAVYYCIVRV
    DDDDDDDDDDDDDDDDD=======M====================================================================================
    -----------------DAKTTQPPSMDCAEGRAANLPCNHSTISGNEYVYWYRQIHSQGPQYIIHGLKNNETNEMASLIITEDRKSSTLILPHATLRDTAVYYCIVRV
    """
    nt_p = 0
    if mn != 0:
        nt_p = mn
        translation_table.append((ref_nt_seq[0:nt_p], SmithWatermanAligner.GAP, SmithWatermanAligner.GAP, "X"))
    for trans_aa, status, real_aa in zip(*backtrack):
        if status == SmithWatermanAligner.MATCH or status == SmithWatermanAligner.MISMATCH:
            translation_table.append((ref_nt_seq[nt_p : nt_p + 3], trans_aa, real_aa, status))
            nt_p += 3
        elif status == SmithWatermanAligner.INS:
            translation_table.append((SmithWatermanAligner.GAP, SmithWatermanAligner.GAP, real_aa, status))
        elif status == SmithWatermanAligner.DEL:
            translation_table.append((ref_nt_seq[nt_p : nt_p + 3], trans_aa, SmithWatermanAligner.GAP, status))
            nt_p += 3
    if nt_p < len(ref_nt_seq):
        translation_table.append((ref_nt_seq[nt_p:], SmithWatermanAligner.GAP, SmithWatermanAligner.GAP, "X"))
    return translation_table


class FullGenerationRecord:
    uuid: str
    trav: str
    traj: str
    trbv: str
    trbj: str
    trac: str
    trbc: str
    initial_cdr3: Optional[str]

    def __init__(
        self,
        *,
        uuid: str,
        trav: str,
        traj: str,
        trbv: str,
        trbj: str,
        trac: str,
        trbc: str,
    ):
        self.uuid = uuid
        self.trav = trav
        self.traj = traj
        self.trbj = trbj
        self.trbv = trbv
        self.trac = trac
        self.trbc = trbc
        self.initial_cdr3 = None

    def set_initial_cdr3(self, initial_cdr3: str):
        self.initial_cdr3 = initial_cdr3

    def as_tuple(self):
        return (
            self.uuid,
            self.trav,
            self.traj,
            self.trbv,
            self.trbj,
            self.trac,
            self.trbc,
            self.initial_cdr3,
        )


class GenerationFailure(RuntimeError):
    _fgr: FullGenerationRecord

    def __init__(self, fgr: FullGenerationRecord):
        self._fgr = fgr  # TODO: NOT FINISHED!

    @property
    def fgr(self) -> FullGenerationRecord:
        return self._fgr


class Cdr3InsertionTable:
    _table: Dict[str, Dict[int, List[List[int]]]]
    """
    {chain: {num_to_insert: [[freq]]}}
    """
    _num_to_insert: Dict[str, Dict[int, int]]
    """
    {chain: {num_to_insert: freq}}
    """

    def __init__(self, cdr3_insertion_table_path: str) -> None:
        with get_reader(cdr3_insertion_table_path) as reader:
            self._table = json.load(reader)
        self._table["A"] = {int(k): v for k, v in self._table["A"].items()}
        self._table["B"] = {int(k): v for k, v in self._table["B"].items()}
        self._num_to_insert = {
            "A": {i: sum(itertools.chain(*self._table["A"][i])) for i in self._table["A"].keys()},
            "B": {i: sum(itertools.chain(*self._table["B"][i])) for i in self._table["B"].keys()},
        }

    def generate_cdr3(self, chain: Literal["A", "B"]) -> TCRTranslationTableType:
        rets = []
        rng = random.SystemRandom()
        cdr3_length = rng.choices(
            list(self._num_to_insert[chain].keys()),
            list(self._num_to_insert[chain].values()),
        )[0]
        cdr3_pos_freq: List[List[int]] = self._table[chain][cdr3_length]
        for i in range(cdr3_length):
            aa = rng.choices(self._table["AANames"], cdr3_pos_freq[i])[0]
            aa_transl_table = TRANSL_TABLES[1]["AA"]
            rets.append(
                (
                    TRANSL_TABLES_NT[
                        rng.choice([index for index in range(len(aa_transl_table)) if aa_transl_table[index] == aa])
                    ],
                    aa,
                    aa,
                    "=",
                )
            )
        return rets


class TCell(AbstractJSONSerializable):
    """
    T-Cell representation.

    .. versionadded:: 0.1.0
    .. versionchanged:: 0.1.1
         Supports non-productive TCR and TRC genes.
    """

    _title: Final[str] = "TCR"

    @staticmethod
    def _dump_versions() -> Optional[Mapping[str, Any]]:
        return {"yasim-sctcr": yasim_sctcr.__version__}

    @staticmethod
    def _dump_metadata() -> Optional[Mapping[str, Any]]:
        return {}  # Disabled

    @staticmethod
    def _validate_versions(versions: Mapping[str, Any]) -> None:
        pass  # Disabled

    @classmethod
    def from_dict(cls, in_dict: Mapping[str, Any]):
        return cls(**in_dict)

    _traj_name: str
    _trbj_name: str
    _trav_name: str
    _trbv_name: str
    _trac_name: str
    _trbc_name: str

    _traj_tt: TCRTranslationTableType
    _trav_tt: TCRTranslationTableType
    _trac_tt: TCRTranslationTableType

    _trbj_tt: TCRTranslationTableType
    _trbv_tt: TCRTranslationTableType
    _trbc_tt: TCRTranslationTableType

    _tra_cdr3_tt: TCRTranslationTableType
    _trb_cdr3_tt: TCRTranslationTableType
    _tcr_uuid: str

    def __init__(
        self,
        *,
        traj_name: str,
        trbj_name: str,
        trav_name: str,
        trbv_name: str,
        trac_name: str,
        trbc_name: str,
        traj_tt: TCRTranslationTableType,
        trav_tt: TCRTranslationTableType,
        trac_tt: TCRTranslationTableType,
        trbj_tt: TCRTranslationTableType,
        trbv_tt: TCRTranslationTableType,
        trbc_tt: TCRTranslationTableType,
        tra_cdr3_tt: TCRTranslationTableType,
        trb_cdr3_tt: TCRTranslationTableType,
        tcr_uuid: str,
    ):
        self._traj_name = traj_name
        self._trbj_name = trbj_name
        self._trav_name = trav_name
        self._trbv_name = trbv_name
        self._trac_name = trac_name
        self._trbc_name = trbc_name
        self._traj_tt = traj_tt
        self._trav_tt = trav_tt
        self._trbj_tt = trbj_tt
        self._trbv_tt = trbv_tt
        self._trac_tt = trac_tt
        self._trbc_tt = trbc_tt
        self._tra_cdr3_tt = tra_cdr3_tt
        self._trb_cdr3_tt = trb_cdr3_tt
        self._tcr_uuid = tcr_uuid

    @classmethod
    def from_gene_names(
        cls,
        *,
        tcr_genelist: Dict[str, List[str]],
        cdr3_deletion_table: Cdr3DeletionTableType,
        cdr3_insertion_table: Cdr3InsertionTable,
        tcr_cache: Dict[str, TCRTranslationTableType],
        usage_bias_tra: Dict[str, int],
        usage_bias_trb: Dict[str, int],
        tcr_uuid: str,
        is_productive: bool,
    ):
        def choose_name_jv(
            chain_type: Literal["a", "b"]
        ) -> Tuple[Tuple[str, TCRTranslationTableType], Tuple[str, TCRTranslationTableType]]:
            usage_bias = usage_bias_tra if chain_type == "a" else usage_bias_trb
            rg = random.SystemRandom()
            j_name, v_name = rg.choices(list(usage_bias.keys()), list(usage_bias.values()), k=1)[0].split(":")
            return (j_name, copy.deepcopy(tcr_cache[j_name])), (v_name, copy.deepcopy(tcr_cache[v_name]))

        def choose_name(tcr_type: str) -> Tuple[str, TCRTranslationTableType]:
            rg = random.SystemRandom()
            while True:
                tcr_name = rg.choice(tcr_genelist[tcr_type])
                if tcr_name in tcr_cache:
                    return tcr_name, copy.deepcopy(tcr_cache[tcr_name])

        def clip_aa(
            tr_cdr3_tt: TCRTranslationTableType,
            trv_tt: TCRTranslationTableType,
            trj_tt: TCRTranslationTableType,
            j_term: str,
        ) -> Tuple[TCRTranslationTableType, TCRTranslationTableType, TCRTranslationTableType]:
            trv_tt_real_aa = "".join(list(zip(*trv_tt))[2])
            trj_tt_real_aa = "".join(list(zip(*trj_tt))[2])
            c_idx = trv_tt_real_aa[::-1].find("C")
            if c_idx == -1 or c_idx > 0.5 * len(trv_tt_real_aa):
                raise GenerationFailure(fgr)
            f_idx = trj_tt_real_aa.find(j_term)
            if f_idx == -1 or f_idx > 0.5 * len(trj_tt_real_aa):
                raise GenerationFailure(fgr)
            tr_cdr3_tt = trv_tt[-c_idx - 1 : -1] + tr_cdr3_tt + trj_tt[0 : f_idx + 1]
            trv_tt = trv_tt[0 : -c_idx - 1]
            trj_tt = trj_tt[f_idx + 2 :]
            if len(trv_tt) * len(trj_tt) == 0 or not (10 < len(tr_cdr3_tt) < 20):
                raise GenerationFailure(fgr)
            return tr_cdr3_tt, trv_tt, trj_tt

        def clip_ltr(tt: TCRTranslationTableType) -> TCRTranslationTableType:
            """
            TODO: This function clips 5' UTRs only. May also needs to clip 3' UTRs.
            """
            i = 0
            for i in range(len(tt)):
                tt_a = tt[i]
                if "-" in tt_a[0] or tt_a[1] in ("*", "-") or tt_a[2] in ("*", "-"):
                    continue
                else:
                    break
            return tt[i:]

        def productive_nonproductive_adjustment(tt: TCRTranslationTableType) -> TCRTranslationTableType:
            ret_tt = []
            if is_productive:
                for tt_a in tt:
                    if "-" in tt_a[0] or tt_a[1] in ("*", "-") or tt_a[2] in ("*", "-"):
                        continue
                    ret_tt.append(tt_a)
            else:
                for tt_a in tt:
                    if "-" in tt_a[0] or tt_a[1] == "-" or tt_a[2] == "-":
                        tt_a[1] = "*"
                        tt_a[2] = "*"
                        tt_a[0] = tt_a[0].replace("-", "")
                    if tt_a[0] == "":
                        continue
                    ret_tt.append(tt_a)
            return ret_tt

        (traj_name, traj_tt), (trav_name, trav_tt) = choose_name_jv("a")
        (trbj_name, trbj_tt), (trbv_name, trbv_tt) = choose_name_jv("b")
        trbc_name, trbc_tt = choose_name("trbc_names")
        trac_name, trac_tt = choose_name("trac_names")
        fgr = FullGenerationRecord(
            uuid=tcr_uuid,
            traj=traj_name,
            trav=trav_name,
            trac=trac_name,
            trbj=trbj_name,
            trbv=trbv_name,
            trbc=trbc_name,
        )
        rng = random.SystemRandom()

        chosen_deletion: Dict[str, int] = {
            k: rng.choices(
                population=list(cdr3_deletion_table[k].keys()),
                weights=list(cdr3_deletion_table[k].values()),
            )[0]
            for k in cdr3_deletion_table.keys()
        }

        if functools.reduce(operator.mul, chosen_deletion.values()) > 625:
            raise GenerationFailure(fgr)

        def clip_nt(
            tr_tt: TCRTranslationTableType,
            gene_name: str,
            pos: int,
        ) -> None:
            """
            Delete terminal untranslated NTs and generated deletions

            :param pos: Position to clip from. Use 0 for 5' and -1 for 3'.
            :param gene_name: Name of the gene.
            :param tr_tt: TranslationTable object for the gene.
            """
            num_clip = 0
            try:
                while tr_tt[pos][-1] == "X":
                    _ = tr_tt.pop(pos)
                    num_clip += 1
                for _ in range(chosen_deletion[gene_name]):
                    _ = tr_tt.pop(pos)
                    num_clip += 1
            except IndexError:
                raise GenerationFailure(fgr)
            if len(tr_tt) == 0:
                raise GenerationFailure(fgr)

        clip_nt(trav_tt, "trav", -1)
        clip_nt(trbv_tt, "trbv", -1)
        clip_nt(traj_tt, "traj", 0)
        clip_nt(trbj_tt, "trbj", 0)

        tra_cdr3_tt = cdr3_insertion_table.generate_cdr3("A")
        trb_cdr3_tt = cdr3_insertion_table.generate_cdr3("B")

        tra_cdr3_tt, trav_tt, traj_tt = clip_aa(
            tra_cdr3_tt,
            trav_tt,
            traj_tt,
            "W" if traj_name.upper() in {"TRAJ33", "TRAJ38", "TRAJ55"} else "F",
        )
        trb_cdr3_tt, trbv_tt, trbj_tt = clip_aa(trb_cdr3_tt, trbv_tt, trbj_tt, "F")
        trac_tt = trac_tt[:15]
        trbc_tt = trbc_tt[:15]
        tra_cdr3_tt = productive_nonproductive_adjustment(clip_ltr(tra_cdr3_tt))
        trav_tt = productive_nonproductive_adjustment(clip_ltr(trav_tt))
        traj_tt = productive_nonproductive_adjustment(clip_ltr(traj_tt))
        trac_tt = productive_nonproductive_adjustment(clip_ltr(trac_tt))
        trb_cdr3_tt = productive_nonproductive_adjustment(clip_ltr(trb_cdr3_tt))
        trbv_tt = productive_nonproductive_adjustment(clip_ltr(trbv_tt))
        trbj_tt = productive_nonproductive_adjustment(clip_ltr(trbj_tt))
        trbc_tt = productive_nonproductive_adjustment(clip_ltr(trbc_tt))

        retc = cls(
            tcr_uuid=tcr_uuid,
            trav_tt=trav_tt,
            traj_tt=traj_tt,
            tra_cdr3_tt=tra_cdr3_tt,
            trbv_tt=trbv_tt,
            trbj_tt=trbj_tt,
            trb_cdr3_tt=trb_cdr3_tt,
            traj_name=traj_name,
            trav_name=trav_name,
            trbj_name=trbj_name,
            trbv_name=trbv_name,
            trac_name=trac_name,
            trbc_name=trbc_name,
            trac_tt=trac_tt,
            trbc_tt=trbc_tt,
        )
        if is_productive and ("*" in retc.alpha_aa or "*" in retc.beta_aa):
            raise GenerationFailure(fgr)
        if not is_productive and not ("*" in retc.alpha_aa or "*" in retc.beta_aa):
            raise GenerationFailure(fgr)
        return retc

    def to_nt_fasta_record(self) -> str:
        return "\n".join(
            (
                f">{self._tcr_uuid}_A",
                self.alpha_nt,
                f">{self._tcr_uuid}_B",
                self.beta_nt,
            )
        )

    def to_aa_fasta_record(self) -> str:
        return "\n".join(
            (
                f">{self._tcr_uuid}_A",
                self.alpha_aa,
                f">{self._tcr_uuid}_B",
                self.beta_aa,
            )
        )

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "tcr_uuid": self._tcr_uuid,
            "traj_name": self._traj_name,
            "trbj_name": self._trbj_name,
            "trav_name": self._trav_name,
            "trbv_name": self._trbv_name,
            "trac_name": self._trac_name,
            "trbc_name": self._trbc_name,
            "traj_tt": self._traj_tt,
            "trav_tt": self._trav_tt,
            "trac_tt": self._trac_tt,
            "trbj_tt": self._trbj_tt,
            "trbv_tt": self._trbv_tt,
            "trbc_tt": self._trbc_tt,
            "tra_cdr3_tt": self._tra_cdr3_tt,
            "trb_cdr3_tt": self._trb_cdr3_tt,
        }

    @property
    def tcr_uuid(self):
        return self._tcr_uuid

    @property
    def alpha_names(self) -> Tuple[str, str, str]:
        return self._trav_name, self._traj_name, self._trac_name

    @property
    def beta_names(self) -> Tuple[str, str, str]:
        return self._trbv_name, self._trbj_name, self._trbc_name

    @property
    def alpha_nt(self) -> str:
        return "".join(
            itertools.chain(
                list(zip(*self._trav_tt))[0],
                list(zip(*self._tra_cdr3_tt))[0],
                list(zip(*self._traj_tt))[0],
                list(zip(*self._trac_tt))[0],
            )
        )

    @property
    def beta_nt(self) -> str:
        return "".join(
            itertools.chain(
                list(zip(*self._trbv_tt))[0],
                list(zip(*self._trb_cdr3_tt))[0],
                list(zip(*self._trbj_tt))[0],
                list(zip(*self._trbc_tt))[0],
            )
        )

    @property
    def alpha_aa(self) -> str:
        return "".join(
            itertools.chain(
                list(zip(*self._trav_tt))[1],
                list(zip(*self._tra_cdr3_tt))[1],
                list(zip(*self._traj_tt))[1],
                list(zip(*self._trac_tt))[1],
            )
        )

    @property
    def beta_aa(self) -> str:
        return "".join(
            itertools.chain(
                list(zip(*self._trbv_tt))[1],
                list(zip(*self._trb_cdr3_tt))[1],
                list(zip(*self._trbj_tt))[1],
                list(zip(*self._trbc_tt))[1],
            )
        )

    @property
    def cdr3_aa(self) -> Tuple[str, str]:
        return (
            "".join(list(zip(*self._tra_cdr3_tt))[1]),
            "".join(list(zip(*self._trb_cdr3_tt))[1]),
        )

    @property
    def cdr3_nt(self) -> Tuple[str, str]:
        return (
            "".join(list(zip(*self._tra_cdr3_tt))[0]),
            "".join(list(zip(*self._trb_cdr3_tt))[0]),
        )
