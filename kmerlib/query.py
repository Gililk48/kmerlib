# kmerlib/query.py
from pathlib import Path
from typing import Iterable, Iterator, Tuple, Union
import json, struct, numpy as np

HDR = "<8sIIQQ"
HDR_LEN = struct.calcsize(HDR)

class K30Index:
    def __init__(self, prefix: Union[str, Path]):
        p = Path(prefix)
        self._index = p if p.suffix == ".bin" else p.with_suffix(".index.bin")
        self._base  = self._index.with_suffix("").with_suffix("")          # strip .index.bin
        self._load_meta()
        self._mmap_arrays()

    # ── helpers ─────────────────────────────────────────
    def _load_meta(self):
        with open(self._base.with_suffix(".bkidx.json")) as fh:
            self.meta = json.load(fh)
        self.is_reads = self.meta["reads"]

    def _mmap_arrays(self):
        with open(self._index, "rb") as fp:
            magic, ver, _, nK, nO = struct.unpack(HDR, fp.read(HDR_LEN))
            assert magic.startswith(b"KIDX") and ver == 1
        self._keys = np.memmap(self._index, "<u8", mode="r",
                               offset=HDR_LEN, shape=(nK,))
        off_off    = HDR_LEN + nK*8
        self._offs = np.memmap(self._index, "<u8", mode="r",
                               offset=off_off, shape=(nK+1,))
        id_off     = off_off + (nK+1)*8
        self._ids  = np.memmap(self._index, "<u4", mode="r",
                               offset=id_off, shape=(nO,))
        pos_off    = id_off + nO*4
        self._pos  = np.memmap(self._index, "<u4", mode="r",
                               offset=pos_off, shape=(nO,))

    def _loc(self, k:int):
        i = int(np.searchsorted(self._keys, k))
        return i if i < len(self._keys) and self._keys[i]==k else None

    # ── public API ─────────────────────────────────────
    def present(self, k:int) -> bool:
        return self._loc(k) is not None

    def count(self, k:int) -> int:
        i = self._loc(k)
        return 0 if i is None else int(self._offs[i+1] - self._offs[i])

    def count_many(self, keys: Iterable[int]):
        # order-preserving, handles duplicates
        return np.fromiter((self.count(int(k)) for k in keys), dtype="<u8")

    def positions(self, k:int) -> Iterator[Tuple[int,int]]:
        i = self._loc(k)
        if i is None:
            return iter(())
        lo, hi = map(int, (self._offs[i], self._offs[i+1]))
        return ((int(self._ids[j]), int(self._pos[j])) for j in range(lo, hi))
