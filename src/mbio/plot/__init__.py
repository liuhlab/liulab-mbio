"""Draw a sequence record as a map through molviz.

`draw_map` is the way in: it takes a record, or a file to read one from, hands it to molviz with
the options chosen, and returns molviz's `Drawing`, whose `write` writes the format a file's
suffix names. `converter` writes a record in molviz's own format, as `molviz_record`, and `cli`
holds the `plot` verbs.
"""

from mbio.plot.converter import molviz_record
from mbio.plot.drawing import draw_map

__all__ = ["draw_map", "molviz_record"]
