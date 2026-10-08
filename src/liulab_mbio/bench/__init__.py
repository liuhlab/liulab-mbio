"""The bench any pipeline shares: reactions, programs, validation and protocol steps.

Every public name in its modules imports from here as well, `goldengate` and `coverage` excepted.
A module that cites a source keeps it as its own `REFERENCES`; `REFERENCES` here gathers the ones
re-exported.

- `amounts`: the weight of DNA, picomoles from nanograms, and what to pipette.
- `reactions`: one reaction table, filled to volume, and what it refuses.
- `pcr`: the PCR and colony PCR reactions and programs.
- `gels`: the ladder and agarose percentage a range of bands takes.
- `goldengate`: NEB's Golden Gate enzymes, assembly reaction and cycling, which the Golden Gate
  cloning pipeline and a library build both run. The other methods' numbers stay with their own
  pipeline, because only this chemistry has two callers. Imported by module: its `REFERENCES`
  and `SOURCES` are one chemistry's, not this package's.
- `validation`: the colony PCR that reads an assembly's junctions, and the Sanger reads.
- `coverage`: how many constructs a set of part lists yields, the colonies a round takes before
  it loses one, and the marks a pooled library's representation read is held to. Imported by
  module: its `REFERENCES` are one sizing rule's, not this package's.
- `readback`: what one well's read came to, where a construct is read back one well at a time,
  and how often one picked colony is clean.
- `inactivation`: an enzyme's heat inactivation.
- `materials`: what a material brings with it -- its own program, what it puts in the tube, and
  the rules it carries. A number attaches to the thing that changes it.
- `plates`: a plate at any format, where each thing sits in it, and the moves between wells.
- `prices`: a banded price record the user holds, the charge it gives a quantity, and the bill.
- `phenotype`: what a clone is expected to show, read off the product's own features.
- `oligos`: the primer order sheet.
- `steps`: the protocol steps any pipeline reuses, built from plain facts, and the smaller
  pieces both pipelines shape the same way. Its two references are cited only by a protocol
  that runs the step they belong to.

Nothing here imports `liulab_mbio.cloning`.
"""

from liulab_mbio.bench import amounts, gels, pcr
from liulab_mbio.bench.amounts import (
    DNA_VOLUME_UL,
    Amount,
    dna_amount,
    molecular_weight,
    to_nanograms,
    to_pmol,
)
from liulab_mbio.bench.gels import LADDER_1_KB_PLUS, LADDER_100_BP, agarose_percent, choose_ladder
from liulab_mbio.bench.inactivation import heat_inactivation
from liulab_mbio.bench.materials import Electroporation, electroporation, material
from liulab_mbio.bench.oligos import (
    SHEET_COLUMNS,
    OrderedOligo,
    oligo_row,
    ordered_row,
    primer_sheet,
)
from liulab_mbio.bench.pcr import (
    COLONY_HOLD_CELSIUS,
    COLONY_LYSIS_SECONDS,
    COLONY_PCR_MASTER_MIX,
    COLONY_PCR_VOLUME_UL,
    DNTP_STOCK_MM,
    PRIMER_STOCK_UM,
    colony_pcr_master_mix_component,
    colony_pcr_program,
    colony_pcr_reaction,
    pcr_program,
    pcr_reaction,
)
from liulab_mbio.bench.phenotype import (
    MEDIUM,
    SELECTION,
    Phenotype,
    read_phenotype,
    selection_marker,
)
from liulab_mbio.bench.plates import (
    PrimerPlates,
    compact,
    plate,
    pool,
    primer_plates,
    seat,
    wells_of,
)
from liulab_mbio.bench.prices import Item, PriceRecord, PriceRow, bill, read_prices
from liulab_mbio.bench.reactions import WATER, dna_components, fits, reaction_table
from liulab_mbio.bench.readback import (
    CLEAN_COLONY_CURVE,
    WellVerdict,
    clean_colony_chance,
    identity_check,
    reformat,
)
from liulab_mbio.bench.steps import (
    ASSEMBLY_UL,
    CELLS_UL,
    COLONY_PCR_TITLE,
    DAM_SITE,
    DPNI_CELSIUS,
    DPNI_REFERENCE,
    DPNI_SECONDS,
    DPNI_UNITS,
    HEAT_SHOCK_CELSIUS,
    HEAT_SHOCK_SECONDS,
    ICE_SECONDS,
    IPTG_UM,
    NEB_TRANSFORMATION,
    OUTGROWTH_CELSIUS,
    OUTGROWTH_SECONDS,
    OUTGROWTH_UL,
    PLATE_DILUTION,
    PLATE_REFERENCE,
    PLATE_UL,
    PRIMER_PLATE_EQUIPMENT,
    PRIMER_PLATE_STORAGE,
    RECOVER_SECONDS,
    SEQUENCING_TITLE,
    THAW_SECONDS,
    WORKING_PLATE_SINGLE_USE,
    XGAL_UG_ML,
    Transformation,
    badges,
    card,
    catalogued,
    cleanup_step,
    colony_pcr_step,
    dam_sites,
    dpni_step,
    enzyme_material,
    gel_step,
    listed,
    pcr_step,
    pcr_title,
    phenotype_sentences,
    primer_plate_protocol,
    primer_plate_steps,
    quantify_step,
    sequencing_step,
    transform_step,
)
from liulab_mbio.bench.validation import (
    COLONY_ALLOWANCE,
    COLONY_FLANK,
    CORRECT_CLONE,
    EMPTY_CLONE,
    JUNCTION_OFFSET,
    REVERSE_FLANK,
    REVERSED_CLONE,
    SANGER_ALLOWANCE,
    SANGER_FLANK,
    Clone,
    ColonyCheck,
    SangerRead,
    colony_pcr_check,
    sanger_primers,
)
from liulab_mbio.protocol.model import Reference

#: Every source the modules cite, in the order a protocol lists them.
REFERENCES: tuple[Reference, ...] = (*pcr.REFERENCES, *amounts.REFERENCES, *gels.REFERENCES)

__all__ = [
    "ASSEMBLY_UL",
    "CELLS_UL",
    "CLEAN_COLONY_CURVE",
    "COLONY_ALLOWANCE",
    "COLONY_FLANK",
    "COLONY_HOLD_CELSIUS",
    "COLONY_LYSIS_SECONDS",
    "COLONY_PCR_MASTER_MIX",
    "COLONY_PCR_TITLE",
    "COLONY_PCR_VOLUME_UL",
    "CORRECT_CLONE",
    "DAM_SITE",
    "DNA_VOLUME_UL",
    "DNTP_STOCK_MM",
    "DPNI_CELSIUS",
    "DPNI_REFERENCE",
    "DPNI_SECONDS",
    "DPNI_UNITS",
    "EMPTY_CLONE",
    "HEAT_SHOCK_CELSIUS",
    "HEAT_SHOCK_SECONDS",
    "ICE_SECONDS",
    "IPTG_UM",
    "JUNCTION_OFFSET",
    "LADDER_1_KB_PLUS",
    "LADDER_100_BP",
    "MEDIUM",
    "NEB_TRANSFORMATION",
    "OUTGROWTH_CELSIUS",
    "OUTGROWTH_SECONDS",
    "OUTGROWTH_UL",
    "PLATE_DILUTION",
    "PLATE_REFERENCE",
    "PLATE_UL",
    "PRIMER_PLATE_EQUIPMENT",
    "PRIMER_PLATE_STORAGE",
    "PRIMER_STOCK_UM",
    "RECOVER_SECONDS",
    "REFERENCES",
    "REVERSED_CLONE",
    "REVERSE_FLANK",
    "SANGER_ALLOWANCE",
    "SANGER_FLANK",
    "SELECTION",
    "SEQUENCING_TITLE",
    "SHEET_COLUMNS",
    "THAW_SECONDS",
    "WATER",
    "WORKING_PLATE_SINGLE_USE",
    "XGAL_UG_ML",
    "Amount",
    "Clone",
    "ColonyCheck",
    "Electroporation",
    "Item",
    "OrderedOligo",
    "Phenotype",
    "PriceRecord",
    "PriceRow",
    "PrimerPlates",
    "SangerRead",
    "Transformation",
    "WellVerdict",
    "agarose_percent",
    "badges",
    "bill",
    "card",
    "catalogued",
    "choose_ladder",
    "clean_colony_chance",
    "cleanup_step",
    "colony_pcr_check",
    "colony_pcr_master_mix_component",
    "colony_pcr_program",
    "colony_pcr_reaction",
    "colony_pcr_step",
    "compact",
    "dam_sites",
    "dna_amount",
    "dna_components",
    "dpni_step",
    "electroporation",
    "enzyme_material",
    "fits",
    "gel_step",
    "heat_inactivation",
    "identity_check",
    "listed",
    "material",
    "molecular_weight",
    "oligo_row",
    "ordered_row",
    "pcr_program",
    "pcr_reaction",
    "pcr_step",
    "pcr_title",
    "phenotype_sentences",
    "plate",
    "pool",
    "primer_plate_protocol",
    "primer_plate_steps",
    "primer_plates",
    "primer_sheet",
    "quantify_step",
    "reaction_table",
    "read_phenotype",
    "read_prices",
    "reformat",
    "sanger_primers",
    "seat",
    "selection_marker",
    "sequencing_step",
    "to_nanograms",
    "to_pmol",
    "transform_step",
    "wells_of",
]
