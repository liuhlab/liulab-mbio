"""The bench any cloning pipeline shares: reactions, programs, validation and protocol steps.

Every public name in its modules imports from here as well. A module that cites a source keeps
it as its own `REFERENCES`; `REFERENCES` here gathers them.

- `amounts`: the weight of DNA, picomoles from nanograms, and what to pipette.
- `reactions`: one reaction table, filled to volume, and what it refuses.
- `pcr`: the PCR and colony PCR reactions and programs.
- `gels`: the ladder and agarose percentage a range of bands takes.
- `validation`: the colony PCR that reads an assembly's junctions, and the Sanger reads.
- `inactivation`: an enzyme's heat inactivation.
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
from liulab_mbio.bench.oligos import SHEET_COLUMNS, oligo_row, primer_sheet
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
from liulab_mbio.bench.phenotype import SELECTION, Phenotype, read_phenotype
from liulab_mbio.bench.reactions import WATER, dna_components, fits, reaction_table
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
    RECOVER_SECONDS,
    SEQUENCING_TITLE,
    THAW_SECONDS,
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
    "NEB_TRANSFORMATION",
    "OUTGROWTH_CELSIUS",
    "OUTGROWTH_SECONDS",
    "OUTGROWTH_UL",
    "PLATE_DILUTION",
    "PLATE_REFERENCE",
    "PLATE_UL",
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
    "XGAL_UG_ML",
    "Amount",
    "Clone",
    "ColonyCheck",
    "Phenotype",
    "SangerRead",
    "Transformation",
    "agarose_percent",
    "badges",
    "card",
    "catalogued",
    "choose_ladder",
    "cleanup_step",
    "colony_pcr_check",
    "colony_pcr_master_mix_component",
    "colony_pcr_program",
    "colony_pcr_reaction",
    "colony_pcr_step",
    "dam_sites",
    "dna_amount",
    "dna_components",
    "dpni_step",
    "enzyme_material",
    "fits",
    "gel_step",
    "heat_inactivation",
    "listed",
    "molecular_weight",
    "oligo_row",
    "pcr_program",
    "pcr_reaction",
    "pcr_step",
    "pcr_title",
    "phenotype_sentences",
    "primer_sheet",
    "quantify_step",
    "reaction_table",
    "read_phenotype",
    "sanger_primers",
    "sequencing_step",
    "to_nanograms",
    "to_pmol",
    "transform_step",
]
