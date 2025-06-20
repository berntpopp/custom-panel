"""SNP panel scrapers for extracting rsIDs from various sources."""

from .base_snp_scraper import BaseSNPScraper
from .parse_ampliseq import AmpliseqParser
from .parse_eurogentest import EurogentestParser
from .parse_idt_amplicon import IDTAmpliconParser
from .parse_idt_hybridization import IDTHybridizationParser
from .parse_nimagen_hest import NimaGenHESTParser
from .parse_nimagen_hist import NimaGenHISTParser
from .parse_pengelly import PengellyParser

__all__ = [
    "BaseSNPScraper",
    "PengellyParser",
    "EurogentestParser",
    "IDTAmpliconParser",
    "IDTHybridizationParser",
    "AmpliseqParser",
    "NimaGenHESTParser",
    "NimaGenHISTParser",
]
