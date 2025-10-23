from pydantic import BaseModel, Field

NON_HOM_CHROM = ['X', 'Y', 'MT', 'M']
CHROM_ORDER = list(map(str, range(1, 23))) + NON_HOM_CHROM


class Coordinates(BaseModel):
    """
    A representation of genomic coordinates
    """

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def string_format(self) -> str:
        """
        forms a string representation: chr-pos-ref-alt
        """
        return f'{self.chrom}-{self.pos}-{self.ref}-{self.alt}'

    def __lt__(self, other) -> bool:
        """
        enables positional sorting
        """
        # this will return False for same chrom and position
        if self.chrom == other.chrom:
            return self.pos < other.pos
        # otherwise take the relative index from sorted chromosomes list
        if self.chrom in CHROM_ORDER and other.chrom in CHROM_ORDER:
            return CHROM_ORDER.index(self.chrom) < CHROM_ORDER.index(other.chrom)
        # if self is on a canonical chromosome, sort before HLA/Decoy etc.
        return self.chrom in CHROM_ORDER


class VariantAf(BaseModel):
    coordinates: Coordinates = Field(repr=True)
    info: dict[str, str | int | float] = Field(default_factory=dict)
    het_samples: set[str] = Field(default_factory=set, exclude=True)
    hom_samples: set[str] = Field(default_factory=set, exclude=True)
    phased: dict = Field(default_factory=dict, exclude=True)
    transcript_consequences: list[dict[str, str | float | int]]

    def __str__(self):
        return repr(self)

    def __lt__(self, other):
        return self.coordinates < other.coordinates

    def __eq__(self, other):
        return self.coordinates == other.coordinates

    def __hash__(self):
        return hash(self.coordinates)


class ReportableVariant(BaseModel):
    """
    A variant passing MOI tests, to be reported
    """

    sample_id: str
    var_data: VariantAf
    reasons: set[str] = Field(default_factory=set)
    gene: str = Field(default_factory=str)
    support_vars: set[str] = Field(default_factory=set)

    def __eq__(self, other):
        """
        makes reported variants comparable
        """
        return self.sample_id == other.sample_id and self.var_data.coordinates == other.var_data.coordinates

    def __hash__(self):
        return hash((self.sample_id, self.var_data.coordinates))

    def __lt__(self, other):
        return self.var_data.coordinates < other.var_data.coordinates


class ResultsAf(BaseModel):
    """
    A representation of a result set
    """

    variants: dict[str, list[ReportableVariant]] = Field(default_factory=dict)
