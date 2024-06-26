#!/usr/bin/env python

from typing import List

from sqlalchemy import Table, Column, Integer, String, ForeignKey, Float, Boolean
from sqlalchemy.orm import relationship, Mapped, mapped_column
from sqlalchemy.orm import declarative_base

from sqlalchemy.ext.declarative import DeferredReflection
from sqlalchemy.ext.declarative import declared_attr

Base = declarative_base(cls=DeferredReflection)

# helper table for many2many relationship gene <-> phenotype
resfinder_phenotype_association_table = Table(
        "resfinder_sequence_phenotype",
        Base.metadata,
        Column("id", Integer(), primary_key=True, autoincrement=True),
        Column("sequence_id", ForeignKey("resfinder_sequence.id", ondelete="CASCADE"), nullable=False),
        Column("phenotype_id", ForeignKey("phenotype.id", ondelete="CASCADE"), nullable=False),
)

amrfinder_phenotype_association_table = Table(
        "amrfinder_sequence_phenotype",
        Base.metadata,
        Column("id", Integer(), primary_key=True, autoincrement=True),
        Column("sequence_id", ForeignKey("amrfinder_sequence.id", ondelete="CASCADE"), nullable=False),
        Column("phenotype_id", ForeignKey("phenotype.id", ondelete="CASCADE"), nullable=False),
)

# helper table for many2many relationship point mutation <-> phenotype
pointfinder_phenotype_association_table = Table(
        "resfinder_point_phenotype",
        Base.metadata,
        Column("id", Integer(), primary_key=True, autoincrement=True),
        Column("pointfinder_result_id", ForeignKey("pointfinder_result.id", ondelete="CASCADE"), nullable=False),
        Column("phenotype_id", ForeignKey("phenotype.id", ondelete="CASCADE"), nullable=False),
)

amrfinder_point_phenotype_association_table = Table(
        "amrfinder_point_phenotype",
        Base.metadata,
        Column("amrfinder_point_result_id", ForeignKey("amrfinder_point_result.id", ondelete="CASCADE"), nullable=False),
        Column("phenotype_id", ForeignKey("phenotype.id", ondelete="CASCADE"), nullable=False),
)


class ResfinderSequence(Base):
    __tablename__ = "resfinder_sequence"

    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String(200), nullable=False)
    short_name: Mapped[str] = mapped_column(String(50), nullable=True)
    subseq_numbering: Mapped[str] = mapped_column(String(50), nullable=True)
    accession: Mapped[str] = mapped_column(String(50), nullable=False)
    main_numbering: Mapped[str] = mapped_column(String(50), nullable=True)
    internal_numbering: Mapped[str] = mapped_column(String(100), nullable=True)
    crc32_hash: Mapped[str] = mapped_column(String(10), nullable=False)
    sequence: Mapped[str] = mapped_column(String(10000), nullable=False)

    # Relationships
    results: Mapped[List["ResfinderResult"]] = relationship(back_populates="stored_sequence")
    phenotypes: Mapped[List["Phenotype"]] = relationship(
        secondary=resfinder_phenotype_association_table, back_populates="resfinder_sequences",
    )


class AmrfinderSequence(Base):
    __tablename__ = "amrfinder_sequence"
    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String(200), nullable=False)
    short_name: Mapped[str] = mapped_column(String(50), nullable=True)
    accession: Mapped[str] = mapped_column(String(50), nullable=False)
    activity_type: Mapped[str] = mapped_column(String(50), nullable=True)
    is_core: Mapped[bool] = mapped_column(Boolean(), nullable=True)
    long_name: Mapped[str] = mapped_column(String(1000), nullable=False)
    internal_numbering: Mapped[str] = mapped_column(String(100), nullable=True)
    crc32_hash: Mapped[str] = mapped_column(String(10), nullable=False)
    sequence: Mapped[str] = mapped_column(String(10000), nullable=False)

    # Relationships
    results: Mapped[List["AmrfinderResult"]] = relationship(back_populates="stored_sequence")
    phenotypes: Mapped[List["Phenotype"]] = relationship(
        secondary=amrfinder_phenotype_association_table, back_populates="amrfinder_sequences",
    )
    

class Phenotype(Base):
    __tablename__ = "phenotype"

    id: Mapped[int] = mapped_column(primary_key=True)
    resfinder_sequences: Mapped[List["ResfinderSequence"]] = relationship(
        secondary=resfinder_phenotype_association_table, back_populates="phenotypes",
    )
    amrfinder_sequences:  Mapped[List["AmrfinderSequence"]] = relationship(
        secondary=amrfinder_phenotype_association_table, back_populates="phenotypes",
    )
    pointmutations: Mapped[List["PointfinderResult"]] = relationship(
        secondary=pointfinder_phenotype_association_table,  back_populates="phenotypes",
    )
    amrfinder_pointmutations: Mapped[List["AmrfinderPointResult"]] = relationship(
        secondary=amrfinder_point_phenotype_association_table, back_populates="phenotypes",
    )
    class_name: Mapped[str] = mapped_column(String(50), nullable=True)
    phenotype: Mapped[str] = mapped_column(String(50), nullable=False)
    invitroresults: Mapped[List["InVitroResult"]] = relationship(back_populates="phenotype_associated")


class ResfinderResult(Base):
    __tablename__ = "resfinder_result"

    # Columns in Table (physically)
    id: Mapped[int] = mapped_column(primary_key=True)
    identity: Mapped[float] = mapped_column(Float(), nullable=False)
    coverage: Mapped[float] = mapped_column(Float(), nullable=False)
    ref_pos_start: Mapped[int] = mapped_column(Integer(), nullable=True)
    ref_pos_end: Mapped[int] = mapped_column(Integer(), nullable=True)
    qc_issues: Mapped[str] = mapped_column(String(1000), nullable=True)
    orientation: Mapped[str] = mapped_column(String(1), nullable=True)

    # Foreign Keys
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id", ondelete="CASCADE"), nullable=False)
    sequence_id:  Mapped[int] = mapped_column(ForeignKey("resfinder_sequence.id", ondelete="CASCADE"), nullable=False)
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id", ondelete="CASCADE"), nullable=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    sample_associated: Mapped["Sample"] = relationship(back_populates="resfinderresults")
    stored_sequence: Mapped[ResfinderSequence] = relationship(back_populates="results")
    contig_associated: Mapped["Contig"] = relationship(back_populates="resfinderresults")
    version_associated: Mapped["ToolVersion"] = relationship()


class Sample(Base):
    # the sample is the core connecting table which can receive a name for displaying purposes
    # external id can be optionally used if you have a separate db managing your samples
    # in this case the name can be left empty
    # sample may conntect results and contigs separately
    # e.g. sample phenotypes needs to be accessed via sample.resfinderresults.stored_sequence.phenotypes
    __tablename__ = "sample"

    id: Mapped[int] = mapped_column(primary_key=True)
    external_id: Mapped[int] = mapped_column(Integer(), nullable=True)
    name: Mapped[str] = mapped_column(String(500), nullable=True)

    # Relationships
    stored_contigs: Mapped[List["Contig"]] = relationship(back_populates="sample_associated")
    resfinderresults: Mapped[List[ResfinderResult]] = relationship(back_populates="sample_associated")
    pointfinderresults: Mapped[List["PointfinderResult"]] = relationship(back_populates="sample_associated")
    invitroresults: Mapped["InVitroResult"] = relationship(back_populates="sample_associated")
    speciesfinderresults: Mapped[List["SpeciesfinderResult"]] = relationship(back_populates="sample_associated")
    amrfinderresults: Mapped[List["AmrfinderResult"]] = relationship(back_populates="sample_associated")
    amrfinderpointresults: Mapped[List["AmrfinderPointResult"]] = relationship(back_populates="sample_associated")
    mlstresults: Mapped[List["MlstResult"]] = relationship(back_populates="sample_associated")


class Contig(Base):
    __tablename__ = "contig"

    id: Mapped[int] = mapped_column(primary_key=True)
    name:  Mapped[str] = mapped_column(String(500), nullable=False)
    length: Mapped[int] = mapped_column(Integer(), nullable=True)

    # Foreign keys
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id", ondelete="CASCADE"), nullable=False)
    
    # Relationships
    sample_associated: Mapped["Sample"] = relationship(back_populates="stored_contigs")
    resfinderresults: Mapped[List[ResfinderResult]] = relationship(back_populates="contig_associated")
    isescanresults: Mapped[List["ISEScanResult"]] = relationship(back_populates="contig_associated")
    baktaresults: Mapped[List["BaktaResult"]] = relationship(back_populates="contig_associated")
    mobtyperresults: Mapped["MobTyperResult"] = relationship(back_populates="contig_associated")
    plasmidfinderresults: Mapped[List["PlasmidfinderResult"]] = relationship(back_populates="contig_associated")
    phipsyresults:  Mapped[List["PhispyResults"]] = relationship(back_populates="contig_associated")
    amrfinderresults: Mapped[List["AmrfinderResult"]] = relationship(back_populates="contig_associated")
    amrfinderpointresults: Mapped[List["AmrfinderPointResult"]] = relationship(back_populates="contig_associated")


class PointfinderResult(Base):
    __tablename__ = "pointfinder_result"

    id: Mapped[int] = mapped_column(primary_key=True)
    mutation: Mapped[str] = mapped_column(String(50), nullable=False)
    nuc_change: Mapped[str] = mapped_column(String(50), nullable=False)

    # Foreign keys
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id", ondelete="CASCADE"), nullable=False)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    sample_associated: Mapped["Sample"] = relationship(back_populates="pointfinderresults")
    version_associated: Mapped["ToolVersion"] = relationship()
    phenotypes: Mapped[List["Phenotype"]] = relationship(
        secondary=pointfinder_phenotype_association_table, back_populates="pointmutations",
    )


class AmrfinderResultMixin(object):
    """
    Blueprint for AmrfinderResult-Table (splitted into Point and Sequence)
    not an actual table
    """

    # Columns in Table (physically)
    id: Mapped[int] = mapped_column(primary_key=True)
    identity: Mapped[float] = mapped_column(Float(), nullable=False)
    coverage: Mapped[float] = mapped_column(Float(), nullable=False)
    ref_pos_start: Mapped[int] = mapped_column(Integer(), nullable=True)
    ref_pos_end: Mapped[int] = mapped_column(Integer(), nullable=True)
    qc_issues: Mapped[str] = mapped_column(String(1000), nullable=True)
    orientation: Mapped[str] = mapped_column(String(1), nullable=True)
    method: Mapped[str] = mapped_column(String(30), nullable=True)

    # Foreign Keys
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id", ondelete="CASCADE"), nullable=False)
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id", ondelete="CASCADE"), nullable=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships (only declareable for non-back-populating)
    @declared_attr
    def version_associated(self) -> Mapped["ToolVersion"]:
        return relationship()


class AmrfinderResult(AmrfinderResultMixin, Base):
    __tablename__ = "amrfinder_result"

    # Foreign Keys
    sequence_id:  Mapped[int] = mapped_column(ForeignKey("amrfinder_sequence.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    stored_sequence: Mapped["AmrfinderSequence"] = relationship(back_populates="results")
    sample_associated: Mapped["Sample"] = relationship(back_populates="amrfinderresults")
    contig_associated: Mapped["Contig"] = relationship(back_populates="amrfinderresults")
    #version_associated: Mapped["ToolVersion"] = relationship()


class AmrfinderPointResult(AmrfinderResultMixin, Base):
    __tablename__ = "amrfinder_point_result"

    # Columns in Table (physically)
    mutation: Mapped[str] = mapped_column(String(100))

    # Relationships:
    phenotypes: Mapped[List["Phenotype"]] = relationship(
        secondary=amrfinder_point_phenotype_association_table, back_populates="amrfinder_pointmutations",
    )
    sample_associated: Mapped["Sample"] = relationship(back_populates="amrfinderpointresults")
    contig_associated: Mapped["Contig"] = relationship(back_populates="amrfinderpointresults")
    #version_associated: Mapped["ToolVersion"] = relationship()


class ISEScanResult(Base):
    __tablename__ = "isescan_result"

    id: Mapped[int] = mapped_column(primary_key=True)
    family: Mapped[str] = mapped_column(String(50), nullable=False)
    cluster: Mapped[str] = mapped_column(String(50), nullable=False)
    is_start_pos: Mapped[int] = mapped_column(Integer(), nullable=False)
    is_end_pos: Mapped[int] = mapped_column(Integer(), nullable=False)
    is_copy_number: Mapped[int] = mapped_column(Integer(), nullable=False)
    # inverted repeats coordinates
    ir_start_pos1: Mapped[int] = mapped_column(Integer(), nullable=True)
    ir_end_pos1: Mapped[int] = mapped_column(Integer(), nullable=False)
    ir_start_pos2: Mapped[int] = mapped_column(Integer(), nullable=True)
    ir_end_pos2: Mapped[int] = mapped_column(Integer(), nullable=True)
    # inverted repeats scoring
    score: Mapped[int] = mapped_column(Integer(), nullable=True)
    irId: Mapped[int] = mapped_column(Integer(), nullable=True)
    irLen: Mapped[int] = mapped_column(Integer(), nullable=True)
    nGaps: Mapped[int] = mapped_column(Integer(), nullable=True)
    # transposase orf
    orf_start_pos: Mapped[int] = mapped_column(Integer(), nullable=True)
    orf_end_pos: Mapped[int] = mapped_column(Integer(), nullable=True)
    orientation: Mapped[str] = mapped_column(String(1), nullable=True)
    e_value: Mapped[float] = mapped_column(Float(), nullable=True)
    complete: Mapped[bool] = mapped_column(Boolean(), nullable=True) # "type" (c or p)
    ov: Mapped[str] =  mapped_column(Integer(), nullable=True) # numeric value for counting
    tir: Mapped[str] = mapped_column(String(500), nullable=True) #terminal inverted repeat


    # Foreign keys:
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id", ondelete="CASCADE"), nullable=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships:
    contig_associated: Mapped["Contig"] = relationship(back_populates="isescanresults")
    version_associated: Mapped["ToolVersion"] = relationship()


class BaktaResult(Base):
    __tablename__ = "bakta_result"

    # Columns in Table (physically)
    id: Mapped[int] = mapped_column(primary_key=True)
    ref_pos_start: Mapped[int] = mapped_column(Integer(), nullable=False)
    ref_pos_end: Mapped[int] = mapped_column(Integer(), nullable=False)
    orientation: Mapped[str] = mapped_column(String(1), nullable=False)
    gene_name: Mapped[str] = mapped_column(String(20), nullable=True)
    product_type: Mapped[str] = mapped_column(String(20), nullable=True)
    product: Mapped[str] = mapped_column(String(1000), nullable=True)
    db_xref: Mapped[str] = mapped_column(String(1000), nullable=True)

    # Foreign Keys
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id", ondelete="CASCADE"), nullable=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    contig_associated: Mapped["Contig"] = relationship(back_populates="baktaresults")
    version_associated: Mapped["ToolVersion"] = relationship()


class MobTyperResult(Base):
    __tablename__ = "mobtyper_result"
    # Columns in Table (physically)
    id: Mapped[int] = mapped_column(primary_key=True)
    gc_content: Mapped[float] = mapped_column(Float(), nullable=True)
    rep_type: Mapped[str] = mapped_column(String(100), nullable=True)
    rep_type_accession: Mapped[str] = mapped_column(String(1000), nullable=True)

    relaxase_type: Mapped[str] = mapped_column(String(100), nullable=True)
    relaxase_type_accession: Mapped[str] = mapped_column(String(1000), nullable=True)

    mpf_type: Mapped[str] = mapped_column(String(100), nullable=True)
    mpf_type_accession:  Mapped[str] = mapped_column(String(1000), nullable=True)

    orit_type: Mapped[str] = mapped_column(String(100), nullable=True)
    orit_accession:  Mapped[str] = mapped_column(String(1000), nullable=True)

    predicted_mobility: Mapped[str] = mapped_column(String(100), nullable=True)
    mash_nearest_neighbor: Mapped[str] = mapped_column(String(100), nullable=True)
    mash_neighbor_distance: Mapped[int] = mapped_column(Integer(), nullable=True)
    mash_neighbor_identification: Mapped[str] = mapped_column(String(200), nullable=True)
    primary_cluster_id: Mapped[str] = mapped_column(String(100), nullable=True)
    secondary_cluster_id: Mapped[str] = mapped_column(String(100), nullable=True)

    predicted_host_range_overall_rank: Mapped[str] = mapped_column(String(100), nullable=True)
    predicted_host_range_overall_name: Mapped[str] = mapped_column(String(100), nullable=True)
    observed_host_range_ncbi_rank: Mapped[str] = mapped_column(String(100), nullable=True)
    observed_host_range_ncbi_name: Mapped[str] = mapped_column(String(100), nullable=True)
    reported_host_range_lit_rank: Mapped[str] = mapped_column(String(100), nullable=True)
    reported_host_range_lit_name: Mapped[str] = mapped_column(String(100), nullable=True)
    associated_pmid: Mapped[str] = mapped_column(String(100), nullable=True)

    # Foreign Keys
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id", ondelete="CASCADE"), nullable=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    contig_associated: Mapped["Contig"] = relationship(back_populates="mobtyperresults")
    version_associated: Mapped["ToolVersion"] = relationship()


class InVitroResult(Base):
    # result for minimal inhibition concentration
    __tablename__ = "invitro_result"
    id: Mapped[int] = mapped_column(primary_key=True)
    mic: Mapped[float] = mapped_column(Float(), nullable=True)
    category: Mapped[str] = mapped_column(String(10), nullable=True)

    # Foreign keys
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id", ondelete="CASCADE"), nullable=False)
    phenotype_id: Mapped[int] = mapped_column(ForeignKey("phenotype.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    sample_associated: Mapped["Sample"] = relationship(back_populates="invitroresults") 
    phenotype_associated: Mapped["Phenotype"] = relationship(back_populates="invitroresults")


class PlasmidfinderResult(Base):
    __tablename__ = "plasmidfinder_result"
    #plasmidfinder output columns are:
    #Database, Plasmid, Identity, Query / Template length, Contig, Position in contig, Note, Accession number
    id: Mapped[int] = mapped_column(primary_key=True)
    database_name: Mapped[str] = mapped_column(String(100), nullable=True)
    plasmid: Mapped[str] = mapped_column(String(100), nullable=True)
    identity: Mapped[float] = mapped_column(Float(), nullable=False) 
    query_length: Mapped[float] = mapped_column(Float(), nullable=False) 
    template_length: Mapped[float] = mapped_column(Float(), nullable=False) 
    ref_pos_start: Mapped[int] = mapped_column(Integer(), nullable=True)
    ref_pos_end: Mapped[int] = mapped_column(Integer(), nullable=True)
    note: Mapped[str] = mapped_column(String(1000), nullable=True)
    accession_number: Mapped[str] = mapped_column(String(1000), nullable=True)

    # Foreign Keys
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id", ondelete="CASCADE"), nullable=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    contig_associated: Mapped["Contig"] = relationship(back_populates="plasmidfinderresults")
    version_associated: Mapped["ToolVersion"] = relationship()
    
    
class PhispyResults(Base):
    __tablename__ = "phispy_result"
    #phispy output columns are:
    # Prophage number; The contig upon which the prophage resides; The start location of the prophage; The stop location of the prophage If we can detect the att sites; 
    #the additional columns are: 
    # start of attL; end of attL; start of attR; end of attR; sequence of attL; sequence of attR; The explanation of why this att site was chosen for this prophage.
    
    id: Mapped[int] = mapped_column(primary_key=True)
    prophage_number: Mapped[str] = mapped_column(String(20), nullable=True)
    start: Mapped[int] = mapped_column(Integer(), nullable=True)
    stop: Mapped[int] = mapped_column(Integer(), nullable=True)
    start_attL: Mapped[int] = mapped_column(Integer(), nullable=True)
    end_attL: Mapped[int] = mapped_column(Integer(), nullable=True)
    start_attR: Mapped[int] = mapped_column(Integer(), nullable=True)
    end_attR: Mapped[int] = mapped_column(Integer(), nullable=True)
    sequence_attL: Mapped[str] = mapped_column(String(200), nullable=True)
    sequence_attR: Mapped[str] = mapped_column(String(200), nullable=True)
    description: Mapped[str] = mapped_column(String(200), nullable=True)
        
    
    # Foreign Keys
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id", ondelete="CASCADE"), nullable=True)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    contig_associated: Mapped["Contig"] = relationship(back_populates="phipsyresults")
    version_associated: Mapped["ToolVersion"] = relationship()


class SpeciesfinderResult(Base):
    __tablename__ = "speciesfinder_result"
    id: Mapped[int] = mapped_column(primary_key=True)
    template: Mapped[str] = mapped_column(String(500), nullable=True)
    species: Mapped[str] = mapped_column(String(200), nullable=True)
    match_id: Mapped[str] = mapped_column(String(200), nullable=True)
    database_name: Mapped[str] = mapped_column(String(200), nullable=True)
    confidence_of_result: Mapped[str] = mapped_column(String(20), nullable=True)
    file_format: Mapped[str] = mapped_column(String(20), nullable=True)
    method: Mapped[str] = mapped_column(String(20), nullable=True)
    
    # Foreign Keys
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id", ondelete="CASCADE"), nullable=False)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    sample_associated: Mapped["Sample"] = relationship(back_populates="speciesfinderresults")
    version_associated: Mapped["ToolVersion"] = relationship()


class MlstResult(Base):
    __tablename__ = "mlst_result"
    id: Mapped[int] = mapped_column(primary_key=True)
    scheme_name: Mapped[str] = mapped_column(String(100), nullable=True)
    sequence_type: Mapped[int] = mapped_column(Integer(), nullable=True)
    allele_types: Mapped[str] = mapped_column(String(500), nullable=True)

    # Foreign Keys
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id", ondelete="CASCADE"), nullable=False)
    version_id: Mapped[int] = mapped_column(ForeignKey("tool_version.id", ondelete="CASCADE"), nullable=False)

    # Relationships
    sample_associated: Mapped["Sample"] = relationship(back_populates="mlstresults")
    version_associated: Mapped["ToolVersion"] = relationship()


class ToolVersion(Base):
    # Table which contains information about the version of db and tool
    __tablename__ = "tool_version"
    id: Mapped[int] = mapped_column(primary_key=True)
    tool_name: Mapped[str] = mapped_column(String(100), nullable=True)
    tool_version: Mapped[str] = mapped_column(String(100), nullable=True)
    input_type: Mapped[str] = mapped_column(String(10), nullable=True)
    db_version: Mapped[str] = mapped_column(String(100), nullable=True)

    # Relationships:
    # TODO write all associated back-populations, this does not seem to
    # be neccessary at this point
