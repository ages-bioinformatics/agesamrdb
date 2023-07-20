#!/usr/bin/env python

from typing import List

from sqlalchemy import Table, Column, Integer, String, ForeignKey, Float, Boolean
from sqlalchemy.orm import relationship, Mapped, mapped_column
from sqlalchemy.orm import declarative_base

from sqlalchemy.ext.declarative import DeferredReflection

Base = declarative_base(cls=DeferredReflection)

phenotype_association_table = Table(
        "sequencephenotype",
        Base.metadata,
        Column("id", Integer(), primary_key=True, autoincrement=True),
        Column("sequence_id", ForeignKey("resfinder_sequence.id"), nullable=False),
        Column("phenotype_id", ForeignKey("phenotype.id"), nullable=False),
)

class ResfinderSequence(Base):
        __tablename__ = "resfinder_sequence"

        id: Mapped[int] = mapped_column(primary_key=True)
        phenotypes: Mapped[List["Phenotype"]] = relationship(
            secondary=phenotype_association_table, back_populates="sequences",
        )
        results: Mapped[List["ResfinderResult"]] = relationship(back_populates="stored_sequence")
        sequence: Mapped[str] = mapped_column(String(10000), nullable=False)
        name: Mapped[str] = mapped_column(String(200), nullable=False)
        short_name: Mapped[str] = mapped_column(String(50), nullable=True)
        subseq_numbering: Mapped[str] = mapped_column(String(50), nullable=True)
        accession: Mapped[str] = mapped_column(String(50), nullable=False)
        main_numbering: Mapped[str] = mapped_column(String(50), nullable=True)
        internal_numbering: Mapped[str] = mapped_column(String(50), nullable=True)
        crc32_hash: Mapped[str] = mapped_column(String(10), nullable=False)


class Phenotype(Base):
        __tablename__ = "phenotype"

        id: Mapped[int] = mapped_column(primary_key=True)
        sequences: Mapped[List[ResfinderSequence]] = relationship(
            secondary=phenotype_association_table, back_populates="phenotypes",
        )
        class_name: Mapped[str] = mapped_column(String(50), nullable=True)
        phenotype: Mapped[str] = mapped_column(String(50), nullable=False)


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
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id"), nullable=False)
    sequence_id:  Mapped[int] = mapped_column(ForeignKey("resfinder_sequence.id"), nullable=False)
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id"), nullable=True)

    # Relationships
    sample_associated: Mapped["Sample"] = relationship(back_populates="resfinderresults")
    stored_sequence: Mapped[ResfinderSequence] = relationship(back_populates="results")
    contig_associated: Mapped["Contig"] = relationship(back_populates="resfinderresults")


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


class Contig(Base):
    __tablename__ = "contig"

    id: Mapped[int] = mapped_column(primary_key=True)
    name:  Mapped[str] = mapped_column(String(500), nullable=False)
    length: Mapped[int] = mapped_column(Integer(), nullable=True)

    # Foreign keys
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id"), nullable=False)
    
    # Relationships
    sample_associated: Mapped["Sample"] = relationship(back_populates="stored_contigs")
    resfinderresults: Mapped[List[ResfinderResult]] = relationship(back_populates="contig_associated")
    isescanresults: Mapped[List["ISEScanResult"]] = relationship(back_populates="contig_associated")
    baktaresults: Mapped[List["BaktaResult"]] = relationship(back_populates="contig_associated")


class PointfinderResult(Base):
    __tablename__ = "pointfinder_result"

    id: Mapped[int] = mapped_column(primary_key=True)
    mutation: Mapped[str] = mapped_column(String(50), nullable=False)
    nuc_change: Mapped[str] = mapped_column(String(50), nullable=False)
    phenotype: Mapped[str] = mapped_column(String(50), nullable=False)
    input_type: Mapped[str] = mapped_column(String(10), nullable=False)

    # Foreign keys
    sample_id: Mapped[int] = mapped_column(ForeignKey("sample.id"), nullable=False)

    # Relationships
    sample_associated: Mapped["Sample"] = relationship(back_populates="pointfinderresults")


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
    tir: Mapped[str] = mapped_column(String(100), nullable=True) #terminal inverted repeat


    # Foreign keys:
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id"), nullable=True)

    # Relationships:
    contig_associated: Mapped["Contig"] = relationship(back_populates="isescanresults")


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
    contig_id: Mapped[int] = mapped_column(ForeignKey("contig.id"), nullable=True)

    # Relationships
    contig_associated: Mapped["Contig"] = relationship(back_populates="baktaresults")


