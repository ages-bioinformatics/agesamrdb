#!/usr/bin/env python

from typing import List

from sqlalchemy import Table, Column, Integer, String, ForeignKey, Float
from sqlalchemy.orm import relationship, Mapped, mapped_column
from sqlalchemy.orm import declarative_base

from sqlalchemy.ext.declarative import DeferredReflection

Base = declarative_base(cls=DeferredReflection)

phenotype_association_table = Table(
        "sequencephenotype",
        Base.metadata,
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
