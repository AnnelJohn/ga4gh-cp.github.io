#!/usr/bin/env -S uv run --quiet --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "phenopackets>=2.0.2.post5",
#     "protobuf"
# ]
# ///
"""
Create phenopacket Family objects from example variant CSV files.
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any

from google.protobuf.json_format import MessageToJson
from phenopackets import *

# Constants
METADATA = MetaData(
    phenopacket_schema_version='2.0',
    created=datetime.now(),
    created_by="csv_to_phenopackets.py",
    resources=[
        Resource(
            id="hp",
            name="Human Phenotype Ontology",
            url="http://purl.obolibrary.org/obo/hp.owl",
            version="2026-01-08",
            namespace_prefix="HP",
            iri_prefix="http://purl.obolibrary.org/obo/HP_",
        ),
        Resource(
            id="mondo",
            name="Mondo Disease Ontology",
            url="http://purl.obolibrary.org/obo/mondo/mondo-international.owl",
            version="2026-02-03",
            namespace_prefix="MONDO",
            iri_prefix="http://purl.obolibrary.org/obo/MONDO_"
        ),
        Resource(
            id="geno",
            name="Genotype Ontology",
            url="http://purl.obolibrary.org/obo/geno.owl",
            version="2025-07-25",
            namespace_prefix="GENO",
            iri_prefix="http://purl.obolibrary.org/obo/GENO_"
        )
    ]
)


# Data structures
@dataclass
class FamilyMember:
    """Holds phenopacket and metadata for a family member"""
    family_id: str
    individual_id: str
    role: str  # FATHER, MOTHER, PROBAND, SIBLING
    sex: int  # Sex enum value
    affected: str  # AFFECTED, UNAFFECTED, MISSING
    phenopacket: Phenopacket


# Helper functions
def parse_hpo_terms(hpo_string: str) -> List[PhenotypicFeature]:
    """
    Parse HPO terms from string.
    Format: HP:0001250 (Seizure); HP:0001263 (Global developmental delay)
    """
    phenotypic_features = []

    if not hpo_string:
        return phenotypic_features

    hpo_terms = [term.strip() for term in re.split(r'[;|]', hpo_string) if term.strip()]
    for hpo_term in hpo_terms:
        match = re.match(r'(HP:\d+)\s*\(([^)]+)\)', hpo_term)
        if match:
            hpo_id = match.group(1)
            hpo_label = match.group(2)
            phenotypic_features.append(
                PhenotypicFeature(type=OntologyClass(id=hpo_id, label=hpo_label))
            )

    return phenotypic_features


def build_variation_descriptor(chrom: str, pos: str, ref: str, alt: str,
                               gene: str, transcript: str, hgvsc: str,
                               hgvsp: str, zygosity: str) -> VariationDescriptor:
    """Build a GenomicInterpretation from variant data."""
    return VariationDescriptor(
        id=f'{chrom}-{pos}-{ref}-{alt}',
        gene_context=GeneDescriptor(value_id='', symbol=gene),
        expressions=[
            Expression(syntax='hgvs', value=f'{transcript}:{hgvsc}'),
            Expression(syntax='hgvs', value=hgvsp)
        ],
        vcf_record=VcfRecord(
            genome_assembly='GRCh38',
            chrom=chrom,
            pos=int(pos),
            ref=ref,
            alt=alt
        ),
        allelic_state=OntologyClass(
            id='GENO:0000136',
            label='homozygous'
        ) if zygosity == 'Homozygous' else OntologyClass(
            id='GENO:0000135',
            label='heterozygous'
        )
    )


def get_interpretation_status(affected_status: str, num_variants: int) -> str:
    if affected_status == 'AFFECTED':
        return 'CONTRIBUTORY' if num_variants == 2 else 'CAUSATIVE'
    return 'CANDIDATE'


def parse_variants_from_row(row: Dict[str, Any]) -> List[VariationDescriptor]:
    """
    Parse variant data from row and build Interpretation.
    Returns None if no variant data present.
    """
    variant_descriptors = []
    # Check for first variant
    if row.get('Chrom-1') and row['Chrom-1'].strip():
        vd1 = build_variation_descriptor(
            chrom=row['Chrom-1'],
            pos=row['Pos-1'],
            ref=row['Ref-1'],
            alt=row['Alt-1'],
            gene=row['Gene'],
            transcript=row['Transcript'],
            hgvsc=row['HGVSc-1'],
            hgvsp=row['HGVSp-1'],
            zygosity=row['Zygosity']
        )
        variant_descriptors.append(vd1)

    # Check for second variant (compound het)
    if row.get('Chrom-2') and row['Chrom-2'].strip():
        vd2 = build_variation_descriptor(
            chrom=row['Chrom-2'],
            pos=row['Pos-2'],
            ref=row['Ref-2'],
            alt=row['Alt-2'],
            gene=row['Gene'],
            transcript=row['Transcript'],
            hgvsc=row['HGVSc-2'],
            hgvsp=row['HGVSp-2'],
            zygosity=row['Zygosity']
        )
        variant_descriptors.append(vd2)

    return variant_descriptors


def parse_variant_interpretations(row: Dict[str, Any], individual_id: str, affected_status: str) -> List[Interpretation]:
    variation_descriptors = parse_variants_from_row(row)
    if len(variation_descriptors) == 0:
        return []

    interpretation_status = get_interpretation_status(affected_status, len(variation_descriptors))

    genomic_interpretations = []
    for variation_descriptor in variation_descriptors:
        genomic_interpretation = GenomicInterpretation(
            subject_or_biosample_id=individual_id,  # Will be set when we know the individual_id
            interpretation_status=GenomicInterpretation.InterpretationStatus.Value(interpretation_status),
            variant_interpretation=VariantInterpretation(variation_descriptor=variation_descriptor)
        )
        genomic_interpretations.append(genomic_interpretation)

    interpretation = Interpretation(
        id=f'{individual_id}',
        progress_status=Interpretation.ProgressStatus.Value('SOLVED' if affected_status == 'AFFECTED' else 'COMPLETED'),
        diagnosis=Diagnosis(
            # Strictly the diagnosis requires a disease but the input doesn't contain one, so we're using a stub here
            disease=OntologyClass(id='MONDO:0003847', label='hereditary disease'),
            genomic_interpretations=genomic_interpretations
        )
    )

    return [interpretation]


def parse_row_to_family_member(row: Dict[str, Any], id_key: str) -> FamilyMember:
    """Parse a CSV row directly into a FamilyMember with Phenopacket."""
    # Extract basic info
    individual_id = row[id_key]

    # Parse family ID and role
    match = re.match(r'^(.+)_(MOTHER|FATHER|PROBAND|SIBLING)', individual_id)
    if not match:
        raise ValueError(f"Invalid ID format: {individual_id}")

    family_id = match.group(1)
    role = match.group(2)

    # Map sex and affected status
    sex = row.get('SEX', 'UNKNOWN_SEX').upper()
    affected = row.get('AFFECTED STATUS', 'MISSING').upper()

    # Build phenotypic features
    hpo_string = row.get('HPO PRESENT', '')
    phenotypic_features = parse_hpo_terms(hpo_string)

    # Build time element (age of onset)
    time_element = None
    age_onset = row.get('AGE OF ONSET (yrs)', '')
    if age_onset:
        time_element = TimeElement(age=Age(iso8601duration=f"P{age_onset}Y"))

    # Build interpretations
    interpretations = parse_variant_interpretations(row, individual_id, affected)

    # Build phenopacket
    phenopacket = Phenopacket(
        id=f"{family_id}_{role}",
        subject=Individual(
            id=individual_id,
            sex=Sex.Value(sex),
            time_at_last_encounter=time_element
        ),
        phenotypic_features=phenotypic_features,
        interpretations=interpretations,
        meta_data=METADATA if len(phenotypic_features) > 0  else MetaData(
            phenopacket_schema_version='2.0',
            created=datetime.now(),
            created_by="csv_to_phenopackets.py"
        )
    )

    return FamilyMember(
        family_id=family_id,
        individual_id=individual_id,
        role=role,
        sex=sex,
        affected='MISSING' if affected == 'UNKNOWN' else affected,
        phenopacket=phenopacket
    )


def read_family_members_from_csv(input_file: Path) -> List[FamilyMember]:
    """Read CSV and convert each row to a FamilyMember with Phenopacket."""
    family_members = []

    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)

        # ID column
        id_key = 'ID'

        for row in reader:
            try:
                member = parse_row_to_family_member(row, id_key)
                family_members.append(member)
            except ValueError as e:
                print(f"Warning: Skipping row - {e}", file=sys.stderr)
                continue

    return family_members


def build_pedigree(family_id: str, members: List[FamilyMember]) -> Pedigree:
    """Build pedigree from family members."""
    # Find parent IDs
    father_id = next((m.individual_id for m in members if m.role == 'FATHER'), '0')
    mother_id = next((m.individual_id for m in members if m.role == 'MOTHER'), '0')

    pedigree_persons = []
    for member in members:
        if member.role in ['FATHER', 'MOTHER']:
            paternal_id = '0'
            maternal_id = '0'
        else:
            paternal_id = father_id
            maternal_id = mother_id

        pedigree_persons.append(Pedigree.Person(
            family_id=family_id,
            individual_id=member.individual_id,
            paternal_id=paternal_id,
            maternal_id=maternal_id,
            sex=member.sex,
            affected_status=Pedigree.Person.AffectedStatus.Value(member.affected)
        ))

    return Pedigree(persons=pedigree_persons)


def build_family_from_members(family_id: str, members: List[FamilyMember]) -> Family:
    """Build a Family phenopacket from FamilyMember objects."""
    # Separate proband and relatives
    proband = next((m.phenopacket for m in members if m.role == 'PROBAND'), None)
    relatives = [m.phenopacket for m in members if m.role != 'PROBAND']

    # Build pedigree
    pedigree = build_pedigree(family_id, members)

    return Family(
        id=family_id,
        proband=proband,
        relatives=relatives,
        pedigree=pedigree,
        meta_data=METADATA
    )


def write_family_to_file(family: Family, output_directory: str):
    """Write Family phenopacket to JSON file."""
    output_path = os.path.join(output_directory, f'{family.id}_PROBAND.json')
    with open(output_path, 'w') as f:
        family_json = MessageToJson(family)
        f.write(family_json)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Create phenopacket Family objects from variant CSV files.'
    )
    parser.add_argument('input_file', help='Input CSV file', type=Path)
    parser.add_argument('output_directory', help='Output directory for JSON files', type=Path)

    args = parser.parse_args()

    # Validate input
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output_directory, exist_ok=True)

    # Read and parse CSV rows into phenopackets
    print(f"Reading samples from {args.input_file}...")
    family_members = read_family_members_from_csv(args.input_file)
    print(f"Parsed {len(family_members)} family members")

    # Group by family
    families = defaultdict(list)
    for member in family_members:
        families[member.family_id].append(member)

    print(f"Found {len(families)} families")

    # Create and write family phenopackets
    print("Creating phenopacket files...")
    for family_id, members in families.items():
        family = build_family_from_members(family_id, members)
        write_family_to_file(family, args.output_directory)
        print(f"  Created {family_id}_PROBAND.json")

    print(f"\nCompleted! {len(families)} files created in {args.output_directory}/")


if __name__ == '__main__':
    main()
