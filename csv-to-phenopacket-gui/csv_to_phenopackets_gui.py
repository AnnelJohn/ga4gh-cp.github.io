#!/usr/bin/env python3
"""
Professional GUI wrapper for CSV → Phenopackets Family conversion
Designed for GA4GH / standards community sharing
"""

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import sys
from pathlib import Path
from datetime import datetime
import re
import csv
from collections import defaultdict
from google.protobuf.json_format import MessageToJson

# Phenopackets v2 import (official style)
import phenopackets.schema.v2 as pps2

# ────────────────────────────────────────────────
#  Constants & Helper functions
# ────────────────────────────────────────────────

METADATA = pps2.MetaData()
METADATA.phenopacket_schema_version = "2.0"
METADATA.created.FromDatetime(datetime.now())
METADATA.created_by = "csv_to_phenopacket_gui.py"

hp = pps2.Resource()
hp.id = "hp"
hp.name = "Human Phenotype Ontology"
hp.url = "http://purl.obolibrary.org/obo/hp.owl"
hp.version = "2026-01-08"
hp.namespace_prefix = "HP"
hp.iri_prefix = "http://purl.obolibrary.org/obo/HP_"

mondo = pps2.Resource()
mondo.id = "mondo"
mondo.name = "Mondo Disease Ontology"
mondo.url = "http://purl.obolibrary.org/obo/mondo/mondo-international.owl"
mondo.version = "2026-02-03"
mondo.namespace_prefix = "MONDO"
mondo.iri_prefix = "http://purl.obolibrary.org/obo/MONDO_"

geno = pps2.Resource()
geno.id = "geno"
geno.name = "Genotype Ontology"
geno.url = "http://purl.obolibrary.org/obo/geno.owl"
geno.version = "2025-07-25"
geno.namespace_prefix = "GENO"
geno.iri_prefix = "http://purl.obolibrary.org/obo/GENO_"

METADATA.resources.extend([hp, mondo, geno])

class FamilyMember:
    def __init__(self, family_id, individual_id, role, sex, affected, phenopacket):
        self.family_id = family_id
        self.individual_id = individual_id
        self.role = role
        self.sex = sex
        self.affected = affected
        self.phenopacket = phenopacket

def parse_hpo_terms(hpo_string: str):
    phenotypic_features = []
    if not hpo_string:
        return phenotypic_features
    hpo_terms = [term.strip() for term in re.split(r'[;|]', hpo_string) if term.strip()]
    for hpo_term in hpo_terms:
        match = re.match(r'(HP:\d+)\s*\(([^)]+)\)', hpo_term)
        if match:
            hpo_id, hpo_label = match.groups()
            oc = pps2.OntologyClass()
            oc.id = hpo_id
            oc.label = hpo_label
            pf = pps2.PhenotypicFeature()
            pf.type.CopyFrom(oc)
            phenotypic_features.append(pf)
    return phenotypic_features

def build_variation_descriptor(chrom, pos, ref, alt, gene, transcript, hgvsc, hgvsp, zygosity):
    vd = pps2.VariationDescriptor()
    vd.id = f'{chrom}-{pos}-{ref}-{alt}'
    gene_desc = pps2.GeneDescriptor()
    gene_desc.symbol = gene
    vd.gene_context.CopyFrom(gene_desc)
    expr1 = pps2.Expression()
    expr1.syntax = "hgvs"
    expr1.value = f'{transcript}:{hgvsc}'
    vd.expressions.append(expr1)
    expr2 = pps2.Expression()
    expr2.syntax = "hgvs"
    expr2.value = hgvsp
    vd.expressions.append(expr2)
    vcf = pps2.VcfRecord()
    vcf.genome_assembly = "GRCh38"
    vcf.chrom = chrom
    vcf.pos = int(pos)
    vcf.ref = ref
    vcf.alt = alt
    vd.vcf_record.CopyFrom(vcf)
    allelic = pps2.OntologyClass()
    if zygosity == 'Homozygous':
        allelic.id = "GENO:0000136"
        allelic.label = "homozygous"
    else:
        allelic.id = "GENO:0000135"
        allelic.label = "heterozygous"
    vd.allelic_state.CopyFrom(allelic)
    return vd

def get_interpretation_status(affected_status: str, num_variants: int) -> str:
    if affected_status == 'AFFECTED':
        return 'CONTRIBUTORY' if num_variants == 2 else 'CAUSATIVE'
    return 'CANDIDATE'

def parse_variants_from_row(row):
    vds = []
    if row.get('Chrom-1') and row['Chrom-1'].strip():
        vd1 = build_variation_descriptor(
            row['Chrom-1'], row['Pos-1'], row['Ref-1'], row['Alt-1'],
            row['Gene'], row['Transcript'], row['HGVSc-1'], row['HGVSp-1'], row['Zygosity']
        )
        vds.append(vd1)
    if row.get('Chrom-2') and row['Chrom-2'].strip():
        vd2 = build_variation_descriptor(
            row['Chrom-2'], row['Pos-2'], row['Ref-2'], row['Alt-2'],
            row['Gene'], row['Transcript'], row['HGVSc-2'], row['HGVSp-2'], row['Zygosity']
        )
        vds.append(vd2)
    return vds

def parse_variant_interpretations(row, individual_id, affected_status):
    variation_descriptors = parse_variants_from_row(row)
    if not variation_descriptors:
        return []
    interpretation_status_str = get_interpretation_status(affected_status, len(variation_descriptors))
    genomic_interpretations = []
    for vd in variation_descriptors:
        gi = pps2.GenomicInterpretation()
        gi.subject_or_biosample_id = individual_id
        gi.interpretation_status = getattr(pps2.GenomicInterpretation.InterpretationStatus, interpretation_status_str)
        vi = pps2.VariantInterpretation()
        vi.variation_descriptor.CopyFrom(vd)
        gi.variant_interpretation.CopyFrom(vi)
        genomic_interpretations.append(gi)
    interp = pps2.Interpretation()
    interp.id = individual_id
    interp.progress_status = pps2.Interpretation.ProgressStatus.SOLVED if affected_status == 'AFFECTED' else pps2.Interpretation.ProgressStatus.COMPLETED
    disease = pps2.OntologyClass()
    disease.id = "MONDO:0003847"
    disease.label = "hereditary disease"
    diag = pps2.Diagnosis()
    diag.disease.CopyFrom(disease)
    diag.genomic_interpretations.extend(genomic_interpretations)
    interp.diagnosis.CopyFrom(diag)
    return [interp]

def parse_row_to_family_member(row, id_key='ID'):
    individual_id = row[id_key]
    match = re.match(r'^(.+)_(MOTHER|FATHER|PROBAND|SIBLING)', individual_id)
    if not match:
        raise ValueError(f"Invalid ID format: {individual_id}")
    family_id, role = match.groups()
    sex_str = row.get('SEX', 'UNKNOWN_SEX').upper()
    affected_str = row.get('AFFECTED STATUS', 'MISSING').upper()
    sex_enum = getattr(pps2.Sex, sex_str, pps2.Sex.UNKNOWN_SEX)
    affected = 'MISSING' if affected_str in ('UNKNOWN', 'MISSING') else affected_str
    hpo_string = row.get('HPO PRESENT', '')
    phenotypic_features = parse_hpo_terms(hpo_string)
    time_element = None
    age_onset = row.get('AGE OF ONSET (yrs)', '').strip()
    if age_onset and age_onset.replace('.', '', 1).isdigit():
        age = pps2.Age()
        age.iso8601duration = f"P{age_onset}Y"
        time_element = pps2.TimeElement()
        time_element.age.CopyFrom(age)
    interpretations = parse_variant_interpretations(row, individual_id, affected_str)
    subject = pps2.Individual()
    subject.id = individual_id
    subject.sex = sex_enum
    if time_element:
        subject.time_at_last_encounter.CopyFrom(time_element)
    pp = pps2.Phenopacket()
    pp.id = f"{family_id}_{role}"
    pp.subject.CopyFrom(subject)
    pp.phenotypic_features.extend(phenotypic_features)
    pp.interpretations.extend(interpretations)
    if phenotypic_features:
        pp.meta_data.CopyFrom(METADATA)
    else:
        minimal_md = pps2.MetaData()
        minimal_md.phenopacket_schema_version = "2.0"
        minimal_md.created.FromDatetime(datetime.now())
        minimal_md.created_by = "csv_to_phenopacket_gui.py"
        pp.meta_data.CopyFrom(minimal_md)
    return FamilyMember(family_id, individual_id, role, sex_enum, affected, pp)

def read_family_members_from_csv(input_path: Path):
    family_members = []
    with open(input_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                member = parse_row_to_family_member(row)
                family_members.append(member)
            except ValueError as e:
                print(f"Warning: Skipping row - {e}", file=sys.stderr)
    return family_members

def build_pedigree(family_id: str, members):
    father_id = next((m.individual_id for m in members if m.role == 'FATHER'), '0')
    mother_id = next((m.individual_id for m in members if m.role == 'MOTHER'), '0')
    persons = []
    for m in members:
        paternal_id = '0' if m.role in ['FATHER', 'MOTHER'] else father_id
        maternal_id = '0' if m.role in ['FATHER', 'MOTHER'] else mother_id
        p = pps2.Pedigree.Person()
        p.family_id = family_id
        p.individual_id = m.individual_id
        p.paternal_id = paternal_id
        p.maternal_id = maternal_id
        p.sex = m.sex
        p.affected_status = getattr(pps2.Pedigree.Person.AffectedStatus, m.affected)
        persons.append(p)
    pedigree = pps2.Pedigree()
    pedigree.persons.extend(persons)
    return pedigree

def build_family_from_members(family_id: str, members):
    proband = next((m.phenopacket for m in members if m.role == 'PROBAND'), None)
    relatives = [m.phenopacket for m in members if m.role != 'PROBAND']
    pedigree = build_pedigree(family_id, members)
    family = pps2.Family()
    family.id = family_id
    if proband:
        family.proband.CopyFrom(proband)
    family.relatives.extend(relatives)
    family.pedigree.CopyFrom(pedigree)
    family.meta_data.CopyFrom(METADATA)
    return family

def convert_csv_to_phenopackets(input_path: Path, output_dir: Path):
    os.makedirs(output_dir, exist_ok=True)
    family_members = read_family_members_from_csv(input_path)
    families = defaultdict(list)
    for member in family_members:
        families[member.family_id].append(member)
    for family_id, members in families.items():
        family = build_family_from_members(family_id, members)
        output_path = output_dir / f"{family_id}_PROBAND.json"
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(MessageToJson(family, indent=2))
    return len(families)

# ────────────────────────────────────────────────
#                Professional GUI
# ────────────────────────────────────────────────

class ConverterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("CSV to Phenopackets Converter")
        self.root.geometry("680x420")
        self.root.resizable(False, False)
        self.root.configure(bg="#f8f9fa")

        # Modern theme & styles
        style = ttk.Style()
        style.theme_use('clam')

        style.configure("TLabel", font=("Segoe UI", 10), foreground="#212529", background="#f8f9fa")
        style.configure("TEntry", font=("Segoe UI", 10), fieldbackground="#ffffff", foreground="#212529")
        style.configure("Main.TFrame", background="#f8f9fa")

        style.configure("Accent.TButton",
                        font=("Segoe UI", 12, "bold"),
                        foreground="white",
                        background="#005f99",
                        padding=12,
                        borderwidth=0)
        style.map("Accent.TButton",
                  background=[('active', '#004080'), ('disabled', '#cccccc')],
                  foreground=[('active', 'white'), ('disabled', '#666666')])

        # Main container
        main_frame = ttk.Frame(root, padding="40 30 40 20", style="Main.TFrame")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Header
        ttk.Label(main_frame,
                  text="CSV to Phenopackets Family Converter",
                  font=("Segoe UI", 16, "bold"),
                  foreground="#005f99").pack(pady=(0, 30))

        # ── Input section ──
        input_frame = ttk.Frame(main_frame)
        input_frame.pack(fill=tk.X, pady=(0, 20))

        ttk.Label(input_frame, text="Input CSV file:").pack(side=tk.LEFT, padx=(0, 12))

        self.input_path = tk.StringVar()
        ttk.Entry(input_frame, textvariable=self.input_path, width=60).pack(
            side=tk.LEFT, expand=True, fill=tk.X
        )

        ttk.Button(input_frame, text="Browse", command=self.browse_input, width=10).pack(
            side=tk.LEFT, padx=(12, 0)
        )

        # ── Output section ──
        output_frame = ttk.Frame(main_frame)
        output_frame.pack(fill=tk.X, pady=(0, 30))

        ttk.Label(output_frame, text="Output folder:").pack(side=tk.LEFT, padx=(0, 12))

        self.output_path = tk.StringVar(value=str(Path.home() / "Desktop" / "phenopackets_output"))
        ttk.Entry(output_frame, textvariable=self.output_path, width=60).pack(
            side=tk.LEFT, expand=True, fill=tk.X
        )

        ttk.Button(output_frame, text="Browse", command=self.browse_output, width=10).pack(
            side=tk.LEFT, padx=(12, 0)
        )

        # Convert button
        self.convert_btn = ttk.Button(main_frame,
                                      text="Convert to Phenopackets",
                                      command=self.run_conversion,
                                      style="Accent.TButton")
        self.convert_btn.pack(pady=25)

        # Status
        self.status_label = ttk.Label(main_frame, text="", foreground="#006600", font=("Segoe UI", 10))
        self.status_label.pack(pady=(10, 0))

        # Footer
        ttk.Label(main_frame,
                  text="GA4GH Phenopackets Schema v2.0 • Demo Conversion Tool",
                  font=("Segoe UI", 9), foreground="#6c757d").pack(side=tk.BOTTOM, pady=(30, 0))

    def browse_input(self):
        path = filedialog.askopenfilename(title="Select CSV file", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if path:
            self.input_path.set(path)

    def browse_output(self):
        path = filedialog.askdirectory(title="Select output folder")
        if path:
            self.output_path.set(path)

    def run_conversion(self):
        input_str = self.input_path.get().strip()
        output_str = self.output_path.get().strip()

        if not input_str or not Path(input_str).is_file():
            messagebox.showerror("Error", "Please select a valid CSV file.")
            return
        if not output_str:
            messagebox.showerror("Error", "Please select an output folder.")
            return

        self.convert_btn.state(['disabled'])
        self.status_label.config(text="Converting... please wait", foreground="#005f99")
        self.root.update_idletasks()

        try:
            input_p = Path(input_str)
            output_p = Path(output_str)
            num_families = convert_csv_to_phenopackets(input_p, output_p)
            msg = f"Conversion successful!\nCreated {num_families} Phenopacket family file(s) in:\n{output_p}"
            messagebox.showinfo("Success", msg)
            self.status_label.config(text="Ready", foreground="#006600")
            self._open_folder(output_p)
        except Exception as e:
            import traceback
            print(traceback.format_exc())
            messagebox.showerror("Conversion Failed", f"An error occurred:\n{str(e)}")
            self.status_label.config(text="Conversion failed – see console", foreground="#c62828")
        finally:
            self.convert_btn.state(['!disabled'])

    def _open_folder(self, path: Path):
        try:
            if sys.platform == "win32":
                os.startfile(path)
            elif sys.platform == "darwin":
                os.system(f"open '{path}'")
            else:
                os.system(f"xdg-open '{path}'")
        except:
            pass

def main():
    root = tk.Tk()
    app = ConverterApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
