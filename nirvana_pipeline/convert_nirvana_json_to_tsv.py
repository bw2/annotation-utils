"""This script converts an ExpansionHunter variant catalog to a GangSTR spec. Variant catalog entries
that include adjacent repeats (eg. with LocusStructure like "(A)*(ACG)*") are split into multiple GangSTR
specs - one per repeat.
"""

import argparse
import gzip
import pandas as pd
import re
import requests
import simplejson as json
import tqdm

from pprint import pprint

from bw2_annotation_utils.spliceai_scores import get_spliceai_scores_from_api

LOF_CONSEQUENCES = {
	'splice_acceptor_variant',
	'splice_donor_variant',
	'stop_gained',
	'frameshift_variant',
}

MISSENSE_CONSEQUENCES = {
	'inframe_insertion',
	'inframe_deletion',
	"stop_lost",
	"initiator_codon_variant",
	"start_lost",
	"protein_altering_variant",
	"missense_variant",
	"splice_region_variant",
}

SYNONYMOUS_CONSEQUENCES = {
	'synonymous_variant',
	'stop_retained_variant',
}

CONSEQUENCES_TO_IGNORE = {
	"downstream_gene_variant",
	"upstream_gene_variant",
}

def parse_nirvana_json(path, call_spliceai_api=False, verbose=False):
	print(f"Parsing {path}" +  (" and calling SpliceAI-lookup to add SpliceAI scores" if call_spliceai_api else ""))

	open_func = gzip.open if path.endswith(".gz") else open
	with open_func(path, "rt") as f:
		json_dict = json.load(f)

	header = json_dict.get("header", {})
	sample_ids = header["samples"]

	positions = json_dict['positions']
	if verbose:
		positions = tqdm.tqdm(positions, unit=" variants")

	for position in positions:
		for i, variant in enumerate(position['variants']):
			transcripts = variant.get('transcripts', [])
			if any(t.get('bioType') == "protein_coding" for t in transcripts):
				transcripts = [t for t in transcripts if t.get('bioType') == "protein_coding"]

			consequences = sorted({c for transcript in transcripts for c in transcript['consequence']})
			lof_consequences = sorted({c for c in consequences if c in LOF_CONSEQUENCES})
			missense_consequences = sorted({c for c in consequences if c in MISSENSE_CONSEQUENCES})
			synonymous_consequences = sorted({c for c in consequences if c in SYNONYMOUS_CONSEQUENCES})
			gene_names = list(sorted({t['hgnc'] for t in transcripts if 'hgnc' in t}))
			gene_name = gene_names[0] if gene_names else ""

			clinvar_significances = sorted({sig for c in variant.get('clinvar', []) for sig in c['significance']})
			spliceai_gain_score = spliceai_loss_score = None
			if 'spliceAI' in variant:
				spliceai_gain_score = max(
					max(t.get('donorGainScore', 0) for t in variant['spliceAI']),
					max(t.get('acceptorGainScore', 0) for t in variant['spliceAI'])
				)
				spliceai_loss_score = max(
					max(t.get('donorLossScore', 0) for t in variant['spliceAI']),
					max(t.get('acceptorLossScore', 0) for t in variant['spliceAI'])
				)

			chrom = position['chromosome']
			pos = position['position']
			ref = position['refAllele']
			alt = position['altAlleles'][i]

			record = {
				'chrom': position['chromosome'].replace("chr", "").upper(),
				'pos': pos,
				'ref': ref,
				'alt': alt,
				#'genotype': position['samples'][0]['genotype'],
				'variant_type': variant['variantType'],
				'consequences': ', '.join([c for c in consequences if c not in CONSEQUENCES_TO_IGNORE]),
				'synonymous_consequences': ', '.join(synonymous_consequences),
				'missense_consequences': ', '.join(missense_consequences),
				'lof_consequences': ', '.join(lof_consequences),
				'is_synonymous': len(synonymous_consequences) > 0,
				'is_missense': len(missense_consequences) > 0,
				'is_lof': len(lof_consequences) > 0,
				'gene': gene_name,
				'clinvar': ','.join(clinvar_significances),
				'spliceAI_gain': spliceai_gain_score,
				'spliceAI_loss': spliceai_loss_score,
				'1kg_af': variant.get('oneKg', {}).get('allAf', 0),
				'gnomad_af': variant.get('gnomad', {}).get('allAf', 0),
				'topmed_af': variant.get('topmed', {}).get('allAf', 0),
			}

			if call_spliceai_api:
				print(f"Calling SpliceAI-lookup for {chrom}-{pos}-{ref}-{alt}")
				record['spliceAI_gain'], record['spliceAI_loss'] = get_spliceai_scores_from_api(
					chrom.replace("chr", "").upper(), pos, ref, alt)

			for sample_id, sample in zip(sample_ids, position['samples']):
				record[sample_id] = sample['genotype']

			yield record




def main():
	p = argparse.ArgumentParser()
	p.add_argument("-o", "--output-file", help="bed file output path")
	p.add_argument("-s", "--call-spliceai-api", action="store_true", help="Use the SpliceAI-lookup API to get spliceAI scores")
	p.add_argument("--verbose", action="store_true")
	p.add_argument("nirvana_json", help="Path of Nirvana JSON file")
	args = p.parse_args()

	if not args.output_file:
		args.output_file = re.sub(".json(.gz)?$", "", args.nirvana_json) + ".tsv.gz"

	output_rows = parse_nirvana_json(args.nirvana_json, call_spliceai_api=args.call_spliceai_api, verbose=args.verbose)
	df = pd.DataFrame(output_rows)
	df.to_csv(args.output_file, sep="\t", index=False)

	print(f"Wrote {len(df)} rows to {args.output_file}")

if __name__ == "__main__":
	main()
