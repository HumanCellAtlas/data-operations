# This script interrogates a project metadata tsv by row
# This currently validates against hard-coded schema rules, rather than a live schema version
# However, ontologies are validated against the live OLS API
# usage is $ python3 metadatatsv_validator.py <file_name.tsv>

# Takes a long time, usually times out due to 1733 IDs for 'cell_suspension.timecourse.unit.ontology'
# Single\ cell\ RNAseq\ characterization\ of\ cell\ types\ produced\ over\ time\ in\ an\ in\ vitro\ model\ of\ human\ inhibitory\ interneuron\ differentiation.\ 2019-10-24\ 18.58.tsv
# https://github.com/HumanCellAtlas/hca-data-wrangling/issues/225

import sys
import json
import re
import requests

methods_10x = [
	"10X 5' v2 sequencing",
	"10x 3' v3 sequencing",
	"10X 3' v2 sequencing",
	"10X 3' v1 sequencing"
	]

methods_other = [
	'Smart-seq',
	'Smart-seq2',
	'inDrop',
	'Drop-seq',
	'DroNc-seq',
	'CITE-seq',
	'MARS-seq'
	]

boolean = ['True', 'False']

# a list of acceptable values for each field
enums = {
	'embedded_in_matrigel': boolean,
	'paired_end': boolean,
	'cold_perfused': boolean,
	'corresponding_contributor': boolean,
	'drop_uniformity': boolean,
	'ncbi_taxon_id': [
		'9606',
		'10090'
		],
	'genus_species.ontology': [
		'NCBITaxon:9606',
		'NCBITaxon:10090'
		],
	'cell_viability_result': [
		'pass',
		'fail'
		],
	'mycoplasma_testing_results': [
		'pass',
		'fail'
		],
	'is_living': [
		'yes',
		'no',
		'unknown',
		'not applicable'
		],
	'multiplexed': [
		'yes',
		'no'
		],
	'normothermic_regional_perfusion': [
		'yes',
		'no',
		'unknown'
		],
	'overlapping_tiles': [
		'yes',
		'no'
		],
	'pluripotency_vector_removed': [
		'yes',
		'no',
		'unknown'
		],
	'mycoplasma_testing_method': [
		'Direct DNA stain',
        'Indirect DNA stain',
        'Broth and agar culture',
        'PCR',
        'Nested PCR',
        'ELISA',
        'Autoradiography',
        'Immunostaining',
        'Cell-based assay',
        'Microbiological assay'
		],
	'feeder_layer_type': [
		'feeder-free',
        'feeder-dependent, JK1 feeder cells',
        'feeder-dependent, SNL 76/7 feeder cells',
        'feeder-dependent, human marrow stromal cells',
        'feeder-dependent, bovine embryonic fibroblast cells',
        'feeder-dependent, mouse embryonic fibroblast cells',
        'feeder-dependent, mouse fibroblast STO cells',
        'feeder-dependent, mouse bone marrow stromal cells',
        'feeder-dependent, mouse yolk sac-derived endothelial cells',
        'feeder-dependent, human foreskin fibroblast cells',
        'feeder-dependent, human newborn fibroblast cells',
        'feeder-dependent, human fetal lung fibroblast cells',
        'feeder-dependent, human uterine endometrial cells',
        'feeder-dependent, human breast parenchymal cells',
        'feeder-dependent, human embryonic fibroblast cells',
        'feeder-dependent, human adipose stromal cells',
        'feeder-dependent, human amniotic epithelial cells',
        'feeder-dependent, human placental fibroblast cells',
        'feeder-dependent, human umbilical cord stromal cells',
        'feeder-dependent, human fetal muscle cells',
        'feeder-dependent, human fetal skin cells',
        'feeder-dependent, human fetal liver stromal cells',
        'feeder-dependent, human fallopian tubal epithelial cells',
        'feeder-dependent, human amniotic mesenchymal cells'
		],
	'nutritional_state': [
		'normal',
		'fasting',
		'feeding tube removed'
		],
	'storage_method': [
		'ambient temperature',
        'cut slide',
        'fresh',
        'frozen at -70C',
        'frozen at -80C',
        'frozen at -150C',
        'frozen in liquid nitrogen',
        'frozen in vapor phase',
        'paraffin block',
        'RNAlater at 4C',
        'RNAlater at 25C',
        'RNAlater at -20C'
		],
	'preservation_method': [
		'cryopreservation in liquid nitrogen (dead tissue)',
        'cryopreservation in dry ice (dead tissue)',
        'cryopreservation of live cells in liquid nitrogen',
        'cryopreservation, other',
        'formalin fixed, unbuffered',
        'formalin fixed, buffered',
        'formalin fixed and paraffin embedded',
        'hypothermic preservation media at 2-8C',
        'fresh'
		],
	'autolysis_score': [
		'none',
		'mild',
		'moderate'
		],
	'ischemic_temperature': [
		'warm',
		'cold'
		],
	'barcode_read': [
		'Read 1',
        'Read 2',
        'i7 Index',
        'i5 Index'
		],
	'well_quality': [
		'OK',
        'control, 2-cell well',
        'control, empty well',
        'low quality cell'
		],
	'input_type': [
		'biomaterial',
		'file'
		],
	'output_type': [
		'biomaterial',
		'file'
		],
	'sex': [
		'male',
		'female',
		'mixed',
		'unknown'
		],
	'reference_type': [
		'genome sequence',
        'transcriptome sequence',
        'annotation reference',
        'transcriptome index',
        'genome sequence index'
		],
	'assembly_type': [
		'primary assembly',
        'complete assembly',
        'patch assembly'
		],
	'read_index': [
		'read1',
        'read2',
        'index1',
        'index2',
        'single-end, non-indexed'
		],
	'nucleic_acid_source': [
        'bulk cell',
        'single cell',
        'single nucleus',
        'bulk nuclei',
        'mitochondria'
		],
	'end_bias': [
        '3 prime tag',
        '3 prime end bias',
        '5 prime tag',
        '5 prime end bias',
        'full length'
		],
	'primer': [
		'poly-dT',
		'random'
		],
	'strand': [
		'first',
		'second',
		'unstranded',
		'not provided'
		],
	'analysis_run_type': [
		'run',
		'copy-forward'
		],
	'cell_line.type': [
        'primary',
        'immortalized',
        'stem cell',
        'stem cell-derived',
        'induced pluripotent',
        'synthetic'
		],
	'ipsc_induction_protocol.method': [
		'lentivirus',
        'sendai virus',
        'Gun particle',
        'piggyBac transposon',
        'miRNA viral',
        'adenovirus',
        'cre-loxP',
        'plasmid',
        'retroviral'
		],
	'analysis_run_type': [
		'run'
		],
	'sequence_file.file_core.format': [
		'fastq'
		],
	'analysis_file.file_core.format': [
        'bam',
        'bai',
        'csv',
        'npy',
        'npz',
        'zarr'
		],
	'supplementary_file.file_core.format': [
		'docx',
        'pdf'
		],
	'cell_viability_method': [
		'Green calcein-AM stain and Red ethidium homodimer-1 stain',
		'Visual inspection',
		'Trypan blue stain',
		'Manual haemocytometer (C-chip)',
		'Flow cytometry staining'
		],
	'library_construction_method.ontology_label':
		methods_10x + methods_other
	}

time_pattern = '\d{4}-(0[1-9]|1[0-2])(-(0[1-9]|[1-2]\d|3[0-1])(T([0-1]\d|2[0-4]):[0-6]\d:[0-6]\d(.\d{3})?Z)?)?'
integer_pattern = '\d+'
float_pattern = '\d+(.\d+)?'
range_float_int_pattern = '\d+.?\d*-?\d*.?\d*' # THIS ISNT GREAT

regex = {
	'protocol_core.document': '.*.pdf',
	'array_express_accessions': 'E-[A-Z]{4}-[0-9]{4,5}',
	'insdc_project_accessions': 'PRJ(NA|EB)\d{6}',
	'insdc_study_accessions': '[ES]RP\d{6}',
	'insdc_sample_accession': '[ES]RS\d{7}',
	'insdc_experiment_accession': '[ES]RX\d{7}',
	'geo_series_accessions': 'GSE\d{5,6}',
	#'biostudies_accessions': ,
	'biosamples_accession': 'SAM((EA\d{7})|(N\d{8}))',
	'protocols_io_doi': '10.17504/protocols.io.\w{7}',
	'publication_doi': '10.\d{4,5}/.+',
	'pmid': '[0-9]{7,8}',
	'collection_time': time_pattern,
	'start_time': time_pattern,
	'date_established': time_pattern,
	'expiry_date': time_pattern,
	'cell_size': range_float_int_pattern,
	'height': range_float_int_pattern,
	'weight': range_float_int_pattern,
	'value': integer_pattern,
	'read_length': integer_pattern,
	'passage_number': integer_pattern,
	'barcode_length': integer_pattern,
	'barcode_offset': integer_pattern,
	'nominal_length': integer_pattern,
	'spike_in_dilution': integer_pattern,
	'lane_index': integer_pattern,
	'ischemic_time': integer_pattern,
	'hardy_scale': '0|1|2|3|4',
	'pooled_channels': integer_pattern,
	'days_on_ventilator': float_pattern,
	'maximum_size': float_pattern,
	'minimum_size': float_pattern,
	'age': float_pattern,
	'gestational_age': float_pattern,
	'organism_age': float_pattern,
	'percent_cell_viability': '1?\d{1,2}(.\d{1})?'
	}

species_dev_db = {'9606':'HsapDv','10090':'EFO'}

# if metadata exists for the key, then metadata should exist for the corresponding value, and vice versa
# the way this is checked, a property can only appear as a key once and a value once
dependencies = {
	'text': 'ontology',
	'ontology': 'ontology_label',
	'storage_time': 'storage_time_unit.ontology',
	'age': 'age_unit.ontology',
	'organism_age': 'organism_age_unit.ontology',
	'gestational_age': 'gestational_age_unit.ontology',
	'height': 'height_unit.ontology',
	'weight': 'weight_unit.ontology',
	'cell_size': 'cell_size_unit.ontology',
	'timecourse.value': 'timecourse.unit.ontology',
	'cell_type.ontology': 'model_organ.ontology'
	}

value_dependencies = {
	'test1': {
		'if_statement':
			('library_preparation_protocol.library_construction_method.ontology_label',methods_10x),
		'then_statement': 
			('sequencing_protocol.paired_end','False')
			}
		}

# only one of each pair should be present
negative_dependencies = {
	'gestational_age': 'organism_age',
	'organ': 'model_organ',
	'organ_parts': 'model_organ_parts'
	}

# acceptable ontology databases for each proprerty
ontology_regex = {
	'genus_species': ['NCBITaxon'],
	'development_stage': ['EFO','HsapDv'],
	'organ': ['UBERON'],
	'model_organ': ['UBERON'],
	'organ_parts': ['UBERON','CL','FMA'],
	'model_organ_parts': ['UBERON','CL','FMA'],
	'cell_type': ['CL'],
	'selected_cell_types': ['CL'],
	'diseases': ['MONDO'],
	'disease': ['MONDO'],
	'strain': ['EFO'],
	'method': ['EFO'],
	'instrument_manufacturer_model': ['EFO'],
	'library_construction_method': ['EFO'],
	'library_preamplification_method': ['EFO','OBI'],
	'cdna_library_amplification_method': ['EFO','OBI'],
	'input_nucleic_acid_molecule': ['CHEBI','OBI','EFO'],
	'content_description': ['data'],
	'unit': ['UO'],
	'cell_size_unit': ['UO'],
	'gestational_age_unit': ['UO'],
	'weight_unit': ['UO'],
	'height_unit': ['UO'],
	'organism_age_unit': ['UO'],
	'storage_time_unit': ['UO'],
	'age_unit': ['UO'],
	'ethnicity': ['HANCESTRO']
	}

child_parent = {
	'organ_parts': 'organ',
	'model_organ_parts': 'model_organ',
	'cell_type': 'model_organ'
	#'selected_cell_types'
	}

ont_iri = {
	'data': (		'edam',		'http://edamontology.org/'),
	'EFO': (		'efo',		'http://www.ebi.ac.uk/efo/'),
	'CHEBI': (		'chebi',	'http://purl.obolibrary.org/obo/'),
	'CL': (			'cl',		'http://purl.obolibrary.org/obo/'),
	'FMA': (		'fma',		'http://purl.obolibrary.org/obo/'),
	'HANCESTRO': (	'hancestro','http://purl.obolibrary.org/obo/'),
	'HsapDv': (		'hsapdv',	'http://purl.obolibrary.org/obo/'),
	'MONDO': (		'mondo',	'http://purl.obolibrary.org/obo/'),
	'NCBITaxon': (	'ncbitaxon','http://purl.obolibrary.org/obo/'),
	'OBI': (		'obi',		'http://purl.obolibrary.org/obo/'),
	'PATO': (		'pato',		'http://purl.obolibrary.org/obo/'),
	'UBERON': (		'uberon',	'http://purl.obolibrary.org/obo/'),
	'UO': (			'uo',		'http://purl.obolibrary.org/obo/')
	}


def formatError(base_object,row,error_message):
	error_fields = [
		project,
		'row:' + str(row),
		'bundle:' + bundle_metadata[row],
		base_object + ':' + full_metadata[base_object + '.provenance.document_id'][row],
		error_message + '\n'
		]
	return('\t'.join(error_fields))


def lookupOntology(term):
	namespace = term.split(':')[0]
	if namespace in ont_iri.keys():
		ontology = ont_iri[namespace][0]
		iri_base = ont_iri[namespace][1].replace(':','%253A').replace('/','%252F')
		iri = iri_base + term.replace(':','_')
		URL = 'http://www.ebi.ac.uk/ols/api/ontologies/' + ontology + '/terms/' + iri
		if requests.get(URL).status_code == 200:
			HEADERS = {'accept': 'application/json'}
			response_json_dict = requests.get(URL, headers=HEADERS).json()	
			return(response_json_dict['label'])
		else:
			return('not a valid term')
	else:
		return('not a valid ontology DB')


def lookupAncestors(term):
	namespace = term.split(':')[0]
	if namespace in ont_iri.keys():
		ontology = ont_iri[namespace][0]
		ancestors = []
		URL = 'http://www.ebi.ac.uk/ols/api/ontologies/' + ontology + '/hierarchicalAncestors?id=' + term
		if requests.get(URL).status_code == 200:
			HEADERS = {'accept': 'application/json'}
			response_json_dict = requests.get(URL, headers=HEADERS).json()	
			for term in response_json_dict['_embedded']['terms']:
				ancestor_term = term['obo_id']
				if ancestor_term:
					ancestors.append(ancestor_term)
		else:
			return('not a valid term')
		return(ancestors)
	else:
		return('not a valid ontology DB')


def checkOntologies(base_object,field,base_field,two_level_module_field):
	errors = ''
	ontology_dict = {} # store looked-up ID/names to prevent duplicate API calls
	ancestor_dict = {}
	ontology_label_field = base_field + '.ontology_label'
	ontology_object = two_level_module_field.split('.')[0]
	expected_namespaces = ontology_regex[ontology_object]
	for row in full_metadata[field].keys():
		term_name = full_metadata[ontology_label_field][row]
		if '||' in full_metadata[field][row]:
			for term_id in full_metadata[field][row].split('||'):
				#check rgex of ontology DB namespace - THIS MAY BE TOO RESTRICTIVE
				metadata_namespace = term_id.split(':')[0]
				if two_level_module_field in ['disease.ontology','diseases.ontology'] and term_id == 'PATO:0000461':
					error_message = (field + ' contains normal [PATO:0000461] but also other disease(s)')
					errors += formatError(base_object,row,error_message)
				elif metadata_namespace not in expected_namespaces:
				 	error_message = (field + ':' + term_id + ' but expected one of ' + str(expected_namespaces))
				 	errors += formatError(base_object,row,error_message)
				#check ID matches label
				if term_id in ontology_dict.keys():
					ontology_term_name = ontology_dict[term_id]
				else:
					ontology_term_name = lookupOntology(term_id)
					ontology_dict[term_id] = ontology_term_name
				if ontology_term_name not in term_name.split('||'):
					error_message = (field + ':' + term_id + ' is ' + ontology_term_name + ' but input as ' + term_name)
					errors += formatError(base_object,row,error_message)
		elif full_metadata[field][row] != '':
			#check rgex of ontology DB namespace - THIS MAY BE TOO RESTRICTIVE
			metadata_namespace = full_metadata[field][row].split(':')[0]
			term_id = full_metadata[field][row]
			if metadata_namespace not in expected_namespaces and not (two_level_module_field in ['disease.ontology','diseases.ontology'] and term_id == 'PATO:0000461'):
			 	error_message = (field + ':' + full_metadata[field][row] + ' but expected one of ' + str(expected_namespaces))
			 	errors += formatError(base_object,row,error_message)
			#check ID matches label
			if term_id in ontology_dict.keys():
				ontology_term_name = ontology_dict[term_id]
			else:
				ontology_term_name = lookupOntology(term_id)
				ontology_dict[term_id] = ontology_term_name
			if ontology_term_name != term_name:
				error_message = (field + ':' + term_id + ' is ' + ontology_term_name + ' but input as ' + term_name)
				errors += formatError(base_object,row,error_message)
			if ontology_object in child_parent.keys():
				parent_term_field = field.replace(ontology_object,child_parent[ontology_object])
				parent_term_id = full_metadata[parent_term_field][row]
				if term_id in ancestor_dict.keys():
					ontology_parents = ancestor_dict[term_id]
				else:
					ontology_parents = lookupAncestors(term_id)
					ancestor_dict[term_id] = ontology_parents
				if parent_term_id != '' and parent_term_id not in ontology_parents:
					error_message = (field + ':' + term_id + ' is not a child of ' + parent_term_field + ':' + parent_term_id)
					errors += formatError(base_object,row,error_message)
	return(errors)


def checkSpecies(base_object,row,species_fields,dev_stage_fields):
	errors = ''
	species_values = {}
	id_list = []
	for field in species_fields:
		species_id = full_metadata[field][row].split(':')[-1]
		if species_id != '':
			species_values[field] = species_id
			id_list.append(species_id)
	if all(v == id_list[0] for v in id_list):
		species = id_list[0]
		for field in dev_stage_fields:
			term_namespace = full_metadata[field][row].split(':')[0]
			if term_namespace != species_dev_db[species]:
				error_message = (field + ':' + full_metadata[field][row] + ' but species is ' + species)
				errors += formatError(base_object,row,error_message)
	else:
		error_message = ('multiple species in ' + str(species_values))
		errors += formatError(base_object,row,error_message)
	return(errors)


def checkDependencies(base_object,field,base_field,module_field):
	errors = ''
	if module_field in dependencies.keys():
		depend_on = base_field + '.' + dependencies[module_field]
		if depend_on not in full_metadata.keys():
			for row in bundle_metadata.keys():
				if full_metadata[field][row] != '':
					value = full_metadata[field][row]
					error_message = (field + ':' + value + ' entered but not ' + depend_on)
					errors += formatError(base_object,row,error_message)
		else:
			for row in bundle_metadata.keys():
				if (full_metadata[field][row] != '' and full_metadata[depend_on][row] == ''):
					value = full_metadata[field][row]
					error_message = (field + ':' + value + ' entered but not ' + depend_on)
					errors += formatError(base_object,row,error_message)
				elif (full_metadata[field][row] == '' and full_metadata[depend_on][row] != ''):
					value = full_metadata[depend_on][row]
					error_message = (depend_on + ':' + value + ' entered but not ' + field)
					errors += formatError(base_object,row,error_message)
	if module_field in dependencies.values():
		for key, pair in dependencies.items():
			if pair == module_field:
				k = key
		depend_on = base_field + '.' + k
		if depend_on not in full_metadata.keys():
			for row in bundle_metadata.keys():
				if full_metadata[field][row] != '':
					value = full_metadata[field][row]
					error_message = (field + ':' + value + ' entered but not ' + depend_on)
					errors += formatError(base_object,row,error_message)
		else:
			for row in bundle_metadata.keys():
				if (full_metadata[field][row] != '' and full_metadata[depend_on][row] == ''):
					value = full_metadata[field][row]
					error_message = (field + ':' + value + ' entered but not ' + depend_on)
					errors += formatError(base_object,row,error_message)
				elif (full_metadata[field][row] == '' and full_metadata[depend_on][row] != ''):
					value = full_metadata[depend_on][row]
					error_message = (depend_on + ':' + value + ' entered but not ' + field)
					errors += formatError(base_object,row,error_message)
	return(errors)


def checkValuedDependencies(base_object,row):
	errors = ''
	for test in value_dependencies.keys():
		if_field = value_dependencies[test]['if_statement'][0]
		if_value = value_dependencies[test]['if_statement'][1]
		if full_metadata[if_field][row] in if_value:
			then_field = value_dependencies[test]['then_statement'][0]
			then_value = value_dependencies[test]['then_statement'][1]
			if full_metadata[then_field][row] != then_value:
				error_message =  (if_field + ':' + full_metadata[if_field][row] + ' but ' + then_field + ':' + full_metadata[then_field][row] + ' should be ' + then_value)
				errors += formatError(base_object,row,error_message)
	return(errors)


def checkNegativeDependencies(base_object,field,base_field,module_field):
	errors = ''
	anti_dependent = base_field + '.' + negative_dependencies[module_field]
	if anti_dependent in full_metadata.keys():
		for row in full_metadata[field].keys():
			if full_metadata[field][row] != '' and full_metadata[anti_dependent][row] != '':
				error_message =  ('both ' + field + ':' + full_metadata[field][row] + ' and ' + anti_dependent + ':' + full_metadata[anti_dependent][row] + ' are entered')
				errors += formatError(base_object,row,error_message)
	return(errors)


def checkRegex(base_object,field,module_field):
	errors = ''
	pattern = '^' + regex[module_field] + '$'
	for row in full_metadata[field].keys():
		if '||' in full_metadata[field][row]:
			for value in full_metadata[field][row].split('||'):
				z = re.match(pattern, value)
			if not z:
				error_message =  (field + ':' + value + ' does not match pattern')
				errors += formatError(base_object,row,error_message)
		elif full_metadata[field][row] != '':
			z = re.match(pattern, full_metadata[field][row])
			if not z:
				error_message =  (field + ':' + full_metadata[field][row] + ' does not match pattern')
				errors += formatError(base_object,row,error_message)
	return(errors)


def checkEnums(base_object,field,module_field):
	errors = ''
	for row in full_metadata[field].keys():
		if '||' in full_metadata[field][row]:
			for value in full_metadata[field][row].split('||'):
				if value.strip('"') not in enums[module_field] and value != '':
					error_message =  (field + ':' + value.strip('"') + ' not in ' + str(enums[module_field]))
					errors += formatError(base_object,row,error_message)
		elif full_metadata[field][row].strip('"') not in enums[module_field] and full_metadata[field][row] != '':
			error_message =  (field + ':' + full_metadata[field][row].strip('"') + ' not in ' + str(enums[module_field]))
			errors += formatError(base_object,row,error_message)
	return(errors)


def checkLeadingCharacter(base_object,field):
	errors = ''
	for row in full_metadata[field].keys():
		for value in full_metadata[field][row].split('||'):
			if value.startswith((' ','" ')):
				error_message =  (field + ':' + full_metadata[field][row] + ' has a leading whitespace')
				errors += formatError(base_object,row,error_message)
			elif value.startswith(('=','"=')):
				error_message =  (field + ':' + full_metadata[field][row]  + ' is an excel formula')
				errors += formatError(base_object,row,error_message)
			if value.endswith((' ','" ')):
				error_message =  (field + ':' + full_metadata[field][row]  + ' has a trailing whitespace')
				errors += formatError(base_object,row,error_message)
	return(errors)


def transformMetadata(file):
	metadata_dictionary = {}
	bundle_index = {}
	header_index = {}

	f = open(file)
	lines = f.read().split('\n')[:-1]
	count = 0
	for i, line in enumerate(lines):
		if i == 0:
			column_names = line.split('\t')
			for header in column_names:
				header_index[column_names.index(header)] = header
				metadata_dictionary[header] = {}
		else:
			count += 1
			data = line.split('\t')
			if len(data) != len(metadata_dictionary.keys()):
				sys.exit('Parsing issue with the field:' + header_index[len(data)-1])
			for column_number in range(0,len(metadata_dictionary.keys())):
				datum = data[column_number]
				header_name = header_index[column_number]
				metadata_dictionary[header_name][count] = datum
				if header_name == 'bundle_uuid':
					bundle_index[count] = datum
	f.close()
	return(metadata_dictionary,bundle_index)


file = sys.argv[1]
metadata = transformMetadata(file)
full_metadata = metadata[0]
bundle_metadata = metadata[1]
project = full_metadata['project.provenance.document_id'][1]
species_fields = []
dev_stage_fields = []
errors = ''

for field in full_metadata.keys():
	base_object = field.split('.')[0]
	property_sections = {}
	base_field = '.'.join(field.split('.')[:-1])
	module_field = field.split('.')[-1]
	property_sections[module_field] = base_field
	if len(field.split('.')) > 1:
		two_level_module_field = '.'.join(field.split('.')[-2:])
		double_base_field = '.'.join(field.split('.')[:-2])
		property_sections[two_level_module_field] = double_base_field
		if len(field.split('.')) > 2:
			three_level_module_field = '.'.join(field.split('.')[-3:])
			three_base_field = '.'.join(field.split('.')[:-3])
			property_sections[three_level_module_field] = three_base_field

	errors += checkLeadingCharacter(base_object,field)

	if module_field == 'ontology' and (base_field + '.ontology_label') in full_metadata.keys():
		errors += checkOntologies(base_object,field,base_field,two_level_module_field)

	for partial_field in property_sections.keys():
		if partial_field in ['ncbi_taxon_id','genus_species.ontology']:
			species_fields.append(field)
		if partial_field == 'development_stage.ontology':
			dev_stage_fields.append(field)
		if partial_field in enums.keys():
			errors += checkEnums(base_object,field,partial_field)
		if partial_field in regex.keys():
			errors += checkRegex(base_object,field,partial_field)
		if partial_field in dependencies.keys() or partial_field in dependencies.values():
			errors += checkDependencies(base_object,field,property_sections[partial_field],partial_field)
		if partial_field in negative_dependencies.keys():
			errors += checkNegativeDependencies(base_object,field,property_sections[partial_field],partial_field)

for row in bundle_metadata.keys():
	errors += checkSpecies(base_object,row,species_fields,dev_stage_fields)
	errors += checkValuedDependencies(base_object,row)

if errors == '':
	print('no validation errors found')
else:
	results_file = open(project + '_errors.txt','w')
	results_file.write(errors)
	results_file.close()
	print(project + '_errors.txt written')
