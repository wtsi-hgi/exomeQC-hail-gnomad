import hail as hl

SUBPOPS = {'NFE': ['BGR', 'EST', 'NWE', 'SEU', 'SWE', 'ONF'],
           'EAS': ['KOR', 'JPN', 'OEA']
           }
GENOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
EXOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
EXAC_POPS = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]

pop_names = {
    'oth': 'Other',
    'afr': 'African-American/African',
    'ami': 'Amish',
    'amr': 'Latino',
    'eas': 'East Asian',
    'fin': 'Finnish',
    'eur': 'European',
    'nfe': 'Non-Finnish European',
    'sas': 'South Asian',
    'mde': 'Middle Eastern',
    'asj': 'Ashkenazi Jewish',
    'uniform': 'Uniform',
    'sas_non_consang': 'South Asian (F < 0.05)',
    'consanguineous': 'South Asian (F > 0.05)',
    'exac': 'ExAC',
    'bgr': 'Bulgarian (Eastern European)',
    'deu': 'German',
    'est': 'Estonian',
    'esp': 'Spanish',
    'gbr': 'British',
    'nwe': 'North-Western European',
    'seu': 'Southern European',
    'ita': 'Italian',
    'swe': 'Swedish',
    'chn': 'Chinese',
    'kor': 'Korean',
    'hkg': 'Hong Kong',
    'sgp': 'Singaporean',
    'twn': 'Taiwanese',
    'jpn': 'Japanese',
    'oea': 'Other East Asian',
    'oeu': 'Other European',
    'onf': 'Other Non-Finnish European',
    'unk': 'Unknown'
}

REFERENCES = {'GRCh37', 'GRCh38'}
