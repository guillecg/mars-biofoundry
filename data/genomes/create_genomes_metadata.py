import pandas as pd


metadata_df = pd.DataFrame.from_records([
    {
        "Species": "Tessaracoccus lapidicaptus IPBSL-7",
        "Code": "tel",
        "Database": "IMG",
        "ID": "2791354959",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/MBQD00000000.1/",
        "Download": "https://www.ncbi.nlm.nih.gov/Traces/wgs/MBQD01",
        "Protein annotation file": "MBQD01P.1.fsa_aa"
    },
    {
        "Species": "Tessaracoccus sp. T2.5-30",
        "Code": "tez",
        "Database": "IMG",
        "ID": "2751185744",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/CP003629.1",
        "Download": "https://www.ncbi.nlm.nih.gov/nuccore/CP003629.1",
        "Protein annotation file": "CP019229.1.faa"
    },
    {
        "Species": "Desulfosporosinus meridiei DEEP",
        "Code": "dmi",
        "Database": "IMG",
        "ID": "2721755100",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/CP003629.1",
        "Download": "https://www.ncbi.nlm.nih.gov/nuccore/CP003629.1",
        "Protein annotation file": "CP003629.1.faa"
    },
    {
        "Species": "Brevundimonas sp. T2.26MG-97",
        "Code": "bme",
        "Database": "NCBI",
        "ID": "NZ_UXHF01000001.1",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/NZ_UXHF01000001.1",
        "Download": "https://www.ncbi.nlm.nih.gov/nuccore/NZ_UXHF01000001.1",
        "Protein annotation file": "NZ_UXHF01000001.1.faa"
    },
    {
        "Species": "Rhizobium sp. T2.30D-1.1",
        "Code": "rhi1",
        "Database": "NCBI",
        "ID": "NZ_UEYP01000001.1",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/NZ_UEYP01000001.1",
        "Download": "https://www.ncbi.nlm.nih.gov/nuccore/NZ_UEYP01000001.1",
        "Protein annotation file": "NZ_UEYP01000001.1.faa"
    },
    {
        "Species": "Rhizobium sp. T2.26MG-112.2",
        "Code": "rhi2",
        "Database": "NCBI",
        "ID": "NZ_UEYQ01000001.1",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/NZ_UEYQ01000001.1",
        "Download": "https://www.ncbi.nlm.nih.gov/nuccore/NZ_UEYQ01000001.1",
        "Protein annotation file": "NZ_UEYQ01000001.1.faa"
    },
    {
        "Species": "Rhodoplanes sp. T2.26MG-98",
        "Code": "rho",
        "Database": "NCBI",
        "ID": "NZ_UWOC01000001.1",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/NZ_UWOC01000001.1",
        "Download": "https://www.ncbi.nlm.nih.gov/nuccore/NZ_UWOC01000001.1",
        "Protein annotation file": "NZ_UWOC01000001.1.faa"
    },
    {
        "Species": "Pseudomonas sp. T2.31D-1",
        "Code": "pse",
        "Database": "ENA",
        "ID": "CAJFAG010000000.1",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/CAJFAG010000000",
        "Download": "https://www.ncbi.nlm.nih.gov/Traces/wgs/CAJFAG01",
        "Protein annotation file": "CAJFAG01.1.fsa_aa"
    },
    {
        "Species": "Shewanella sp. T2.3D-1.1",
        "Code": "shw",
        "Database": "ENA",
        "ID": "CACVBT0200000010",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/CACVBT000000000.3/",
        "Download": "https://www.ncbi.nlm.nih.gov/Traces/wgs/CACVBT03",
        "Protein annotation file": "CACVBT03.1.fsa_aa"
    },
    {
        "Species": "Cyanobacteria IPBSL",
        "Code": "cya",
        "Database": "MG-RAST",
        "ID": "mgp83581",
        "Record": "https://www.mg-rast.org/mgmain.html?mgpage=project&project=mgp83581",
        "Download": "https://www.mg-rast.org/mgmain.html?mgpage=download&metagenome=mgm4729322.3#annotationDownloads",
        "Protein annotation file": None
    },
    {
        "Species": "Acidovorax BoFeN1",
        "Code": "aci",
        "Database": "NCBI",
        "ID": "QOZT00000000.1",
        "Record": "https://www.ncbi.nlm.nih.gov/nuccore/QOZT00000000.1/",
        "Download": "https://www.ncbi.nlm.nih.gov/Traces/wgs/QOZT01",
        "Protein annotation file": "QOZT01.1.fsa_aa"
    }
])


metadata_df.to_csv(
    "genome-metadata.csv",
    index=False,
    header=True,
    sep=",",
    mode="w"
)
