"""Reviewed gene-level equivalence, without inferring primer-region identity."""


INDEXED_GENE_EQUIVALENCES = {
    "16s_1": "16s",
    "16s_2": "16s",
    "18s_1": "18s",
    "18s_2": "18s",
    "28s_1": "28s",
    "28s_2": "28s",
    "co1_1": "co1",
    "co1_2": "co1",
    "co2_1": "co2",
    "co2_2": "co2",
}


REVIEWED_PANEL_ITS_SIGNATURES = {
    "insecta": {
        "its_1": (
            "TACACACCGCCCGTCGCTACTA",
            "ACTCGCCGTTACTRRGG",
        ),
        "its_2": (
            "GTAGGTGAACCTGCAGAAGGATCA",
            "ACTCGCCGTTACTRRGG",
        ),
    },
    "cnidaria": {
        "its_1": (
            "TACACACCGCCCGTCGCTACTA",
            "ACTCGCCGTTACTRRGG",
        ),
        "its_2": (
            "GTAGGTGAACCTGCAGAAGGATCA",
            "ACTCGCCGTTACTRRGG",
        ),
    },
    "c_elegans": {
        "its_1": (
            "TACACACCGCCCGTCGCTATCC",
            "ACTCGCCGTTACTAAGG",
        ),
        "its_2": (
            "GTAGGTGAACCTGCAGCTGGATCA",
            "ACTCGCCGTTACTAAGG",
        ),
    },
}


def target_name(entry):
    gene = entry.get("gene", "")
    region = entry.get("region")
    index = entry.get("index")
    name = f"{gene}-{region}" if region is not None else gene
    return f"{name}_{index}" if index is not None else name


def _reviewed_panel_its_equivalence(panel_data):
    panel_name = str((panel_data or {}).get("name", "")).casefold()
    expected = REVIEWED_PANEL_ITS_SIGNATURES.get(panel_name)
    if expected is None:
        return False
    observed = {}
    for primer in panel_data.get("primers", []):
        name = target_name(primer).casefold()
        if name in expected:
            if name in observed:
                return False
            observed[name] = (
                primer.get("forward_seq"),
                primer.get("reverse_seq"),
            )
    return observed == expected


def target_logical_genes(panel_data):
    reviewed_its = _reviewed_panel_its_equivalence(panel_data)
    targets = {}
    for section in ("primers", "references"):
        for entry in (panel_data or {}).get(section, []):
            name = target_name(entry)
            normalized = name.casefold()
            if reviewed_its and normalized in {"its_1", "its_2"}:
                targets[name] = "its_rdna_cluster"
            else:
                targets[name] = INDEXED_GENE_EQUIVALENCES.get(normalized, normalized)
    return targets


def logical_gene_name(target, panel_data=None, target_mapping=None):
    normalized = (target or "").casefold()
    if target_mapping:
        for name, logical_gene in target_mapping.items():
            if str(name).casefold() == normalized:
                return str(logical_gene).casefold()
    if panel_data is not None:
        mapping = target_logical_genes(panel_data)
        for name, logical_gene in mapping.items():
            if name.casefold() == normalized:
                return logical_gene
    return INDEXED_GENE_EQUIVALENCES.get(normalized, normalized)


def reference_target_names(panel_data, reference_gene_names):
    mapping = target_logical_genes(panel_data)
    reference_logical_genes = {
        logical_gene_name(gene, panel_data=panel_data) for gene in reference_gene_names
    }
    return {
        target
        for target, logical_gene in mapping.items()
        if logical_gene in reference_logical_genes
    }
