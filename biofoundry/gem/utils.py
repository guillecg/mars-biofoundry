def get_gene_counts(path: str) -> int:
    """
    Retrieve the number of genes in an annotated genome file.

    Parameters
    ----------
    path : str
        The path to the annotated genome.

    Returns
    -------
    gene_counts : int
        The number of genes in the annotated genome.

    Examples
    --------
    None

    """

    with open(path, mode="r") as filein:
        gene_counts = filein.read().count(">")

    return gene_counts
