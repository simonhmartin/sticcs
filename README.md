# sticcs

`sticcs` is a method for inferring the series of genealogies along the genome, also called an Ancestral Recombination Graph (ARG) or [tree sequence](https://tskit.dev/tutorials/what_is.html). Unlike some other methods, `sticcs` **does not require phased haplotypes**, and it can work on **any ploidy level**.

The input for `sticcs` is polarised genotype data. This means you need to know the ancestral allele at each site, or you need an appropriate outgroup(s) to allow inference of the derived allele.

The method is described in [this paper](https://doi.org/10.1093/genetics/iyaf181).

### Installation

`sticcs` requires [`cyvcf2`](https://github.com/brentp/cyvcf2) and [`numpy`](https://numpy.org/). If these are not already installed, they will be downloaded and installed when you run the install command.

If you would like to export tree sequence objects from [`tskit`](https://tskit.dev/tskit/docs/stable/introduction.html), you will also need to install `tskit` yourself before running `sticcs`. Then install `sticcs`:

```bash
git clone https://github.com/simonhmartin/sticcs.git

cd sticcs

pip install -e .
```

### Command line tool

The command line tool takes as input a modified vcf file that contains a `DC` field, giving the count of derived alleles for each individual at each site.

You can make this from your standard vcf by running:
```bash
sticcs prep -i <input vcf> -o <output vcf>  --outgroup <outgroup sample ID>
```

If your vcf file already has the ancestral allele (provided in the `AA` field in the `INFO` section), then you do not need to specifiy outrgoups for polarising.

Now you can run the main command to make the tree sequence:

```bash
sticcs ts -i <input vcf> -o <output prefix> --output_format tskit
```
This will make a treesequence file that can be loaded and analysed using [tskit](https://tskit.dev/tskit/docs/stable/introduction.html). The default for `--output_format` is `newick`, which makes a file of newick trees and a separate file giving the chromosme coordinates of each tree interval.

### Python API

Classes and functions from `sticcs` can be used by importing `sticcs` in your python script. Full documention is not yet available, but some example can be seen in the [`twisst2 code`](https://github.com/simonhmartin/twisst2/blob/main/twisst2/twisst2.py).
