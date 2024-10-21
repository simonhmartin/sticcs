# sticcs

`sticcs` is a method for inferring the series of genealogies along the genome, also called a [tree sequence](https://tskit.dev/tutorials/what_is.html). Unlike some other methods, `sticcs` **does not require phased haplotypes**, and it can work on **any ploidy level**.

The input for `sticcs` is polarised genotype data. This means you need to know the ancestral allele at each site, or you need an appropriate outgroup(s) to allow inference of the derived allele.

A publication is in preparation, but you can check out the [How it works](#How_it_works) section below to learn about the apporach.

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

### How it works

The approach will be described in a publication, but for now here is an overview.

`sticcs` stands for Sequential Tree Inference by Collecting Compatible Sites. It uses a simple heuristic approach and works in two steps:
1. Find collections of compatible variant sites
2. Build the tree for each collection using Perfect Phylogeny

To understand the approach, it makes sense to first look at the tree building step.

#### Tree building using Perfect Phylogeny

Assuming an infinite-sites mutation model (no recurrent mutation), each variant site tells us about a 'node' in the tree. All the individuals that carry the derived allele must share a common ancestor at that node. Given a set of biallelic variants that are **all mutually compatible**, there is always at least one bifurcating tree compatible with the set of variants. This concept is called Perfect Phylogeny.

##### Haploid Perfect Phylogeny

Say we have six haploid individuals and 4 SNPs:

|        |  `A`  |  `B`  |  `C`  |  `D`  |  `E`  |  `F`  |
| :----: | ----- | ----- | ----- | ----- | ----- | ----- |
| `SNP1` |**`1`**|**`1`**|  `0`  |  `0`  |  `0`  |  `0`  |
| `SNP2` |  `0`  |  `0`  |  `0`  |**`1`**|**`1`**|  `0`  |
| `SNP3` |**`1`**|**`1`**|**`1`**|  `0`  |  `0`  |  `0`  |
| `SNP4` |**`1`**|**`1`**|**`1`**|**`1`**|**`1`**|  `0`  |

The SNPs are polarised, so that 0 indicates the ancestral state and **`1`** indicates the derived state.

`SNP1` tells us that the tree includes a clade containing individuals `A` and `B`. `SNP2` telss us that the tree includes a clade containing individuals `D` and `E`, and so on. Using this logic, the inferred tree for this cluster of sites is:

`(((A,B),C),(D,E)),F);`

If we were lacking `SNP1`, for example, we would still be able to infer the tree, but it would contain a *polytomy* for individuals `A`, `B` and `C`, but there is still only ever one tree compatible with our set of SNPs:

`((A,B,C),(D,E)),F);`


##### Diploid Perfect Phylogeny

Now, if instead of six phased haplotypes we have three unphased diploid genotypes.:

|        |   `A`   |   `B`   |   `C`   |
| :----: | ------- | ------- | ------- |
| `SNP1` |`1/1`|`0/0`|  `0/0`  |
| `SNP2` |  `0/0`  |  `0/1`  |  `0/1`  |
| `SNP3` |`1/1`|`0/1`|`0/0`|
| `SNP4` |`1/1`|`1/1`|`0/1`|

We can write these out in terms of the number of derived alleles they each carry:

|        |  `A`  |  `B`  |  `C`  |
| :----: | ----- | ----- | ----- |
| `SNP1` |`2`|`0`|  `0`  |
| `SNP2` |  `0`  |  `1`  |  `1`  |
| `SNP3` |`2`|`1`|`0`|
| `SNP4` |`2`|`2`|`1`|

Once again we can use the same logic to infer the tree.

`SNP1` tells us that the tree includes a clade containing *both* tips from individual `A`.
`SNP2` tells us that the tree includes a clade containing one tip from individual `B` and one from `C`.
`SNP3` tells us that the *other* tip from `B` groups with the two tips from `A`.
`SNP4` tells us that one of the tips from `C` is basal to all the rest.

So our inferred tree is:

`(((A1,A2),B2),(B1,C1)),C2);`

Here the two tips from each individual are **arbitrary**. We don't know which of the two haplotypes from individual B is grouping with the the two from individual A, but we know that one of them is.

You might have noticed that some SNPs in this diploid example are compatible with more than one possible clade. For example, considering `SNP2` in isolation, it might tell us about a clade of `(B1,C1)`, or of `(B2, C1)` etc. In this case, we arbitrarily decided on `(B1,C1)`. This means that when we got to `SNP3`, we only have `B2` left to create the clade `((A1,A2),B2)`. `sticcs` considers the tips from a single individual as interchangeable.

The `sticcs` algorithm works through SNP patterns in order from the smallest to the largest clades (i.e. working from leaves to root). Once all SNPs have been considered, it runs through all SNP patterns again, mopping up any remaining available tips into clades, until there are no more loose tips that can be mopped up.

When building a tree sequence (from left-to-right along the chromosme), `sticcs` retains any compatible clades that were present in the previous tree. Thus, although the decisions to cluster say tips `B1` and `C1` was arbitrary the first time `SNP2` was encountered, this information is retained if `SNP2` is also present next tree. In this way, `sticcs` does capture some information about short range haplotype phase, but it is not a replacement for a phasing tool.

But how exactly does `sticcs` decide which SNPs are included in which local genealogy. This takes us back to the initial step of collecting compatible sites.

#### Collecting of compatible sites

Our goal is to infer genealogies along the chromosome, and the recombination points at which they change. To do this using perfect phylogeny, we want to make collections of compatible variants. However, this is not as simple as just dividing the genome into blocks of compatible sites, for two reasons.

First, compatibility is assesed in a pairwise manner (see next section). If for example we have three adjacent SNPs, and SNP2 is compatible with both SNP1 and SNP3, but SNP1 and SNP3 are not compatible with each other, we need to decide whether the recombination point is between SNP1 and SNP2 or between SNP2 and SNP3. Secondly, adjacent genealogies tend to be very similar - usually only differing by one recombinatyion event and therefore one branch on the tree. This means that, for a focal genealogy with a given interval, there will usually be variant sites outside of that interval that are informative about the focal genealogy, and we want to include those in our inference.

`sticcs` addresses these challenges by using a heuristic clustering algorithm with the following steps:
* First each variant is labeled in terms of its 'pattern': a list giving the number of derived alleles observed in each individual
* For each variant site *i* on the chromosome, we compile two sets of compatible patterns. The left-compatible set is the set of patterns for all varaint sites from 1 to *i* that are compatible with *i*, and where a site carrying that pattern is not separated from *i* by a variant with an incompatible pattern. The right-compatible set is similar, but for all variants between *i* and the end of the chromosome.
* Cluster *i* is initialised as the set of patterns shared by both the left- and right-compatible sets for site *i*
* Additional patterns from both the left- and right-compatible sets are then added to cluster *i*, starting from the pattern associated with the nearest variant site to site *i* on the chromosome and moving outwards in both directions, provided this pattern is compatible with all existing patterns in cluster *i*.
* Any adjacent groups of sites that have identical clusters are then merged into a single cluster, and each cluster has an interval defined by the first and last of the sites that were merged.
* Chromosome intervals are then extended to the midpoints between the clusters (and to chromosome start and end for terminal clusters) such that there are no gaps.


##### Defining compatible pairs of variants

Whether any pair of variants is compatible with the same genealogy can be tested using the [Four-genete test](https://en.wikipedia.org/wiki/Four-gamete_test). If the dataset contains all four pairs of alleles at the two sites:
`--0--0--`
`--0--1--`
`--1--0--`
`--1--1--`

These sites are not compartible with the same genealogy. For polarised data, in which `1` represents the derived state, we know that the ancestral haplotype `--0--0--` does (or did) exist, so we only need to observe the latter three derived haplotypes to know that the sites are incompatible.

##### Testing for compatibility when using diploids and polyploids

For unphased diploid or polyploid genotypes, we have less information, but we can still test for compatibility under the most parsimonius case. For example:

If we saw these three diploid genotypes in three individuals:

`--0/1--0/0--`
`--0/0--0/1--`
`--0/1--0/1--`

We can infer that individual 1 carries haplotypes `--0--0--` and `--1--0--`, while individual 2 carries `--0--0--` and `--0--1--`. Thus, the most parsimonius haplotypes for individual 3 are `--1--0--` and `--0--1--`. Hence, we conclude here that only two of the three possible derived haplotypes are present:
`--0--1--`
`--1--0--`

So these two sites are deemed compatible.

As a general rule that applies to any ploidy level, we can say that all three derived haplotypes must be present if:
* There is at least one individual at which the number of derived alleles at the first site exceeds that at the second (implies `--1--0--` exists)
* There is at least one individual at which the number of derived alleles at the second site exceeds that at the first (implies `--0--1--` exists)
* There is at least one individual at which the sum of derived alleles across the two sites exceeds the ploidy (implies `--1--1--` exists)

