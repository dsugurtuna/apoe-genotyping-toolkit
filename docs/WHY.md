# Why it's built this way

## The problem

APOE e2/e3/e4 status comes from two SNPs, rs429358 and rs7412, and is one of the most requested genotypes in dementia research. The mapping is simple, but getting it wrong is easy: an assumed allele, a flipped strand or a mis-parsed file silently turns e4 carriers into non-carriers.

## Design choices

**Why read the counted allele from the PLINK column name?** Because PLINK `--recode A` counts its A1 allele, which is usually the minor allele (C at rs429358, T at rs7412) but not always. The column is named `rs429358_C` or `rs429358_T`, so the code can check instead of assume.

**Why refuse A/G calls instead of flipping them?** Because both SNPs are C/T, and A/G usually means the data are on the opposite strand. Flipping is correct only if you know the strand; guessing wrong swaps e3 and e4. Stopping with a clear message is safer than a plausible wrong answer.

**Why accept genotype strings and 0/1/2 counts in CSV?** Because both turn up in practice, and the example data shipped with the repository used counts while the reader only understood strings. The counted alleles are stated in the docstring and README.

**Why call C/T at both SNPs e2/e4?** Because without phase the double heterozygote could be e2/e4 or e1/e3, and e1 is extremely rare. That is the usual convention, and the code says so.

**Why a Wilson interval for feasibility?** Because target genotypes such as e4/e4 are uncommon, and the textbook normal interval can go below zero or be far too narrow for small proportions. The interval is only meaningful when the cohort is treated as a sample of a wider population; for the cohort itself, the count is exact.

**Why is recall selection deterministic?** Because a recall list may need to be regenerated and audited. Taking candidates in file order within each age band gives the same answer every time; a random seed is on the roadmap.

**Why is mypy strict in CI?** Because it would have caught the bug that made `apoe-toolkit feasibility` crash on every call: the CLI passed keyword arguments the function did not have.

## Questions worth asking

**"Would you use this output to tell someone their APOE status?"**
No. It is a research genotype from array or imputed data. Returning APOE status to individuals needs an accredited test, genetic counselling and consent for that purpose. This tool supports counting and selecting participants for studies that have those arrangements.

**"How do you know the calls are right?"**
The unit tests pin the lookup table and the allele handling. In real use you would check call rates for both SNPs, compare genotype frequencies with published frequencies for the population, and, where available, compare with directly genotyped samples. rs429358 is missing from some arrays and can impute poorly, which is worth checking first.

**"What if the stratification cannot fill every age band?"**
It selects what is available and reports the shortfall; it does not borrow from neighbouring bands, because that would change the age profile the study asked for. The summary shows available versus selected for each arm.

## What's next

- Read VCF directly as well as PLINK output.
- Add an optional random seed to recall selection.
- Report call rates for rs429358 and rs7412 with every run.
