#!/bin/bash

# to concatenate and index fa files

cat mouse.fa Pchabaudi_ens.fa > mousePchabaudi.fa
cat mouse.fa Pyoelii_ens.fa > mousePyoelii.fa
cat monkey.fa Pcoatneyi_ens.fa > monkeyPcoatneyi.fa
cat monkey.fa Pvivax_ens.fa > monkeyPvivax.fa
cat mouse.fa Pberghei_ens.fa > mousePberghei.fa

cat human.fa Pfalciparum_ens.fa > humanPfalciparum.fa
cat human.fa Pberghei_ens.fa > humanPberghei.fa
cat human.fa Pvivax_ens.fa > humanPvivax.fa


# index

bwa index mousePchabaudi.fa
bwa index mousePyoelii.fa
bwa index monkeyPcoatneyi.fa
bwa index monkeyPvivax.fa
bwa index mousePberghei.fa

bwa index humanPfalciparum.fa
bwa index humanPberghei.fa
bwa index humanPvivax.fa

