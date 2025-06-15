Tool name : Identification and Phasing de novo DNA variation

![image](https://github.com/user-attachments/assets/f165de0f-93d5-498c-b98c-686bb7d79c81)



Intro:

de novo mutations are genetic changes that occur for the first time in an individual and are not
inherited from their parents, potentially contributing to genetic disorders or diversity. Profiling de
novo mutations is vital in clinical contexts, understanding population dynamics and also sheds
light on developmental biology, DNA replication and maintenance. Delineating the
parent-of-origin (PoO) of de novo variants additionally reveals parental specific effects on
germline mutation rates. DeMinTR (pronounced dementor) is a tool to identify and phase de
novo mutations from WGS datasets for STRs, SNVs and Indels. Addressing STR variations
requires additional methods owing to their unique length polymorphism. The first step DeMinTR
identifies de novo TR by considering genotype likelihoods as input and outputs estimated
posterior probability of mutation occurring at each TR loci in child. Secondly, Phasing of de
novo variants is based on two methods, allele sharing and read tracing. In allele sharing we
infer the PoO of the de novo allele where the inherited allele has only a single possibility of
origin. In read tracing, PoO is determined by building a haplotype with proximal phased
heterozygous SNV using read/read-pairs. We phased de novo variants identified from 33 Trio
samples of the 1000 genome project. DeMinTR results overlap more than 99% with existing tool
MonSTR for STRs and adding up the Read tracing method we are able to phase 20% additional
de novo variants. DeMinTR determines the PoO using two methods which uncovers the most
cases where existing tools fail to phase. MonSTR limits to only allele sharing methods which
leads to drop in many cases. Unfazed, a tool overlooks nearby heterozygous SNVs while
phasing. Therefore we see the requirement of developing a methodology which improves the
sensitivity of results. We are also further focusing on phasing SVs.
