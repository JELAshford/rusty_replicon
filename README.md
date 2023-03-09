<img src="src/silly_little_logo.png" alt="drawing" width="200"/>

# Rusty Replicon: Simulating Replication Timing from DNA Annotations

An implementation of "Replicon" (Gindin et. al. 2014) as a Rust programming and bioinformatics exercise, based on the <ins>[original](https://doi.org/10.1002/msb.134859)</ins> <ins>[works](https://doi.org/10.3389/fgene.2014.00378)</ins>. Follows the same strategy: predict replication timing (RT) for a full genome using only a probability of initiation at each position in the genome and the number of replication machineries as an input. 

## Implementation Details
Each cell simulation uses a very space-efficient representation of the replication-state, which is independant of the genome length. The representation stores the replicated state as a series of alternating replicated (`R`) and unreplicated (`U`) runs, much like a [run-length encoding (RLE)](https://en.wikipedia.org/wiki/Run-length_encoding). This allows for a fixed size represntation driven by the number of replication machineries (`M`) of size = `(M * 2) + 3`.

The alternating states in the representation allow the bi-directional progress of the fork to be done by unreplicated regions "giving" their bases to replicating regions on either side. Merging replicating regions with a shared adjacent empty unreplicated region conserves the number of replication forks implicitly, and requires only one edge-condition for 5' chromosome end. I have some ideas for other versions of this encoding that take up more space but may allow for fewer values being moved around such as: 
- extending the array to remove the 5' edge case,
- doubling the length of the array and starting in the middle, with heuristics to choose which way to shift values to minimise work,

but these are not high-priority optimisations. 

So far, the basic simulation is completed and I'm generating an initiation probability landscape (IPLS) for the MCF-7 data from [the Hansen paper](https://doi.org/10.1073/pnas.0912402107).

## To-Do
- [x] Basic replication tracking on a fixed size genome, with uniform IPLS
- [ ] Process ENCODE data: 
  - [x] Replication Timing for MCF-7
  - [ ] DNAse-seq for MCF-7
- [ ] Read DNAse IPLS into Rust
- [ ] Multi-cell asynchronous run
- [ ] Flow-gate sorting simulation
- [ ] Comparison to reference RT data
