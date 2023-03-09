use std::time::Instant;

use anyhow::Result;
use rand::{
    distributions::Distribution,
    distributions::{Uniform, WeightedIndex},
    prelude::*,
};
use rand_chacha::ChaCha8Rng;

#[derive(Debug, PartialEq)]
enum CellState {
    GPhase,
    SPhase,
}

#[derive(Debug, Default, Clone)]
pub struct Genome {
    genome_length: usize,
    num_origins: usize,
    replication_state: Vec<usize>,
}
impl Genome {
    fn new(genome_length: usize, num_origins: usize) -> Self {
        let mut start_vec: Vec<usize> = vec![0; (num_origins * 2) + 3];
        start_vec[1] = genome_length;
        Genome {
            genome_length,
            num_origins,
            replication_state: start_vec,
        }
    }
    fn is_replicated(&self, position: usize) -> bool {
        // Handle out of bounds
        if position >= self.genome_length {
            panic!(
                "Index {} is too large, cannot index beyond genome length {}",
                position, self.genome_length
            )
        }
        // Identify correct insertion location
        let mut check_index: usize = 0;
        let mut cumsum: usize = 0;
        for (ind, value) in self.replication_state.iter().enumerate() {
            check_index = ind;
            cumsum += value;
            if position < cumsum {
                break;
            }
        }
        // All even indexes are replicated ranges
        check_index % 2 == 0
    }
    fn is_fully_replicated(&self) -> bool {
        // genome is fully replicated if there's no positions in unreplicated (odd) storage indexes
        for (ind, val) in self.replication_state.iter().enumerate() {
            if (ind % 2 != 0) & (*val != 0) {
                return false;
            }
        }
        true
    }
    fn random_unreplicated_point(&self, rng_obj: &mut ChaCha8Rng) -> Result<usize, anyhow::Error> {
        // Generate a random point from each unreplicated range, and store the region lengths
        let mut cumsum: usize = 0;
        let mut random_vals: Vec<usize> = vec![];
        let mut region_lengths: Vec<usize> = vec![];
        for (ind, length) in self.replication_state.iter().enumerate() {
            let ilength = *length as isize;
            if (ind % 2 != 0) & (ilength != 0) {
                region_lengths.push(*length);
                random_vals.push(Uniform::new(cumsum, cumsum + length).sample(rng_obj));
            }
            cumsum += length;
        }
        // Choose one of these random points, weighted by region length
        let dist = WeightedIndex::new(&region_lengths)?;
        Ok(random_vals[dist.sample(rng_obj)])
    }
    fn assign_replicators(&mut self, num_replicators: &usize, rng_obj: &mut ChaCha8Rng) {
        // If there are unassigned replicators, assign them
        for _rep in 0..*num_replicators {
            // Sample random position to replicate
            let mut new_initiation_pos: isize = -1;
            while new_initiation_pos < 0 {
                let sampling_out = self.random_unreplicated_point(rng_obj);
                match sampling_out {
                    Ok(sample_pos) => {
                        let assign_prob: f64 = rng_obj.gen();
                        if assign_prob > 0.9 {
                            new_initiation_pos = sample_pos as isize;
                        }
                    }
                    Err(_err) => return,
                };
            }
            let position: usize = new_initiation_pos as usize;

            // Identify insertion location
            let mut insert_index: usize = 0;
            let mut cumsum: usize = 0;
            for (ind, length) in self.replication_state.iter().enumerate() {
                insert_index = ind;
                cumsum += length;
                if position < cumsum {
                    break;
                }
            }
            // Get current bin state and work out adjacent values
            let current_length = self.replication_state[insert_index];
            let left_count = position + current_length - cumsum;
            let right_count = (cumsum - 1) - position;
            // Move all values forward 2 positions until 2 after current
            for index in ((insert_index + 2)..self.replication_state.len()).rev() {
                self.replication_state[index] = self.replication_state[index - 2];
            }
            // Insert the new values
            self.replication_state[insert_index + 2] = right_count;
            self.replication_state[insert_index + 1] = 1;
            self.replication_state[insert_index] = left_count;
        }
    }
    fn replicate_and_merge(&mut self, step_size: usize) -> usize {
        let num_entries = self.replication_state.len();
        let mut num_merged: usize = 0;

        for index in (1..(num_entries - 1)).step_by(2).rev() {
            // At each unreplicated region, give one of the values to
            // adjacent occupied replication regions
            let left_occupied = self.replication_state[index - 1] > 0;
            let right_occupied = self.replication_state[index + 1] > 0;

            if self.replication_state[index] > 0 {
                if left_occupied {
                    let move_amount = self.replication_state[index].min(step_size);
                    self.replication_state[index - 1] += move_amount;
                    self.replication_state[index] -= move_amount;
                }
                if (right_occupied) && (self.replication_state[index] > 0) {
                    let move_amount = self.replication_state[index].min(step_size);
                    self.replication_state[index + 1] += move_amount;
                    self.replication_state[index] -= move_amount;
                }
            }

            // Merge if now 0 and both neighbours are occupied
            if (self.replication_state[index] == 0) && left_occupied && right_occupied {
                // Update left by addding right, then shift all rest
                self.replication_state[index - 1] += self.replication_state[index + 1];
                for step_index in index..(num_entries - 2) {
                    self.replication_state[step_index] = self.replication_state[step_index + 2];
                }
                self.replication_state[&num_entries - 2] = 0;
                self.replication_state[&num_entries - 1] = 0;

                // Count the merge
                num_merged += 1;
            }
        }
        // Edge case for merging genome start
        if (self.replication_state[0] == 0) && (self.replication_state[1] == 0) {
            for step_index in 0..(num_entries - 2) {
                self.replication_state[step_index] = self.replication_state[step_index + 2];
            }
            self.replication_state[&num_entries - 2] = 0;
            self.replication_state[&num_entries - 1] = 0;
        }
        num_merged
    }
}

#[derive(Debug)]
struct SingleChromCell {
    num_replicators: usize,
    replication_rate: usize,
    cell_state: CellState,
    genome: Genome,
}
impl SingleChromCell {
    fn new(genome_length: usize, replication_rate: usize, num_replicators: usize) -> Self {
        SingleChromCell {
            num_replicators,
            replication_rate,
            cell_state: CellState::GPhase,
            genome: Genome::new(genome_length, num_replicators),
        }
    }
    fn run_replication(mut self, g_phase_prob: f64) {
        let mut rng = ChaCha8Rng::seed_from_u64(1701);

        // Loop until enters G-phase
        let mut num_warmup_iters: isize = 0;
        while self.cell_state == CellState::GPhase {
            let test_val: f64 = rng.gen();
            if test_val > g_phase_prob {
                self.cell_state = CellState::SPhase;
            }
            num_warmup_iters += 1;
        }
        println!("Entered S phase after {num_warmup_iters:?} warmups!");

        // Replication run
        let now = Instant::now();
        let mut unassigned_replicators = self.num_replicators;
        let mut num_iterations: usize = 0;
        while !self.genome.is_fully_replicated() {
            // Assign unassigned replicators
            self.genome
                .assign_replicators(&unassigned_replicators, &mut rng);
            unassigned_replicators = 0;
            // Carry out replication and merge steps
            let num_merged = self.genome.replicate_and_merge(self.replication_rate);
            unassigned_replicators += num_merged;
            // Metric store
            num_iterations += 1;
        }

        println!(
            "Converged in {} iterations to: {:?}",
            &num_iterations, &self.genome.replication_state
        );
        println!("Time taken: {:.2?}", now.elapsed());
    }
}

fn main() -> Result<()> {
    // Create a prototype genome
    let chrom_size: usize = 500_000_000;
    let num_replicators: usize = chrom_size / 1_600_000;
    let cell = SingleChromCell::new(chrom_size, 50, num_replicators);

    // Basic checking
    println!("{:}", cell.genome.is_replicated(100_000));
    println!("{:}", cell.genome.is_fully_replicated());

    // Run replication
    cell.run_replication(0.9);

    // We done!
    Ok(())
}

#[cfg(test)]
mod tests {}
