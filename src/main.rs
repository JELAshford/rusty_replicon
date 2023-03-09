use rand::{
    distributions::Distribution,
    distributions::{Uniform, WeightedIndex},
    prelude::*,
};
use rand_chacha::ChaCha8Rng;
use std::time::Instant;

#[derive(Debug, Default, Clone, PartialEq)]
enum CellState {
    #[default]
    GPhase,
    SPhase,
}

#[derive(Debug, Default, Clone)]
struct Cell {
    genome_length: usize,
    num_replicators: usize,
    unassigned_replicators: usize,
    cell_state: CellState,
    replication_rate: usize,
    replication_state: Vec<usize>,
}
impl Cell {
    fn new(genome_length: usize, num_replicators: usize, replication_rate: usize) -> Self {
        let mut start_vec: Vec<usize> = vec![0; (num_replicators * 2) + 3];
        start_vec[1] = genome_length;
        Cell {
            genome_length,
            num_replicators,
            unassigned_replicators: num_replicators,
            cell_state: CellState::GPhase,
            replication_rate,
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
    fn assign_replicators(&mut self, rng_obj: &mut ChaCha8Rng) {
        // If there are unassigned replicators, assign them
        while self.unassigned_replicators > 0 {
            // Setup RNG
            let mut cumsum: usize = 0;
            let mut random_ranges: Vec<Uniform<usize>> = Vec::with_capacity(self.num_replicators);
            let mut region_lengths: Vec<usize> = Vec::with_capacity(self.num_replicators);
            for (ind, length) in self.replication_state.iter().enumerate() {
                let ilength = *length as isize;
                if (ind % 2 != 0) & (ilength != 0) {
                    let unif_range = Uniform::new(cumsum, cumsum + length);
                    region_lengths.push(*length);
                    random_ranges.push(unif_range);
                }
                cumsum += length;
            }
            // Choose a randomly generated point, weighted between ranges by range length
            let mut new_initiation_pos: isize = -1;
            match WeightedIndex::new(&region_lengths) {
                Ok(valid_dist) => {
                    while new_initiation_pos < 0 {
                        let samp_range = random_ranges[valid_dist.sample(rng_obj)];
                        let sample_pos = samp_range.sample(rng_obj);
                        if rng_obj.gen::<f64>() > 0.9 {
                            new_initiation_pos = sample_pos as isize;
                        };
                    }
                }
                Err(_err) => return, // no more places to choose
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

            // Update number of repliactors
            self.unassigned_replicators -= 1;
        }
    }
    fn replicate_and_merge(&mut self) {
        let num_entries = self.replication_state.len();

        for index in (1..(num_entries - 1)).step_by(2).rev() {
            // At each unreplicated region, give one of the values to
            // adjacent occupied replication regions
            let left_occupied = self.replication_state[index - 1] > 0;
            let right_occupied = self.replication_state[index + 1] > 0;

            if self.replication_state[index] > 0 {
                if left_occupied {
                    let move_amount = self.replication_state[index].min(self.replication_rate);
                    self.replication_state[index - 1] += move_amount;
                    self.replication_state[index] -= move_amount;
                }
                if (right_occupied) && (self.replication_state[index] > 0) {
                    let move_amount = self.replication_state[index].min(self.replication_rate);
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
                self.unassigned_replicators += 1;
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
    }
    fn full_replication(&mut self, g_phase_prob: f64) {
        let mut rng = ChaCha8Rng::seed_from_u64(1701);

        // Loop until enters G-phase
        let mut num_warmup_iters: isize = 0;
        while self.cell_state == CellState::GPhase {
            if rng.gen::<f64>() > g_phase_prob {
                self.cell_state = CellState::SPhase;
            }
            num_warmup_iters += 1;
        }
        println!("Entered S phase after {num_warmup_iters:?} warmups!");

        // Replication run
        let now = Instant::now();
        let mut num_iterations: usize = 0;
        while !self.is_fully_replicated() {
            self.assign_replicators(&mut rng);
            self.replicate_and_merge();
            num_iterations += 1;
        }
        println!("Time taken: {:.2?}", now.elapsed());
        println!(
            "Converged in {} iterations to: {:?}",
            &num_iterations, &self.replication_state
        );
    }
}

fn main() {
    // Create a prototype genome
    let chrom_size: usize = 500_000_000;
    let num_replicators: usize = chrom_size / 1_600_000;
    let mut cell = Cell::new(chrom_size, num_replicators, 50);

    // Basic checking
    println!("{:}", cell.is_replicated(100_000));
    println!("{:}", cell.is_fully_replicated());

    // Run replication
    cell.full_replication(0.9);
}

#[cfg(test)]
mod tests {}
