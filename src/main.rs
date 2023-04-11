use rand::{
    distributions::Distribution,
    distributions::{Uniform, WeightedIndex},
    prelude::*,
};
use rand_chacha::ChaCha8Rng;
use range_set_blaze::prelude::*;
use std::{ops::RangeInclusive, thread, time::Duration, time::Instant};

#[derive(Debug, Default, Clone, PartialEq)]
enum CellState {
    #[default]
    GPhase,
    SPhase,
}

#[derive(Debug, Default, Clone)]
struct Cell {
    genome_range: RangeSetBlaze<usize>,
    genome_length: usize,
    max_forks: usize,
    cell_state: CellState,
    replication_rate: usize,
    replication_state: RangeSetBlaze<usize>,
}
impl Cell {
    fn new(genome_length: usize, num_replicators: usize, replication_rate: usize) -> Self {
        Cell {
            genome_range: RangeSetBlaze::from_iter([0..=genome_length]),
            genome_length,
            max_forks: num_replicators,
            cell_state: CellState::GPhase,
            replication_rate,
            replication_state: RangeSetBlaze::<usize>::new(),
        }
    }
    fn is_replicated(&self, position: usize) -> bool {
        self.replication_state.contains(position)
    }
    fn is_fully_replicated(&self) -> bool {
        let total_length: usize = self
            .replication_state
            .ranges()
            .map(|r| r.end() - r.start())
            .sum();
        total_length == self.genome_length
    }
    fn new_random_position(&self, rng_obj: &mut ChaCha8Rng) -> isize {
        let compliment_ranges = self.genome_range.clone() - self.replication_state.clone();
        let range_uniforms: Vec<Uniform<usize>> = compliment_ranges
            .ranges()
            .map(|r| Uniform::new(r.start(), r.end() + 1))
            .collect();
        let region_lengths: Vec<usize> = compliment_ranges
            .ranges()
            .map(|r| r.end() + 1 - r.start())
            .collect();
        let mut new_initiation_pos: isize = -1;
        match WeightedIndex::new(&region_lengths) {
            Ok(valid_dist) => {
                while new_initiation_pos < 0 {
                    let samp_range = range_uniforms[valid_dist.sample(rng_obj)];
                    if rng_obj.gen::<f64>() > 0.9 {
                        new_initiation_pos = samp_range.sample(rng_obj) as isize;
                    };
                }
                new_initiation_pos
            }
            Err(_err) => new_initiation_pos,
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
            // Work out how many forks are current assigned, based on replication_state
            let mut assigned_forks: usize = self.replication_state.ranges().len();
            if assigned_forks > 1 {
                if self.replication_state.contains(0) {
                    assigned_forks -= 1;
                }
                if self.replication_state.contains(self.genome_length) {
                    assigned_forks -= 1;
                }
            }

            // Assign all the necessary un-assigned forks
            // while assigned_forks < (self.max_forks - 1) {
            for _ind in 0..(self.max_forks - assigned_forks - 1) {
                // Get new random position
                let sample_out = self.new_random_position(&mut rng);
                if sample_out < 0 {
                    break;
                }
                let new_pos = sample_out as usize;
                // Assign two new forks
                self.replication_state |= RangeSetBlaze::from_iter([
                    new_pos..=(new_pos + self.replication_rate),
                    (new_pos - self.replication_rate)..=new_pos,
                ]);
                assigned_forks += 2
            }

            // Move all forks
            self.replication_state |= self
                .replication_state
                .ranges()
                .map(|range| {
                    let lower = (range.start() - self.replication_rate).max(0);
                    let upper = (range.end() + self.replication_rate).min(self.genome_length);
                    lower..=upper
                })
                .collect::<RangeSetBlaze<usize>>();
            // println!("{:?}", self.replication_state);
            // println!(" ");
            // thread::sleep(Duration::from_millis(200));
            // Update iterations
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
    // let chrom_size: usize = 100_000;
    // let num_replicators: usize = 12;
    let mut cell = Cell::new(chrom_size, num_replicators, 50);

    // Basic checking
    println!("{:}", cell.is_replicated(100_000));
    println!("{:}", cell.is_fully_replicated());

    // Run replication
    cell.full_replication(0.9);
}

#[cfg(test)]
mod tests {}
