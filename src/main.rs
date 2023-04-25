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
    genome_range: RangeSetBlaze<isize>,
    genome_length: isize,
    max_forks: isize,
    cell_state: CellState,
    replication_rate: isize,
    replication_state: RangeSetBlaze<isize>,
}
impl Cell {
    fn new(genome_length: isize, num_replicators: isize, replication_rate: isize) -> Self {
        Cell {
            genome_range: RangeSetBlaze::from_iter([0..=genome_length]),
            genome_length,
            max_forks: num_replicators,
            cell_state: CellState::GPhase,
            replication_rate,
            replication_state: RangeSetBlaze::<isize>::new(),
        }
    }
    fn is_replicated(&self, position: isize) -> bool {
        self.replication_state.contains(position)
    }
    fn is_fully_replicated(&self) -> bool {
        (self.replication_state.len() as isize) == (self.genome_length + 1)
    }
    fn new_random_range(&self, rng_obj: &mut ChaCha8Rng) -> Option<RangeInclusive<isize>> {
        let compliment_ranges = self.genome_range.clone() - self.replication_state.clone();
        let num_samples = (self.max_forks * 3) as usize;
        let mut range_uniforms: Vec<Uniform<isize>> = Vec::with_capacity(num_samples);
        let mut region_lengths: Vec<isize> = Vec::with_capacity(num_samples);
        for range in compliment_ranges.ranges() {
            range_uniforms.push(Uniform::new(range.start(), range.end() + 1));
            region_lengths.push(range.end() + 1 - range.start());
        }

        match WeightedIndex::new(&region_lengths) {
            Ok(valid_dist) => {
                let mut new_pos: isize = -1;
                while new_pos < 0 {
                    let samp_range = range_uniforms[valid_dist.sample(rng_obj)];
                    if rng_obj.gen::<f64>() > 0.9 {
                        new_pos = samp_range.sample(rng_obj);
                    };
                }
                Some((new_pos - self.replication_rate)..=(new_pos + self.replication_rate))
            }
            Err(_err) => None,
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
        let mut num_iterations: isize = 0;
        while !self.is_fully_replicated() {
            // Work out how many forks are current assigned, based on replication_state
            let mut assigned_forks: isize = self.replication_state.ranges().len() as isize;
            if assigned_forks > 1 {
                assigned_forks += (self.replication_state.contains(0) as isize)
                    + (self.replication_state.contains(self.genome_length) as isize);
            }

            // Assign all the necessary un-assigned forks
            self.replication_state |= (0..(self.max_forks - assigned_forks - 1))
                .filter_map(|_| self.new_random_range(&mut rng))
                .collect::<RangeSetBlaze<isize>>();

            // Move all forks
            self.replication_state = self
                .replication_state
                .ranges()
                .map(|range| {
                    let lower = (range.start() - self.replication_rate).max(0);
                    let upper = (range.end() + self.replication_rate).min(self.genome_length);
                    lower..=upper
                })
                .collect::<RangeSetBlaze<isize>>();

            if num_iterations > 20_000 {
                println!("{:?}", self.replication_state);
                println!(" ");
                thread::sleep(Duration::from_millis(20));
            }
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
    let chrom_size: isize = 500_000_000;
    let num_replicators: isize = chrom_size / 1_600_000;
    // let chrom_size: usize = 10_000;
    // let num_replicators: usize = 10;
    let mut cell = Cell::new(chrom_size, num_replicators, 50);

    // Basic checking
    println!("{:}", cell.is_replicated(100_000));
    println!("{:}", cell.is_fully_replicated());

    // Run replication
    cell.full_replication(0.9);
}

#[cfg(test)]
mod tests {}
