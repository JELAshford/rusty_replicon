use rand::{
    distributions::Distribution,
    distributions::{Uniform, WeightedIndex},
    prelude::*,
};
use rand_chacha::ChaCha8Rng;
use std::time::Instant;
use std::{thread, time};

#[derive(Debug, Clone)]
struct Fork {
    origin: isize,
    position: isize,
    step: isize,
    to_delete: bool,
}

struct Cell {
    genome_length: usize,
    num_replicators: usize,
    step_size: isize,
    fork_state: Vec<Fork>,
}
impl Cell {
    fn new(genome_length: usize, num_replicators: usize, step_size: isize) -> Self {
        Cell {
            genome_length,
            num_replicators,
            step_size,
            fork_state: Vec::with_capacity(num_replicators),
        }
    }
    fn fully_replicated(&mut self) -> bool {
        // If state is empty, not replicated and can't be summed
        if self.fork_state.is_empty() {
            return false;
        }
        // Otherwise, sum up and check against genome length
        self.genome_length
            == self
                .fork_state
                .iter()
                .map(|fork| (fork.origin - fork.position).unsigned_abs())
                .sum()
    }
    fn insert_fork_pair(&mut self, insertion_position: isize) {
        // calcualte insertion index of this value
        let mut insertion_index: usize = 0;
        if !self.fork_state.is_empty() {
            insertion_index = self.fork_state.len();
            for (ind, fork) in self.fork_state.iter().enumerate() {
                if insertion_position <= fork.position {
                    insertion_index = ind;
                    break;
                }
            }
        }
        // Do the insertion
        self.fork_state.insert(
            insertion_index,
            Fork {
                origin: insertion_position,
                position: (insertion_position + self.step_size)
                    .max(0)
                    .min(self.genome_length as isize),
                step: self.step_size,
                to_delete: false,
            },
        );
        self.fork_state.insert(
            insertion_index,
            Fork {
                origin: insertion_position,
                position: (insertion_position - self.step_size)
                    .max(0)
                    .min(self.genome_length as isize),
                step: -self.step_size,
                to_delete: false,
            },
        );
    }
    fn step_forks(&mut self) {
        for fork in self.fork_state.iter_mut() {
            fork.position += fork.step;
            fork.position = fork.position.max(0).min(self.genome_length as isize)
        }
    }
    fn merge_forks(&mut self) {
        // Iterate over the forks, check if overlapped
        let mut end_point_reached = false;
        while !end_point_reached {
            for ind in (0..self.fork_state.len() - 3).step_by(2) {
                if self.fork_state[ind + 1].position >= self.fork_state[ind + 2].position {
                    self.fork_state[ind].origin = self.fork_state[ind + 1].position;
                    self.fork_state[ind + 3].origin = self.fork_state[ind + 1].position;
                    self.fork_state[ind].position = self.fork_state[ind]
                        .position
                        .min(self.fork_state[ind + 2].position);
                    self.fork_state[ind + 3].position = self.fork_state[ind + 1]
                        .position
                        .max(self.fork_state[ind + 3].position);
                    self.fork_state.remove(ind + 1);
                    self.fork_state.remove(ind + 2);
                    break;
                }
            }
            end_point_reached = true;
        }
    }
    fn replenish_forks(&mut self, rng_obj: &mut ChaCha8Rng) {
        // Keep all forks in use by sampling a new position
        while self.fork_state.len() < self.num_replicators {
            // Grab the unreplicated ranges from the fork state
            let mut unreplicated_ranges: Vec<Uniform<isize>> =
                Vec::with_capacity(self.num_replicators);
            let mut unreplicated_range_lengths: Vec<usize> =
                Vec::with_capacity(self.num_replicators);
            if self.fork_state.is_empty() {
                unreplicated_ranges.push(Uniform::new(0, self.genome_length as isize));
                unreplicated_range_lengths.push(self.genome_length);
            } else {
                let mut cumsum: isize = 0;
                for fork in self.fork_state.iter() {
                    // if fork.step < 0 && fork.position != 0 {
                    //     println!("{:?}", &cumsum);
                    //     println!("{:?}", &fork.position);
                    //     let range_len: usize = (fork.position - cumsum) as usize;
                    //     if range_len > 1 {
                    //         println!("Inner loop");
                    //         println!("{:?}", &unreplicated_ranges);
                    //         println!("{:?}", &cumsum);
                    //         println!("{:?}", &fork.position);
                    //         unreplicated_ranges.push(Uniform::new(cumsum, fork.position));
                    //         unreplicated_range_lengths.push(range_len);
                    //     }
                    // }
                    let mut left_val = fork.origin;
                    let mut right_val = fork.position;
                    if fork.step < 0 {
                        left_val = fork.position;
                        right_val = fork.origin;
                    }
                    if left_val > cumsum {
                        let range_len: usize = (left_val - cumsum) as usize;
                        if range_len > 1 {
                            unreplicated_ranges.push(Uniform::new(cumsum, left_val));
                            unreplicated_range_lengths.push(range_len);
                        }
                    }
                    cumsum = right_val;
                }
                if cumsum != self.genome_length as isize {
                    unreplicated_ranges.push(Uniform::new(cumsum, self.genome_length as isize));
                    unreplicated_range_lengths.push(self.genome_length - cumsum as usize);
                }
            }
            println!("{:?}", &unreplicated_ranges);
            // Sample from these ranges weighted by their length
            let mut new_position: isize = -1;
            match WeightedIndex::new(&unreplicated_range_lengths) {
                Ok(valid_dist) => {
                    while new_position < 0 {
                        let samp_range = unreplicated_ranges[valid_dist.sample(rng_obj)];
                        let sample_pos = samp_range.sample(rng_obj);
                        if rng_obj.gen::<f64>() > 0.9 {
                            new_position = sample_pos;
                        };
                    }
                }
                Err(_err) => return, // no more places to choose
            } // Insert forks at this position
            println!("{:?}", &new_position);
            self.insert_fork_pair(new_position);
        }
    }
}

fn main() {
    let mut rng = ChaCha8Rng::seed_from_u64(700);
    let mut cell: Cell = Cell::new(100, 4, 2);
    let mut iteration = 0;
    let now = Instant::now();
    while !cell.fully_replicated() {
        println!("{iteration:?}");
        println!("{:?}", &cell.fork_state);
        cell.replenish_forks(&mut rng);
        println!("{:?}", &cell.fork_state);
        cell.step_forks();
        println!("{:?}", &cell.fork_state);
        cell.merge_forks();
        println!("{:?}", &cell.fork_state);
        thread::sleep(time::Duration::from_millis(100));
        iteration += 1;
    }
    println!("Time taken: {:.2?}", now.elapsed());
    println!("Finished in {iteration:?} iterations.");
    println!("Final state: {:?}", &cell.fork_state);
}

#[cfg(test)]
mod tests {}
