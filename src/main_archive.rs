use anyhow::Result;
use std::collections::VecDeque;

#[derive(Debug, Default, Clone, Copy)]
struct Range {
    replicated: bool,
    length: isize,
}

#[derive(Debug, Default)]
pub struct Chromosome {
    genome_length: isize,
    replication_state: VecDeque<Range>,
}
impl Chromosome {
    fn new(genome_length: isize, num_origins: usize) -> Self {
        let mut start_vec = VecDeque::with_capacity((num_origins * 2) + 1);
        start_vec.push_front(Range {
            replicated: false,
            length: genome_length,
        });
        Chromosome {
            genome_length,
            replication_state: start_vec,
        }
    }
    fn replicate(mut self, position: isize) -> Self {
        // Handle out of bounds
        if position < 0 {
            panic!("Index {position} is too small, cannot use negative indexes")
        }
        if position >= self.genome_length {
            panic!(
                "Index {} is too large, cannot index beyond genome length {}",
                position, self.genome_length
            )
        }

        // Identify correct insertion location
        let mut insert_index: usize = 0;
        let mut cumsum: isize = 0;
        for (ind, range) in self.replication_state.iter().enumerate() {
            insert_index = ind;
            cumsum += range.length;
            if position < cumsum {
                break;
            }
        }

        // Get current bin state
        let current_length: isize = self.replication_state[insert_index].length;

        // If already replicated, leave!
        if self.replication_state[insert_index].replicated {
            return self;
        };

        // Otherwise, do a full split:
        let left_count = position - ((cumsum - 1) - current_length) - 1; // = position - (cumsum - current_length - 1) - 1 = position - cumsum + current_length + 1 - 1 = positions - cumsum + current_length
        let right_count = (cumsum - 1) - position;
        self.replication_state.remove(insert_index);
        self.replication_state.insert(
            insert_index,
            Range {
                replicated: true,
                length: 1,
            },
        );
        if right_count > 0 {
            self.replication_state.insert(
                (insert_index + 1).min(self.replication_state.len()),
                Range {
                    replicated: false,
                    length: right_count,
                },
            )
        }
        if left_count > 0 {
            self.replication_state.insert(
                insert_index.max(0),
                Range {
                    replicated: false,
                    length: left_count,
                },
            )
        }

        // Merge width adjacent
        let new_data_loc = insert_index + (left_count > 0) as usize;
        if new_data_loc < (self.replication_state.len() - 1)
            && self.replication_state[new_data_loc].replicated
                == self.replication_state[new_data_loc + 1].replicated
        {
            self.replication_state[new_data_loc].length +=
                self.replication_state[new_data_loc].length;
            self.replication_state.remove(new_data_loc + 1);
        };
        if new_data_loc > 0
            && self.replication_state[new_data_loc].replicated
                == self.replication_state[new_data_loc - 1].replicated
        {
            self.replication_state[new_data_loc].length += 1;
            self.replication_state.remove(new_data_loc - 1);
        };

        // Bring it all back!
        self
    }

    fn is_replicated(&self, position: isize) -> bool {
        // Handle out of bounds
        if position < 0 {
            panic!("Index {position} is too small, cannot use negative indexes")
        }
        if position >= self.genome_length {
            panic!(
                "Index {} is too large, cannot index beyond genome length {}",
                position, self.genome_length
            )
        }

        // Identify correct insertion location
        let mut check_index: usize = 0;
        let mut cumsum: isize = 0;
        for (ind, range) in self.replication_state.iter().enumerate() {
            check_index = ind;
            cumsum += range.length;
            if position < cumsum {
                break;
            }
        }
        // Return replication state
        self.replication_state[check_index].replicated
    }
}

fn main() -> Result<()> {
    // Create a prototype genome
    let mut chromosome = Chromosome::new(500, 10);
    chromosome = chromosome.replicate(100);
    println!("{:?}", &chromosome);
    println!("{:}", chromosome.is_replicated(100));
    // We done!
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::Chromosome;
    #[test]
    fn replicate_single_pos() {
        let mut chromosome = Chromosome::new(500, 10);
        chromosome = chromosome.replicate(100);
        assert!(chromosome.is_replicated(100));
    }

    #[test]
    #[should_panic]
    fn out_of_bounds_low() {
        let chromosome = Chromosome::new(500, 10);
        chromosome.replicate(-2);
    }

    #[test]
    #[should_panic]
    fn out_of_bounds_high() {
        let chromosome = Chromosome::new(500, 10);
        chromosome.replicate(510);
    }
}
