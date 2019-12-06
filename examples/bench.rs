use nalgebra as na;
use rand::{rngs::StdRng, Rng, SeedableRng as _};

use quickhull::ConvexHull;

fn main() {
    let mut rng = StdRng::seed_from_u64(1234567890);
    let size = 1000.0;
    let points: Vec<_> = (0..5_000_000)
        .map(|_| {
            na::Point3::new(
                rng.gen_range(-size, size),
                rng.gen_range(-size, size),
                rng.gen_range(-size, size),
            )
        })
        .collect();

    for _ in 0..10 {
        ConvexHull::from_points(&points).unwrap();
    }
}
