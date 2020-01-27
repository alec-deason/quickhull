use nalgebra as na;
use rand::Rng;

use quickhull::ConvexHull;

fn main() {
    let mut rng = rand::thread_rng();
    let size = 1000.0;
    let points: Vec<_> = (0..5_000)
        .map(|_| {
            na::Point3::new(
                rng.gen_range(-size, size),
                rng.gen_range(-size, size),
                rng.gen_range(-size, size),
            )
        })
        .collect();

    let mut hull = ConvexHull::new();
    let mut count = 0;
    for _ in 0..5000 {
        count+= hull.from_points(&points).map(|p| p.len()).unwrap_or(0);
    }
    println!("Emitted a total of {} vertices", count);
}
