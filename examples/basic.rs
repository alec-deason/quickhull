use nalgebra as na;
use rand::Rng;

use quickhull::from_points;

fn main() {
    let mut rng = rand::thread_rng();
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

    let hull = from_points(&points).unwrap();
    println!("Hull has {} triangles", hull.len());
}
