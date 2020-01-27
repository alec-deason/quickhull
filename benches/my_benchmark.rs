use criterion::{black_box, criterion_group, criterion_main, Criterion};

use nalgebra as na;
use rand::{rngs::StdRng, Rng, SeedableRng as _};
use quickhull::ConvexHull;


fn criterion_benchmark(c: &mut Criterion) {
    let size = 1000.0;
    let mut rng = StdRng::seed_from_u64(1234567892);
    let points: Vec<_> = (0..500)
    .map(|_| {
        na::Point3::new(
            rng.gen_range(-size, size),
            rng.gen_range(-size, size),
            rng.gen_range(-size, size),
        )
    })
    .collect();
    let mut hull = ConvexHull::new();
    c.bench_function("from_points", |b| b.iter(|| {
        hull.from_points(black_box(&points));
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
