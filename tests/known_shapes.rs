use std::collections::HashSet;

use nalgebra as na;
use ordered_float::NotNan;

use quickhull::ConvexHull;

macro_rules! test_point_clouds {
    ($($name:ident: ($points:expr, $expected:expr),)*) => {
        $(
            #[test]
            fn $name() {
                let points = $points;
                let hull = ConvexHull::from_points(&points).unwrap();
                let vertices = hull.vertices();
                assert!(vertices.len() == 12);

                let expected:HashSet<_> = $expected.iter().map(|p| na::Point3::new(
                    NotNan::new(p.x).unwrap(),
                    NotNan::new(p.y).unwrap(),
                    NotNan::new(p.z).unwrap(),
                )).collect();
                let output:HashSet<_> = vertices.iter().map(|(p, _)| na::Point3::new(
                    NotNan::new(p.x).unwrap(),
                    NotNan::new(p.y).unwrap(),
                    NotNan::new(p.z).unwrap(),
                )).collect();
                let diff:HashSet<_> = expected.difference(&output).collect();

                assert!(diff.is_empty());
            }
        )*
    }
}

test_point_clouds! {
    hollow_tetrahedron: (vec![
        na::Point3::new(0.0, 0.0, 1.0),
        na::Point3::new(0.0, 1.0, 0.0),
        na::Point3::new(-1.0, 0.0, 0.0),
        na::Point3::new(1.0, 0.0, 0.0),
    ],
    vec![
        na::Point3::new(0.0, 0.0, 1.0),
        na::Point3::new(0.0, 1.0, 0.0),
        na::Point3::new(-1.0, 0.0, 0.0),
        na::Point3::new(1.0, 0.0, 0.0),
    ]),

    non_hollow_tetrahedron: (vec![
        na::Point3::new(0.0, 0.0, 1.0),
        na::Point3::new(0.0, 1.0, 0.0),
        na::Point3::new(-1.0, 0.0, 0.0),
        na::Point3::new(1.0, 0.0, 0.0),

        na::Point3::new(0.1, 0.1, 0.1),
    ],
    vec![
        na::Point3::new(0.0, 0.0, 1.0),
        na::Point3::new(0.0, 1.0, 0.0),
        na::Point3::new(-1.0, 0.0, 0.0),
        na::Point3::new(1.0, 0.0, 0.0),
    ]),
}
