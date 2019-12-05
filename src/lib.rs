use nalgebra as na;
use ordered_float::NotNan;

#[cfg(feature = "genmesh")]
mod genmesh;

type Point = na::Point3<f64>;
#[derive(Clone)]
pub struct Plane {
    a: Point,
    b: Point,
    c: Point,
    normal: na::Vector3<f64>,
}

impl Plane {
    fn new(a: Point, b: Point, c: Point) -> Self {
        let point1 = a - b;
        let point2 = b - c;
        let normal = point1.cross(&point2).normalize();

        Self { a, b, c, normal }
    }
}

pub struct ConvexHull {
    triangles: Vec<Plane>,
}

impl ConvexHull {
    pub fn quickhull(points: &[Point]) -> Result<Self, ()> {
        if points.len() < 4 {
            return Err(());
        };
        let plane = initial_triangle(points);

        let triangles: Vec<Plane> = quick_hull(&plane, points, false)
            .into_iter()
            .chain(quick_hull(&plane, points, true))
            .collect();

        Ok(Self { triangles })
    }

    #[cfg(feature = "genmesh")]
    pub fn mesh_generator(&self) -> crate::genmesh::ConvexHullMeshGenerator<'_> {
        crate::genmesh::ConvexHullMeshGenerator::new(self)
    }
}

fn quick_hull(plane: &Plane, points: &[Point], invert: bool) -> Vec<Plane> {
    let mut max_dist = 0.0;
    let mut max_point = None;
    let live_points: Vec<_> = points
        .iter()
        .filter(|p| {
            let mut d = signed_distance_to_plane(plane, p);
            if invert {
                d *= -1.0;
            }
            if d > 10.0f64.powf(-10.0) {
                if d > max_dist {
                    max_dist = d;
                    max_point = Some(*p);
                }
                true
            } else {
                false
            }
        })
        .cloned()
        .collect();

    if let Some(max_point) = max_point {
        let max_point = *max_point;
        let possible_internal_points = [plane.a, plane.b, plane.c, max_point];
        let mut plane_b = Plane::new(plane.a, plane.b, max_point);
        set_correct_normal(&mut plane_b, &possible_internal_points);
        let mut plane_c = Plane::new(plane.b, plane.c, max_point);
        set_correct_normal(&mut plane_c, &possible_internal_points);
        let mut plane_d = Plane::new(plane.c, plane.a, max_point);
        set_correct_normal(&mut plane_d, &possible_internal_points);
        quick_hull(&plane_b, &live_points, false)
            .into_iter()
            .chain(quick_hull(&plane_c, &live_points, false))
            .chain(quick_hull(&plane_d, &live_points, false))
            .collect()
    } else {
        vec![plane.clone()]
    }
}

fn initial_triangle(points: &[Point]) -> Plane {
    //FIXME: be less lazy and do this in a single pass
    let x_max = *points
        .iter()
        .max_by_key(|p| NotNan::new(p.x).unwrap())
        .unwrap();
    let y_max = *points
        .iter()
        .max_by_key(|p| NotNan::new(p.y).unwrap())
        .unwrap();
    let z_max = *points
        .iter()
        .max_by_key(|p| NotNan::new(p.z).unwrap())
        .unwrap();
    let x_min = *points
        .iter()
        .min_by_key(|p| NotNan::new(p.x).unwrap())
        .unwrap();
    let y_min = *points
        .iter()
        .min_by_key(|p| NotNan::new(p.y).unwrap())
        .unwrap();
    let z_min = *points
        .iter()
        .min_by_key(|p| NotNan::new(p.z).unwrap())
        .unwrap();
    let extreams = [x_max, y_max, z_max, x_min, y_min, z_min];

    let mut max_dist = 0.0;
    let mut line = [x_max, x_max];
    for (i, a) in extreams.iter().enumerate() {
        for b in &extreams[i + 1..] {
            let d = na::distance(a, b);
            if d > max_dist {
                max_dist = d;
                line = [*a, *b];
            }
        }
    }
    let third_point = *points
        .iter()
        .max_by_key(|p| NotNan::new(distance_to_line(&line, p).unwrap_or(0.0)).unwrap())
        .unwrap();

    Plane::new(line[0], line[1], third_point)
}

fn set_correct_normal(plane: &mut Plane, points: &[Point]) {
    for p in points.iter() {
        let d = plane.normal.dot(&(p - plane.a));
        if d > 10.0f64.powf(-10.0) {
            plane.normal *= -1.0;
            return;
        }
    }
}

fn distance_to_line(line: &[Point; 2], point: &Point) -> Result<f64, ()> {
    let vec1 = point - line[0];
    let vec2 = point - line[1];
    let vec4 = vec1.cross(&vec2);

    if vec2.magnitude() <= 0.0 {
        Err(())
    } else {
        Ok(vec4.magnitude() / vec2.magnitude())
    }
}

fn signed_distance_to_plane(plane: &Plane, point: &Point) -> f64 {
    plane.normal.dot(&(point - plane.a))
}
