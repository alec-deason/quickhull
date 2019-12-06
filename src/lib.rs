use nalgebra as na;
use ordered_float::NotNan;

#[cfg(feature = "genmesh")]
mod genmesh;

type Point = na::Point3<f64>;
#[derive(Clone)]
pub struct Triangle {
    a: Point,
    b: Point,
    c: Point,
    normal: na::Vector3<f64>,
}

impl Triangle {
    fn new(a: Point, b: Point, c: Point) -> Self {
        let point1 = a - b;
        let point2 = b - c;
        let normal = point1.cross(&point2).normalize();

        Self { a, b, c, normal }
    }
}

pub struct ConvexHull {
    pub triangles: Vec<Triangle>,
}

impl ConvexHull {
    pub fn from_points(points: &[Point]) -> Result<Self, ()> {
        if points.len() < 4 {
            return Err(());
        };
        let plane = initial_triangle(points);

        let triangles: Vec<Triangle> = quick_hull(&plane, points, false)
            .into_iter()
            .chain(quick_hull(&plane, points, true))
            .collect();

        Ok(Self { triangles })
    }

    pub fn vertices(&self) -> Vec<(&Point, &na::Vector3<f64>)> {
        let mut vertices = Vec::with_capacity(self.triangles.len() * 3);
        for triangle in &self.triangles {
            vertices.extend(&[
                (&triangle.a, &triangle.normal),
                (&triangle.b, &triangle.normal),
                (&triangle.c, &triangle.normal),
            ]);
        }
        vertices
    }

    #[cfg(feature = "genmesh")]
    pub fn mesh_generator(&self) -> crate::genmesh::ConvexHullMeshGenerator<'_> {
        crate::genmesh::ConvexHullMeshGenerator::new(self)
    }
}

fn quick_hull(plane: &Triangle, points: &[Point], invert: bool) -> Vec<Triangle> {
    let mut max_dist = 0.0;
    let mut max_point = None;
    let live_points: Vec<_> = points
        .iter()
        .filter(|p| {
            let mut d = plane.normal.dot(&(*p - &plane.a));
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
        let mut plane_b = Triangle::new(plane.a, plane.b, max_point);
        set_correct_normal(&mut plane_b, &possible_internal_points);
        let mut plane_c = Triangle::new(plane.b, plane.c, max_point);
        set_correct_normal(&mut plane_c, &possible_internal_points);
        let mut plane_d = Triangle::new(plane.c, plane.a, max_point);
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

fn initial_triangle(points: &[Point]) -> Triangle {
    let mut x_max = &points[0];
    let mut x_min = &points[0];
    let mut y_max = &points[0];
    let mut y_min = &points[0];
    let mut z_max = &points[0];
    let mut z_min = &points[0];
    for p in points.iter() {
        if p.x > x_max.x {
            x_max = p;
        } else if p.x < x_min.x {
            x_min = p;
        }
        if p.y > y_max.y {
            y_max = p;
        } else if p.y < x_min.y {
            y_min = p;
        }
        if p.z > z_max.z {
            z_max = p;
        } else if p.z < z_min.z {
            z_min = p;
        }
    }

    let extreams = [x_max, y_max, z_max, x_min, y_min, z_min];

    let mut max_dist = 0.0;
    let mut line = [x_max, x_max];
    for (i, a) in extreams.iter().enumerate() {
        for b in &extreams[i + 1..] {
            let d = na::distance(*a, *b);
            if d > max_dist {
                max_dist = d;
                line = [*a, *b];
            }
        }
    }
    let third_point = *points
        .iter()
        .max_by_key(|p| NotNan::new(distance_to_line(&line, p)).unwrap())
        .unwrap();

    Triangle::new(*line[0], *line[1], third_point)
}

fn set_correct_normal(plane: &mut Triangle, points: &[Point]) {
    for p in points.iter() {
        let d = plane.normal.dot(&(p - plane.a));
        if d > 10.0f64.powf(-10.0) {
            plane.normal *= -1.0;
            return;
        }
    }
}

fn distance_to_line(line: &[&Point; 2], point: &Point) -> f64 {
    let vec1 = point - line[0];
    let vec2 = point - line[1];
    let vec4 = vec1.cross(&vec2);

    if vec2.magnitude() <= 0.0 {
        0.0
    } else {
        vec4.magnitude() / vec2.magnitude()
    }
}
