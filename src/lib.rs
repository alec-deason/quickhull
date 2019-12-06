use nalgebra as na;
use ordered_float::NotNan;

#[cfg(feature = "rayon")]
use rayon::prelude::*;

#[cfg(feature = "genmesh")]
mod genmesh;

#[derive(Clone)]
pub struct Triangle<T: na::Scalar> {
    a: na::Point3<T>,
    b: na::Point3<T>,
    c: na::Point3<T>,
    normal: na::Vector3<T>,
}

impl<T: na::RealField> Triangle<T> {
    fn new(a: na::Point3<T>, b: na::Point3<T>, c: na::Point3<T>) -> Self {
        let point1 = a - b;
        let point2 = b - c;
        let normal = point1.cross(&point2).normalize();

        Self { a, b, c, normal }
    }
}

pub struct ConvexHull<T: na::Scalar> {
    pub triangles: Vec<Triangle<T>>,
}

impl<T: na::RealField + num_traits::Float> ConvexHull<T> {
    pub fn from_points(points: &[na::Point3<T>]) -> Result<Self, ()> {
        if points.len() < 4 {
            return Err(());
        };
        let plane = initial_triangle(points);

        let triangles: Vec<Triangle<T>> = quick_hull(&plane, points, false)
            .into_iter()
            .chain(quick_hull(&plane, points, true))
            .collect();

        Ok(Self { triangles })
    }

    pub fn vertices(&self) -> Vec<(&na::Point3<T>, &na::Vector3<T>)> {
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
    pub fn mesh_generator(&self) -> crate::genmesh::ConvexHullMeshGenerator<'_, T> {
        crate::genmesh::ConvexHullMeshGenerator::new(self)
    }
}

fn quick_hull<T: na::RealField + num_traits::Float + num_traits::FromPrimitive>(plane: &Triangle<T>, points: &[na::Point3<T>], invert: bool) -> Vec<Triangle<T>> {
    let threshold = T::from_f64(10.0f64.powf(-10.0)).unwrap();
    let mut max_dist = T::zero();
    let mut max_point = None;
    let live_points: Vec<_> = points
        .iter()
        .filter(|p| {
            let mut d = plane.normal.dot(&(*p - &plane.a));
            if invert {
                d = -d;
            }
            if d > threshold {
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

        #[cfg(not(feature = "rayon"))]
        {
            [plane_b, plane_c, plane_d]
                .iter()
                .map(|plane| quick_hull(plane, &live_points, false))
                .flatten()
                .collect()
        }
        #[cfg(feature = "rayon")]
        {
            [plane_b, plane_c, plane_d]
                .par_iter()
                .map(|plane| quick_hull(plane, &live_points, false))
                .flatten()
                .collect()
        }
    } else {
        vec![plane.clone()]
    }
}

fn initial_triangle<T: na::RealField + num_traits::Float>(points: &[na::Point3<T>]) -> Triangle<T> {
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

    let mut max_dist = T::zero();
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

fn set_correct_normal<T: na::RealField + num_traits::Float>(plane: &mut Triangle<T>, points: &[na::Point3<T>]) {
    let threshold = T::from_f64(10.0f64.powf(-10.0)).unwrap();
    for p in points.iter() {
        let d = plane.normal.dot(&(p - plane.a));
        if d > threshold {
            plane.normal = -plane.normal;
            return;
        }
    }
}

fn distance_to_line<T: na::RealField + num_traits::Zero>(line: &[&na::Point3<T>; 2], point: &na::Point3<T>) -> T {
    let vec1 = point - line[0];
    let vec2 = point - line[1];
    let vec4 = vec1.cross(&vec2);

    if vec2.magnitude() <= T::zero() {
        T::zero()
    } else {
        vec4.magnitude() / vec2.magnitude()
    }
}
