use nalgebra as na;

use genmesh::{
    generators::{IndexedPolygon, SharedVertex},
    Polygon, Triangle, Vertex,
};

use crate::ConvexHull;

pub struct ConvexHullMeshGenerator<'a, T: na::RealField> {
    parent: &'a ConvexHull<T>,
    triangle: usize,
}

impl<'a, T: na::RealField> ConvexHullMeshGenerator<'a, T> {
    pub fn new(parent: &'a ConvexHull<T>) -> Self {
        Self {
            parent,
            triangle: 0,
        }
    }
}

pub trait ToArray32: na::RealField {
    fn to_array_32(vec: &na::Matrix3x1<Self>) -> [f32; 3];
}

impl ToArray32 for f64 {
    fn to_array_32(vec: &na::Matrix3x1<Self>) -> [f32; 3] {
        [vec.x as f32, vec.y as f32, vec.z as f32]
    }
}

impl ToArray32 for f32 {
    fn to_array_32(vec: &na::Matrix3x1<Self>) -> [f32; 3] {
        [vec.x, vec.y, vec.z]
    }
}

impl<T: na::RealField + ToArray32> Iterator for ConvexHullMeshGenerator<'_, T> {
    type Item = Polygon<Vertex>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.triangle >= self.parent.triangles.len() {
            None
        } else {
            let triangle = &self.parent.triangles[self.triangle];
            self.triangle += 1;
            let a = Vertex {
                pos: T::to_array_32(&triangle.a.coords).into(),
                normal: T::to_array_32(&triangle.normal).into(),
            };
            let b = Vertex {
                pos: T::to_array_32(&triangle.b.coords).into(),
                normal: T::to_array_32(&triangle.normal).into(),
            };
            let c = Vertex {
                pos: T::to_array_32(&triangle.c.coords).into(),
                normal: T::to_array_32(&triangle.normal).into(),
            };
            Some(Polygon::PolyTri(Triangle::new(a, b, c)))
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.parent.triangles.len() - self.triangle;
        (remaining, Some(remaining))
    }
}

impl<T: na::RealField + Into<f32>> SharedVertex<Vertex> for ConvexHullMeshGenerator<'_, T> {
    fn shared_vertex(&self, idx: usize) -> Vertex {
        let triangle = idx / 3;
        let vert = triangle % 3;
        let t = &self.parent.triangles[triangle];
        let p = match vert {
            0 => t.a,
            1 => t.b,
            _ => t.c,
        };
        Vertex {
            pos: [p.x.into(), p.y.into(), p.z.into()].into(),
            normal: [t.normal.x.into(), t.normal.y.into(), t.normal.z.into()].into(),
        }
    }

    fn shared_vertex_count(&self) -> usize {
        0
    }
}

impl<T: na::RealField> IndexedPolygon<Triangle<usize>> for ConvexHullMeshGenerator<'_, T> {
    fn indexed_polygon(&self, idx: usize) -> Triangle<usize> {
        let idx = idx * 3;
        Triangle::new(idx, idx + 1, idx + 2)
    }

    fn indexed_polygon_count(&self) -> usize {
        self.parent.triangles.len()
    }
}
