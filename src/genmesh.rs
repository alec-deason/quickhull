use genmesh::{
    generators::{IndexedPolygon, SharedVertex},
    Polygon, Triangle, Vertex,
};

use crate::ConvexHull;

pub struct ConvexHullMeshGenerator<'a> {
    parent: &'a ConvexHull,
    triangle: usize,
}

impl<'a> ConvexHullMeshGenerator<'a> {
    pub fn new(parent: &'a ConvexHull) -> Self {
        Self {
            parent,
            triangle: 0,
        }
    }
}

impl Iterator for ConvexHullMeshGenerator<'_> {
    type Item = Polygon<Vertex>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.triangle >= self.parent.triangles.len() {
            None
        } else {
            let triangle = &self.parent.triangles[self.triangle];
            self.triangle += 1;
            let a = Vertex {
                pos: [
                    triangle.a.x as f32,
                    triangle.a.y as f32,
                    triangle.a.z as f32,
                ]
                .into(),
                normal: [
                    triangle.normal.x as f32,
                    triangle.normal.y as f32,
                    triangle.normal.z as f32,
                ]
                .into(),
            };
            let b = Vertex {
                pos: [
                    triangle.b.x as f32,
                    triangle.b.y as f32,
                    triangle.b.z as f32,
                ]
                .into(),
                normal: [
                    triangle.normal.x as f32,
                    triangle.normal.y as f32,
                    triangle.normal.z as f32,
                ]
                .into(),
            };
            let c = Vertex {
                pos: [
                    triangle.c.x as f32,
                    triangle.c.y as f32,
                    triangle.c.z as f32,
                ]
                .into(),
                normal: [
                    triangle.normal.x as f32,
                    triangle.normal.y as f32,
                    triangle.normal.z as f32,
                ]
                .into(),
            };
            Some(Polygon::PolyTri(Triangle::new(a, b, c)))
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.parent.triangles.len() - self.triangle;
        (remaining, Some(remaining))
    }
}

impl SharedVertex<Vertex> for ConvexHullMeshGenerator<'_> {
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
            pos: [p.x as f32, p.y as f32, p.z as f32].into(),
            normal: [t.normal.x as f32, t.normal.y as f32, t.normal.z as f32].into(),
        }
    }

    fn shared_vertex_count(&self) -> usize {
        0
    }
}

impl IndexedPolygon<Triangle<usize>> for ConvexHullMeshGenerator<'_> {
    fn indexed_polygon(&self, idx: usize) -> Triangle<usize> {
        let idx = idx * 3;
        Triangle::new(idx, idx + 1, idx + 2)
    }

    fn indexed_polygon_count(&self) -> usize {
        self.parent.triangles.len()
    }
}
