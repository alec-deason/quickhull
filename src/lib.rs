use std::collections::{HashSet, HashMap};

use nalgebra as na;

const EPSILON: f32 = 0.0001;
type Point=na::Point3<f32>;
type Vector=na::Vector3<f32>;

pub fn from_points(points: &[Point]) -> Result<Vec<Point>, ()> {
    let mut triangles:HashMap<usize, (Triangle, Vec<Point>)> = {
        let mut cloud: Vec<_> = points.iter().cloned().collect();
        initial_tetrahedron(points).into_iter().enumerate().map(|(idx, t)| {
            let outside_set = extract_outside_set(&mut cloud, &t);
            (idx, (t, outside_set))
        }).collect()
    };
    let mut next_idx = triangles.keys().max().unwrap() + 1;

    loop {
        let mut to_replace = None;
        for (idx, (_, outside_set)) in triangles.iter() {
            if outside_set.len() == 0 {
                continue
            }
            to_replace = Some(*idx);
            break
        }
        if let Some(idx) = to_replace {
            let max_p = {
                let (triangle, outside_set) = triangles.get(&idx).unwrap();
                let mut max_d = std::f32::NEG_INFINITY;
                let mut max_p = outside_set[0];
                for p in outside_set.iter() {
                    let d = triangle.normal.dot(&(p - triangle.a));
                    if d > max_d {
                        max_d = d;
                        max_p = *p;
                    }
                }
                max_p
            };

            let (horizon, to_remove) = calculate_horizon(&mut triangles, idx, &max_p)?;
            let mut outside_set: Vec<Point> = to_remove.iter().map(|idx| triangles.get_mut(idx).unwrap().1.drain(..).collect::<Vec<Point>>()).flatten().collect();

            for idx in to_remove {
                triangles.remove(&idx);
            }
            let mut new_triangles = Vec::with_capacity(horizon.len());
            for (neighbor, edge) in horizon.iter() {
                let idx = next_idx;
                next_idx += 1;
                {
                    let (neighbor_triangle, _) = triangles.get_mut(&neighbor).unwrap();
                    if (na::distance(&neighbor_triangle.a, &edge.0) < 0.00001 && na::distance(&neighbor_triangle.b, &edge.1) < 0.00001) ||
                       (na::distance(&neighbor_triangle.b, &edge.0) < 0.00001 && na::distance(&neighbor_triangle.a, &edge.1) < 0.00001) {
                        neighbor_triangle.a_neighbor = idx;
                    } else  if (na::distance(&neighbor_triangle.b, &edge.0) < 0.00001 && na::distance(&neighbor_triangle.c, &edge.1) < 0.00001) ||
                               (na::distance(&neighbor_triangle.c, &edge.0) < 0.00001 && na::distance(&neighbor_triangle.b, &edge.1) < 0.00001) {
                        neighbor_triangle.b_neighbor = idx;
                    } else if (na::distance(&neighbor_triangle.c, &edge.0) < 0.00001 && na::distance(&neighbor_triangle.a, &edge.1) < 0.00001) ||
                              (na::distance(&neighbor_triangle.a, &edge.0) < 0.00001 && na::distance(&neighbor_triangle.c, &edge.1) < 0.00001) {
                        neighbor_triangle.c_neighbor = idx;
                    } else { panic!("wrong neighbor somehow?"); }
                }
                let mut triangle = Triangle::new(edge.0, edge.1, max_p, *neighbor, 0, 0);
                let mut hull_points = Vec::with_capacity(triangles.len()*3);
                //FIXME: I can certainly do this with an iterator and avoid the bounds checks
                for (t,_) in triangles.values() {
                    hull_points.push(t.a);
                    hull_points.push(t.b);
                    hull_points.push(t.c);
                }
                set_correct_normal(&mut triangle, &hull_points);
                let local_outside_set = extract_outside_set(&mut outside_set, &triangle);
                new_triangles.push((idx, (triangle, local_outside_set)));
            }
            for i in 0..new_triangles.len() {
                let (_, (triangle, _)) = &new_triangles[i];
                let mut b_neighbor = None;
                let mut c_neighbor = None;
                for (j, (other_idx, (other, _))) in new_triangles.iter().enumerate() {
                    if i==j {
                        continue
                    }
                    if na::distance(&triangle.a, &other.a) < 0.00001 || na::distance(&triangle.a, &other.b) < 0.00001 {
                        c_neighbor = Some(*other_idx);
                    }
                    if na::distance(&triangle.b, &other.b) < 0.00001 || na::distance(&triangle.b, &other.a) < 0.00001 {
                        b_neighbor = Some(*other_idx);
                    }
                }
                if let (Some(b_idx), Some(c_idx)) = (b_neighbor, c_neighbor) {
                    let (_, (triangle, _)) = &mut new_triangles[i];
                    triangle.b_neighbor = b_idx;
                    triangle.c_neighbor = c_idx;
                } else {
                    eprintln!("Something: {:#?} {:#?} {:?} {:?}", horizon, triangle, b_neighbor, c_neighbor);
                    eprintln!("FAIL");
                    return Err(());
                }
            }
            triangles.extend(new_triangles);
        } else {
            break
        }
    }

    let mut result = Vec::with_capacity(triangles.len()*3);
    //FIXME: I can certainly do this with an iterator and avoid the bounds checks
    for (_, (t,_)) in triangles {
        result.push(t.a);
        result.push(t.b);
        result.push(t.c);
    }
    eprintln!("SUCCESS");
    Ok(result)
}

fn calculate_horizon(triangles: &mut HashMap<usize, (Triangle, Vec<Point>)>, initial: usize, max_p: &Point) -> Result<(Vec<(usize, (Point, Point))>, Vec<usize>), ()> {
    let mut to_remove = Vec::with_capacity(triangles.len());

    let mut visited = HashSet::with_capacity(triangles.len());
    let mut to_visit = Vec::with_capacity(triangles.len());
    {
        let (triangle, _) = &triangles[&initial];
        to_visit.push((initial, (triangle.a, triangle.b)));
    }
    let mut horizon = Vec::with_capacity(triangles.len());
    while let Some((idx, edge)) = to_visit.pop() {
        let (triangle, _) = triangles.get(&idx).ok_or(())?;
        let d = triangle.normal.dot(&(*max_p - triangle.a));
        if d > EPSILON {
            if !visited.contains(&idx) {
                visited.insert(idx);
            } else {
                continue;
            }
            to_remove.push(idx);
            for (n, edge) in &[
                (triangle.a_neighbor, (triangle.a, triangle.b)),
                (triangle.b_neighbor, (triangle.b, triangle.c)),
                (triangle.c_neighbor, (triangle.c, triangle.a)),
            ] {
                to_visit.push((*n, *edge));
            }
        } else {
            horizon.push((idx, edge));
        }
    }

    to_remove.shrink_to_fit();
    horizon.shrink_to_fit();
    Ok((horizon, to_remove))
}

fn extract_outside_set(points: &mut Vec<Point>, triangle: &Triangle) -> Vec<Point> {
    let mut outside_set = Vec::with_capacity(points.len());
    points.retain(|p| {
        let d = triangle.normal.dot(&(*p - triangle.a));
        if d > EPSILON {
            outside_set.push(*p);
            false
        } else {
            true
        }
    });
    outside_set.shrink_to_fit();
    outside_set
}

#[derive(Clone, Debug)]
struct Triangle {
    a: Point,
    b: Point,
    c: Point,
    a_neighbor: usize,
    b_neighbor: usize,
    c_neighbor: usize,
    normal: Vector,
}

impl Triangle {
    fn new(a: Point, b: Point, c: Point, a_neighbor: usize, b_neighbor: usize, c_neighbor: usize) -> Self {
        let point1 = b - a;
        let point2 = c - a;
        Self {
            a, b, c, a_neighbor, b_neighbor, c_neighbor,
            normal: point1.cross(&point2).normalize(),
        }
    }
}

fn initial_tetrahedron(points: &[Point]) -> Vec<Triangle> {
    let mut extremas = extremas(points);
    assert!(extremas.len() >= 4, "degenerate hull");
    extremas.sort_by_key(|p| (p.z*1000.0) as i32);
    let mut other = extremas.split_off(3);
    let tip = other.pop().unwrap();
    extremas.sort_by_key(|p| (p.x*1000.0)as i32);
    let a = extremas.remove(0);
    let b = extremas.pop().unwrap();
    let c = extremas.pop().unwrap();
    let mut base = Triangle::new(a, b, c, 1, 2, 3);
    set_correct_normal(&mut base, &[tip]);
    let mut a_tri = Triangle::new(a, b, tip, 0, 2, 3);
    set_correct_normal(&mut a_tri, &[c]);
    let mut b_tri = Triangle::new(b, c, tip, 0, 3, 1);
    set_correct_normal(&mut b_tri, &[a]);
    let mut c_tri = Triangle::new(c, a, tip, 0, 1, 2);
    set_correct_normal(&mut c_tri, &[b]);

    vec![base, a_tri, b_tri, c_tri]
}

fn set_correct_normal(triangle: &mut Triangle, points: &[Point]) {
    for point in points {
        let d = triangle.normal.dot(&(point - triangle.a));
        if d > EPSILON {
            triangle.normal *= -1.0;
            break
        }
    }
}

fn extremas(points: &[Point]) -> Vec<Point> {
    let mut min_x = &points[0];
    let mut max_x = &points[0];
    let mut min_y = &points[0];
    let mut max_y = &points[0];
    let mut min_z = &points[0];
    let mut max_z = &points[0];
    for p in points {
        if p.x < min_x.x {
            min_x = p;
        }
        if p.x > max_x.x {
            max_x = p;
        }
        if p.y < min_y.y {
            min_y = p;
        }
        if p.y > max_y.y {
            max_y = p;
        }
        if p.z < min_z.z {
            min_z = p;
        }
        if p.z > max_z.z {
            max_z = p;
        }
    }

    let mut extremas = Vec::with_capacity(6);
    for p in &[min_x, max_x, min_y, max_y, min_z, max_z] {
        if let None = extremas.iter().find(|pp| p==pp) {
            extremas.push(**p);
        }
    }

    extremas.shrink_to_fit();
    extremas
}
