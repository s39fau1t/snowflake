use std::f64::consts::PI;

use macroquad::prelude::*;

const BRANCH_LENGTH: f64 = 1.0;
const EPSILON: f64 = 0.000001;

#[derive(Clone)]
struct Branch {
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
}

impl Branch {
    fn new(x1: f64, y1: f64, x2: f64, y2: f64) -> Self {
        Branch { x1, y1, x2, y2 }
    }
}

fn approx_eq(x: f64, y: f64) -> bool {
    (x - y).abs() < EPSILON
}

fn approx_le(lhs: f64, rhs: f64) -> bool {
    lhs < rhs || approx_eq(lhs, rhs)
}

fn approx_ge(lhs: f64, rhs: f64) -> bool {
    lhs > rhs || approx_eq(lhs, rhs)
}

fn approx_in_range_inclusive(x: f64, low: f64, high: f64) -> bool {
    assert!(low <= high);

    approx_ge(x, low) && approx_le(x, high)
}

//TODO: check distance from line to point
fn point_on_branch(x: f64, y: f64, b: &Branch) -> bool {
    let l = b.x1.min(b.x2);
    let r = b.x1.max(b.x2);
    let o = b.y1.min(b.y2);
    let t = b.y1.max(b.y2);

    let lhs = (b.x2 - b.x1) * (y - b.y1);
    let rhs = (x - b.x1) * (b.y2 - b.y1);

    approx_eq(lhs, rhs) && approx_in_range_inclusive(x, l, r) && approx_in_range_inclusive(y, o, t)
}

//TODO: clean up
fn colinear_intersecting(b1: &Branch, b2: &Branch) -> bool {
    let l1 = b1.x1.min(b1.x2);
    let r1 = b1.x1.max(b1.x2);
    let l2 = b2.x1.min(b2.x2);
    let r2 = b2.x1.max(b2.x2);
    let o1 = b1.y1.min(b1.y2);
    let t1 = b1.y1.max(b1.y2);
    let o2 = b2.y1.min(b2.y2);
    let t2 = b2.y1.max(b2.y2);

    approx_le(l1, r2) && approx_le(l2, r1) && approx_le(o1, t2) && approx_le(o2, t1)
}

fn intersect(b1: &Branch, b2: &Branch) -> Option<(f64, f64)> {
    let denom = (b2.y2 - b2.y1) * (b1.x2 - b1.x1) - (b2.x2 - b2.x1) * (b1.y2 - b1.y1);

    if approx_eq(denom, 0.0) {
        let t_num = (b2.x2 - b2.x1) * (b1.y1 - b2.y1) - (b2.y2 - b2.y1) * (b1.x1 - b2.x1);

        if approx_eq(t_num, 0.0) && colinear_intersecting(&b1, &b2) {
            let denom = (b1.x2 - b1.x1) - (b2.x2 - b2.x1);

            //vertical lines check
            let t = if !approx_eq(denom, 0.0) {
                (b2.x1 - b1.x1) / denom
            } else {
                (b2.y1 - b1.y1) / ((b1.y2 - b1.y1) - (b2.y2 - b2.y1))
            };

            //this works because all lines are the same length
            Some((t, t))
        } else {
            None
        }
    } else {
        let t = ((b2.x2 - b2.x1) * (b1.y1 - b2.y1) - (b2.y2 - b2.y1) * (b1.x1 - b2.x1)) / denom;
        let u = ((b1.x2 - b1.x1) * (b1.y1 - b2.y1) - (b1.y2 - b1.y1) * (b1.x1 - b2.x1)) / denom;

        if approx_in_range_inclusive(t, 0.0, 1.0) && approx_in_range_inclusive(u, 0.0, 1.0) {
            Some((t, u))
        } else {
            None
        }
    }
}

fn resolve(grew: &mut [Branch], all_branches: &[Branch]) -> Vec<Branch> {
    let mut intersections = Vec::new();

    for i in 0..grew.len() {
        for j in (i + 1)..grew.len() {
            let [b1, b2] = grew.get_disjoint_mut([i, j]).unwrap();

            if approx_eq(b1.x1, b2.x1) && approx_eq(b1.y1, b2.y1) {
                continue;
            }

            if let Some((t, u)) = intersect(b1, b2) {
                let ix = b1.x1 + t * (b1.x2 - b1.x1);
                let iy = b1.y1 + t * (b1.y2 - b1.y1);

                intersections.push((i, Some(j), ix, iy, t.max(u)));
            }
        }
    }

    for (i, b1) in grew.iter_mut().enumerate() {
        for b2 in all_branches {
            if approx_eq(b2.x2, b1.x1) && approx_eq(b2.y2, b1.y1) {
                continue;
            }

            if let Some((t, u)) = intersect(b1, b2) {
                let ix = b1.x1 + t * (b1.x2 - b1.x1);
                let iy = b1.y1 + t * (b1.y2 - b1.y1);

                intersections.push((i, None, ix, iy, t.max(u)));
            }
        }
    }

    let mut has_intersected = vec![false; grew.len()];

    intersections
        .sort_unstable_by(|(_, _, _, _, ta), (_, _, _, _, tb)| ta.partial_cmp(tb).unwrap());

    for (i, j, ix, iy, _) in intersections {
        if j.is_some() {
            let j = j.unwrap();
            let [b1, b2] = grew.get_disjoint_mut([i, j]).unwrap();

            if point_on_branch(ix, iy, b1) && point_on_branch(ix, iy, b2) {
                b1.x2 = ix;
                b2.x2 = ix;
                b1.y2 = iy;
                b2.y2 = iy;

                has_intersected[i] = true;
                has_intersected[j] = true;
            }
        } else {
            let b = &mut grew[i];

            if point_on_branch(ix, iy, b) {
                b.x2 = ix;
                b.y2 = iy;

                has_intersected[i] = true;
            }
        }
    }

    has_intersected
        .iter()
        .zip(grew)
        .filter_map(|(intersect, b)| if !intersect { Some(b.clone()) } else { None })
        .collect()
}

fn grow(
    branches: Vec<Branch>,
    all_branches: &mut Vec<Branch>,
    max_iters: i32,
    i: i32,
) -> Vec<Branch> {
    if i == max_iters {
        return branches;
    }

    println!(
        "----------{}----------; {}, {}",
        i + 1,
        branches.len(),
        all_branches.len()
    );

    let mut grew = Vec::new();

    for branch in branches {
        let angle = (branch.y2 - branch.y1).atan2(branch.x2 - branch.x1);
        let left = angle + PI / 4.0;
        let right = angle - PI / 4.0;

        grew.push(Branch::new(
            branch.x2,
            branch.y2,
            branch.x2 + left.cos() * BRANCH_LENGTH,
            branch.y2 + left.sin() * BRANCH_LENGTH,
        ));
        grew.push(Branch::new(
            branch.x2,
            branch.y2,
            branch.x2 + right.cos() * BRANCH_LENGTH,
            branch.y2 + right.sin() * BRANCH_LENGTH,
        ));
    }

    let to_grow = resolve(&mut grew, all_branches);

    all_branches.append(&mut grew);

    grow(to_grow, all_branches, max_iters, i + 1)
}

#[macroquad::main("Snowflake")]
async fn main() {
    let mut branches = vec![
        Branch::new(0.0, 0.0, BRANCH_LENGTH, 0.0),
        Branch::new(0.0, 0.0, -BRANCH_LENGTH, 0.0),
        Branch::new(0.0, 0.0, 0.0, BRANCH_LENGTH),
        Branch::new(0.0, 0.0, 0.0, -BRANCH_LENGTH),
    ];

    println!("generating...");

    let mut to_grow = grow(branches.clone(), &mut branches, 200, 0);

    println!("final: {}", branches.len());

    let mut transx = 0.0;
    let mut transy = 0.0;
    let mut scale = 1.0;

    loop {
        if is_key_down(KeyCode::Equal) {
            scale *= 1.02;
        }
        if is_key_down(KeyCode::Minus) {
            scale /= 1.02;
        }
        if is_key_down(KeyCode::Left) {
            transx += 6.0 / scale;
        }
        if is_key_down(KeyCode::Right) {
            transx -= 6.0 / scale;
        }
        if is_key_down(KeyCode::Up) {
            transy -= 6.0 / scale;
        }
        if is_key_down(KeyCode::Down) {
            transy += 6.0 / scale;
        }

        if is_mouse_button_pressed(MouseButton::Left) {
            to_grow = grow(to_grow, &mut branches, 1, 0);
        }

        for branch in &branches {
            draw_line(
                ((branch.x1 + transx) * scale) as f32 + screen_width() / 2.0,
                ((branch.y1 + transy) * scale) as f32 * -1.0 + screen_height() / 2.0,
                ((branch.x2 + transx) * scale) as f32 + screen_width() / 2.0,
                ((branch.y2 + transy) * scale) as f32 * -1.0 + screen_height() / 2.0,
                1.0,
                WHITE,
            );
        }

        next_frame().await
    }
}
