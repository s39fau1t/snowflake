use std::f64::consts::PI;

use macroquad::prelude::*;

const BRANCH_LENGTH: f64 = 20.0;
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

fn colinear_intersecting(b1: &Branch, b2: &Branch) -> bool {
    let l1 = b1.x1.min(b1.x2);
    let r1 = b1.x1.max(b1.x2);
    let l2 = b2.x1.min(b2.x2);
    let r2 = b2.x1.max(b2.x2);
    let o1 = b1.y1.min(b1.y2);
    let t1 = b1.y1.max(b1.y2);
    let o2 = b2.y1.min(b2.y2);
    let t2 = b2.y1.max(b2.y2);

    l1 <= r2 && l2 <= r1 && o1 <= t2 && o2 <= t1
}

fn grow(branches: &[Branch], all_branches: &[Branch], max_iters: i32, iters: i32) -> Vec<Branch> {
    if iters == max_iters {
        return Vec::new();
    }

    println!(
        "----------{}----------; {}, {}",
        iters + 1,
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

    for i in 0..grew.len() {
        let (b1, others) = grew[i..].split_first_mut().unwrap();

        for b2 in others {
            let denom = (b2.y2 - b2.y1) * (b1.x2 - b1.x1) - (b2.x2 - b2.x1) * (b1.y2 - b1.y1);

            if approx_eq(denom, 0.0) {
                let t_num = (b2.x2 - b2.x1) * (b1.y1 - b2.y1) - (b2.y2 - b2.y1) * (b1.x1 - b2.x1);
                if approx_eq(t_num, 0.0) && colinear_intersecting(b1, b2) {
                    let mx = (b1.x1 + b2.x1) / 2.0;
                    let my = (b1.y1 + b2.y1) / 2.0;

                    b1.x2 = mx;
                    b2.x2 = mx;
                    b1.y2 = my;
                    b2.y2 = my;
                }
            } else {
                let t =
                    ((b2.x2 - b2.x1) * (b1.y1 - b2.y1) - (b2.y2 - b2.y1) * (b1.x1 - b2.x1)) / denom;
                let u =
                    ((b1.x2 - b1.x1) * (b1.y1 - b2.y1) - (b1.y2 - b1.y1) * (b1.x1 - b2.x1)) / denom;

                if EPSILON < t && t < 1.0 - EPSILON && EPSILON < u && u < 1.0 - EPSILON {
                    let ix = b1.x1 + t * (b1.x2 - b1.x1);
                    let iy = b1.y1 + t * (b1.y2 - b1.y1);

                    b1.x2 = ix;
                    b2.x2 = ix;
                    b1.y2 = iy;
                    b2.y2 = iy;
                }
            }
        }
    }

    for b1 in &mut grew {
        for b2 in all_branches {
            let denom = (b2.y2 - b2.y1) * (b1.x2 - b1.x1) - (b2.x2 - b2.x1) * (b1.y2 - b1.y1);

            if approx_eq(denom, 0.0) {
                let t_num = (b2.x2 - b2.x1) * (b1.y1 - b2.y1) - (b2.y2 - b2.y1) * (b1.x1 - b2.x1);
                if approx_eq(t_num, 0.0) && colinear_intersecting(b1, b2) {
                    let mx = (b1.x1 + b2.x1) / 2.0;
                    let my = (b1.y1 + b2.y1) / 2.0;

                    b1.x2 = mx;
                    b1.y2 = my;
                }
            } else {
                let t =
                    ((b2.x2 - b2.x1) * (b1.y1 - b2.y1) - (b2.y2 - b2.y1) * (b1.x1 - b2.x1)) / denom;
                let u =
                    ((b1.x2 - b1.x1) * (b1.y1 - b2.y1) - (b1.y2 - b1.y1) * (b1.x1 - b2.x1)) / denom;

                if EPSILON < t
                    && (t < 1.0 || approx_eq(t, 1.0))
                    && EPSILON < u
                    && (u < 1.0 || approx_eq(u, 1.0))
                {
                    let ix = b1.x1 + t * (b1.x2 - b1.x1);
                    let iy = b1.y1 + t * (b1.y2 - b1.y1);

                    b1.x2 = ix;
                    b1.y2 = iy;
                }
            }
        }
    }

    grew.append(&mut grow(
        &grew
            .iter()
            .enumerate()
            .filter(|(i, b1)| {
                !grew.iter().chain(all_branches).enumerate().any(|(j, b2)| {
                    if *i == j {
                        return false;
                    }

                    let denom =
                        (b2.y2 - b2.y1) * (b1.x2 - b1.x1) - (b2.x2 - b2.x1) * (b1.y2 - b1.y1);

                    if approx_eq(denom, 0.0) {
                        let t_num =
                            (b2.x2 - b2.x1) * (b1.y1 - b2.y1) - (b2.y2 - b2.y1) * (b1.x1 - b2.x1);
                        approx_eq(t_num, 0.0) && colinear_intersecting(b1, b2)
                    } else {
                        let t = ((b2.x2 - b2.x1) * (b1.y1 - b2.y1)
                            - (b2.y2 - b2.y1) * (b1.x1 - b2.x1))
                            / denom;
                        let u = ((b1.x2 - b1.x1) * (b1.y1 - b2.y1)
                            - (b1.y2 - b1.y1) * (b1.x1 - b2.x1))
                            / denom;

                        EPSILON < t
                            && (t < 1.0 || approx_eq(1.0, t))
                            && EPSILON < u
                            && (u < 1.0 || approx_eq(1.0, u))
                    }
                })
            })
            .map(|(_, branch)| branch)
            .cloned()
            .collect::<Vec<_>>(),
        &[grew.as_slice(), all_branches].concat(),
        max_iters,
        iters + 1,
    ));

    grew
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

    branches.append(&mut grow(&branches, &branches, 300, 0));

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
