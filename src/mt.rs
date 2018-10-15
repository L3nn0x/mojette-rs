use super::*;

pub fn direct<T: Default + std::clone::Clone + std::ops::BitXorAssign, S: Support<T>>(support: S, projections: &mut Transform<T>) {
    let P = support.width(); // k
    let Q = support.height(); // l

    for proj in projections.iter_mut() {
        let p = proj.p();
        let q = proj.q();

        let size = p.abs() as usize * (Q - 1) + q as usize * (P - 1) + 1;

        proj.bins = vec![T::default(); size];

        let offset = if p < 0 {
            (Q - 1) as isize * p as isize
        } else { 0 };

        for l in 0..Q {
            for k in 0..P {
                proj.bins[(k as isize * q as isize + (l as isize * p as isize - offset)) as usize] ^= support.get_data(l, k).clone();
            }
        }
    }
}

struct Univoc(usize, usize);

pub fn inverse<T: Default + std::clone::Clone + std::ops::BitXorAssign + std::cmp::PartialEq>(support: &mut Support<T>, mut projections: Transform<T>) {
    let P = support.width(); // k
    let Q = support.height(); // l
    let mut offsets = Vec::new();
    let mut unitary_projs = Vec::new();
    let mut dietmar_projs = Vec::new();

    // Build unitary and dietmar projections.
    for proj in projections.iter() {
        let p = proj.p();
        let q = proj.q();
        let size = proj.bins.len();

        let offset = if p < 0 {
            (Q - 1) as isize * p as isize
        } else { 0 };
        offsets.push(offset);
        unitary_projs.push(Projection::<usize>::new(p, q, size));
        dietmar_projs.push(Projection::<usize>::new(p, q, size));
    }

    let mut dietmar = 0;
    for l in 0..Q {
        for k in 0..P {
            for (n, proj) in projections.iter().enumerate() {
                let index = (k as isize* proj.q() as isize + (l as isize * proj.p() as isize - offsets[n])) as usize;
                unitary_projs[n].bins[index] += 1;
                dietmar_projs[n].bins[index] += dietmar;
            }
            dietmar += 1;
        }
    }

    // Find initial univocs bins.
    let mut list = std::collections::VecDeque::new();
    for (n, proj) in projections.iter().enumerate() {
        for i in 0..proj.bins.len() {
            if unitary_projs[n].bins[i] == 1 {
                list.push_back(Univoc(n, i));
            }
        }
    }

    // Reconstruct.
    while let Some(univoc) = list.pop_front() {
        let dietmar = dietmar_projs[univoc.0].bins[univoc.1];
        let bin = projections[univoc.0].bins[univoc.1].clone();

        if unitary_projs[univoc.0].bins[univoc.1] == 1 {
            // update support
            let l = dietmar / P;
            let k = dietmar - l * P;
            support.set_data(l, k, bin.clone());

            // update projections
            for (n, proj) in projections.iter_mut().enumerate() {
                let index = (k as isize * proj.q() as isize + (l as isize * proj.p() as isize - offsets[n])) as usize;
                proj.bins[index] ^= bin.clone();
                dietmar_projs[n].bins[index] -= dietmar;
                unitary_projs[n].bins[index] -= 1;
                if unitary_projs[n].bins[index] == 1 {
                    list.push_back(Univoc(n, index));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

        #[test]
    fn basic_test() {
        let mut support = RectangularSupport::<usize>::new(3, 3);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;
        support.data[3] = 4;
        support.data[4] = 5;
        support.data[5] = 42;
        support.data[6] = 90;
        support.data[7] = 34;
        support.data[8] = 1;

        let mut transform = Transform::<usize>::new();
        transform.push(Projection::<usize>::new(1, 0, 0));
        transform.push(Projection::<usize>::new(1, 1, 0));
        transform.push(Projection::<usize>::new(0, 1, 0));
        transform.push(Projection::<usize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = RectangularSupport::<usize>::new(3, 3);

        inverse(&mut support, transform);

        let res = vec![3, 8, 2, 4, 5, 42, 90, 34, 1];
        
        assert_eq!(res, support.data);
    }

    #[test]
    fn rectangular_test() {
        let mut support = RectangularSupport::<usize>::new(3, 4);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;
        support.data[3] = 4;
        support.data[4] = 5;
        support.data[5] = 42;
        support.data[6] = 90;
        support.data[7] = 34;
        support.data[8] = 1;
        support.data[9] = std::usize::MAX;
        support.data[10] = 23;
        support.data[11] = 99;

        let mut transform = Transform::<usize>::new();
        transform.push(Projection::<usize>::new(1, 0, 0));
        transform.push(Projection::<usize>::new(1, 1, 0));
        transform.push(Projection::<usize>::new(0, 1, 0));
        transform.push(Projection::<usize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = RectangularSupport::<usize>::new(3, 4);

        inverse(&mut support, transform);

        let res = vec![3, 8, 2, 4, 5, 42, 90, 34, 1, std::usize::MAX, 23, 99];
        
        assert_eq!(res, support.data);
    }

    #[test]
    fn rectangular_test2() {
        let mut support = RectangularSupport::<isize>::new(4, 3);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;
        support.data[3] = 4;
        support.data[4] = 5;
        support.data[5] = 42;
        support.data[6] = 90;
        support.data[7] = 34;
        support.data[8] = 1;
        support.data[9] = 55;
        support.data[10] = 23;
        support.data[11] = 99;

        let mut transform = Transform::<isize>::new();
        transform.push(Projection::<isize>::new(1, 0, 0));
        transform.push(Projection::<isize>::new(1, 1, 0));
        transform.push(Projection::<isize>::new(0, 1, 0));
        transform.push(Projection::<isize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = RectangularSupport::<isize>::new(4, 3);

        inverse(&mut support, transform);

        let res = vec![3, 8, 2, 4, 5, 42, 90, 34, 1, 55, 23, 99];
        
        assert_eq!(res, support.data);
    }

    #[test]
    fn rectangular_test3() {
        let mut support = RectangularSupport::<isize>::new(1, 3);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;

        let mut transform = Transform::<isize>::new();
        transform.push(Projection::<isize>::new(1, 0, 0));
        transform.push(Projection::<isize>::new(1, 1, 0));
        transform.push(Projection::<isize>::new(0, 1, 0));
        transform.push(Projection::<isize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = RectangularSupport::<isize>::new(1, 3);

        inverse(&mut support, transform);

        let res = vec![3, 8, 2];
        
        assert_eq!(res, support.data);
    }

    #[test]
    fn rectangular_test_missing_data() {
        let mut support = RectangularSupport::<isize>::new(1, 3);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;

        let mut transform = Transform::<isize>::new();
        transform.push(Projection::<isize>::new(1, 0, 0));
        transform.push(Projection::<isize>::new(1, 1, 0));
        transform.push(Projection::<isize>::new(0, 1, 0));
        transform.push(Projection::<isize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = RectangularSupport::<isize>::new(1, 3);

        transform.pop();
        transform.pop();

        inverse(&mut support, transform);

        let res = vec![3, 8, 2];
        
        assert_eq!(res, support.data);
    }

    #[test]
    fn linear_test() {
        let mut support = LinearSupport::<isize>::new(3);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;

        let mut transform = Transform::<isize>::new();
        transform.push(Projection::<isize>::new(1, 0, 0));
        transform.push(Projection::<isize>::new(1, 1, 0));
        transform.push(Projection::<isize>::new(0, 1, 0));
        transform.push(Projection::<isize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = LinearSupport::<isize>::new(3);

        inverse(&mut support, transform);

        let res = vec![3, 8, 2];
        
        assert_eq!(res, support.data);
    }

    #[test]
    fn linear_test_missing_data() {
        let mut support = LinearSupport::<isize>::new(3);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;

        let mut transform = Transform::<isize>::new();
        transform.push(Projection::<isize>::new(1, 0, 0));
        transform.push(Projection::<isize>::new(1, 1, 0));
        transform.push(Projection::<isize>::new(0, 1, 0));
        transform.push(Projection::<isize>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = LinearSupport::<isize>::new(3);

        transform.pop();
        transform.pop();

        inverse(&mut support, transform);

        let res = vec![3, 8, 2];
        
        assert_eq!(res, support.data);
    }

    #[test]
    fn linear_test2() {
        let mut support = LinearSupport::<u8>::new(3);
        support.data[0] = 3;
        support.data[1] = 8;
        support.data[2] = 2;

        let mut transform = Transform::<u8>::new();
        transform.push(Projection::<u8>::new(1, 0, 0));
        transform.push(Projection::<u8>::new(1, 1, 0));
        transform.push(Projection::<u8>::new(0, 1, 0));
        transform.push(Projection::<u8>::new(-1, 1, 0));

        direct(support, &mut transform);

        let mut support = LinearSupport::<u8>::new(3);

        inverse(&mut support, transform);

        let res = vec![3, 8, 2];
        
        assert_eq!(res, support.data);
    }
}